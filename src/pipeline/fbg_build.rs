//! `haricot-detect`'s command line, and the graph build + layout outputs of `haricot-layout` (the
//! steps after the forests, in memory), whose graph flags the same code parses.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};

use crate::graph;
use crate::layout;
use crate::overlap::Overlap;
use crate::params::{GraphParams, Params};

/// Mean coverage over `nodes`, rounded to the nearest integer — the honest `dp:i:`
/// for a collapsed segment (unitig or contig), matching the weight `unitig_links`
/// already uses for best-overlap (`layout/unitig.rs:203-210`). Empty => 0.
pub fn mean_cov(g: &graph::FuzzyGraph, nodes: &[(u32, graph::Dir)]) -> u32 {
    if nodes.is_empty() {
        return 0;
    }
    let sum: u64 = nodes
        .iter()
        .map(|&(n, _)| g.nodes()[n as usize].cov as u64)
        .sum();
    ((sum as f64 / nodes.len() as f64).round()) as u32
}

/// A `/proc/self/status` size field in GiB (0.0 if unreadable).
pub fn status_gb(field: &str) -> f64 {
    fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find_map(|l| l.strip_prefix(field))
                .and_then(|r| r.split_whitespace().next().map(str::to_owned))
        })
        .and_then(|kb| kb.parse::<f64>().ok())
        .map(|kb| kb / 1024.0 / 1024.0)
        .unwrap_or(0.0)
}

/// Emit one machine-parseable per-stage profiling line to stderr:
/// `PROFILE<TAB>name<TAB>seconds<TAB>RSS=<live>GB<TAB>HWM=<peak-so-far>GB<TAB>extra`.
///
/// `RSS` is the live VmRSS at the stage boundary, so successive lines trace memory
/// growth. `HWM` is VmHWM — the kernel's high-water mark, i.e. the same quantity
/// `/usr/bin/time -v` reports as "Maximum resident set size", but sampled *per stage*.
/// Without it a stage's true peak is invisible: pod5 reading alone allocates ~32 GB
/// on the full D. melanogaster set, which the allocator does not hand back, so the
/// end-of-run maximum cannot say whether the overlap stage exceeded that floor.
pub fn prof(name: &str, since: Instant, extra: &str) {
    eprintln!(
        "PROFILE\t{name}\t{:.3}\tRSS={:.2}GB\tHWM={:.2}GB\t{extra}",
        since.elapsed().as_secs_f64(),
        status_gb("VmRSS:"),
        status_gb("VmHWM:"),
    );
}

/// The graph build's layout options (`--clean`, `--unitig`, `--dump-paths`). The unitig outputs are
/// always laid out by wtdbg2's fragment layout (`--frg-layout`), and the build always breaks ties by
/// the tie rank (`--canonical-ties`).
#[derive(Default)]
pub struct Layout {
    pub clean: bool,
    /// Also emit the unitig graph (unitigs.gfa + unitig_nodes.tsv).
    pub unitig: bool,
    /// Per-contig ordered node paths (unitig_paths.tsv, contig_paths.tsv); implies `unitig`.
    pub dump_paths: bool,
}

/// `haricot-detect`'s command line.
pub struct Opts {
    pub input: Option<PathBuf>,
    pub outdir: PathBuf,
    pub params: Params,
    pub layout: Layout,
    pub save_overlaps: Option<PathBuf>,
}

/// The candidate detection of `haricot assemble` (its first step) as haricot-detect flags: the Hough
/// overlapper with the z-primary gate (z >= 3) and the null gate (z >= 4), >= 4 aligned bins, the partner
/// pre-filter (nlp2n_1 >= 2), the repeat-mask cap max(2000, 17 x the input's mean seed occurrence),
/// the canonical pair frame, bin 192.
pub const DETECT_FLAGS: &str = "--hough --hough-min-conf 0 --hough-min-z 3 --null-gate --null-z 4 \
     --min-aln-bins 4 --max-seed-occ 2000 --seed-occ-rel 17 --prefilter-nlp 2 --symmetric --bin 192";

/// The graph builds of `haricot assemble` (haricot-layout's): bin 192, K = 4 bins, wtdbg2's regovl,
/// the relative edge gate 0.20, cleaning with
/// wtdbg2's branched bubble popper, the unitig outputs and paths, wtdbg2's fragment layout, ties
/// by the tie rank.
pub const GRAPH_FLAGS: &str = "--bin 192 --node-bins 4 --regovl-bins 1 --min-edge-support 0.20 \
     --clean --unitig --dump-paths --branched-bubbles 2 --frg-layout --canonical-ties";

/// Parse `haricot-detect`'s flags (the program's arguments after its name). Unknown `--` flags are an
/// error; a bare word is the input pod5. `--production` stands for [`DETECT_FLAGS`] then
/// [`GRAPH_FLAGS`], parsed in its place (a later flag overrides it): the detection and graph flags of
/// `haricot assemble`, not its pipeline (no forests, no fill). The graph
/// flags leave the detection unchanged (no overlap step reads a graph parameter).
pub fn parse(args: impl Iterator<Item = String>) -> Result<Opts> {
    let mut args = args.flat_map(|a| {
        if a == "--production" {
            let f = DETECT_FLAGS
                .split_whitespace()
                .chain(GRAPH_FLAGS.split_whitespace());
            f.map(String::from).collect()
        } else {
            vec![a]
        }
    });
    let mut input: Option<PathBuf> = None;
    let mut outdir = PathBuf::from("out_fbg");
    let mut params = Params::r10_4_1();
    let mut do_clean = false;
    let mut do_unitig = false;
    let mut dump_paths = false;
    let mut save_overlaps: Option<PathBuf> = None;

    while let Some(a) = args.next() {
        match a.as_str() {
            "--outdir" => outdir = args.next().context("--outdir")?.into(),
            // Save the candidate cache (read lengths, names, overlaps, content hashes) for the
            // later steps.
            "--save-overlaps" => {
                save_overlaps = Some(args.next().context("--save-overlaps")?.into())
            }
            "--bin" => params.overlap.bin_events = args.next().context("--bin")?.parse()?,
            "--min-aln-bins" => {
                params.overlap.min_aln_bins = args.next().context("--min-aln-bins")?.parse()?
            }
            "--seed-occ-rel" => {
                params.overlap.seed_occ_rel = args.next().context("--seed-occ-rel")?.parse()?
            }
            "--max-seed-occ" => {
                params.overlap.max_seed_occ = args.next().context("--max-seed-occ")?.parse()?
            }
            // The tie rank is the graph build's only tie-break (see `build`).
            "--canonical-ties" => {}
            // The Hough-extent overlapper (dp::hough_aln), the only one.
            "--hough" => {}
            // The Hough confidence gate: off (0), the only setting.
            "--hough-min-conf" => {
                let v: f32 = args.next().context("--hough-min-conf")?.parse()?;
                anyhow::ensure!(
                    v == 0.0,
                    "--hough-min-conf {v}: only 0 (no confidence gate)"
                );
            }
            // Hough SIGNIFICANCE floor on the fit (the z-primary gate). See
            // params::OverlapParams::hough_min_z. 0 = off.
            "--hough-min-z" => {
                params.overlap.hough_min_z = args.next().context("--hough-min-z")?.parse()?
            }
            // The null-hypothesis significance gate and its z threshold.
            "--null-gate" => params.overlap.null_gate = true,
            "--null-z" => params.overlap.null_z = args.next().context("--null-z")?.parse()?,
            // The partner pre-filter (`overlap::prefilter`): align a gathered pair only when
            // its `nlp2n_1` reaches this level. Absent = no pre-filter.
            "--prefilter-nlp" => {
                params.overlap.prefilter_nlp =
                    Some(args.next().context("--prefilter-nlp")?.parse()?)
            }
            // The canonical pair frame (params::OverlapParams::symmetric_pairs): each pair
            // aligned with the longer read as the query; the pre-filter scores both frames.
            "--symmetric" => params.overlap.symmetric_pairs = true,
            // The node size K, in bins.
            "--node-bins" => {
                params.graph.node_bins = args.next().context("--node-bins")?.parse()?
            }
            // Relative edge gate: drop edges below f * min(cov(from), cov(to)).
            "--min-edge-support" => {
                params.graph.min_edge_support =
                    args.next().context("--min-edge-support")?.parse()?
            }
            "--regovl-bins" => {
                params.graph.regovl_bins = args.next().context("--regovl-bins")?.parse()?
            }
            // wtdbg2's branched popper, safe then unsafe: the only cleaning popper.
            "--branched-bubbles" => {
                let v: u8 = args.next().context("--branched-bubbles")?.parse()?;
                anyhow::ensure!(v == 2, "--branched-bubbles {v}: only 2 (safe, then unsafe)");
            }
            // Lay the unitig outputs out through wtdbg2's second level (layout::frg).
            "--frg-layout" => params.graph.frg_layout = true,
            "--clean" => do_clean = true,
            // Detection only: haricot-detect stops after the overlap cache is saved.
            "--no-graph" => {}
            // Also emit the unitig graph (unitigs.gfa + unitig_nodes.tsv).
            "--unitig" => do_unitig = true,
            // Also write each unitig's and contig's ordered node path (unitig_paths.tsv,
            // contig_paths.tsv). Implies --unitig.
            "--dump-paths" => {
                dump_paths = true;
                do_unitig = true;
            }
            s if s.starts_with("--") => anyhow::bail!("unknown flag {s}"),
            _ => input = Some(a.into()),
        }
    }
    Ok(Opts {
        input,
        outdir,
        params,
        layout: Layout {
            clean: do_clean,
            unitig: do_unitig,
            dump_paths,
        },
        save_overlaps,
    })
}

/// A graph directory's read tables: read_names.tsv, and read_content.tsv when the reads carry
/// content hashes.
pub fn write_read_tables(outdir: &Path, names: &[String], content: Option<&[u64]>) -> Result<()> {
    // Read-id -> name table, in sketch order (= `read_id` order). The threads dump and
    // the rescue's tables key reads by `read_id` while truth is keyed by read name, so a
    // join needs this map. Written unconditionally: it is ~1 MB, and reconstructing the
    // order any other way is an assumption rather than a fact.
    {
        let path = outdir.join("read_names.tsv");
        let f = fs::File::create(&path).with_context(|| format!("writing {}", path.display()))?;
        let mut w = std::io::BufWriter::new(f);
        for (i, n) in names.iter().enumerate() {
            writeln!(w, "{i}\t{n}")?;
        }
    }
    // Every read's content hash in hash order (`cache::content_table`), for the end bridges' ties.
    if let Some(c) = content {
        let path = outdir.join("read_content.tsv");
        fs::write(&path, crate::pipeline::cache::content_table(names, c))
            .with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

/// The overlap stats line (stderr): the count and the aligned-span distribution, in bins and in
/// events.
pub fn overlap_stats(overlaps: &[Overlap], params: &Params) {
    {
        let be = params.overlap.bin_events;
        let mut ab: Vec<u32> = overlaps.iter().map(|o| o.aligned_bins).collect();
        if ab.is_empty() {
            eprintln!("overlaps: 0  (bin_events={be})");
        } else {
            ab.sort_unstable();
            let med = ab[ab.len() / 2];
            let p90 = ab[ab.len() * 9 / 10];
            eprintln!(
                "overlaps: {}  | aligned_bins med={} p90={}  | aligned_ev med={} (bin_events={be}, \
                 max_anchor_density={}, min_aln_bins={})",
                overlaps.len(),
                med,
                p90,
                med * be,
                params.overlap.max_anchor_density,
                params.overlap.min_aln_bins,
            );
        }
    }
}

/// The reads of an overlap set: every read's length in retained events, its id and, when known,
/// its signal content hash.
#[derive(Clone, Copy)]
pub struct Reads<'a> {
    pub counts: &'a [u32],
    pub names: &'a [String],
    pub content: Option<&'a [u64]>,
}

/// What a build laid out.
pub struct Built {
    pub unitigs: usize,
    pub contigs: usize,
}

/// The graph build, cleaning and layout over `overlaps`, written to `outdir` (fbg.gfa,
/// fbg_nodes.tsv; with `lay.unitig` unitigs.gfa, unitig_nodes.tsv, contigs.gfa, contig_nodes.tsv;
/// with `lay.dump_paths` unitig_paths.tsv and contig_paths.tsv), the Stage-3 dump to `threads_dump`
/// (`graph::build_graph_dumping`).
pub fn build(
    outdir: &Path,
    gp: &GraphParams,
    lay: &Layout,
    reads: &Reads,
    overlaps: &[Overlap],
    threads_dump: Option<&Path>,
) -> Result<Built> {
    let Reads {
        counts,
        names,
        content,
    } = *reads;
    let t = Instant::now();
    // --canonical-ties: equal-coverage candidates commit by their anchor's tie rank (the order of the
    // reads' signal content hashes, readorder::rank_by_content) instead of its number
    // (graph::build::build_graph_annotated_ranked), so the graph depends neither on the read
    // numbering nor on the read names.
    let c = content.context(
        "--canonical-ties needs a cache with read content hashes (re-run the detection with \
         this build)",
    )?;
    let rank = Some(crate::readorder::rank_by_content(c)?);
    let (mut g, annos) = graph::build_graph_dumping(overlaps, gp, rank.as_deref(), threads_dump);
    prof(
        "build_graph",
        t,
        &format!("nodes={} edges={}", g.n_nodes(), g.n_edges()),
    );
    eprintln!("FBG: {} nodes, {} edges", g.n_nodes(), g.n_edges());
    if lay.clean {
        let t = Instant::now();
        graph::clean(&mut g, gp);
        prof(
            "clean",
            t,
            &format!("live_nodes={} edges={}", g.n_live_nodes(), g.n_edges()),
        );
        eprintln!(
            "after clean: {} live nodes, {} edges",
            g.n_live_nodes(),
            g.n_edges()
        );
    }

    let t_write = Instant::now();
    let mut gfa = BufWriter::new(File::create(outdir.join("fbg.gfa"))?);
    layout::write_gfa(&g, &mut gfa)?;
    gfa.flush()?;

    let mut tsv = BufWriter::new(File::create(outdir.join("fbg_nodes.tsv"))?);
    writeln!(
        tsv,
        "node\tanchor_read\tev_beg\tev_end\tcov\tanchor_coord_events"
    )?;
    for (i, n) in g.nodes().iter().enumerate() {
        if n.del {
            continue; // tombstoned by cleaning — not in the GFA either
        }
        let a = annos[i];
        writeln!(
            tsv,
            "n{}\t{}\t{}\t{}\t{}\t{}",
            n.id,
            names[a.anchor_read as usize],
            a.ev_beg,
            a.ev_end,
            n.cov,
            counts[a.anchor_read as usize]
        )?;
    }
    tsv.flush()?;
    prof("write_fbg", t_write, "fbg.gfa + fbg_nodes.tsv");
    eprintln!(
        "wrote {}/fbg.gfa and {}/fbg_nodes.tsv",
        outdir.display(),
        outdir.display()
    );

    let mut built = Built {
        unitigs: 0,
        contigs: 0,
    };
    if lay.unitig {
        let t_unitig = Instant::now();
        // Collapse the (cleaned) graph into unitigs and write the unitig graph: one
        // segment per unitig, one link per FBG edge that is not internal to a unitig
        // (`layout::unitig_links`). Provenance = the anchor read of each unitig's
        // first node, so it can be Bandage-colored the same way as the FBG.
        //
        // Both the collapse and the linking live in `layout::unitig` so they are
        // unit-tested. A "two unitigs link if they share an endpoint node" rule would
        // link all k unitigs at a junction pairwise (a k-clique where the truth is a
        // star) and invent cycles. Verify any change here against `gfatools asm -u`,
        // which must agree on component and cycle counts.
        let unitigs = {
            let (us, st) = layout::frg_contigs(&g, &layout::FrgParams::default());
            eprintln!(
                "FRG layout: {} branching nodes masked, {} fragments, {} read links, {} weak cut, \
                 {} transitive, {} tips -> {} contigs ({} joined)",
                st.masked,
                st.fragments,
                st.links,
                st.weak_cut,
                st.transitive,
                st.tips,
                st.contigs,
                st.joined
            );
            us
        };
        built.unitigs = unitigs.len();
        let ulinks = layout::unitig_links(&g, &unitigs);
        let mut ug = BufWriter::new(File::create(outdir.join("unitigs.gfa"))?);
        layout::write_header(&mut ug)?;
        // dp:i: must be the unitig's MEAN NODE COVERAGE: Bandage renders it as
        // "depth" (the node count there would show a 267-node unitig as 267x).
        for (ui, u) in unitigs.iter().enumerate() {
            let cov = mean_cov(&g, &u.nodes);
            layout::write_segment(&mut ug, &format!("u{ui}"), u.len_events, Some(cov))?;
        }
        for l in &ulinks {
            layout::write_link(
                &mut ug,
                &format!("u{}", l.from),
                l.from_dir.sign(),
                &format!("u{}", l.to),
                l.to_dir.sign(),
                l.overlap_events,
                None,
                None,
            )?;
        }
        ug.flush()?;
        // Leading columns match fbg_nodes.tsv exactly (node, anchor_read, ev_beg,
        // ev_end, cov, anchor_coord_events) so a unitig can be placed and coloured the
        // same way as an FBG node; n_nodes/len_events are appended.
        //
        // CAVEAT: the placement is the unitig's FIRST node — a point, not a span.
        // A unitig's nodes may be keyed on different anchor reads, so there is no
        // single read to interpolate a whole unitig within. Good enough to colour by
        // locus; do NOT compute union coverage from this file.
        let mut utsv = BufWriter::new(File::create(outdir.join("unitig_nodes.tsv"))?);
        writeln!(
            utsv,
            "node\tanchor_read\tev_beg\tev_end\tcov\tanchor_coord_events\tn_nodes\tlen_events"
        )?;
        for (ui, u) in unitigs.iter().enumerate() {
            let rep = u.nodes.first().map(|&(n, _)| n).unwrap_or(0) as usize;
            let a = &annos[rep];
            writeln!(
                utsv,
                "u{ui}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                names[a.anchor_read as usize],
                a.ev_beg,
                a.ev_end,
                g.nodes()[rep].cov,
                counts[a.anchor_read as usize],
                u.nodes.len(),
                u.len_events
            )?;
        }
        utsv.flush()?;
        prof(
            "layout_unitig",
            t_unitig,
            &format!("unitigs={}", unitigs.len()),
        );
        eprintln!(
            "wrote {}/unitigs.gfa ({} unitigs)",
            outdir.display(),
            unitigs.len()
        );

        let t_contig = Instant::now();
        // Best-overlap contig layout: each contig = one segment (length in events),
        // colored by its first node's anchor read. Contigs are the linearised output
        // (the dominant tiling of each cluster → one long contig).
        let contigs = layout::gen_contigs_from(&g, &unitigs);
        built.contigs = contigs.len();
        let mut cg = BufWriter::new(File::create(outdir.join("contigs.gfa"))?);
        layout::write_header(&mut cg)?;
        let mut ctsv = BufWriter::new(File::create(outdir.join("contig_nodes.tsv"))?);
        writeln!(ctsv, "node\tanchor_read\tn_unitigs\tlen_events")?;
        let mut lens: Vec<u32> = Vec::new();
        for (ci, c) in contigs.iter().enumerate() {
            // Mean node coverage over every node of every unitig in the contig — see
            // the unitigs.gfa note above.
            let nodes: Vec<_> = c
                .unitigs
                .iter()
                .flat_map(|&ui| unitigs[ui].nodes.iter().copied())
                .collect();
            let cov = mean_cov(&g, &nodes);
            layout::write_segment(&mut cg, &format!("c{ci}"), c.len_events, Some(cov))?;
            // representative anchor read = first node of the contig's first unitig
            let first_node = unitigs[c.unitigs[0]].nodes[0].0 as usize;
            let read = &names[annos[first_node].anchor_read as usize];
            writeln!(ctsv, "c{ci}\t{read}\t{}\t{}", c.unitigs.len(), c.len_events)?;
            lens.push(c.len_events);
        }
        cg.flush()?;
        ctsv.flush()?;
        lens.sort_unstable_by(|a, b| b.cmp(a));
        let big = lens.iter().filter(|&&l| l >= 5000).count();
        prof(
            "layout_contig",
            t_contig,
            &format!("contigs={}", contigs.len()),
        );
        eprintln!(
            "wrote {}/contigs.gfa ({} contigs, {} ≥5kev; top lens {:?})",
            outdir.display(),
            contigs.len(),
            big,
            &lens[..lens.len().min(6)]
        );

        // Each unitig's and contig's ordered node path, from which a sequence can be stitched
        // (`evaluation/anchor_seq.py`). Row: contig (the unitig or contig name), step,
        // anchor_read, ev_beg, ev_end (in that read's events), dir, advance_events, coord_events.
        //
        // `dir` is the node's orientation along the path ('+' = traversed L end to R
        // end) — a stitcher must know which way to lay the node's events down.
        //
        // `advance_events` is THE column that makes stitching correct: how far this
        // node advances before the next one starts, i.e. `node_len + offset` of the
        // edge to the next step (`clean::step_len`). Consecutive nodes overlap by
        // ~regovl (96 events at bin96), so emitting each node's full ev_beg..ev_end
        // would re-emit the shared bins and inflate every contig by ~384/288 = 1.33x
        // (duplicated sequence). A stitcher should emit only
        // `advance_events` worth of each node, and the full node length for the last
        // step of a contig. Consecutive steps have different anchor reads, so this
        // cannot be recovered downstream by comparing coordinates.
        //
        // `coord_events` is the anchor read's total retained-event count. An exact map
        // from events to bases needs the basecaller's move table
        // (`evaluation/anchor_seq.py`); the uniform rate `ev / coord_events * read_len`
        // is off by up to about one node width.
        if lay.dump_paths {
            // (from,to) -> offset, both directions: a path step may traverse an edge
            // either way.
            let mut off: std::collections::HashMap<(u32, u32), i32> =
                std::collections::HashMap::new();
            for e in g.edges() {
                off.insert((e.from, e.to), e.offset);
                off.insert((e.to, e.from), e.offset);
            }
            // unitig_paths.tsv: one row per step of every unitig (named in the `contig`
            // column). A unitig's node order and orientation are unambiguous, so its path can
            // be stitched into a sequence; it is the path file the evaluation reads.
            //
            // contig_paths.tsv records no unitig orientation: `Contig` stores only
            // `unitigs: Vec<usize>` (`gen_contigs` tracks which END it enters each unitig by
            // and discards it), so a unitig entered at its right end is still written forward
            // there. Its paths are therefore not stitchable into a sequence.
            let mut uf = BufWriter::new(File::create(outdir.join("unitig_paths.tsv"))?);
            writeln!(
                uf,
                "contig\tstep\tanchor_read\tev_beg\tev_end\tdir\tadvance_events\tcoord_events"
            )?;
            for (ui, u) in unitigs.iter().enumerate() {
                for (step, &(nid, dir)) in u.nodes.iter().enumerate() {
                    let a = annos[nid as usize];
                    let len = g.nodes()[nid as usize].len_events as i64;
                    let advance = match u.nodes.get(step + 1) {
                        Some(&(next, _)) => off
                            .get(&(nid, next))
                            .map(|&o| (len + o as i64).max(1))
                            .unwrap_or(len),
                        None => len,
                    };
                    writeln!(
                        uf,
                        "u{ui}\t{step}\t{}\t{}\t{}\t{}\t{}\t{}",
                        names[a.anchor_read as usize],
                        a.ev_beg,
                        a.ev_end,
                        dir.sign(),
                        advance,
                        counts[a.anchor_read as usize]
                    )?;
                }
            }
            uf.flush()?;
            eprintln!("wrote {}/unitig_paths.tsv", outdir.display());

            let mut pf = BufWriter::new(File::create(outdir.join("contig_paths.tsv"))?);
            writeln!(
                pf,
                "contig\tstep\tanchor_read\tev_beg\tev_end\tdir\tadvance_events\tcoord_events"
            )?;
            for (ci, c) in contigs.iter().enumerate() {
                let path: Vec<(u32, crate::graph::Dir)> = c
                    .unitigs
                    .iter()
                    .flat_map(|&ui| unitigs[ui].nodes.clone())
                    .collect();
                for (step, &(nid, dir)) in path.iter().enumerate() {
                    let a = annos[nid as usize];
                    let len = g.nodes()[nid as usize].len_events as i64;
                    // Advance to the next step; the last node contributes its whole
                    // length. A missing edge (unitigs chained by gen_contigs without a
                    // direct link) also falls back to the full length.
                    let advance = match path.get(step + 1) {
                        Some(&(next, _)) => off
                            .get(&(nid, next))
                            .map(|&o| (len + o as i64).max(1))
                            .unwrap_or(len),
                        None => len,
                    };
                    writeln!(
                        pf,
                        "c{ci}\t{step}\t{}\t{}\t{}\t{}\t{}\t{}",
                        names[a.anchor_read as usize],
                        a.ev_beg,
                        a.ev_end,
                        dir.sign(),
                        advance,
                        counts[a.anchor_read as usize]
                    )?;
                }
            }
            pf.flush()?;
            eprintln!("wrote {}/contig_paths.tsv", outdir.display());
        }
    }
    Ok(built)
}
