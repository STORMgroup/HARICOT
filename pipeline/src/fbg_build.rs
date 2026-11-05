//! `haricot-detect`'s command line and its graph build + layout outputs, shared with `haricot-layout` (the
//! steps after the forests, in memory; `#[path]`-included like `cache.rs`): the same flags parsed
//! by the same code, the same build and the same files whether the overlaps come from a cache
//! (`--load-overlaps`), from the detection, or from memory.
#![allow(dead_code)]

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};

use haricot::graph;
use haricot::layout;
use haricot::overlap::Overlap;
use haricot::params::{GraphParams, Params, QuantMethod, SketchMethod};

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

/// The graph build's layout options (`--clean`, `--unzip`, `--unitig`, `--dump-paths`,
/// `--frg-layout` / `--frg-weak` / `--frg-links`, `--canonical-ties`).
#[derive(Default)]
pub struct Layout {
    pub clean: bool,
    pub unzip: bool,
    /// Also emit the unitig graph (unitigs.gfa + unitig_nodes.tsv).
    pub unitig: bool,
    /// Per-contig ordered node paths (unitig_paths.tsv, contig_paths.tsv); implies `unitig`.
    pub dump_paths: bool,
    /// --frg-layout: lay the unitig outputs out through wtdbg2's second level (layout::frg).
    pub frg: Option<layout::FrgParams>,
    /// --canonical-ties: the graph build's tie-break by tie rank, the signal content ([`build`]).
    pub canonical_ties: bool,
}

/// `haricot-detect`'s command line.
pub struct Opts {
    pub input: Option<PathBuf>,
    pub outdir: PathBuf,
    pub params: Params,
    pub layout: Layout,
    /// --no-graph: detection only; stop after the overlap cache is saved.
    pub no_graph: bool,
    pub ovl_min_z: Option<f32>,
    pub save_overlaps: Option<PathBuf>,
    pub load_overlaps: Option<PathBuf>,
    pub max_reads: Option<usize>,
    /// --canonical-order: number the reads by read id, not by their place in the pod5
    /// (haricot::readorder), so the candidates are a function of the read set alone.
    pub canonical_order: bool,
    pub stream_sketch: bool,
    pub stream_chunk_bytes: usize,
}

/// The candidate detection of `haricot assemble` (its first step) as haricot-detect flags: the Hough
/// overlapper with the z-primary gate (z >= 3) and the null gate (z >= 4), >= 4 aligned bins, the partner
/// pre-filter (nlp2n_1 >= 2), the repeat-mask cap max(2000, 17 x the input's mean seed occurrence),
/// the canonical pair frame, bin 192.
pub const DETECT_FLAGS: &str = "--hough --hough-min-conf 0 --hough-min-z 3 --null-gate --null-z 4 \
     --min-aln-bins 4 --max-seed-occ 2000 --seed-occ-rel 17 --prefilter-nlp 2 --symmetric --bin 192";

/// The graph builds of `haricot assemble` (haricot-layout's default `--graph`, before policy7's read
/// threading): bin 192, K = 4 bins, wtdbg2's regovl, the relative edge gate 0.20, cleaning with
/// wtdbg2's branched bubble popper, the unitig outputs and paths, wtdbg2's fragment layout, ties
/// by the tie rank.
pub const GRAPH_FLAGS: &str = "--bin 192 --node-bins 4 --regovl-bins 1 --min-edge-support 0.20 \
     --clean --unitig --dump-paths --branched-bubbles 2 --frg-layout --canonical-ties";

/// Parse `haricot-detect`'s flags (the program's arguments after its name). Unknown `--` flags are an
/// error; a bare word is the input pod5. `--production` stands for [`DETECT_FLAGS`] then
/// [`GRAPH_FLAGS`], parsed in its place (a later flag overrides it): the detection and graph flags of
/// `haricot assemble`, not its pipeline (no forests, no fill, no read threading). The graph
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
    // --no-graph: detection only; stop after the overlap cache is saved.
    let mut no_graph = false;
    let mut do_unitig = false;
    let mut dump_paths = false;
    let mut unzip = false;
    // --frg-layout: lay the unitig outputs out through wtdbg2's second level (layout::frg).
    let mut frg: Option<layout::FrgParams> = None;
    let mut ovl_min_z: Option<f32> = None;
    let mut save_overlaps: Option<PathBuf> = None;
    let mut load_overlaps: Option<PathBuf> = None;
    let mut max_reads: Option<usize> = None;
    // --canonical-order: number the reads by read id, not by their place in the pod5
    // (haricot::readorder), so the candidates are a function of the read set alone.
    let mut canonical_order = false;
    // --canonical-ties: the graph build's tie-break by tie rank, the signal content (see below).
    let mut canonical_ties = false;
    // Stream the read->sketch pass by default (bounds the raw signal's share of
    // peak RSS); --no-stream-sketch restores the load-everything-then-sketch form.
    let mut stream_sketch = true;
    // 4 GiB of signal in flight. 4 GiB was faster than 1 GiB on every dataset measured,
    // and faster than not streaming on two bacterial both-strand sets, because sketching a
    // chunk whose signal is still cache-resident more than repays the lost read-ahead once
    // chunks are large. On a D. melanogaster set it cost 24 s of 10 min and no memory: the
    // seed store dominates that peak, not the signal.
    let mut stream_chunk_bytes: usize = 4 << 30;

    // Deferred: depends on bin_events, which a later arg may still change.
    let mut pending_node_events: Option<u32> = None;
    while let Some(a) = args.next() {
        match a.as_str() {
            "--outdir" => outdir = args.next().context("--outdir")?.into(),
            // Cache the (counts, names, overlaps) so later builds skip sketch+overlap.
            "--save-overlaps" => {
                save_overlaps = Some(args.next().context("--save-overlaps")?.into())
            }
            // Load a cached overlap set and go straight to the build (no pod5 needed).
            "--load-overlaps" => {
                load_overlaps = Some(args.next().context("--load-overlaps")?.into())
            }
            "--dense" => params.sketch.method = SketchMethod::Dense,
            // Re-threshold a cache offline. The anchor work is identical for every
            // gate, so one permissive pass plus this flag replaces one expensive pass
            // per operating point. NaN (synthesized overlaps) always passes.
            "--ovl-min-z" => ovl_min_z = Some(args.next().context("--ovl-min-z")?.parse()?),
            // Alphabet. `--n-buckets` widens `q` to fit, since `Params::validate`
            // requires `n_buckets <= 2^q` and `q` is the e-mer packing width.
            "--n-buckets" => {
                let v: u32 = args.next().context("--n-buckets")?.parse()?;
                params.quant.n_buckets = v;
                let need = u32::BITS - (v.max(1) - 1).leading_zeros();
                if params.quant.q < need {
                    params.quant.q = need;
                }
            }
            "--quant" => {
                let v = args.next().context("--quant")?;
                params.quant.method = match v.as_str() {
                    "piecewise" => QuantMethod::Piecewise,
                    "normal-quantile" => QuantMethod::NormalQuantile,
                    other => anyhow::bail!("unknown --quant {other} (piecewise|normal-quantile)"),
                };
            }
            // Pre-quantization float event merge (z-units); ~0.3 = aggressive. Only
            // bites when overlaps are (re)computed from pod5, not on --load-overlaps.
            "--event-merge" => {
                params.signal.event_merge = args.next().context("--event-merge")?.parse()?
            }
            "--bin" => params.overlap.bin_events = args.next().context("--bin")?.parse()?,
            // Bin-invariance knobs: rescale the acceptance gates when coarsening --bin
            // so the test isolates bin granularity from the per-bin gate thresholds.
            "--max-anchor-density" => {
                params.overlap.max_anchor_density =
                    args.next().context("--max-anchor-density")?.parse()?
            }
            "--min-aln-bins" => {
                params.overlap.min_aln_bins = args.next().context("--min-aln-bins")?.parse()?
            }
            "--max-gap" => params.overlap.max_gap = args.next().context("--max-gap")?.parse()?,
            "--max-similarity" => {
                params.overlap.max_similarity = args.next().context("--max-similarity")?.parse()?
            }
            "--seed-occ-rel" => {
                params.overlap.seed_occ_rel = args.next().context("--seed-occ-rel")?.parse()?
            }
            "--max-seed-occ" => {
                params.overlap.max_seed_occ = args.next().context("--max-seed-occ")?.parse()?
            }
            // Read-count ladder: keep only the first N reads of the pod5 (pod5 order is
            // by read id, not genome position, so a prefix is spatially unbiased — but
            // it IS a lower-depth dataset, exactly as a smaller real run would be).
            "--max-reads" => max_reads = Some(args.next().context("--max-reads")?.parse()?),
            "--canonical-order" => canonical_order = true,
            "--canonical-ties" => canonical_ties = true,
            "--no-stream-sketch" => stream_sketch = false,
            "--stream-chunk-mb" => {
                stream_chunk_bytes =
                    args.next().context("--stream-chunk-mb")?.parse::<usize>()? << 20
            }
            // DP band / gap / penalty knobs (the aligned-span truncation levers).
            "--max-bvar" => params.overlap.max_bvar = args.next().context("--max-bvar")?.parse()?,
            "--max-bgap" => params.overlap.max_bgap = args.next().context("--max-bgap")?.parse()?,
            "--pgap" => params.overlap.pgap = args.next().context("--pgap")?.parse()?,
            "--pvar" => params.overlap.pvar = args.next().context("--pvar")?.parse()?,
            // Use the Hough-extent overlapper (dp::hough_aln) instead of the bin-DP.
            "--hough" => params.overlap.hough = true,
            "--hough-min-conf" => {
                params.overlap.hough_min_conf = args.next().context("--hough-min-conf")?.parse()?
            }
            // Hough SIGNIFICANCE floor, at the same early exit as --hough-min-conf and
            // on the same fit. `--hough-min-conf 0 --hough-min-z X` is the z-primary
            // gate. See params::OverlapParams::hough_min_z. 0 = off.
            "--hough-min-z" => {
                params.overlap.hough_min_z = args.next().context("--hough-min-z")?.parse()?
            }
            // Opt-in null-hypothesis significance gate (overlap-precision lever) and
            // its z threshold. Only bites when overlaps are (re)computed from pod5,
            // not on --load-overlaps (the gate acts during overlap computation).
            "--null-gate" => params.overlap.null_gate = true,
            "--null-z" => params.overlap.null_z = args.next().context("--null-z")?.parse()?,
            // Opt-in partner pre-filter (`overlap::prefilter`): align a gathered pair only
            // when its owner-frame `nlp2n_1` reaches this level. Absent = no pre-filter.
            "--prefilter-nlp" => {
                params.overlap.prefilter_nlp =
                    Some(args.next().context("--prefilter-nlp")?.parse()?)
            }
            // Opt-in canonical pair frame (params::OverlapParams::symmetric_pairs): each pair
            // aligned with the longer read as the query; the pre-filter scores both frames.
            "--symmetric" => params.overlap.symmetric_pairs = true,
            // Slope-beam knobs (see params::BeamCfg).
            "--slope-win-lo" => {
                params.overlap.slope_beam.win_lo = args.next().context("--slope-win-lo")?.parse()?
            }
            "--slope-win-hi" => {
                params.overlap.slope_beam.win_hi = args.next().context("--slope-win-hi")?.parse()?
            }
            "--slope-coarse-step" => {
                params.overlap.slope_beam.coarse_step =
                    args.next().context("--slope-coarse-step")?.parse()?
            }
            "--slope-k" => {
                params.overlap.slope_beam.beam_k = args.next().context("--slope-k")?.parse()?
            }
            "--slope-levels" => {
                params.overlap.slope_beam.levels = args.next().context("--slope-levels")?.parse()?
            }
            // First filter: floor on the level-0 peak (params::BeamCfg::screen_peak).
            "--hough-screen" => {
                params.overlap.slope_beam.screen_peak =
                    args.next().context("--hough-screen")?.parse()?
            }
            // Density-normalised first filter (params::BeamCfg::screen_z).
            "--hough-screen-z" => {
                params.overlap.slope_beam.screen_z =
                    args.next().context("--hough-screen-z")?.parse()?
            }
            // Stage-0b filter against the gate's null (params::BeamCfg::screen_null).
            "--hough-screen-null" => {
                params.overlap.slope_beam.screen_null =
                    args.next().context("--hough-screen-null")?.parse()?
            }
            // Pick the slope by the gate's z among the top-k basins (params::BeamCfg::select_z).
            "--slope-select-z" => params.overlap.slope_beam.select_z = true,
            // Coarsen the node (K-bin) size, in BINS; the dominant de-tangling lever.
            "--node-bins" => {
                params.graph.node_bins = args.next().context("--node-bins")?.parse()?
            }
            // Deprecated: node size is defined in bins. Converted AFTER the whole arg
            // loop — it depends on bin_events, and this is a sequential parser, so
            // converting here would read whatever --bin had been seen SO FAR.
            // `--node-events 360 --bin 96` then silently yielded 1440-event nodes (it
            // divided by the default 24, giving K=15, which --bin 96 later scaled).
            // cli.rs is immune because clap collects every flag before we read any.
            "--node-events" => {
                pending_node_events = Some(args.next().context("--node-events")?.parse()?)
            }
            "--min-node-cov" => {
                params.graph.min_node_cov = args.next().context("--min-node-cov")?.parse()?
            }
            "--min-edge-cov" => {
                params.graph.min_edge_cov = args.next().context("--min-edge-cov")?.parse()?
            }
            // Relative edge gate: drop edges below f * min(cov(from), cov(to)).
            "--min-edge-support" => {
                params.graph.min_edge_support =
                    args.next().context("--min-edge-support")?.parse()?
            }
            "--node-drop" => {
                params.graph.node_drop = args.next().context("--node-drop")?.parse()?
            }
            "--regovl-bins" => {
                params.graph.regovl_bins = args.next().context("--regovl-bins")?.parse()?
            }
            // Cleaning reach, in nodes (defaults 10 / 40).
            "--tip-steps" => {
                params.graph.tip_steps = args.next().context("--tip-steps")?.parse()?
            }
            "--bubble-steps" => {
                params.graph.bubble_steps = args.next().context("--bubble-steps")?.parse()?
            }
            "--branched-bubbles" => {
                params.graph.branched_bubbles =
                    args.next().context("--branched-bubbles")?.parse()?
            }
            "--read-threading" => params.graph.read_threading = true,
            "--unzip" => unzip = true,
            "--frg-layout" => {
                params.graph.frg_layout = true;
                frg.get_or_insert_with(layout::FrgParams::default);
            }
            "--frg-weak" => {
                frg.get_or_insert_with(layout::FrgParams::default).weak =
                    args.next().context("--frg-weak")?.parse()?
            }
            "--frg-links" => {
                frg.get_or_insert_with(layout::FrgParams::default).links =
                    match args.next().context("--frg-links")?.as_str() {
                        "firstlast" => layout::LinkMode::FirstLast,
                        "consecutive" => layout::LinkMode::Consecutive,
                        x => anyhow::bail!("--frg-links firstlast|consecutive, not {x}"),
                    }
            }
            "--thread-min-reads" => {
                params.graph.thread_min_reads =
                    args.next().context("--thread-min-reads")?.parse()?
            }
            "--clean" => do_clean = true,
            "--no-graph" => no_graph = true,
            // Also emit the unitig graph (unitigs.gfa + unitig_nodes.tsv).
            "--unitig" => do_unitig = true,
            // Dump per-contig ordered node paths (contig_paths.tsv) for layout->FASTA
            // stitching (miniasm-style QUAST eval). Implies --unitig.
            "--dump-paths" => {
                dump_paths = true;
                do_unitig = true;
            }
            s if s.starts_with("--") => anyhow::bail!("unknown flag {s}"),
            _ => input = Some(a.into()),
        }
    }
    // Now that every arg is parsed, bin_events is final. Round up to a whole bin
    // like wtdbg2 (`wtdbg.c:638-640`) rather than silently rounding to nearest.
    if let Some(v) = pending_node_events {
        let bin = params.overlap.bin_events.max(1);
        params.graph.node_bins = v.div_ceil(bin).max(1);
        eprintln!(
            "warning: --node-events is deprecated; {v} events / {bin} per bin -> \
             --node-bins {} ({} events)",
            params.graph.node_bins,
            params.graph.node_bins * bin
        );
    }
    Ok(Opts {
        input,
        outdir,
        params,
        layout: Layout {
            clean: do_clean,
            unzip,
            unitig: do_unitig,
            dump_paths,
            frg,
            canonical_ties,
        },
        no_graph,
        ovl_min_z,
        save_overlaps,
        load_overlaps,
        max_reads,
        canonical_order,
        stream_sketch,
        stream_chunk_bytes,
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
        fs::write(&path, crate::cache::content_table(names, c))
            .with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

/// The overlap stats line (stderr).
pub fn overlap_stats(overlaps: &[Overlap], params: &Params) {
    // Overlap stats: count + aligned-span distribution (in bins AND events), so a
    // bin-size sweep can compare span and not just final coverage.
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
    // numbering nor on the read names. Off: ties by read number.
    let rank = if lay.canonical_ties {
        let c = content.context(
            "--canonical-ties needs a cache with read content hashes (re-run the detection with \
             this build)",
        )?;
        Some(haricot::readorder::rank_by_content(c)?)
    } else {
        None
    };
    let (mut g, mut annos) =
        graph::build_graph_dumping(overlaps, gp, rank.as_deref(), threads_dump);
    prof(
        "build_graph",
        t,
        &format!("nodes={} edges={}", g.n_nodes(), g.n_edges()),
    );
    eprintln!("FBG: {} nodes, {} edges", g.n_nodes(), g.n_edges());
    // How much of the branching could read threading actually resolve? Reported before
    // AND after cleaning, because cleaning changes the branch set.
    let thread_stats = |g: &graph::fbg::FuzzyGraph, when: &str| {
        if !gp.read_threading {
            return;
        }
        let s = graph::threading::branch_stats(g, g.transitions(), gp.thread_min_reads);
        eprintln!(
            "THREADING {when}: triples={} branch_ends={} entries={} \
             resolvable={} ({:.1}%) ambiguous={} ({:.1}%) dead={} ({:.1}%)",
            g.transitions().len(),
            s.branch_ends,
            s.entries,
            s.resolvable,
            100.0 * s.resolvable as f64 / s.entries.max(1) as f64,
            s.ambiguous,
            100.0 * s.ambiguous as f64 / s.entries.max(1) as f64,
            s.dead,
            100.0 * s.dead as f64 / s.entries.max(1) as f64,
        );
    };
    thread_stats(&g, "raw");
    if lay.clean {
        let t = Instant::now();
        graph::clean(&mut g, gp);
        thread_stats(&g, "clean");
        // Unzip AFTER cleaning: the edge set has to be final, and transitions are keyed
        // on node ids, which cleaning soft-deletes but never renumbers.
        if lay.unzip {
            let st = graph::threading::unzip(&mut g, &mut annos, gp.thread_min_reads);
            eprintln!(
                "UNZIP: split {} repeat nodes into {} copies ({} consumed); \
                 graph now {} nodes / {} edges",
                st.nodes,
                st.copies,
                st.consumed,
                g.live_nodes().count(),
                g.n_edges()
            );
            thread_stats(&g, "unzip");
        }
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
        let unitigs = match lay.frg {
            Some(fp) => {
                let (us, st) = layout::frg_contigs(&g, &fp);
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
            }
            None => layout::gen_unitigs_threaded(&g, gp.thread_min_reads),
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

        // Per-contig ordered node path for miniasm-style layout->FASTA stitching.
        // Row: contig, step, anchor_read, ev_beg, ev_end (in that read's events),
        // dir, advance_events, coord_events.
        //
        // `dir` is the node's orientation along the path ('+' = traversed L end to R
        // end) — a stitcher must know which way to lay the node's events down.
        //
        // `advance_events` is THE column that makes stitching correct: how far this
        // node advances before the next one starts, i.e. `node_len + offset` of the
        // edge to the next step (`clean::step_len`). Consecutive nodes overlap by
        // ~regovl (96 events at bin96), so emitting each node's full ev_beg..ev_end
        // would re-emit the shared bins and inflate every contig by ~384/288 = 1.33x
        // — which QUAST would score as duplication. A stitcher should emit only
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
            // Emit the UNITIG paths too. A unitig is a maximal non-branching path —
            // its node order and orientation are unambiguous, so stitching it is
            // sound. Prefer unitig_paths.tsv for any QUAST-style eval.
            //
            // contig_paths.tsv is NOT trustworthy: `Contig` stores only
            // `unitigs: Vec<usize>` — `gen_contigs` tracks which END it enters each
            // unitig by and then discards it — so a unitig entered at its right end
            // is still laid down forward here. Fixing that needs an orientation on
            // `Contig`. Symptom: consecutive steps land megabases apart on the
            // reference, and QUAST genome fraction collapses to ~29%.
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
                let path: Vec<(u32, haricot::graph::Dir)> = c
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
