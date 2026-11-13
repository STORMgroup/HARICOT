//! `haricot-layout`: the steps of `haricot assemble` after the overlap forests, in two runs that each
//! hold the candidate overlaps in memory.
//!
//! * `layout`: the depth-adaptive keep cut (the depth of the overlaps the stage-2 forest keeps at
//!   0.995 picks the stage-2 keep level, Algorithm S2 lines 21-24), the graph build of the kept
//!   overlaps (W/pol7, W/thr.tsv), the end bridges (W/q16_mbbc.drop.tsv) and the build with them
//!   (W/eb, W/thr_eb.tsv);
//! * `final`: the contig-end rescue's admission (W/admit.tsv) and the final build into OUT.
//!
//! Between the two, `haricot assemble` runs the contig-end rescue (`haricot-rescue rescue`, into
//! W/q9c) and the rescue forest, which read the W/eb build.
//!
//! The steps run the code in `src/pipeline/`: `keep_cut.rs` (drop lists, depth estimate, keep
//! rule, candidate dump), `fill.rs` (each path straightened and filled), `fbg_build.rs` (graph
//! flags, build and layout files), `endbridge.rs` (end bridges) and `admit.rs` (the rescue's
//! admission). Each step works on
//! the overlaps in memory, which are the values an intermediate cache would hold (a cache's bincode
//! round trip is lossless). Outputs in the work directory W: depth.txt, pol7/, thr.tsv,
//! q16_mbbc.drop.tsv, eb/, thr_eb.tsv and admit.tsv; the final build goes to OUT (with
//! threads.tsv). The intermediate caches no later step reads (s2.bin, f.bin, ff.bin, cand.dump.tsv,
//! keb.bin, q9c.bin) are written only with `--keep-intermediates`.
//!
//! Usage:
//!   haricot-layout layout <ovl.bin> --work W --key KEY --stage2-drop PREFIX --scores S1.scores.tsv
//!       [--keep-intermediates] [--keep-level depth-adaptive]
//!   haricot-layout final <ovl.bin> --work W --out OUT [--keep-intermediates]
//!
//! PREFIX: the stage-2 drop lists are PREFIX0.995.drop.tsv, PREFIX0.990.drop.tsv and
//! PREFIX0.980.drop.tsv (`haricot-forest`'s). Every graph build takes `fbg_build::GRAPH_FLAGS`. The
//! fill and the rescue's admission work in the canonical pair frame, and the end bridges break
//! score ties by the reads' content hashes. Unknown flags are an error. `--keep-level
//! depth-adaptive` (`layout` only) is the default and the only keep rule; its earlier name
//! `policy7` and the flag's earlier name `--s2-drop` are still accepted.
//! `final` reads from W: q16_mbbc.drop.tsv (written by `layout`), acc.tsv,
//! q9c/pool_detail.tsv and q9c/pool_ovl.bin (the rescue's).

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail, ensure};
use haricot::overlap::Overlap;

use haricot::pipeline::{admit, cache, endbridge, fbg_build, fill, keep_cut};

use fbg_build::prof;
use keep_cut::KEEP_RULE;

/// The command line.
#[derive(Default)]
struct Args {
    phase: String,
    cands: String,
    work: String,
    out: Option<String>,
    key: Option<String>,
    /// The stage-2 drop lists' common prefix (`--stage2-drop`).
    stage2_drop: Option<String>,
    scores: Option<String>,
    keep: bool,
}

fn parse_args(argv: &[String]) -> Result<Args> {
    let usage = "usage: haricot-layout layout|final <ovl.bin> --work W [--out OUT] [--key KEY] \
                 [--stage2-drop PREFIX] [--scores S1.scores.tsv] [--keep-intermediates] \
                 [--keep-level depth-adaptive]";
    let mut a = Args::default();
    let mut level = false;
    let mut pos: Vec<String> = Vec::new();
    let mut work = None;
    let mut it = argv.iter().skip(1);
    while let Some(x) = it.next() {
        let mut val = |name: &str| -> Result<String> {
            it.next()
                .cloned()
                .with_context(|| format!("{name} needs a value; {usage}"))
        };
        match x.as_str() {
            "--work" => work = Some(val("--work")?),
            "--out" => a.out = Some(val("--out")?),
            "--key" => a.key = Some(val("--key")?),
            // `--s2-drop`: the flag's earlier name, still accepted
            "--stage2-drop" | "--s2-drop" => a.stage2_drop = Some(val(x.as_str())?),
            "--scores" => a.scores = Some(val("--scores")?),
            "--keep-intermediates" => a.keep = true,
            "--keep-level" => {
                let l = val("--keep-level")?;
                ensure!(
                    l == KEEP_RULE || l == keep_cut::KEEP_RULE_ALIAS,
                    "--keep-level {l}: the keep rule is {KEEP_RULE}; {usage}"
                );
                level = true;
            }
            s if s.starts_with("--") => bail!("unknown option {s}; {usage}"),
            _ => pos.push(x.clone()),
        }
    }
    ensure!(pos.len() == 2, "{usage}");
    a.phase = pos[0].clone();
    a.cands = pos[1].clone();
    a.work = work.with_context(|| format!("--work is required; {usage}"))?;
    if a.phase == "layout" {
        ensure!(
            a.key.is_some() && a.stage2_drop.is_some() && a.scores.is_some() && a.out.is_none(),
            "layout needs --key, --stage2-drop and --scores (no --out); {usage}"
        );
    } else if a.phase == "final" {
        ensure!(
            a.out.is_some() && !level,
            "final needs --out (and takes no --keep-level); {usage}"
        );
    } else {
        bail!("{}: the phase is layout or final; {usage}", a.phase);
    }
    Ok(a)
}

/// A path as UTF-8 text (the programs take their paths as strings).
fn text(p: &Path) -> Result<String> {
    p.to_str().map(String::from).context("a non-UTF-8 path")
}

/// The candidates and what every step shares.
struct Ctx {
    counts: Vec<u32>,
    names: Vec<String>,
    ovl: Vec<Overlap>,
    content: Option<Vec<u64>>,
    work: PathBuf,
    /// Every read's tie rank, for the fill's canonical frame.
    rank: Vec<u32>,
    a: Args,
}

impl Ctx {
    fn load(a: Args) -> Result<Ctx> {
        let t = Instant::now();
        let cache::Cache {
            counts,
            names,
            ovl,
            content,
        } = cache::read(&a.cands)?;
        // distinct read ids
        haricot::readorder::ensure_distinct_ids(names.iter().map(String::as_str))
            .with_context(|| format!("reading {}", a.cands))?;
        // the fill's canonical frame: the tie rank; a cache without content hashes is refused
        let rank = cache::tie_rank(content.as_deref())?;
        prof(
            "load",
            t,
            &format!("overlaps={} reads={}", ovl.len(), names.len()),
        );
        Ok(Ctx {
            counts,
            names,
            ovl,
            content,
            work: PathBuf::from(&a.work),
            rank,
            a,
        })
    }

    fn w(&self, name: &str) -> PathBuf {
        self.work.join(name)
    }

    /// Write `ovl` as a cache at W/`name` (`--keep-intermediates`), with the candidates' reads.
    fn keep_cache(&self, name: &str, ovl: &[Overlap]) -> Result<()> {
        if self.a.keep {
            let p = text(&self.w(name))?;
            cache::write(&p, &self.counts, &self.names, ovl, self.content.as_deref())?;
        }
        Ok(())
    }

    /// Which candidates survive the drop list `text` (read from `list`), with the drop list's report
    /// lines. A non-empty list that matches no pair stops the run; pairs not in the cache are reported.
    fn survivors(&self, list: &str, text: &str) -> Result<Vec<&Overlap>> {
        let idx = keep_cut::name_index(&self.names);
        let (pairs, hit, miss) = keep_cut::parse_drop_pairs(text, &idx);
        eprintln!("drop list {list}: {hit} pairs matched, {miss} not in the cache");
        if !text.is_empty() {
            ensure!(
                hit > 0,
                "the cache matched no pairs of the non-empty drop list {list}"
            );
            if miss > 0 {
                println!("WARNING: {miss} pairs of {list} are not in the cache ({hit} matched)");
            }
        }
        let f = keep_cut::Filters { drop_pairs: pairs };
        let mut t = keep_cut::Tally::default();
        let out: Vec<&Overlap> = self.ovl.iter().filter(|o| f.keep(o, &mut t)).collect();
        eprintln!("{}", t.line(self.ovl.len(), out.len()));
        Ok(out)
    }

    /// [`Ctx::survivors`] of the list at `list`, as owned overlaps.
    fn keep_without(&self, step: &str, list: &str, text: &str) -> Result<Vec<Overlap>> {
        let t = Instant::now();
        let out: Vec<Overlap> = self.survivors(list, text)?.into_iter().cloned().collect();
        prof(step, t, &format!("kept={}", out.len()));
        Ok(out)
    }

    /// One build: the fill (Theil-Sen, canonical frame), then the graph build with `graph`'s
    /// flags into `outdir`, the Stage-3 read placements to `threads`.
    fn build(
        &self,
        step: &str,
        mut ovl: Vec<Overlap>,
        g: &fbg_build::Opts,
        outdir: &Path,
        threads: &Path,
    ) -> Result<()> {
        let reads = fbg_build::Reads {
            counts: &self.counts,
            names: &self.names,
            content: self.content.as_deref(),
        };
        let t = Instant::now();
        let frame = (self.counts.as_slice(), self.rank.as_slice());
        let st = fill::fill_all(&mut ovl, frame);
        eprintln!("{}", fill::report(ovl.len(), &st));
        self.keep_cache("ff.bin", &ovl)?;
        prof(&format!("{step}_fill"), t, "");
        let t = Instant::now();
        std::fs::create_dir_all(outdir)?;
        fbg_build::write_read_tables(outdir, &self.names, self.content.as_deref())?;
        fbg_build::overlap_stats(&ovl, &g.params);
        let b = fbg_build::build(
            outdir,
            &g.params.graph,
            &g.layout,
            &reads,
            &ovl,
            Some(threads),
        )?;
        prof(&format!("{step}_graph"), t, "");
        println!("[{step}] {} unitigs, {} contigs", b.unitigs, b.contigs);
        Ok(())
    }
}

/// haricot-detect's flags for the graph builds: `fbg_build::GRAPH_FLAGS`.
fn graph_opts() -> Result<fbg_build::Opts> {
    fbg_build::parse(fbg_build::GRAPH_FLAGS.split_whitespace().map(String::from))
}

/// `layout`: the depth-adaptive keep cut, its build, the end bridges, their build.
fn layout(c: &Ctx) -> Result<()> {
    let a = &c.a;
    let (key, prefix) = (a.key.as_deref().unwrap(), a.stage2_drop.as_deref().unwrap());
    let read = |p: &str| std::fs::read_to_string(p).with_context(|| format!("reading {p}"));
    // the keep cut: the candidates PREFIX0.995.drop.tsv keeps (s2.bin) and their depth estimate
    // (depth.txt) pick the stage-2 keep level
    let t = Instant::now();
    let l995 = format!("{prefix}0.995.drop.tsv");
    let kept995 = c.survivors(&l995, &read(&l995)?)?;
    let dep = keep_cut::py_fixed1(keep_cut::depth_estimate(&kept995, &c.counts));
    if a.keep {
        let v: Vec<Overlap> = kept995.iter().copied().cloned().collect();
        c.keep_cache("s2.bin", &v)?;
    }
    drop(kept995);
    std::fs::write(c.w("depth.txt"), format!("{dep}\n")).context("writing depth.txt")?;
    let keep = keep_cut::depth_adaptive_level(&dep);
    println!("[{key}] depth {dep} -> keep {keep}");
    prof("keep_depth", t, &format!("depth={dep}"));
    let g = graph_opts()?;

    // the kept overlaps: the candidates the drop list of that level keeps, built into W/pol7
    // (W/thr.tsv)
    let drop_list = format!("{prefix}{keep}.drop.tsv");
    let f = c.keep_without("keep", &drop_list, &read(&drop_list)?)?;
    c.keep_cache("f.bin", &f)?;
    c.build("keep", f, &g, &c.w("pol7"), &c.w("thr.tsv"))?;

    // end bridges (src/pipeline/endbridge.rs) on W/pol7 and W/thr.tsv: --drop DROP --scores S --out-dir W
    // --only q16_mbbc --content-ties, the candidates from memory
    let t = Instant::now();
    if a.keep {
        std::fs::write(
            c.w("cand.dump.tsv"),
            keep_cut::dump_tsv(&c.ovl, &c.counts, &c.names),
        )
        .context("writing cand.dump.tsv")?;
    }
    let argv = vec![
        "haricot-end-bridges".to_string(),
        text(&c.w("pol7"))?,
        text(&c.w("thr.tsv"))?,
        "(candidates in memory)".to_string(),
        "--drop".to_string(),
        drop_list,
        "--scores".to_string(),
        a.scores.clone().unwrap(),
        "--out-dir".to_string(),
        a.work.clone(),
        "--only".to_string(),
        "q16_mbbc".to_string(),
        "--content-ties".to_string(),
    ];
    let ea = endbridge::parse_args(&argv)?;
    let mut nm = endbridge::Names::default();
    let mem = endbridge::Mem {
        ovl: &c.ovl,
        counts: &c.counts,
        names: &c.names,
    };
    let sel = endbridge::run_on(&ea, &mut nm, mem)?;
    for l in &sel.lines {
        eprintln!("{l}");
        if l.contains("contig ends") || l.contains("mbbc>=0.5 ") {
            println!("{l}");
        }
    }
    let arms = endbridge::write_arms(&sel, &nm, ea.only.as_deref(), &ea.out_dir)?;
    let (_, q16) = arms.first().context("no q16_mbbc arm")?;
    prof("end_bridges", t, "");

    // their layout: the candidates W/q16_mbbc.drop.tsv keeps, the build into W/eb (W/thr_eb.tsv)
    let q16p = format!("{}/q16_mbbc.drop.tsv", a.work);
    let f = c.keep_without("eb_keep", &q16p, q16)?;
    c.keep_cache("f.bin", &f)?;
    c.build("eb", f, &g, &c.w("eb"), &c.w("thr_eb.tsv"))
}

/// `final`: the rescue's keep and admission, the final build into OUT (OUT/threads.tsv).
fn final_build(c: &Ctx) -> Result<()> {
    let a = &c.a;
    let w = |n: &str| text(&c.w(n));
    // the candidates W/q16_mbbc.drop.tsv keeps (keb.bin)
    let q16p = w("q16_mbbc.drop.tsv")?;
    let q16 = std::fs::read_to_string(&q16p).with_context(|| format!("reading {q16p}"))?;
    let mut keb = c.keep_without("rescue_keep", &q16p, &q16)?;
    c.keep_cache("keb.bin", &keb)?;
    // the rescue's admission of W/q9c/pool_detail.tsv, W/q9c/pool_ovl.bin into keb.bin as allowed by
    // W/acc.tsv, into q9c.bin, the admitted pairs to W/admit.tsv
    let t = Instant::now();
    let admit_tsv = w("admit.tsv")?;
    let line = admit::admit(
        &w("q9c/pool_detail.tsv")?,
        &w("q9c/pool_ovl.bin")?,
        &w("acc.tsv")?,
        &mut keb,
        c.content.as_deref(),
        Some(admit_tsv.as_str()),
    )?;
    eprintln!("{line}");
    println!("{line}");
    c.keep_cache("q9c.bin", &keb)?;
    prof("rescue_admit", t, "");
    let g = graph_opts()?;
    let out = PathBuf::from(a.out.as_deref().unwrap());
    c.build("final", keb, &g, &out, &out.join("threads.tsv"))
}

fn main() -> Result<()> {
    let t = Instant::now();
    let argv: Vec<String> = std::env::args().collect();
    let a = parse_args(&argv)?;
    let c = Ctx::load(a)?;
    match c.a.phase.as_str() {
        "layout" => layout(&c)?,
        _ => final_build(&c)?,
    }
    prof("TOTAL", t, &c.a.phase);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ovl(a: u32, b: u32) -> Overlap {
        Overlap {
            a,
            b,
            strand: haricot::params::Strand::Forward,
            n_anchors: 9,
            matched: 400,
            cnt: 9,
            aligned_bins: 5,
            bin_events: 192,
            gap: 1,
            slope: 1.0,
            a_bin_span: (0, 4),
            b_bin_span: (2, 6),
            path: vec![(0, 2), (4, 6)],
            z: 6.5,
        }
    }

    /// A non-empty list that matches no pair stops the run, an empty list drops nothing, a matched
    /// pair is dropped in either name order.
    #[test]
    fn drop_lists_are_checked_and_applied() {
        let c = Ctx {
            counts: vec![900, 1200, 700],
            names: ["r0", "r1", "r2"].map(String::from).to_vec(),
            ovl: vec![ovl(0, 1), ovl(1, 2), ovl(0, 2)],
            content: None,
            work: PathBuf::from("."),
            rank: Vec::new(),
            a: Args::default(),
        };
        let e = c.survivors("bogus", "nope_a\tnope_b\n").unwrap_err();
        assert!(e.to_string().contains("matched no pairs"));
        assert_eq!(c.survivors("empty", "").unwrap().len(), 3);
        let kept = c.survivors("one", "r2\tr1\nnope\tr0\n").unwrap();
        let pairs: Vec<(u32, u32)> = kept.iter().map(|o| (o.a, o.b)).collect();
        assert_eq!(pairs, [(0, 1), (0, 2)]);
    }

    #[test]
    fn graph_flags_are_the_builds_only() {
        // the graph flags of `haricot assemble`, spelled out
        let gx = "--bin 192 --node-bins 4 --regovl-bins 1 --min-edge-support 0.20 --clean --unitig \
                  --dump-paths --branched-bubbles 2 --frg-layout --canonical-ties";
        assert_eq!(
            gx.split_whitespace().collect::<Vec<_>>(),
            fbg_build::GRAPH_FLAGS
                .split_whitespace()
                .collect::<Vec<_>>()
        );
        let g = graph_opts().unwrap();
        assert!(g.layout.clean && g.layout.unitig && g.layout.dump_paths);
        assert!(g.params.graph.frg_layout);
        assert_eq!(g.params.graph.node_bins, 4);
        assert!(g.input.is_none() && g.save_overlaps.is_none());
    }

    /// `haricot-detect --production` is the detection and graph flags spelled out, parsed in place;
    /// the detection's parameters (everything but the graph's) are those of the detection flags alone.
    #[test]
    fn production_is_the_detection_and_graph_flags() {
        let det = "--hough --hough-min-conf 0 --hough-min-z 3 --null-gate --null-z 4 --min-aln-bins 4 \
                   --max-seed-occ 2000 --seed-occ-rel 17 --prefilter-nlp 2 --symmetric --bin 192";
        let gx = "--bin 192 --node-bins 4 --regovl-bins 1 --min-edge-support 0.20 --clean --unitig \
                  --dump-paths --branched-bubbles 2 --frg-layout --canonical-ties";
        let parse = |l: String| fbg_build::parse(l.split_whitespace().map(String::from)).unwrap();
        let tail = "--outdir w/p2f --save-overlaps w/ovl.bin --no-graph";
        let p = parse(format!("r.pod5 --production {tail}"));
        let e = parse(format!("r.pod5 {det} {gx} {tail}"));
        let d = parse(format!("r.pod5 {det} {tail}"));
        assert_eq!(format!("{:?}", p.params), format!("{:?}", e.params));
        let (lp, le) = (&p.layout, &e.layout);
        assert!(lp.clean && lp.unitig && lp.dump_paths);
        assert_eq!(
            (lp.clean, lp.unitig, lp.dump_paths),
            (le.clean, le.unitig, le.dump_paths)
        );
        let pp = &p.params;
        let dp = &d.params;
        assert_eq!(
            format!(
                "{:?}",
                (
                    &pp.signal,
                    &pp.segmentation,
                    &pp.quant,
                    &pp.sketch,
                    &pp.overlap
                )
            ),
            format!(
                "{:?}",
                (
                    &dp.signal,
                    &dp.segmentation,
                    &dp.quant,
                    &dp.sketch,
                    &dp.overlap
                )
            )
        );
        assert!(p.input.is_some() && p.save_overlaps.is_some());
        // a flag after it overrides it
        assert_eq!(
            parse("--production --bin 96".into())
                .params
                .overlap
                .bin_events,
            96
        );
    }

    #[test]
    fn arguments() {
        let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        let ok = parse_args(&s(&[
            "haricot-layout",
            "layout",
            "ovl.bin",
            "--work",
            "w",
            "--key",
            "k",
            "--stage2-drop",
            "m/stage2__k__keep",
            "--scores",
            "m/stage1__k.scores.tsv",
            "--keep-level",
            "depth-adaptive",
        ]))
        .unwrap();
        assert!(!ok.keep);
        let a = |l: &str| parse_args(&l.split(' ').map(String::from).collect::<Vec<_>>());
        let lay = "haricot-layout layout c --work w --key k --stage2-drop p --scores s";
        let fin = "haricot-layout final c --work w --out o";
        assert!(a(lay).is_ok());
        assert!(a(&format!("{lay} --keep-intermediates")).unwrap().keep);
        assert!(a(&format!("{lay} --keep-level 0.990")).is_err());
        assert!(a(&format!("{lay} --keep-level")).is_err());
        // the earlier names, still accepted
        assert!(a(&format!("{lay} --keep-level policy7")).is_ok());
        let old = "haricot-layout layout c --work w --key k --s2-drop p --scores s";
        assert_eq!(a(old).unwrap().stage2_drop.as_deref(), Some("p"));
        assert!(a(&format!("{lay} --no-symmetric")).is_err());
        assert!(a(fin).is_ok());
        assert!(a(&format!("{fin} --keep-level depth-adaptive")).is_err());
        assert!(parse_args(&s(&["haricot-layout", "final", "ovl.bin", "--work", "w"])).is_err());
        assert!(parse_args(&s(&["haricot-layout", "other", "ovl.bin", "--work", "w"])).is_err());
    }
}
