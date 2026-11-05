//! The steps of `haricot assemble` after the forests, in two processes, the candidates held in memory:
//!
//! * `layout`: policy7 (the depth of the stage-2 keep-0.995 set picks the keep cut and the read
//!   threading), the build of the kept candidates (W/pol7), the end bridges (W/q16_mbbc.drop.tsv) and
//!   their build (W/eb);
//! * `final`: the contig-end rescue's admission (W/admit.tsv) and the final build into OUT.
//!
//! Between the two, `haricot assemble` runs the rescue's scan (`haricot-rescue rescue`, into W/q9c),
//! features and forest, which read the W/eb build.
//!
//! Every step runs shared code (`src/lo.rs`: haricot-ovl-filter's drop lists, filter, dump and depth;
//! `src/fill.rs`: the fill; `src/fbg_build.rs`: haricot-detect's flags, build and layout files;
//! `src/endbridge.rs`: the end bridges; `src/admit.rs`: haricot-rescue admit), on the values a program
//! would read back from the intermediate caches (a cache's bincode round trip is lossless). The outputs are the work directory's depth.txt,
//! pol7/, thr.tsv, q16_mbbc.drop.tsv, eb/, thr_eb.tsv and admit.tsv, and the final build in OUT (with
//! threads.tsv). The intermediates nothing downstream reads (s2.bin, f.bin, ff.bin, rm.bin,
//! cand.dump.tsv, keb.bin, q9c.bin) are written only with `--keep-intermediates`.
//!
//! Usage:
//!   haricot-layout layout <ovl.bin> --work W --key KEY --s2-drop PREFIX --scores S1.scores.tsv
//!       [--graph "FLAGS"] [--no-symmetric] [--no-content-ties] [--read-model]
//!       [--keep-intermediates] [--keep-level policy7|0.995|0.990|0.980]
//!   haricot-layout final <ovl.bin> --work W --out OUT [--graph "FLAGS"] [--no-symmetric] [--read-model]
//!       [--keep-intermediates]
//!
//! PREFIX: the stage-2 drop lists are PREFIX0.995.drop.tsv, PREFIX0.990.drop.tsv and
//! PREFIX0.980.drop.tsv (haricot-forest's). FLAGS: haricot-detect's flags for every graph build, by default
//! `fbg_build::GRAPH_FLAGS` (an empty FLAGS is refused); policy7's read threading is appended. The
//! canonical pair frame (the fill's and haricot-rescue admit's `--symmetric`) and the content-hash
//! tie rule of the end bridges (`--content-ties`) are the defaults (`--symmetric`
//! and `--content-ties` are accepted); `--no-symmetric` and `--no-content-ties` turn them off.
//! `--read-model`: haricot-ovl-filter's read model before each fill. Unknown flags are an error.
//! `--keep-level` (layout only): `policy7` (the default) lets the depth pick the stage-2 keep level; a
//! level (0.995, 0.990 or 0.980) is used at every depth instead. It replaces the keep cut only: the read
//! threading stays policy7's (by depth), in both phases.
//! Inputs `final` reads from W: q16_mbbc.drop.tsv and depth.txt (written by `layout`), acc.tsv,
//! q9c/pool_detail.tsv and q9c/pool_ovl.bin (the rescue's).
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail, ensure};
use haricot::overlap::Overlap;

#[path = "../admit.rs"]
mod admit;
#[path = "../cache.rs"]
mod cache;
#[path = "../endbridge.rs"]
mod endbridge;
#[path = "../fbg_build.rs"]
mod fbg_build;
#[path = "../fill.rs"]
mod fill;
#[path = "../lo.rs"]
mod lo;

use fbg_build::prof;

/// The command line.
#[derive(Default)]
struct Args {
    phase: String,
    cands: String,
    work: String,
    out: Option<String>,
    key: Option<String>,
    s2_drop: Option<String>,
    scores: Option<String>,
    graph: Option<String>,
    symmetric: bool,
    content_ties: bool,
    read_model: bool,
    keep: bool,
    keep_level: KeepRule,
}

/// The stage-2 keep cut: policy7's depth rule, or one fixed level (`--keep-level`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
enum KeepRule {
    #[default]
    Policy7,
    Fixed(&'static str),
}

impl KeepRule {
    /// `policy7`, or one of the three levels haricot-forest writes drop lists for (as their file names
    /// spell them).
    fn parse(s: &str) -> Result<KeepRule> {
        Ok(match s {
            "policy7" => KeepRule::Policy7,
            "0.995" => KeepRule::Fixed("0.995"),
            "0.990" => KeepRule::Fixed("0.990"),
            "0.980" => KeepRule::Fixed("0.980"),
            _ => bail!("--keep-level {s}: policy7, 0.995, 0.990 or 0.980"),
        })
    }
}

fn parse_args(argv: &[String]) -> Result<Args> {
    let usage = "usage: haricot-layout layout|final <ovl.bin> --work W [--out OUT] [--key KEY] \
                 [--s2-drop PREFIX] [--scores S1.scores.tsv] [--graph FLAGS] [--no-symmetric] \
                 [--no-content-ties] [--read-model] [--keep-intermediates] \
                 [--keep-level policy7|0.995|0.990|0.980]";
    let mut a = Args::default();
    let mut level = None;
    let mut pos: Vec<String> = Vec::new();
    let mut work = None;
    // the canonical frame and the content-hash ties: the defaults, each with its explicit form
    let (mut sym, mut nosym, mut ties, mut noties) = (false, false, false, false);
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
            "--s2-drop" => a.s2_drop = Some(val("--s2-drop")?),
            "--scores" => a.scores = Some(val("--scores")?),
            "--graph" => a.graph = Some(val("--graph")?),
            "--symmetric" => sym = true,
            "--no-symmetric" => nosym = true,
            "--content-ties" => ties = true,
            "--no-content-ties" => noties = true,
            "--read-model" => a.read_model = true,
            "--keep-intermediates" => a.keep = true,
            "--keep-level" => level = Some(KeepRule::parse(&val("--keep-level")?)?),
            s if s.starts_with("--") => bail!("unknown option {s}; {usage}"),
            _ => pos.push(x.clone()),
        }
    }
    ensure!(pos.len() == 2, "{usage}");
    a.phase = pos[0].clone();
    a.cands = pos[1].clone();
    a.work = work.with_context(|| format!("--work is required; {usage}"))?;
    ensure!(
        !(sym && nosym),
        "--symmetric and --no-symmetric together; {usage}"
    );
    ensure!(
        !(ties && noties),
        "--content-ties and --no-content-ties together; {usage}"
    );
    a.symmetric = !nosym;
    a.content_ties = !noties;
    let g = a
        .graph
        .get_or_insert_with(|| fbg_build::GRAPH_FLAGS.to_string());
    ensure!(
        !g.trim().is_empty(),
        "--graph '{g}': empty (leave --graph out for the default graph flags); {usage}"
    );
    if a.phase == "layout" {
        ensure!(
            a.key.is_some() && a.s2_drop.is_some() && a.scores.is_some() && a.out.is_none(),
            "layout needs --key, --s2-drop and --scores (no --out); {usage}"
        );
    } else if a.phase == "final" {
        ensure!(
            a.out.is_some() && !ties && !noties && level.is_none(),
            "final needs --out (and takes no --content-ties, --no-content-ties or --keep-level); \
             {usage}"
        );
    } else {
        bail!("{}: the phase is layout or final; {usage}", a.phase);
    }
    a.keep_level = level.unwrap_or_default();
    Ok(a)
}

/// policy7: the depth estimate as printed (depth.txt) picks the stage-2 keep level and the read
/// threading: below 8.5x 0.995, below 19x 0.990, else 0.980 with read threading. `NaN` (no overlaps)
/// reads as 0: 0.995.
fn policy7(dep: &str) -> (&'static str, &'static [&'static str]) {
    let d: f64 = dep.trim().parse().unwrap_or(0.0);
    if d >= 19.0 {
        ("0.980", THREADING)
    } else if d >= 8.5 {
        ("0.990", NO_THREADING)
    } else {
        ("0.995", NO_THREADING)
    }
}

/// The keep level and the read threading under `rule`: policy7's, or the fixed level with policy7's
/// threading (the fixed level replaces the keep cut only).
fn keep_rule(dep: &str, rule: KeepRule) -> (&'static str, &'static [&'static str]) {
    let (keep, th) = policy7(dep);
    match rule {
        KeepRule::Policy7 => (keep, th),
        KeepRule::Fixed(level) => (level, th),
    }
}

/// policy7's read threading at 19x and above, as haricot-detect flags.
const THREADING: &[&str] = &["--read-threading", "--thread-min-reads", "3"];
const NO_THREADING: &[&str] = &[];

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
    /// Every read's tie rank, for the fill's canonical frame (empty with `--no-symmetric`).
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
        // as haricot-detect --load-overlaps checks every cache it builds from
        haricot::readorder::ensure_distinct_ids(names.iter().map(String::as_str))
            .with_context(|| format!("reading {}", a.cands))?;
        // the fill's canonical frame (--symmetric): the tie rank; a cache without content hashes is refused
        let rank = if a.symmetric {
            cache::tie_rank(content.as_deref())?
        } else {
            Vec::new()
        };
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

    /// `haricot-ovl-filter <ovl.bin> --drop-pairs LIST`: which candidates survive the drop list `text` (read
    /// from `list`), with haricot-ovl-filter's report lines. A non-empty list that matches no pair stops
    /// the run; pairs not in the cache are reported.
    fn survivors(&self, list: &str, text: &str) -> Result<Vec<&Overlap>> {
        let idx = lo::name_index(&self.names);
        let (pairs, hit, miss) = lo::parse_drop_pairs(text, &idx);
        eprintln!("--drop-pairs {list}: {hit} pairs matched, {miss} not in the cache");
        if !text.is_empty() {
            ensure!(
                hit > 0,
                "haricot-ovl-filter matched no pairs of the non-empty drop list {list}"
            );
            if miss > 0 {
                println!("WARNING: {miss} pairs of {list} are not in the cache ({hit} matched)");
            }
        }
        let f = lo::Filters {
            drop_pairs: pairs,
            ..Default::default()
        };
        let mut t = lo::Tally::default();
        let out: Vec<&Overlap> = self
            .ovl
            .iter()
            .filter(|o| f.keep(o, &self.counts, &mut t))
            .collect();
        eprintln!("{}", t.line(self.ovl.len(), out.len()));
        Ok(out)
    }

    /// [`Ctx::survivors`] of the list at `list`, as owned overlaps (the cache haricot-ovl-filter wrote).
    fn keep_without(&self, step: &str, list: &str, text: &str) -> Result<Vec<Overlap>> {
        let t = Instant::now();
        let out: Vec<Overlap> = self.survivors(list, text)?.into_iter().cloned().collect();
        prof(step, t, &format!("kept={}", out.len()));
        Ok(out)
    }

    /// One build: [haricot-ovl-filter's read model,] the fill (`fill::DRIVER` [, canonical frame]),
    /// then `haricot-detect --load-overlaps` with `graph`'s flags into `outdir`, the Stage-3 read
    /// placements to `threads`.
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
        if self.a.read_model {
            // haricot-ovl-filter IN --out rm.bin --read-model --read-model-tsv W/STEP.rm.tsv: no filter
            let t = Instant::now();
            eprintln!("{}", lo::Tally::default().line(ovl.len(), ovl.len()));
            let tsv = text(&self.w(&format!("{step}.rm.tsv")))?;
            let n = ovl.len();
            ovl = lo::read_model_apply(ovl, &self.counts, &self.names, Some(&tsv), None)?;
            println!("[{step}] read model: {n} -> {} overlaps", ovl.len());
            self.keep_cache("rm.bin", &ovl)?;
            prof(&format!("{step}_rm"), t, "");
        }
        let t = Instant::now();
        let sym = self.a.symmetric;
        let frame = sym.then_some((self.counts.as_slice(), self.rank.as_slice()));
        let st = fill::fill_all(&mut ovl, &fill::DRIVER, &self.names, frame)?;
        eprintln!("{}", fill::report(ovl.len(), &st, &fill::DRIVER));
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

/// haricot-detect's flags for the graph builds: FLAGS and the read threading `th`. Only the build's
/// flags: the ones that make haricot-detect read, detect, filter or save candidates, or choose its
/// output directory, are refused (haricot-layout builds from memory into its own directories).
fn graph_opts(flags: &str, th: &[&str]) -> Result<fbg_build::Opts> {
    let args = flags
        .split_whitespace()
        .chain(th.iter().copied())
        .map(String::from);
    let g = fbg_build::parse(args)?;
    ensure!(
        g.input.is_none()
            && g.load_overlaps.is_none()
            && g.save_overlaps.is_none()
            && g.ovl_min_z.is_none()
            && !g.no_graph
            && g.outdir == Path::new("out_fbg"),
        "--graph '{flags}': only graph-build flags (no input, --load-overlaps, --save-overlaps, \
         --ovl-min-z, --no-graph or --outdir)"
    );
    Ok(g)
}

/// `layout`: policy7, its build, the end bridges, their build.
fn layout(c: &Ctx) -> Result<()> {
    let a = &c.a;
    let (key, prefix) = (a.key.as_deref().unwrap(), a.s2_drop.as_deref().unwrap());
    let read = |p: &str| std::fs::read_to_string(p).with_context(|| format!("reading {p}"));
    // policy7: haricot-ovl-filter <ovl.bin> --out s2.bin --drop-pairs PREFIX0.995.drop.tsv; haricot-ovl-filter s2.bin
    // --depth > depth.txt
    let t = Instant::now();
    let l995 = format!("{prefix}0.995.drop.tsv");
    let s2 = c.survivors(&l995, &read(&l995)?)?;
    let dep = lo::py_fixed1(lo::depth_estimate(&s2, &c.counts));
    if a.keep {
        let v: Vec<Overlap> = s2.iter().copied().cloned().collect();
        c.keep_cache("s2.bin", &v)?;
    }
    drop(s2);
    std::fs::write(c.w("depth.txt"), format!("{dep}\n")).context("writing depth.txt")?;
    let (keep, th) = keep_rule(&dep, a.keep_level);
    let thx = if th.is_empty() {
        "(no threading)".to_string()
    } else {
        th.join(" ")
    };
    let fixed = match a.keep_level {
        KeepRule::Policy7 => String::new(),
        KeepRule::Fixed(l) => {
            let p7 = policy7(&dep).0;
            format!(" (--keep-level {l}: policy7 would keep {p7})")
        }
    };
    println!("[{key}] depth {dep} -> keep {keep} {thx}{fixed}");
    prof("pol7_depth", t, &format!("depth={dep}"));
    let g = graph_opts(a.graph.as_deref().unwrap(), th)?;

    // the policy7 layout: haricot-ovl-filter --drop-pairs DROP, the build into W/pol7 (W/thr.tsv)
    let drop_list = format!("{prefix}{keep}.drop.tsv");
    let f = c.keep_without("pol7_keep", &drop_list, &read(&drop_list)?)?;
    c.keep_cache("f.bin", &f)?;
    c.build("pol7", f, &g, &c.w("pol7"), &c.w("thr.tsv"))?;

    // end bridges (src/endbridge.rs) on W/pol7 and W/thr.tsv: --drop DROP --scores S --out-dir W
    // --only q16_mbbc [--content-ties], the candidates from memory instead of a haricot-ovl-filter --dump
    let t = Instant::now();
    if a.keep {
        std::fs::write(
            c.w("cand.dump.tsv"),
            lo::dump_tsv(&c.ovl, &c.counts, &c.names),
        )
        .context("writing cand.dump.tsv")?;
    }
    let mut argv = vec![
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
    ];
    if a.content_ties {
        argv.push("--content-ties".to_string());
    }
    let ea = endbridge::parse_args(&argv)?;
    let mut nm = endbridge::Names::default();
    let mem = endbridge::Mem {
        ovl: &c.ovl,
        counts: &c.counts,
        names: &c.names,
    };
    let sel = endbridge::run_on(&ea, &mut nm, Some(mem))?;
    for l in &sel.lines {
        eprintln!("{l}");
        if l.contains("contig ends") || l.contains("mbbc>=0.5 ") {
            println!("{l}");
        }
    }
    let arms = endbridge::write_arms(&sel, &nm, ea.only.as_deref(), &ea.out_dir)?;
    let (_, q16) = arms.first().context("no q16_mbbc arm")?;
    prof("haricot-end-bridges", t, "");

    // their layout: haricot-ovl-filter --drop-pairs W/q16_mbbc.drop.tsv, the build into W/eb (W/thr_eb.tsv)
    let q16p = format!("{}/q16_mbbc.drop.tsv", a.work);
    let f = c.keep_without("eb_keep", &q16p, q16)?;
    c.keep_cache("f.bin", &f)?;
    c.build("eb", f, &g, &c.w("eb"), &c.w("thr_eb.tsv"))
}

/// `final`: the rescue's keep and admission, the final build into OUT (OUT/threads.tsv).
fn final_build(c: &Ctx) -> Result<()> {
    let a = &c.a;
    let w = |n: &str| text(&c.w(n));
    // haricot-ovl-filter <ovl.bin> --out keb.bin --drop-pairs W/q16_mbbc.drop.tsv
    let q16p = w("q16_mbbc.drop.tsv")?;
    let q16 = std::fs::read_to_string(&q16p).with_context(|| format!("reading {q16p}"))?;
    let mut keb = c.keep_without("q9c_keep", &q16p, &q16)?;
    c.keep_cache("keb.bin", &keb)?;
    // haricot-rescue admit W/q9c/pool_detail.tsv W/q9c/pool_ovl.bin keb.bin --allow W/acc.tsv
    // [--symmetric] --out q9c.bin --out-pairs W/admit.tsv
    let t = Instant::now();
    let admit_tsv = w("admit.tsv")?;
    let line = admit::admit(
        &w("q9c/pool_detail.tsv")?,
        &w("q9c/pool_ovl.bin")?,
        &w("acc.tsv")?,
        &mut keb,
        c.content.as_deref(),
        a.symmetric,
        Some(admit_tsv.as_str()),
    )?;
    eprintln!("{line}");
    println!("{line}");
    c.keep_cache("q9c.bin", &keb)?;
    prof("q9c_admit", t, "");
    // policy7's read threading, from the depth `layout` wrote
    let dep = std::fs::read_to_string(c.w("depth.txt")).context("reading depth.txt")?;
    let (_, th) = policy7(&dep);
    let g = graph_opts(a.graph.as_deref().unwrap(), th)?;
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

    /// The rule on the printed depth.
    #[test]
    fn policy7_reads_the_printed_depth() {
        assert_eq!(policy7("8.2").0, "0.995");
        assert_eq!(policy7("8.5\n").0, "0.990");
        assert_eq!(policy7("18.8").0, "0.990");
        assert_eq!(policy7("19.0"), ("0.980", THREADING));
        assert_eq!(policy7("NaN").0, "0.995");
        assert!(policy7("12.2").1.is_empty());
    }

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
        let g = graph_opts(gx, THREADING).unwrap();
        assert!(g.layout.clean && g.layout.unitig && g.layout.dump_paths);
        assert!(g.layout.frg.is_some() && g.layout.canonical_ties);
        assert!(g.params.graph.read_threading);
        assert_eq!(g.params.graph.node_bins, 4);
        assert_eq!(g.params.graph.thread_min_reads, 3);
        assert!(graph_opts("--clean --load-overlaps x.bin", &[]).is_err());
        assert!(graph_opts("--clean --outdir d", &[]).is_err());
        assert!(graph_opts("--clean reads.pod5", &[]).is_err());
        assert!(graph_opts("--clean --no-such-flag", &[]).is_err());
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
        assert!(lp.clean && lp.unitig && lp.dump_paths && lp.canonical_ties && !lp.unzip);
        assert_eq!(
            (
                lp.clean,
                lp.unitig,
                lp.dump_paths,
                lp.canonical_ties,
                lp.unzip
            ),
            (
                le.clean,
                le.unitig,
                le.dump_paths,
                le.canonical_ties,
                le.unzip
            )
        );
        assert_eq!(format!("{:?}", lp.frg), format!("{:?}", le.frg));
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
        assert!(p.no_graph && p.input.is_some() && p.save_overlaps.is_some());
        assert_eq!(
            (
                p.canonical_order,
                p.stream_sketch,
                p.stream_chunk_bytes,
                p.max_reads,
                p.ovl_min_z
            ),
            (
                d.canonical_order,
                d.stream_sketch,
                d.stream_chunk_bytes,
                d.max_reads,
                d.ovl_min_z
            )
        );
        // a flag after it overrides it
        assert_eq!(
            parse("--production --bin 96".into())
                .params
                .overlap
                .bin_events,
            96
        );
    }

    /// A fixed level replaces policy7's keep cut at every depth; the threading stays policy7's.
    #[test]
    fn keep_rule_fixes_the_level_only() {
        for d in ["NaN", "5.0", "8.5", "12.2", "18.8", "19.0", "31.4\n"] {
            assert_eq!(keep_rule(d, KeepRule::Policy7), policy7(d));
            for l in ["0.995", "0.990", "0.980"] {
                let (k, th) = keep_rule(d, KeepRule::Fixed(l));
                assert_eq!(k, l);
                assert_eq!(th, policy7(d).1);
            }
        }
        assert_eq!(
            keep_rule("25.0", KeepRule::Fixed("0.995")),
            ("0.995", THREADING)
        );
        assert_eq!(keep_rule("5.0", KeepRule::Fixed("0.980")).1, NO_THREADING);
    }

    #[test]
    fn keep_level_parses_the_three_levels_and_policy7() {
        assert_eq!(KeepRule::parse("policy7").unwrap(), KeepRule::Policy7);
        for l in ["0.995", "0.990", "0.980"] {
            assert_eq!(KeepRule::parse(l).unwrap(), KeepRule::Fixed(l));
        }
        for bad in ["0.99", "0.98", "0.95", "1", "", "Policy7"] {
            assert!(KeepRule::parse(bad).is_err(), "{bad}");
        }
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
            "--s2-drop",
            "m/s2x__k__keep",
            "--scores",
            "m/x__k.scores.tsv",
            "--graph",
            "--clean",
            "--symmetric",
        ]))
        .unwrap();
        assert!(ok.symmetric && ok.content_ties && !ok.keep);
        assert_eq!(ok.graph.as_deref(), Some("--clean"));
        assert_eq!(ok.keep_level, KeepRule::Policy7);
        let a = |l: &str| parse_args(&l.split(' ').map(String::from).collect::<Vec<_>>());
        let lay = "haricot-layout layout c --work w --key k --s2-drop p --scores s --graph --clean";
        let fin = "haricot-layout final c --work w --out o --graph --clean";
        // the defaults: the default graph flags, the canonical frame, the content-hash ties
        let d = a("haricot-layout layout c --work w --key k --s2-drop p --scores s").unwrap();
        assert_eq!(d.graph.as_deref(), Some(fbg_build::GRAPH_FLAGS));
        assert!(d.symmetric && d.content_ties);
        let o = a(&format!("{lay} --no-symmetric --no-content-ties")).unwrap();
        assert!(!o.symmetric && !o.content_ties);
        assert!(a(&format!("{lay} --symmetric --no-symmetric")).is_err());
        assert!(a(&format!("{lay} --content-ties --no-content-ties")).is_err());
        assert!(a(&format!("{fin} --no-content-ties")).is_err());
        assert!(a(&format!("{fin} --content-ties")).is_err());
        assert!(!a(&format!("{fin} --no-symmetric")).unwrap().symmetric);
        assert!(
            parse_args(&s(&[
                "haricot-layout",
                "final",
                "c",
                "--work",
                "w",
                "--out",
                "o",
                "--graph",
                " "
            ]))
            .is_err()
        );
        let fixed = a(&format!("{lay} --keep-level 0.990")).unwrap();
        assert_eq!(fixed.keep_level, KeepRule::Fixed("0.990"));
        let p7 = a(&format!("{lay} --keep-level policy7")).unwrap();
        assert_eq!(p7.keep_level, KeepRule::Policy7);
        assert!(a(&format!("{lay} --keep-level 0.99")).is_err());
        assert!(a(&format!("{lay} --keep-level")).is_err());
        assert!(a(fin).is_ok());
        assert!(a(&format!("{fin} --keep-level 0.990")).is_err());
        assert!(parse_args(&s(&["haricot-layout", "final", "ovl.bin", "--work", "w"])).is_err());
        assert!(parse_args(&s(&["haricot-layout", "other", "ovl.bin", "--work", "w"])).is_err());
    }
}
