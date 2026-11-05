//! Command-line interface.
//!
//! - `assemble` — the full assembler (reads -> layout), as run for the manuscript ([`crate::assemble`])
//! - `baseline` — the plain wtdbg2 / RawHash2 baseline (bin-DP overlaps, FBG, cleaning) as one GFA
//! - `events`, `sketch`, `overlap` — the library's stages one at a time, for inspection
//! - `info` — print the parameter preset

use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};

use crate::params::Params;
use crate::{pipeline, signal};

#[derive(Parser, Debug)]
#[command(
    name = "haricot",
    version,
    about = "HARICOT: genome assembly from raw nanopore signal (R10.4.1, forward-strand reads)",
    after_help = "Assemble a POD5 file:\n  haricot assemble reads.pod5 -o out\n\nRun `haricot assemble --help` for its options."
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Increase logging verbosity (-v info, -vv debug, -vvv trace).
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    verbose: u8,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Segment raw signals into events and report per-read stats.
    Events {
        #[command(flatten)]
        io: InputArgs,
        #[command(flatten)]
        params: ParamArgs,
    },
    /// Sketch reads into e-mer seeds and report seed counts.
    Sketch {
        #[command(flatten)]
        io: InputArgs,
        #[command(flatten)]
        params: ParamArgs,
        /// Print every seed as `read\thash\tpos`.
        #[arg(long)]
        dump: bool,
    },
    /// Find candidate read overlaps (the all-vs-all KBM analog).
    Overlap {
        #[command(flatten)]
        io: InputArgs,
        #[command(flatten)]
        params: ParamArgs,
    },
    /// Assemble a POD5 file: the full pipeline, as run for the manuscript.
    ///
    /// Writes OUT/unitigs.gfa (the layout: unitigs as paths of read segments), OUT/contigs.gfa, the path and node
    /// tables, OUT/threads.tsv (read placements) and OUT/time.tsv (wall time and peak memory per step); step logs
    /// in OUT/work/.
    #[command(after_help = "Examples:\n  haricot assemble reads.pod5 -o asm\n  \
        haricot assemble reads.pod5 -o asm -t 16 --models path/to/forests\n\n\
        Input: R10.4.1 reads, 5 or 4 kHz. Strands are not merged: reads of the two strands assemble separately.\n\
        The run log states whether the forests read are the published ones (SHA-256).")]
    Assemble {
        /// Input reads (.pod5).
        reads: PathBuf,
        /// Output directory (must not exist, or be empty, or hold an earlier run with --force).
        #[arg(short, long)]
        out: PathBuf,
        /// Directory with the forests (stage1.*, stage2.*, rescue.*).
        #[arg(short, long, default_value = crate::assemble::DEFAULT_MODELS)]
        models: PathBuf,
        /// Threads per step [default: all available cores].
        #[arg(short, long)]
        threads: Option<usize>,
        /// Keep the intermediate files (candidates, features, forest scores) in OUT/work.
        #[arg(long)]
        keep_work: bool,
        /// Replace OUT if it holds an earlier run.
        #[arg(long)]
        force: bool,
        /// Directory of the pipeline programs [default: the directory of this executable].
        #[arg(long, hide = true)]
        bin_dir: Option<PathBuf>,
    },
    /// Baseline (not HARICOT's pipeline): bin-DP overlaps, fuzzy Bruijn graph, cleaning, as one GFA.
    Baseline {
        #[command(flatten)]
        io: InputArgs,
        #[command(flatten)]
        params: ParamArgs,
        /// Output GFA path (default: stdout).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Emit the raw FBG before any graph cleaning (for inspection in Bandage).
        #[arg(long)]
        raw: bool,
        /// Use the OLC read string-graph backend (nodes = reads) instead of the FBG.
        #[arg(long)]
        olc: bool,
        /// With --olc, drop contained reads and emit only the dovetail backbone
        /// (the minimal "one continuous component per cluster" view).
        #[arg(long = "olc-backbone")]
        olc_backbone: bool,
    },
    /// Print the active parameter preset and exit.
    Info {
        #[command(flatten)]
        params: ParamArgs,
    },
}

#[derive(Args, Debug)]
struct InputArgs {
    /// Input raw-signal file (.pod5; fast5 is not supported).
    input: PathBuf,
}

/// Overrides applied on top of the R10.4.1 preset.
#[derive(Args, Debug)]
struct ParamArgs {
    /// Events per e-mer.
    #[arg(long)]
    e: Option<u32>,
    /// Minimizer window (0 = dense sketch).
    #[arg(long)]
    w: Option<u32>,
    /// Nominal k-mer length (affects span).
    #[arg(long)]
    k: Option<u32>,

    /// Sig-diff filter (-1 disables).
    #[arg(long)]
    diff: Option<i32>,
    /// Seed method: dense|minimizer|open-syncmer|closed-syncmer|strobemer.
    #[arg(long)]
    method: Option<String>,
    /// s-mer length (events) for syncmers/strobemers.
    #[arg(long)]
    s: Option<u32>,
    /// Quantization buckets. `q` is raised to fit if it is too small.
    #[arg(long = "n-buckets")]
    n_buckets: Option<u32>,
    /// Quantization scheme: piecewise|normal-quantile.
    ///
    /// The default `piecewise` is RawHash2's 1-6-1 map. `normal-quantile` is
    /// equal-frequency and data-free, so it realizes close to `log2(n_buckets)` bits per
    /// event where piecewise realizes ~2.53 of 3. The pipeline uses `piecewise` at
    /// `n_buckets = 8`.
    #[arg(long = "quant")]
    quant: Option<String>,
    /// Minimum matched length (events) to report an overlap (wtdbg2 min_mat).
    #[arg(long = "min-matches")]
    min_matches: Option<u32>,
    /// Drop seeds occurring more than this many times (GLOBAL absolute cap).
    #[arg(long = "max-seed-occ")]
    max_seed_occ: Option<u32>,
    /// Scale term of the repeat mask: raise --max-seed-occ to this many times the input's
    /// mean seed occurrence (instances / distinct hashes). 0 = off (default).
    #[arg(long = "seed-occ-rel")]
    seed_occ_rel: Option<f32>,
    /// Bin size in events for the bin-DP (wtdbg2 KBM_BSIZE = 256 bp).
    #[arg(long = "bin-events")]
    bin_events: Option<u32>,
    /// Minimum matched/aligned ratio to accept an overlap (wtdbg2 min_sim).
    #[arg(long = "min-similarity")]
    min_similarity: Option<f32>,
    /// Minimum aligned span in bins (wtdbg2 min_aln).
    #[arg(long = "min-aln-bins")]
    min_aln_bins: Option<u32>,
    /// Max gap-bin fraction of a Hough alignment (wtdbg2 `max_gap`, 0.6).
    #[arg(long = "max-gap")]
    max_gap: Option<f32>,
    /// Use the Hough-extent overlapper instead of the banded bin-DP.
    #[arg(long = "hough")]
    hough: bool,
    /// Reject overlaps whose diagonal fails a null-hypothesis significance test (not in
    /// wtdbg2). Off by default.
    #[arg(long = "null-gate")]
    null_gate: bool,
    /// Minimum diagonal significance z-score for `--null-gate` (default 6.0).
    #[arg(long = "null-z")]
    null_z: Option<f32>,
    /// Partner pre-filter: align a gathered pair only when its diagonal statistic
    /// `nlp2n_1` reaches this level (2 = p <= 1e-2 under the uniform null). Off by default.
    /// See `overlap::prefilter`.
    #[arg(long = "prefilter-nlp")]
    prefilter_nlp: Option<f64>,
    /// Hough confidence floor (`OverlapParams::hough_min_conf`, default 6.0): the
    /// winning diagonal's peak / the median peak across slopes.
    #[arg(long = "hough-min-conf")]
    hough_min_conf: Option<f32>,
    /// Hough significance floor `z` (`OverlapParams::hough_min_z`, default 0 = off),
    /// applied at the same early exit as `--hough-min-conf` and on the same fit.
    /// Set `--hough-min-conf 0 --hough-min-z X` for a z-primary gate.
    #[arg(long = "hough-min-z")]
    hough_min_z: Option<f32>,
    /// Slope beam: low/high edges of the plausible slope window (default 0.7/1.4).
    #[arg(long = "slope-win-lo")]
    slope_win_lo: Option<f32>,
    #[arg(long = "slope-win-hi")]
    slope_win_hi: Option<f32>,
    /// Slope beam: coarsest slope step within the window (default 0.2).
    #[arg(long = "slope-coarse-step")]
    slope_coarse_step: Option<f32>,
    /// Slope beam: top-k distinct basins carried per level (default 3).
    #[arg(long = "slope-k")]
    slope_k: Option<u32>,
    /// Slope beam: coarse→fine levels (default 2; more = coarser hierarchical start).
    #[arg(long = "slope-levels")]
    slope_levels: Option<u32>,
    /// Slope beam: first-filter floor on the level-0 peak (`BeamCfg::screen_peak`,
    /// default 0 = off). Pairs below it are rejected before refinement.
    #[arg(long = "hough-screen")]
    hough_screen: Option<u32>,
    /// Slope beam: density-normalised first filter on the level-0 peak
    /// (`BeamCfg::screen_z`, default 0 = off).
    #[arg(long = "hough-screen-z")]
    hough_screen_z: Option<f32>,
    /// Slope beam: stage-0b filter, the gate's null `z` at the best level-0 slope
    /// (`BeamCfg::screen_null`, default 0 = off).
    #[arg(long = "hough-screen-null")]
    hough_screen_null: Option<f32>,
    /// Slope beam: pick the slope among the last level's top-k basins by the gate's null
    /// `z` instead of the raw count (`BeamCfg::select_z`, default off).
    #[arg(long = "slope-select-z")]
    slope_select_z: bool,
    /// Node size K, in DP bins: a node spans `K * --bin-events` events. The main
    /// de-tangling lever (coarsen for long reads). wtdbg2's equivalent is 4.
    #[arg(long = "node-bins")]
    node_bins: Option<u32>,
    /// DEPRECATED alias for `--node-bins`, in events: converted to
    /// `ceil(node_events / bin_events)` bins with a warning. Node size is defined in
    /// bins (see `--node-bins`); an event count that is not a whole number of bins
    /// cannot be honored exactly.
    #[arg(long = "node-events", conflicts_with = "node_bins")]
    node_events: Option<u32>,
    /// Region-overlap conflict threshold, in bins: two regions on a read conflict
    /// once they overlap by MORE than this. wtdbg2's value is 1 bin. `0` makes any
    /// single-bin overlap a conflict (i.e. strictly harsher than wtdbg2).
    #[arg(long = "regovl-bins")]
    regovl_bins: Option<u32>,
    /// Minimum supporting reads to commit a node (raise to ~3 to drop
    /// low-coverage off-consensus spurs).
    #[arg(long = "min-node-cov")]
    min_node_cov: Option<u32>,
    /// Minimum supporting reads to keep an edge.
    #[arg(long = "min-edge-cov")]
    min_edge_cov: Option<u32>,
    /// Relative edge-support gate: drop an edge whose coverage is below
    /// `f * min(cov(from), cov(to))`. 0 = off (default). Depth-aware where
    /// `--min-edge-cov` is not, since the same coverage means very different
    /// things at 4x and 40x. The pipeline uses 0.20.
    #[arg(long = "min-edge-support")]
    min_edge_support: Option<f32>,
    /// Record read-threaded node triples (path-level evidence; see `graph::threading`).
    #[arg(long = "read-threading")]
    read_threading: bool,
    /// Reads that must thread a node triple before a layout traverses it (default 1).
    #[arg(long = "thread-min-reads")]
    thread_min_reads: Option<u32>,
}

impl ParamArgs {
    /// Apply overrides to the R10.4.1 preset and validate.
    fn resolve(&self) -> Result<Params> {
        let mut p = Params::r10_4_1();
        if let Some(v) = self.e {
            p.sketch.e = v;
        }
        if let Some(v) = self.w {
            p.sketch.w = v;
        }
        if let Some(v) = self.k {
            p.sketch.k = v;
        }
        if let Some(v) = self.diff {
            p.sketch.diff = v;
        }
        if let Some(v) = &self.method {
            p.sketch.method =
                crate::emer::select::SketchMethod::parse(v).map_err(|e| anyhow::anyhow!(e))?;
        }
        if let Some(v) = self.s {
            p.sketch.s = v;
        }
        if let Some(v) = self.n_buckets {
            p.quant.n_buckets = v;
            // `validate` requires `n_buckets <= 2^q`, and `q` is the packing width, so
            // asking for more buckets without widening the packing is always a mistake
            // rather than a choice.
            let need = u32::BITS - (v.max(1) - 1).leading_zeros();
            if p.quant.q < need {
                p.quant.q = need;
            }
        }
        if let Some(v) = &self.quant {
            p.quant.method = match v.as_str() {
                "piecewise" => crate::params::QuantMethod::Piecewise,
                "normal-quantile" => crate::params::QuantMethod::NormalQuantile,
                other => anyhow::bail!("unknown --quant {other} (piecewise|normal-quantile)"),
            };
        }
        if let Some(v) = self.min_matches {
            p.overlap.min_matches = v;
        }
        if let Some(v) = self.seed_occ_rel {
            p.overlap.seed_occ_rel = v;
        }
        if let Some(v) = self.max_seed_occ {
            p.overlap.max_seed_occ = v;
        }
        if let Some(v) = self.bin_events {
            p.overlap.bin_events = v;
        }
        if let Some(v) = self.min_similarity {
            p.overlap.min_similarity = v;
        }
        if let Some(v) = self.max_gap {
            p.overlap.max_gap = v;
        }
        if let Some(v) = self.min_aln_bins {
            p.overlap.min_aln_bins = v;
        }
        if self.hough {
            p.overlap.hough = true;
        }
        if self.null_gate {
            p.overlap.null_gate = true;
        }
        if let Some(v) = self.hough_min_z {
            p.overlap.hough_min_z = v;
        }
        if let Some(v) = self.hough_min_conf {
            p.overlap.hough_min_conf = v;
        }
        if let Some(v) = self.null_z {
            p.overlap.null_z = v;
        }
        if let Some(v) = self.prefilter_nlp {
            p.overlap.prefilter_nlp = Some(v);
        }
        if let Some(v) = self.slope_win_lo {
            p.overlap.slope_beam.win_lo = v;
        }
        if let Some(v) = self.slope_win_hi {
            p.overlap.slope_beam.win_hi = v;
        }
        if let Some(v) = self.slope_coarse_step {
            p.overlap.slope_beam.coarse_step = v;
        }
        if let Some(v) = self.slope_k {
            p.overlap.slope_beam.beam_k = v;
        }
        if let Some(v) = self.slope_levels {
            p.overlap.slope_beam.levels = v;
        }
        if let Some(v) = self.hough_screen {
            p.overlap.slope_beam.screen_peak = v;
        }
        if let Some(v) = self.hough_screen_z {
            p.overlap.slope_beam.screen_z = v;
        }
        if let Some(v) = self.hough_screen_null {
            p.overlap.slope_beam.screen_null = v;
        }
        if self.slope_select_z {
            p.overlap.slope_beam.select_z = true;
        }
        if let Some(v) = self.node_bins {
            p.graph.node_bins = v;
        }
        // Deprecated `--node-events`: convert to bins the way wtdbg2 does — round UP
        // to a whole bin and say so (`wtdbg.c:638-640` prints "** Adjust -j to %d").
        // Rounding to nearest instead would silently let a node's declared length
        // differ from the ground it covers.
        if let Some(v) = self.node_events {
            let bin = p.overlap.bin_events.max(1); // --bin-events is applied above
            let bins = v.div_ceil(bin).max(1);
            let adj = bins * bin;
            log::warn!(
                "--node-events is deprecated; use --node-bins. \
                 {v} events / {bin} events-per-bin -> K={bins} bins ({adj} events){}",
                if adj == v {
                    String::new()
                } else {
                    format!(" ** adjusted from {v} to {adj}")
                }
            );
            p.graph.node_bins = bins;
        }
        if let Some(v) = self.regovl_bins {
            p.graph.regovl_bins = v;
        }
        if let Some(v) = self.min_node_cov {
            p.graph.min_node_cov = v;
        }
        if let Some(v) = self.min_edge_cov {
            p.graph.min_edge_cov = v;
        }
        if let Some(v) = self.min_edge_support {
            p.graph.min_edge_support = v;
        }
        if self.read_threading {
            p.graph.read_threading = true;
        }
        if let Some(v) = self.thread_min_reads {
            p.graph.thread_min_reads = v;
        }
        p.validate()
            .map_err(|e| anyhow::anyhow!("invalid parameters: {e}"))?;
        Ok(p)
    }
}

/// Parse args, set up logging, and dispatch. Called by `main`.
pub fn run() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose);

    match cli.command {
        Command::Events { io, params } => cmd_events(&io, &params),
        Command::Sketch { io, params, dump } => cmd_sketch(&io, &params, dump),
        Command::Overlap { io, params } => cmd_overlap(&io, &params),
        Command::Assemble {
            reads,
            out,
            models,
            threads,
            keep_work,
            force,
            bin_dir,
        } => crate::assemble::run(&crate::assemble::Options {
            reads,
            out,
            models,
            threads,
            keep_work,
            force,
            bin_dir,
        }),
        Command::Baseline {
            io,
            params,
            output,
            raw,
            olc,
            olc_backbone,
        } => cmd_baseline(&io, &params, output.as_deref(), raw, olc, olc_backbone),
        Command::Info { params } => {
            println!("{:#?}", params.resolve()?);
            Ok(())
        }
    }
}

fn init_logging(verbose: u8) {
    let mut builder = env_logger::Builder::from_default_env();
    if std::env::var_os("RUST_LOG").is_none() {
        let level = match verbose {
            0 => log::LevelFilter::Warn,
            1 => log::LevelFilter::Info,
            2 => log::LevelFilter::Debug,
            _ => log::LevelFilter::Trace,
        };
        builder.filter_level(level);
    }
    let _ = builder.try_init();
}

fn cmd_events(io: &InputArgs, params: &ParamArgs) -> Result<()> {
    let params = params.resolve()?;
    let reads = signal::read_signals(&io.input)
        .with_context(|| format!("reading {}", io.input.display()))?;
    println!("read_name\tn_samples\tn_events");
    for r in &reads {
        let events = signal::signal_to_events(&r.signal_pa, &params);
        println!("{}\t{}\t{}", r.id, r.signal_pa.len(), events.len());
    }
    log::info!("processed {} reads", reads.len());
    Ok(())
}

fn cmd_sketch(io: &InputArgs, params: &ParamArgs, dump: bool) -> Result<()> {
    let params = params.resolve()?;
    let reads = signal::read_signals(&io.input)
        .with_context(|| format!("reading {}", io.input.display()))?;

    let sketches = pipeline::sketch_reads(&reads, &params);

    if dump {
        // `pos` is the coordinate (merged-stream index); `event_pos` is the raw
        // segmented-event index it came from. See emer::retained_quantized.
        // A seed does not store its raw index: it is read from the retained stream at `pos`.
        println!("read_name\thash\tpos\tevent_pos\tstrand");
        for (r, s) in reads.iter().zip(&sketches) {
            let events = signal::signal_to_events(&r.signal_pa, &params);
            let retained = crate::emer::retained_quantized(&events, &params);
            for seed in &s.seeds {
                println!(
                    "{}\t{}\t{}\t{}\t{:?}",
                    s.name,
                    seed.hash(),
                    seed.pos,
                    retained[seed.pos as usize].event_pos,
                    seed.strand()
                );
            }
        }
    } else {
        // n_events = raw segmented; n_coord_events = after the sig-diff merge (the
        // space coordinates live in). Their ratio is the coordinate compression.
        println!("read_name\tn_events\tn_coord_events\tn_seeds");
        for s in &sketches {
            println!(
                "{}\t{}\t{}\t{}",
                s.name,
                s.n_events,
                s.n_coord_events,
                s.seeds.len()
            );
        }
    }
    Ok(())
}

fn cmd_overlap(io: &InputArgs, params: &ParamArgs) -> Result<()> {
    let params = params.resolve()?;
    let reads = signal::read_signals(&io.input)
        .with_context(|| format!("reading {}", io.input.display()))?;
    let sketches = pipeline::sketch_reads(&reads, &params);
    let overlaps = pipeline::overlaps(&sketches, &params);

    println!("read_a\tread_b\tn_anchors\tmatched\taligned_bins\tsimilarity\tslope\tstrand");
    for o in &overlaps {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:?}",
            sketches[o.a as usize].name,
            sketches[o.b as usize].name,
            o.n_anchors,
            o.matched,
            o.aligned_bins,
            o.similarity(),
            o.slope,
            o.strand,
        );
    }
    log::info!("found {} candidate overlaps", overlaps.len());
    Ok(())
}

fn cmd_baseline(
    io: &InputArgs,
    params: &ParamArgs,
    output: Option<&std::path::Path>,
    raw: bool,
    olc: bool,
    olc_backbone: bool,
) -> Result<()> {
    let params = params.resolve()?;
    let t = std::time::Instant::now();
    let reads = signal::read_signals(&io.input)
        .with_context(|| format!("reading {}", io.input.display()))?;
    log::info!(
        "read {} reads ({:.1} s)",
        reads.len(),
        t.elapsed().as_secs_f64()
    );

    // OLC read string-graph backend (nodes = reads). Emits a GFA directly here
    // because its segments are reads, not FuzzyGraph nodes.
    if olc {
        let sg = crate::layout::StringGraphParams::default();
        let (g, sketches) = pipeline::build_string_graph_layout(&reads, &params, &sg);
        let names: Vec<String> = sketches.iter().map(|s| s.name.clone()).collect();
        let lens: Vec<u32> = sketches.iter().map(|s| s.n_coord_events).collect();
        let n_contained = g.contained.iter().filter(|&&c| c).count();
        log::info!(
            "OLC string graph: {} reads ({} contained), {} backbone edges",
            names.len(),
            n_contained,
            g.edges.len()
        );
        let write = |mut w: &mut dyn std::io::Write| -> std::io::Result<()> {
            if olc_backbone {
                write_backbone_gfa(&g, &names, &lens, w)
            } else {
                // `&mut &mut dyn Write` is a sized `Write`, so the generic writer
                // can take it (a bare `&mut dyn Write` is unsized for the W param).
                crate::layout::write_string_graph_gfa(&g, &names, &lens, &mut w)
            }
        };
        match output {
            Some(path) => {
                let mut f = std::fs::File::create(path)
                    .with_context(|| format!("creating {}", path.display()))?;
                write(&mut f)?;
                eprintln!("wrote OLC GFA to {}", path.display());
            }
            None => {
                let stdout = std::io::stdout();
                let mut lock = stdout.lock();
                write(&mut lock)?;
                lock.flush()?;
            }
        }
        return Ok(());
    }

    let graph = if raw {
        pipeline::build_raw_graph(&reads, &params)
    } else {
        pipeline::assemble(&reads, &params)
    };
    log::info!(
        "assembled graph: {} nodes ({} live), {} edges",
        graph.n_nodes(),
        graph.n_live_nodes(),
        graph.n_edges()
    );

    match output {
        Some(path) => {
            let mut f = std::fs::File::create(path)
                .with_context(|| format!("creating {}", path.display()))?;
            crate::layout::write_gfa(&graph, &mut f)?;
            eprintln!("wrote GFA to {}", path.display());
        }
        None => {
            let stdout = std::io::stdout();
            let mut lock = stdout.lock();
            crate::layout::write_gfa(&graph, &mut lock)?;
            lock.flush()?;
        }
    }
    Ok(())
}

/// Backbone-only string-graph GFA: drop contained reads (transitively implied) and
/// their attachments, leaving just the dovetail backbone.
fn write_backbone_gfa(
    g: &crate::layout::StringGraph,
    names: &[String],
    lens: &[u32],
    out: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    writeln!(out, "H\tVN:Z:1.0")?;
    for (rid, name) in names.iter().enumerate() {
        if !g.contained[rid] {
            writeln!(out, "S\t{name}\t*\tLN:i:{}", lens[rid])?;
        }
    }
    for e in &g.edges {
        writeln!(
            out,
            "L\t{}\t{}\t{}\t{}\t{}M",
            names[e.from as usize],
            e.from_dir.sign(),
            names[e.to as usize],
            e.to_dir.sign(),
            e.overlap_events
        )?;
    }
    Ok(())
}
