//! End-to-end orchestration tying the stages together.
//!
//! ```text
//!   reads (pod5) ─▶ signal_to_events ─▶ sketch (e-mer seeds)
//!                                          │
//!                                          ▼
//!                                   overlaps (KBM analog)
//!                                          │
//!                                          ▼
//!                                build_graph ─▶ clean ─▶ GFA
//! ```
//! [`build_string_graph_layout`] offers an OLC string graph as an alternative
//! layout backend to the FBG path above.

use rayon::prelude::*;

use crate::emer::{self, Seed};
use crate::graph::{self, FuzzyGraph};
use crate::overlap::{self, Overlap};
use crate::params::{Params, Strand};
use crate::signal::{self, RawRead};

/// Per-read sketch result with bookkeeping for reporting.
#[derive(Debug, Clone)]
pub struct ReadSketch {
    pub read_id: u32,
    pub name: String,
    pub n_samples: u32,
    /// Raw segmented events, before the sig-diff merge. Reporting only — this is
    /// **not** the space coordinates live in.
    pub n_events: u32,
    /// Read length in **coordinate space**: the number of events retained by the
    /// sig-diff merge (see [`crate::emer::retained_quantized`]). Bins, anchors and
    /// K-bin windows are all indexed in this space, so this — not `n_events` — is
    /// the read length the graph build must be given.
    pub n_coord_events: u32,
    /// The read's signal content hash ([`crate::readorder::content_hash`] of its pA samples): the key
    /// of the tie rank ([`crate::readorder::rank_by_content`]), a function of the signal alone.
    pub content: u64,
    pub seeds: Vec<Seed>,
}

/// Run signal → events → seeds for every read. Each read is independent, so this
/// is a `par_iter`; the indexed collect preserves input order, so `read_id == i`
/// and results are identical to the serial version.
pub fn sketch_reads(reads: &[RawRead], params: &Params) -> Vec<ReadSketch> {
    sketch_reads_from(reads, 0, params)
}

/// [`sketch_reads`] for a chunk that starts at read id `id_base`.
///
/// Split out for [`sketch_pod5_streaming`]: `read_id` must keep counting across chunk
/// boundaries, since every downstream coordinate is keyed on it.
pub fn sketch_reads_from(reads: &[RawRead], id_base: u32, params: &Params) -> Vec<ReadSketch> {
    reads
        .par_iter()
        .enumerate()
        .map(|(i, r)| sketch_one(id_base + i as u32, &r.id, &r.signal_pa, params))
        .collect()
}

/// [`sketch_reads_from`] over ADC reads ([`signal::AdcRead`]): each read's pA signal is made by the
/// worker that sketches it, so only the chunk's ADC counts are held, not its `f32` signal.
pub fn sketch_adc_from(
    reads: &[signal::AdcRead],
    id_base: u32,
    params: &Params,
) -> Vec<ReadSketch> {
    reads
        .par_iter()
        .enumerate()
        .map(|(i, r)| sketch_one(id_base + i as u32, &r.id, &r.signal_pa(), params))
        .collect()
}

/// One read's signal → events → seeds, and its signal content hash. Both reading paths
/// ([`sketch_reads_from`]'s `f32` reads and [`sketch_adc_from`]'s ADC reads) come here with the
/// read's pA samples, the same bits either way (`signal::AdcRead::signal_pa` is `RawRead`'s
/// conversion), so the hash is the same whichever path read the pod5.
fn sketch_one(id: u32, name: &str, signal_pa: &[f32], params: &Params) -> ReadSketch {
    let events = signal::signal_to_events(signal_pa, params);
    let retained = emer::retained_quantized(&events, params);
    let seeds = emer::sketch_from_retained(&retained, id, Strand::Forward, params);
    // the tie rank's key, while the worker holds the signal (one pass over the samples)
    let content = crate::readorder::content_hash(signal_pa);
    ReadSketch {
        read_id: id,
        name: name.to_string(),
        n_samples: signal_pa.len() as u32,
        n_events: events.len() as u32,
        n_coord_events: retained.len() as u32,
        content,
        seeds,
    }
}

/// Read a pod5 and sketch it **without ever holding the whole signal in memory**.
///
/// The raw signal dominates peak RSS — 32.65 GB of D. melanogaster fwd's 49.40 GB and
/// 7.97 GB of E. coli fwd-30x's 11.30 GB — and the reason is not that it is retained
/// (callers already drop it after sketching) but that *loading is a separate pass*: the
/// whole signal and the whole seed store are both live for the duration of the sketch.
/// Reading a bounded chunk, sketching it, and dropping it before reading the next makes
/// the signal's contribution a constant instead of a term proportional to the dataset,
/// leaving the seed store as the floor.
///
/// Chunks are bounded in **bytes of signal**, not reads: read lengths vary ~10x across
/// these datasets, so a read count is not a memory bound. The next chunk is decoded while
/// this one is sketched ([`signal::for_each_raw_chunk`]: two chunks of `chunk_bytes / 2`).
///
/// Produces exactly what `read_signals_limit` + [`sketch_reads`] produce — same reads,
/// same order, same `read_id`s — because the pod5 iterator is consumed once, in order.
pub fn sketch_pod5_streaming(
    path: &std::path::Path,
    limit: Option<usize>,
    chunk_bytes: usize,
    params: &Params,
) -> anyhow::Result<Vec<ReadSketch>> {
    let mut out: Vec<ReadSketch> = Vec::new();
    if path.extension().and_then(|e| e.to_str()) == Some("pod5") {
        // the chunks as ADC counts: the same reads and chunks in half the memory
        signal::for_each_adc_chunk(path, limit, chunk_bytes, |chunk| {
            out.extend(sketch_adc_from(&chunk, out.len() as u32, params));
            Ok(())
        })?;
    } else {
        signal::for_each_raw_chunk(path, limit, chunk_bytes, |chunk| {
            out.extend(sketch_reads_from(&chunk, out.len() as u32, params));
            Ok(())
        })?;
    }
    Ok(out)
}

/// Find candidate overlaps among already-sketched reads.
///
/// Indexes the sketches' seed slices in place. Cloning them into a
/// fresh `Vec<Vec<Seed>>` first would double the seed store for the whole overlap
/// stage: 13 GB on a 411 M-seed D. melanogaster set.
///
/// # Panics
/// With the canonical pair frame on, when two reads share a signal content hash
/// ([`crate::readorder::rank_by_content`]); a caller with untrusted input checks that first.
pub fn overlaps(sketches: &[ReadSketch], params: &Params) -> Vec<Overlap> {
    let reads: Vec<&[Seed]> = sketches.iter().map(|s| s.seeds.as_slice()).collect();
    let index = overlap::SeedIndex::build(reads.iter().copied(), &params.overlap);
    // The reads' lengths in retained events: the opt-in pre-filter's coordinate space
    // (unused when it is off).
    let lens: Vec<u32> = sketches.iter().map(|s| s.n_coord_events).collect();
    // The opt-in canonical pair frame also needs every read's tie rank (by signal content).
    let rank: Option<Vec<u32>> = params.overlap.symmetric_pairs.then(|| {
        let hashes: Vec<u64> = sketches.iter().map(|s| s.content).collect();
        crate::readorder::rank_by_content(&hashes).unwrap_or_else(|e| panic!("{e}"))
    });
    overlap::find_overlap_candidates_frame(
        &reads,
        &index,
        &params.overlap,
        params.sketch.e,
        Some(&lens),
        rank.as_deref(),
    )
}

/// Build the raw fuzzy-Bruijn graph (nodes + edges) from overlaps, WITHOUT any
/// cleaning — an inspectable snapshot of FBG construction (`--raw` / Bandage).
pub fn build_raw_graph(reads: &[RawRead], params: &Params) -> FuzzyGraph {
    let t = std::time::Instant::now();
    let sketches = sketch_reads(reads, params);
    let n_seeds: usize = sketches.iter().map(|s| s.seeds.len()).sum();
    log::info!(
        "sketch: {} reads, {n_seeds} seeds ({:.1} s)",
        sketches.len(),
        t.elapsed().as_secs_f64()
    );
    let t = std::time::Instant::now();
    let ovls = overlaps(&sketches, params);
    log::info!(
        "overlaps: {} ({:.1} s)",
        ovls.len(),
        t.elapsed().as_secs_f64()
    );
    let t = std::time::Instant::now();
    let g = graph::build_graph(&ovls, &params.graph);
    log::info!(
        "graph: {} nodes, {} edges ({:.1} s)",
        g.n_nodes(),
        g.n_edges(),
        t.elapsed().as_secs_f64()
    );
    g
}

/// Assembly pipeline: build the FBG, then clean it. Stops at the cleaned graph —
/// the unitig/contig collapse ([`crate::layout::unitig`]) is not run here.
pub fn assemble(reads: &[RawRead], params: &Params) -> FuzzyGraph {
    let mut g = build_raw_graph(reads, params);
    let t = std::time::Instant::now();
    graph::clean(&mut g, &params.graph);
    log::info!("clean: {:.1} s", t.elapsed().as_secs_f64());
    g
}

/// Build an OLC read **string graph** (alternative layout backend): sketch →
/// overlap → classify/dovetail/transitive-reduce. Returns the reduced graph plus
/// the sketches (for read names + event lengths).
pub fn build_string_graph_layout(
    reads: &[RawRead],
    params: &Params,
    sg: &crate::layout::StringGraphParams,
) -> (crate::layout::StringGraph, Vec<ReadSketch>) {
    let sketches = sketch_reads(reads, params);
    let ovls = overlaps(&sketches, params);
    let len: Vec<u32> = sketches.iter().map(|s| s.n_coord_events).collect();
    let g = crate::layout::build_string_graph(sketches.len(), &ovls, &len, sg);
    (g, sketches)
}
