//! Pipeline parameters.
//!
//! Defaults follow the RawHash2 ONT **R10.4.1** preset for the signal-processing
//! stages (segmentation, normalization). The concrete values are ported from
//! RawHash2:
//!   - `apply_r10_chemistry`        → `extern/RawHash2/src/main.c:138`
//!   - `apply_r10_sensitive_quant`  → `extern/RawHash2/src/main.c:158`
//!   - non-overridden indexing defaults → `extern/RawHash2/src/roptions.c:7`
//!
//! The overlap and graph defaults are wtdbg2's where it has an equivalent. The pipeline
//! programs set their own values on top (`haricot-detect --production`). Keeping these in
//! one struct makes a parameter change easy: clone [`Params`], tweak a field, re-run.

pub use crate::emer::select::SketchMethod;

/// Strand of a read / seed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Encoding used in seed records (0 = forward, 1 = reverse), matching
    /// RawHash2's low bit of the `y` field.
    pub fn bit(self) -> u32 {
        match self {
            Strand::Forward => 0,
            Strand::Reverse => 1,
        }
    }
}

/// Signal-segmentation (event detection) parameters — the dual t-test detector.
///
/// Ported from `detect_events` / `gen_peaks` in `extern/RawHash2/src/revent.c`.
#[derive(Debug, Clone, Copy)]
pub struct SegmentationParams {
    /// Short detector window length (samples).
    pub window_length1: u32,
    /// Long detector window length (samples).
    pub window_length2: u32,
    /// t-statistic threshold for the short detector.
    pub threshold1: f32,
    /// t-statistic threshold for the long detector.
    pub threshold2: f32,
    /// Minimum prominence for a peak to count.
    pub peak_height: f32,
    /// Skip segments shorter than this (in samples). 0 = no minimum.
    pub min_segment_length: u32,
    /// Skip segments longer than this (in samples).
    pub max_segment_length: u32,
}

/// Quantization scheme — how a normalized event value maps to a bucket.
///
/// [`QuantMethod::Piecewise`] is the RawHash2 design (a fixed piecewise-linear
/// map). Measured on the debug pair, it leaves buckets 0 & 7 nearly empty
/// (entropy 2.52 / 3.0 bits). The *equal-frequency* alternatives split the event
/// distribution so every bucket carries ~equal probability mass, pushing entropy
/// toward `log2(n_buckets)`.
///
/// The bucket boundaries are a *shared codebook*: they must be identical for
/// every read, because cross-read matching is hash equality over bucket sequences
/// ([`crate::emer::sketch::emer_hash`]) — different boundaries map the same signal
/// to different buckets. The schemes that satisfy this by construction (events are
/// z-normalized) are:
///   - [`QuantMethod::NormalQuantile`]: fixed standard-normal quantile cut points
///     (`Φ⁻¹(i/n_buckets)`), data-free and exactly equal-frequency *if* the
///     events were truly `N(0,1)`;
///   - [`QuantMethod::Codebook`]: explicit boundaries fitted once on a *reference*
///     event distribution (e.g. the pooled event means of a read set) and then
///     frozen — see [`QuantParams::boundaries`].
///
/// A *per-read-adaptive* quantile quantizer fits different boundaries per read and
/// so does not satisfy it; there is no distinct variant for it — it is realized by
/// populating [`QuantMethod::Codebook`] with per-read boundaries.
/// A measurement of that control on the debug pair found
/// a mild penalty rather than a collapse: the two reads are from the same organism
/// and run, so their event-mean distributions have near-identical shape and the
/// per-read boundaries nearly coincide with the pooled ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantMethod {
    /// RawHash2's fixed 1-6-1 piecewise-linear map (`dynamic_quantize`).
    Piecewise,
    /// Equal-probability bins of a standard normal — `Φ⁻¹(i/n_buckets)` cut points.
    /// Data-free; exact equal-frequency under the z-normalization assumption.
    NormalQuantile,
    /// Bucket from explicit ascending boundaries in [`QuantParams::boundaries`].
    /// `signal` falls in bucket `b` iff `bounds[b-1] <= signal < bounds[b]`.
    Codebook,
}

/// Maximum number of buckets the boundary codebook can hold (`n_buckets <= 32`).
pub const MAX_BUCKETS: usize = 32;

/// A frozen quantile codebook: `len` ascending cut points (so `len + 1` buckets).
/// Inline fixed-size storage keeps [`QuantParams`] `Copy`.
#[derive(Debug, Clone, Copy)]
pub struct Boundaries {
    pub edges: [f32; MAX_BUCKETS - 1],
    pub len: usize,
}

impl Boundaries {
    /// Empty codebook (no edges → single bucket).
    pub fn empty() -> Self {
        Boundaries {
            edges: [0.0; MAX_BUCKETS - 1],
            len: 0,
        }
    }

    /// Build from a slice of ascending cut points (truncated to `MAX_BUCKETS-1`).
    pub fn from_edges(src: &[f32]) -> Self {
        let mut b = Boundaries::empty();
        b.len = src.len().min(MAX_BUCKETS - 1);
        b.edges[..b.len].copy_from_slice(&src[..b.len]);
        b
    }

    pub fn as_slice(&self) -> &[f32] {
        &self.edges[..self.len]
    }
}

impl PartialEq for Boundaries {
    fn eq(&self, other: &Self) -> bool {
        self.as_slice() == other.as_slice()
    }
}

/// Quantization parameters — bucketing of normalized events into a small alphabet.
///
/// Default ([`QuantMethod::Piecewise`]) is the faithful port of `dynamic_quantize`
/// (`extern/RawHash2/src/rsketch.c:18`). The equal-frequency variants are
/// alternatives, not used by `haricot assemble` (see [`QuantMethod`]).
#[derive(Debug, Clone, Copy)]
pub struct QuantParams {
    /// Number of quantization buckets. Must satisfy `1 <= n_buckets <= 2^q`.
    pub n_buckets: u32,
    /// Bits used per quantized event when packing an e-mer.
    pub q: u32,
    /// Lower edge of the "fine" (high-resolution) region. (Piecewise only.)
    pub fine_min: f32,
    /// Upper edge of the "fine" region. (Piecewise only.)
    pub fine_max: f32,
    /// Fraction of buckets devoted to the fine region. (Piecewise only.)
    pub fine_range: f32,
    /// Which quantization scheme to apply.
    pub method: QuantMethod,
    /// Frozen ascending cut points for [`QuantMethod::Codebook`]. Ignored otherwise.
    pub boundaries: Boundaries,
}

/// E-mer / sketching parameters.
#[derive(Debug, Clone, Copy)]
pub struct SketchParams {
    /// Number of quantized events concatenated into one e-mer.
    pub e: u32,
    /// Minimizer window size (in e-mers). `0` disables minimizing (dense sketch).
    pub w: u32,
    /// The pore model's k-mer length (RawHash2's R10.4.1 `k = 9`). Informational: a seed
    /// nominally spans `k + e - 1` bases, but nothing in the pipeline reads it.
    pub k: u32,
    /// Quantized sig-diff filter: skip an event if `|quant - prev_quant| <= diff`.
    /// `-1` disables the filter; `0` collapses runs of identical quantized values.
    pub diff: i32,
    /// Seed-selection scheme (alternatives to the default are not from a verified preset).
    /// Default [`SketchMethod::Minimizer`] is RawHash2's scheme; the indel-robust
    /// alternatives are in `emer::select`.
    pub method: SketchMethod,
    /// s-mer length (events) for syncmers/strobemers. Ignored by dense/minimizer.
    pub s: u32,
}

/// Signal acquisition parameters (chemistry-dependent).
#[derive(Debug, Clone, Copy)]
pub struct SignalParams {
    /// Sampling rate of the sequencer (Hz).
    pub sample_rate: u32,
    /// Translocation speed (bases per second).
    pub bp_per_sec: u32,
    /// Keep only pA values within this inclusive range (RawHash2 uses 30..=200).
    /// Set to `None` to disable the pA band-pass filter.
    pub pa_filter: Option<(f32, f32)>,
    /// Pre-quantization float event merge: collapse consecutive z-normalized event
    /// means that stay within this threshold (z-units) of the running mean into a
    /// single event. Intended to damp the variable-count segmentation of a sustained
    /// signal level (homopolymer-like), so overlapping reads produce more comparable
    /// event streams. `0.0` (default) disables it, which is the RawHash2-faithful path
    /// the `hash_compat` goldens are pinned to. A typical aggressive value is ~0.3.
    pub event_merge: f32,
}

impl SignalParams {
    /// Samples per base = `sample_rate / bp_per_sec` (12.5 for R10.4.1).
    pub fn sample_per_base(&self) -> f32 {
        self.sample_rate as f32 / self.bp_per_sec as f32
    }
}

/// How [`crate::overlap::dp::hough_aln`] searches for the overlap diagonal's slope: a
/// **prior-windowed beam**. The peak is sought only within a plausible slope window
/// `[win_lo, win_hi]` (the FBG cannot use far-off-slope diagonals anyway); within it a
/// coarse→fine beam keeps the top-`k` *distinct* basins per level and refines each, so the
/// coarse level only has to *nominate* the true basin into a shortlist while the fine level
/// *ranks* it by true apex. The `hough_min_conf` background is still sampled over the full
/// `[0.2, 5.0]` range, coarsely.
///
/// It replaced a 241-slope scan over `[0.2, 5.0]` that was slower and less specific:
/// measured over ~1.1M labeled pairs, the out-of-window slopes were won by 82% (E. coli) /
/// 84% (S. aureus) of random far pairs and 0.0% of true overlaps.
///
/// Defaults give a top-3 two-level `[0.7, 1.4]` search; raise `beam_k` for more robustness
/// to competing diagonals and `levels` for a coarser hierarchical start.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamCfg {
    /// Low edge of the plausible slope window; the peak is only sought at `m >= win_lo`.
    pub win_lo: f32,
    /// High edge of the plausible slope window; the peak is only sought at `m <= win_hi`.
    pub win_hi: f32,
    /// Coarsest slope step within the window (level 0). The schedule refines
    /// geometrically from here down to ~`0.02` over `levels` levels.
    pub coarse_step: f32,
    /// Beam width: top-`k` **distinct** basins carried from each level to the next.
    /// `1` = classic single-winner coarse→fine (loses the true basin to a competitor).
    pub beam_k: u32,
    /// Number of coarse→fine levels (`>= 1`). `2` = coarse then a `0.02` refine; more
    /// levels start coarser and prune the beam gradually (a scale-space hierarchy).
    pub levels: u32,
    /// First filter (`--hough-screen`, `0` = off): reject a pair whose densest `±tol`
    /// offset window over the level-0 slopes holds fewer than this many anchors, before
    /// any refinement. On 4 kb reads the level-0 peak IS the final Hough count for a true
    /// overlap (a coarse step of 0.2 moves a 3,000-event diagonal less than one window),
    /// so the floor removes random pairs whose refinement and null scan would be thrown
    /// away. Measured at floor 15: 0.009-0.074% of accepted true pairs lost, 37-73% of the
    /// alignment time saved on 4 kb / H. pylori sets, 9-18% of accepted FALSE pairs removed
    /// as a side effect.
    pub screen_peak: u32,
    /// First filter, density-normalised (`--hough-screen-z`, `0` = off): reject when the
    /// level-0 peak standardised against a uniform cloud, `(pk0 - lam) / sqrt(max(lam, 1))`
    /// with `lam = n * 2*tol / (a_ext + b_ext)`, is under this. Unlike the raw count it
    /// transfers to D. melanogaster, whose sparse clouds give true overlaps lower counts:
    /// at 5.0, 0.007-0.022% of accepted true pairs lost on all six datasets measured (the
    /// count floor 15 loses 0.62% on D. melanogaster), 51-53% of the alignment time saved on
    /// 4 kb sets, 34% on H. pylori 12x, 15% on D. melanogaster.
    pub screen_z: f32,
    /// Stage 0b, against the gate's own null (`--hough-screen-null`, `0` = off): after level
    /// 0, score the best coarse slope exactly as the final gate scores the fine one,
    /// `(pk0 - e_max) / sqrt(max(e, 1))` ([`crate::overlap::dp::diagonal_significance_parts`]),
    /// and reject when under this. Past `screen_z`, 99-100% of the pairs still refined are
    /// rejected by that gate anyway; this asks its question before the refinement.
    pub screen_null: f32,
    /// Choose the slope by the gate's own statistic (`--slope-select-z`, default off).
    /// Off, the winner is the basin with the largest raw densest-window count, and `z` is
    /// then taken at that one slope. A dense off-diagonal repeat block spreads its offsets
    /// over `a_ext + m*b_ext`, so its count rises as `m` falls and can win at `win_lo`,
    /// where `z` is ~0: the pair is rejected although a `z` > 20 diagonal lies inside the
    /// window, and the verdict depends on which read has the lower id (measured: 2.8-4.5% of
    /// the true pairs of at least 1 kb at bacterial 6x). On, the last level keeps its
    /// top-`beam_k` distinct basins, refits each like the winner, and keeps the one with
    /// the largest null `z` ([`crate::overlap::dp::diagonal_significance_parts`]); ties keep
    /// the higher count.
    pub select_z: bool,
}

impl Default for BeamCfg {
    /// Prior-windowed beam: seek the peak only in `[0.7, 1.4]`, top-3 basins, coarse
    /// `0.2` refined to `0.02` over two levels.
    fn default() -> Self {
        BeamCfg {
            win_lo: 0.7,
            win_hi: 1.4,
            coarse_step: 0.2,
            beam_k: 3,
            levels: 2,
            screen_peak: 0,
            screen_z: 0.0,
            screen_null: 0.0,
            select_z: false,
        }
    }
}

/// wtdbg2-style overlap (KBM) parameters. The bin-DP knobs are ported faithfully
/// from `kbm.h` (`init_kbmpar`, `kbm.h:269`); the only domain change is that bins
/// and lengths are counted in **events**, not bp (see [`crate::overlap::dp`]).
#[derive(Debug, Clone, Copy)]
pub struct OverlapParams {
    /// Bin length, in events (wtdbg2 `KBM_BSIZE` = 256 bp).
    pub bin_events: u32,
    /// Minimum matched length to accept an overlap, in events (wtdbg2 `min_mat`,
    /// 200 bp).
    pub min_matches: u32,
    /// Minimum matched/aligned ratio (wtdbg2 `min_sim`, 0.05).
    pub min_similarity: f32,
    /// Maximum matched/aligned ratio — an upper bound that rejects repeat-saturation
    /// artifacts. A collinear overlap covers a fraction of its aligned span with
    /// e-mer matches (similarity ~0.3-0.95 on the cluster fixture); a short repeat
    /// shared across loci instead piles up more anchors than the span, so its matched
    /// event length exceeds the aligned span (similarity > 1). On that fixture the
    /// cross-cluster false-positive overlaps fall in the latter regime. Not from
    /// wtdbg2 (which masks repeats via seed-occurrence `kmax`); this is the raw-signal
    /// analog of that repeat guard. `f32::INFINITY` disables it.
    pub max_similarity: f32,
    /// Maximum anchor density = `n_anchors / aligned_bins` (shared seeds per aligned
    /// bin). Measured on the cluster fixture: within-cluster overlaps carry a modest
    /// cloud per aligned bin (p90≈1347, max≈4013), while cross-cluster pairs — a repeat
    /// shared across loci dumping a large anchor cloud onto a short diagonal — sit at
    /// ≥2042, up to 4532. The `max_similarity` (saturation) guard catches only the
    /// densest of these; the density guard rejected all 5 cross-cluster false positives
    /// on that fixture, including ones where the DP aligned only a short weak path.
    /// Like `max_similarity`, this is the raw-signal analog of wtdbg2's `kmax`
    /// seed-occurrence repeat mask. `f32::INFINITY` disables it.
    pub max_anchor_density: f32,
    /// DEVIATION from wtdbg2: use the global Hough-extent overlapper
    /// ([`crate::overlap::dp::hough_aln`]) instead of the banded bin-DP. It detects
    /// the diagonal with a 1-D Hough scan and takes the whole on-line extent, so it
    /// recovers fuller alignments (denser FBG tiling) and rejects diffuse
    /// repeat/cross-cluster clouds on Hough significance. `false` = the bin-DP.
    pub hough: bool,
    /// Hough significance (peak/median) an overlap must clear to be accepted, when
    /// `hough` is true. Higher → stricter (fewer cross-cluster false positives).
    pub hough_min_conf: f32,
    /// Hough *significance* floor (`z`), applied at the same early-exit point as
    /// `hough_min_conf` and to the same fitted diagonal. `0` = off (the baseline).
    ///
    /// `conf = peak / median(peak over the scanned slopes)` is a RATIO: its numerator
    /// scales with the overlap and its denominator with the read length, so the same
    /// absolute evidence scores lower in longer reads. `z = (obs - e_max)/sqrt(e)` is
    /// a significance and does not have that property. Measured on a D. melanogaster
    /// set (89,162 true overlaps >= 2 kb): `conf` alone rejects **86.8%** of
    /// all true overlaps -- 98% of the 2-6 kb bands and still 76% of 8-12 kb -- while
    /// `null_z`, `extent`, `matched`, `max_similarity` and `max_anchor_density` reject
    /// ZERO, and only ONE true pair fails to become a candidate at all. The pairs
    /// `conf` discards carry a median `z` of 15.8 (8-12 kb) to 18.8 (>20 kb) on
    /// 161-680 anchors, against the default `null_z` floor of 6.
    ///
    /// Gating here rather than after the DP is deliberate: `hough_aln` returns on a
    /// failed `conf` before building any alignment, so moving the decision keeps that
    /// early exit and the cost profile with it.
    pub hough_min_z: f32,
    /// Drop seeds occurring more than this many times (wtdbg2 `kmax` repeat mask).
    pub max_seed_occ: u32,
    /// Scale term of the repeat mask: raise the cutoff to this many times the input's MEAN
    /// seed occurrence (seed instances / distinct hashes, measured on the index) when that
    /// exceeds [`Self::max_seed_occ`]. `0.0` = off (the default: the fixed cap alone).
    ///
    /// A hash's occurrence grows with genome size x depth, so one absolute cap keeps a
    /// different share of the seeds on every input: 2000 keeps 53% of the seed instances
    /// on E. coli 30x and dmel, but 26% on A. thaliana chr4 30x and 19.5% on chr1 30x,
    /// where true overlaps are left with a fifth of their anchors. In units of the mean
    /// the cutoff is scale-free. At `17.0` (about E. coli 30x's 2000 / 116) it binds only
    /// on inputs larger than every training set, which keep the fixed cap exactly.
    /// See `overlap::kbm::effective_seed_cap`.
    pub seed_occ_rel: f32,
    /// Per-gap-bin score (wtdbg2 `pgap`, −7).
    pub pgap: i32,
    /// Per-bin diagonal-deviation penalty (wtdbg2 `pvar`, −21).
    pub pvar: i32,
    /// DP band half-width in bins (wtdbg2 `max_bvar`, 4).
    pub max_bvar: u32,
    /// Max consecutive gap bins on the diagonal (wtdbg2 `max_bgap`, 4).
    pub max_bgap: u32,
    /// Minimum aligned span in bins (wtdbg2 `min_aln`, 8 = 2048/256).
    pub min_aln_bins: u32,
    /// Max gap-bin fraction of the alignment (wtdbg2 `max_gap`, 0.6).
    pub max_gap: f32,
    /// Max span-length variation (wtdbg2 `aln_var`, 0.25).
    pub aln_var: f32,
    /// DEVIATION from wtdbg2 — **null-hypothesis significance gate** (opt-in). After
    /// an overlap's diagonal is found, require its on-diagonal anchor concentration
    /// to clear a hash-frequency-aware random-seed null: accept only if the Poisson
    /// z-score ([`crate::overlap::dp::diagonal_significance`]) is `>= null_z`. It
    /// keeps strong real diagonals (high anchor count, tight offset) and rejects
    /// diffuse / repeat clouds, augmenting the fixed `max_similarity` /
    /// `max_anchor_density` guards rather than replacing them. `false` = baseline
    /// (byte-identical overlaps); this is the general overlap-precision lever, not a
    /// fix for any one dataset.
    pub null_gate: bool,
    /// Minimum diagonal significance z-score to accept an overlap when `null_gate`
    /// is set. Higher = stricter (more precision, less recall). No effect when
    /// `null_gate` is false.
    pub null_z: f32,
    /// The Hough overlapper's slope search (see [`BeamCfg`]). No effect when `hough` is
    /// false.
    pub slope_beam: BeamCfg,
    /// DEVIATION from wtdbg2 — **opt-in partner pre-filter** (`--prefilter-nlp L`). Before
    /// `align_pair`, every gathered pair is scored by its owner with the diagonal statistic
    /// `nlp2n_1` (512-event slope-1 bins, two-bin windows, every key after
    /// the repeat mask, Poisson tail after the window multiplicity; see
    /// [`crate::overlap::prefilter`]) and aligned only when the score reaches `L`. `L` is a
    /// stated chance level (2 = p <= 1e-2 per pair under the uniform null), not a fitted
    /// threshold. `None` (the default) = no pre-filter.
    /// `haricot-detect --production` sets 2.
    pub prefilter_nlp: Option<f64>,
    /// The canonical pair frame (opt-in; `haricot-detect --production` sets it): every pair is
    /// aligned with the LONGER read as the query (`readorder::query_first`; equal lengths: the
    /// smaller tie rank) and the result is stored in the owner's orientation (transposed exactly
    /// when the owner is
    /// not the query); the pre-filter nominates on the larger of its two frames' scores. A pair's
    /// detection is then a function of the pair, not of which read owns it. Needs the reads'
    /// lengths and tie ranks (`overlap::find_overlap_candidates_frame`). `false` = each pair
    /// aligned in its owner's frame.
    pub symmetric_pairs: bool,
}

/// wtdbg2-style fuzzy-Bruijn-graph + layout parameters. Adapted from `wtdbg.c`.
#[derive(Debug, Clone, Copy)]
pub struct GraphParams {
    /// Node (K-bin) size, in **DP bins**: a node spans `node_bins * bin_events`
    /// events. This is K, and it is the source of truth — node size is *defined* as
    /// a whole number of bins, never re-derived by dividing an event count by
    /// `bin_events`.
    ///
    /// wtdbg2 is the same underneath: `--node-len` is a bp-valued flag, but it is
    /// forced to a multiple of the bin size and converted once at startup
    /// (`wtdbg.c:638-640`, `:960`: `g->reglen = reglen / KBM_BIN_SIZE`), after which
    /// `reglen` *is* K in bins (default 1024/256 = 4). Deriving K per build instead
    /// silently reinterprets it whenever `bin_events` changes — e.g. a 72-event node
    /// at `bin_events=96` rounds to K=1, a 96-event node declared as 72.
    pub node_bins: u32,
    /// Region overlap in **bins**: two regions on the same read conflict once they
    /// overlap by more than this (wtdbg2 `regovl`, `wtdbg.h:1000`). wtdbg2's default
    /// is 256 bp at a 256 bp bin, i.e. **1 bin**; expressing it in bins keeps that
    /// invariant at every bin size instead of letting it floor to 0.
    pub regovl_bins: u32,
    /// Minimum node coverage (supporting reads). wtdbg2 `min_node_cov`.
    pub min_node_cov: u32,
    /// Drop a candidate node if more than this fraction of its regions conflict
    /// (wtdbg2 `node_drop` / `node_max_conflict`, 0.25).
    pub node_drop: f32,
    /// Mask nodes with raw coverage above this (wtdbg2 `max_node_cov`, 200).
    pub max_node_cov: u32,
    /// Minimum edge coverage (supporting reads).
    pub min_edge_cov: u32,
    /// **Relative** edge-support gate: drop an edge whose coverage is below
    /// `min_edge_support · min(cov(from), cov(to))`. `0.0` = off.
    ///
    /// A true chain edge is threaded by most of the reads covering both of its
    /// endpoints; a spurious junction is threaded by a handful. The absolute
    /// [`Self::min_edge_cov`] cannot express that, because it has no idea what the
    /// local depth is: coverage 2 is strong evidence at depth 4 and nearly none at
    /// depth 40. This is the same reasoning wtdbg2 applies by deriving
    /// `--edge-min` from depth (`wtdbg.c:923-932`), made local and continuous.
    ///
    /// Degenerate below ~10x depth, where `min(cov)` is small enough that the
    /// product rounds to the absolute floor anyway.
    pub min_edge_support: f32,
    /// Record read-threaded node triples (`graph::threading`). The FBG builds a path
    /// per read and then reduces it to consecutive pairs; the longer context is thrown
    /// away exactly where it is the only thing that can resolve a repeat. Default
    /// `false`, so nothing is recorded and the graph is unchanged.
    pub read_threading: bool,
    /// Reads that must thread a node triple before a layout will traverse it.
    ///
    /// `1` (the default) trusts a single molecule, which may be chimeric or misplaced.
    /// This is corroboration, not a tuned score: it counts independent observations, so
    /// it does not carry a dataset-specific scale the way a threshold on a statistic
    /// does. wtdbg2's `--edge-min` is the same idea one level down.
    pub thread_min_reads: u32,
    /// Record every read's node placements (as `read_threading` does, without its triples) for
    /// the fragment layout ([`crate::layout::frg_contigs`], wtdbg2's second layout level).
    /// Default `false`: nothing is recorded and every graph and layout is unchanged.
    pub frg_layout: bool,
    /// Max edge reach in events — how far back on a read to connect regions
    /// (wtdbg2 `max_edge_span`, 1024 bins).
    pub max_edge_span_events: u32,
    /// Tip-trimming length, in nodes.
    pub tip_steps: u32,
    /// Bubble-popping search length, in nodes.
    pub bubble_steps: u32,
    /// Bubble popper: `0` = the simple-bubble port [`crate::graph::clean::pop_bubbles`]
    /// (the default), `1` = wtdbg2's own popper in its safe mode, `2` = safe passes, then
    /// its unsafe mode once they converge, as `wtdbg.c:1080-1110` runs them. wtdbg2's
    /// popper (`pop_bubble_core_graph`, `wtdbg-graph.h:2665`) resolves any region with one
    /// entrance and one exit within [`Self::bubble_steps`] -- branched arms, braids --
    /// where the simple port only pops two unbranched chains. See
    /// `graph::clean::pop_branched_bubbles`.
    pub branched_bubbles: u8,
}

/// Full pipeline configuration.
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub signal: SignalParams,
    pub segmentation: SegmentationParams,
    pub quant: QuantParams,
    pub sketch: SketchParams,
    pub overlap: OverlapParams,
    pub graph: GraphParams,
}

impl Default for Params {
    /// The ONT R10.4.1 preset.
    fn default() -> Self {
        Params::r10_4_1()
    }
}

impl Params {
    /// ONT R10.4.1 preset.
    pub fn r10_4_1() -> Self {
        Params {
            signal: SignalParams {
                // apply_r10_chemistry: main.c:152-155
                sample_rate: 5000,
                bp_per_sec: 400,
                // pA band-pass filter from rsig.c read paths.
                pa_filter: Some((30.0, 200.0)),
                event_merge: 0.0,
            },
            segmentation: SegmentationParams {
                // apply_r10_chemistry: main.c:142-144
                window_length1: 3,
                window_length2: 6,
                threshold1: 5.5,
                threshold2: 3.5,
                peak_height: 0.2,
                // roptions.c:135-136
                min_segment_length: 0,
                max_segment_length: 500,
            },
            quant: QuantParams {
                // apply_r10_sensitive_quant: main.c:160-164
                n_buckets: 8,
                q: 3,
                fine_min: -2.0,
                fine_max: 2.0,
                fine_range: 0.75,
                method: QuantMethod::Piecewise,
                boundaries: Boundaries::empty(),
            },
            sketch: SketchParams {
                // roptions.c:7,10 + apply_r10_chemistry k=9 (main.c:140)
                e: 8,
                w: 2,
                k: 9,
                diff: 0,
                method: SketchMethod::Minimizer,
                s: 4,
            },
            overlap: OverlapParams {
                // DEVIATION: wtdbg2's KBM bin is 256 bp (`KBM_BSIZE`, kbm.h:37). Bins
                // here are much finer because the fixture overlaps are short
                // (~457 bp ≈ ~550 events).
                bin_events: 24,
                // wtdbg2 kbm.h:269 (`init_kbmpar`) defaults, in event units, except
                // min_matches: wtdbg2's `min_mat` is 200 bp; 20 events is calibrated
                // against the multi-read null. min_similarity keeps
                // wtdbg2's `min_sim`.
                min_matches: 20,
                min_similarity: 0.05,
                // Within-cluster overlaps top out near sim≈0.95; repeat-saturation
                // FPs sit at sim≈1.45+. 1.2 splits the gap with margin. See
                // OverlapParams::max_similarity.
                max_similarity: 1.2,
                // Within-cluster density p90≈1347 (max≈4013); cross-cluster repeat
                // FPs ≥2042. 2000 rejected all 5 full-fixture cross FPs, dropping only
                // ~3.4% within (tolerable given overlap redundancy). See
                // OverlapParams::max_anchor_density.
                max_anchor_density: 2000.0,
                // hough_min_conf calibrated on overlap_clusters_long: chimeric
                // cross-cluster diagonal confidences top out at ≈5.5, so 6.0 gave ~0
                // cross-cluster FP at ~0.84 within-cluster recall.
                hough: false,
                hough_min_conf: 6.0,
                // Off by default: the baseline gate is `hough_min_conf`.
                hough_min_z: 0.0,
                max_seed_occ: 1000,
                seed_occ_rel: 0.0,
                pgap: -7,
                pvar: -21,
                max_bvar: 4,
                max_bgap: 4,
                min_aln_bins: 8,
                max_gap: 0.6,
                aln_var: 0.25,
                // Opt-in null-hypothesis gate, OFF by default so the baseline overlap
                // set is byte-identical. When enabled, z=6 (~6σ over the random-seed
                // null) is calibrated on two fixtures: on braid_window it removes both
                // basecall-confirmed false-positive pairs (precision 0.999→1.000) at a
                // 0.16-pt recall cost; on an S. aureus D7 cross-strand set it rejects all 3
                // diffuse `spurious` overlaps while keeping all 140 real fold-back
                // diagonals (whose z ≥ 27) and 6 true controls. The safe band is wide
                // (spurious ones top out at z≈5.4, real diagonals start at z≈27.6).
                null_gate: false,
                null_z: 6.0,
                slope_beam: BeamCfg::default(),
                // Opt-in partner pre-filter: off, so the candidate set is unchanged.
                prefilter_nlp: None,
                // Opt-in canonical pair frame: off, so every pair is aligned in its owner's frame.
                symmetric_pairs: false,
            },
            graph: GraphParams {
                node_bins: 3,   // K: a node is 3 DP bins (= 72 events at bin_events 24)
                regovl_bins: 1, // wtdbg2 regovl = 1 bin
                min_node_cov: 2,
                node_drop: 0.25,
                max_node_cov: 200,
                // 1 = keep every threaded edge, so a single read's skip edge becomes a
                // permanent graph edge. wtdbg2 instead auto-sets --edge-min from depth
                // (`wtdbg.c:923-932`: dep<=40 -> 2, >=80 -> 4, else 3). UNRESOLVED: the
                // measurement that justified 1 (cov>=2 shatters the graph) predates the
                // edge-orientation fix and does not reproduce; settle it with a
                // QUAST/NGA50 eval rather than structural proxies. See TODO.md.
                min_edge_cov: 1,
                // Off by default: the cleaned graph stays byte-identical to the
                // baseline until a run opts in. Target under evaluation is 0.20.
                min_edge_support: 0.0,
                read_threading: false,
                thread_min_reads: 1,
                frg_layout: false,
                // 0 = bridge consecutive committed regions on a read regardless
                // of gap (a read spans its own DNA). See graph::build Stage 4. A
                // 1024-event cap was measured to fragment each cluster into 8-14
                // components, but on the pre-fix graph — not re-checked since.
                max_edge_span_events: 0,
                tip_steps: 10,
                bubble_steps: 40,
                branched_bubbles: 0,
            },
        }
    }

    /// Validate cross-field invariants asserted in the RawHash2 sketch code.
    ///
    /// (`e * q` must fit the packed e-mer; `n_buckets` must fit in `q` bits.)
    pub fn validate(&self) -> Result<(), String> {
        let (e, q) = (self.sketch.e, self.quant.q);
        if e * q > 64 {
            return Err(format!("e * q = {e} * {q} = {} exceeds 64 bits", e * q));
        }
        if self.quant.n_buckets < 1 || self.quant.n_buckets > (1u32 << q) {
            return Err(format!(
                "n_buckets ({}) must be in 1..=2^q ({})",
                self.quant.n_buckets,
                1u32 << q
            ));
        }
        if self.overlap.bin_events < 1 {
            return Err("bin_events must be >= 1".to_string());
        }
        if self.graph.node_bins < 1 {
            return Err("node_bins (K) must be >= 1".to_string());
        }
        if let Some(l) = self.overlap.prefilter_nlp
            && !l.is_finite()
        {
            return Err(format!("prefilter_nlp ({l}) must be finite"));
        }
        Ok(())
    }

    /// Node size in events — the *derived* quantity, `node_bins * bin_events`. Use
    /// this (or [`crate::graph::fbg::Node::len_events`]) rather than re-deriving a
    /// node size anywhere; see [`GraphParams::node_bins`].
    pub fn node_events(&self) -> u32 {
        self.graph.node_bins * self.overlap.bin_events
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r10_preset_is_valid() {
        assert!(Params::r10_4_1().validate().is_ok());
    }

    #[test]
    fn sample_per_base_is_12_5() {
        assert_eq!(Params::r10_4_1().signal.sample_per_base(), 12.5);
    }

    #[test]
    fn rejects_overflowing_emer() {
        let mut p = Params::r10_4_1();
        p.sketch.e = 30;
        p.quant.q = 3; // 90 bits > 64
        assert!(p.validate().is_err());
    }
}
