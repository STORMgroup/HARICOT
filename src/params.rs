//! Pipeline parameters.
//!
//! [`Params::r10_4_1`] is the base configuration. Its signal-processing values (segmentation,
//! normalization, quantization, sketching) are the RawHash2 ONT **R10.4.1** preset, ported from
//! RawHash2:
//!   - `apply_r10_chemistry`        → `extern/RawHash2/src/main.c:138`
//!   - `apply_r10_sensitive_quant`  → `extern/RawHash2/src/main.c:158`
//!   - non-overridden indexing defaults → `extern/RawHash2/src/roptions.c:7`
//!
//! Its overlap and graph values are wtdbg2's where it has an equivalent, and the features that are
//! not in wtdbg2 are off. `haricot assemble` keeps the signal-processing values and sets its overlap
//! and graph values on top: the flags `fbg_build::DETECT_FLAGS` and `fbg_build::GRAPH_FLAGS`
//! (`haricot-detect --production`), among them bin 192, K = 4 bins, the null gate, the partner
//! pre-filter, the relative repeat mask, the canonical pair frame, the relative edge gate and
//! wtdbg2's fragment layout.

/// Strand of a read / seed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Strand {
    /// The read as sequenced; every seed and overlap of the pipeline.
    Forward,
    /// The reverse complement: kept in the seed encoding; the pipeline sketches forward only
    /// (signals are not reverse-complemented).
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

/// Quantization parameters — bucketing of normalized events into a small alphabet, by the
/// faithful port of `dynamic_quantize` (`extern/RawHash2/src/rsketch.c:18`; [`crate::emer::quantize`]).
#[derive(Debug, Clone, Copy)]
pub struct QuantParams {
    /// Number of quantization buckets. Must satisfy `1 <= n_buckets <= 2^q`.
    pub n_buckets: u32,
    /// Bits used per quantized event when packing an e-mer.
    pub q: u32,
    /// Lower edge of the "fine" (high-resolution) region.
    pub fine_min: f32,
    /// Upper edge of the "fine" region.
    pub fine_max: f32,
    /// Fraction of buckets devoted to the fine region.
    pub fine_range: f32,
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
}

/// How [`crate::overlap::dp::hough_aln`] searches for the overlap diagonal's slope: a
/// **prior-windowed beam**. The peak is sought only within a plausible slope window
/// `[win_lo, win_hi]` (the FBG cannot use far-off-slope diagonals anyway); within it a
/// coarse→fine beam keeps the top-`k` *distinct* basins per level and refines each, so the
/// coarse level only has to *nominate* the true basin into a shortlist while the fine level
/// *ranks* it by true apex.
///
/// The window excludes slopes that only chance alignments reach: over ~1.1M labeled pairs
/// (E. coli, S. aureus), a scan over `[0.2, 5.0]` found the best slope outside `[0.7, 1.4]` for
/// 82-84% of random far pairs and for no true overlap.
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
        }
    }
}

/// wtdbg2-style overlap (KBM) parameters. The floors are wtdbg2's `kbm.h` ones
/// (`init_kbmpar`, `kbm.h:269`); the only domain change is that bins
/// and lengths are counted in **events**, not bp (see [`crate::overlap::dp`]).
#[derive(Debug, Clone, Copy)]
pub struct OverlapParams {
    /// Bin length, in events (wtdbg2 `KBM_BSIZE` = 256 bp).
    pub bin_events: u32,
    /// Minimum matched length to accept an overlap, in events (wtdbg2 `min_mat`,
    /// 200 bp).
    pub min_matches: u32,
    /// Maximum matched/aligned ratio — an upper bound that rejects repeat-saturation
    /// artifacts. A collinear overlap covers a fraction of its aligned span with
    /// e-mer matches (similarity ~0.3-0.95 on a fixture of clustered reads, not distributed); a
    /// short repeat shared across loci instead piles up more anchors than the span, so its matched
    /// event length exceeds the aligned span (similarity > 1). On that fixture the
    /// cross-cluster false-positive overlaps fall in the latter regime. Not from
    /// wtdbg2 (which masks repeats via seed-occurrence `kmax`); this is the raw-signal
    /// analog of that repeat guard. `f32::INFINITY` disables it.
    pub max_similarity: f32,
    /// Maximum anchor density = `n_anchors / aligned_bins` (shared seeds per aligned
    /// bin). Measured on that fixture: within-cluster overlaps carry a modest
    /// cloud per aligned bin (p90≈1347, max≈4013), while cross-cluster pairs — a repeat
    /// shared across loci dumping a large anchor cloud onto a short diagonal — sit at
    /// ≥2042, up to 4532. The `max_similarity` (saturation) guard catches only the
    /// densest of these; the density guard rejected all 5 cross-cluster false positives
    /// on that fixture, including ones where the DP aligned only a short weak path.
    /// Like `max_similarity`, this is the raw-signal analog of wtdbg2's `kmax`
    /// seed-occurrence repeat mask. `f32::INFINITY` disables it.
    pub max_anchor_density: f32,
    /// Hough *significance* floor (`z`) of [`crate::overlap::dp::hough_aln`], applied to
    /// the fitted diagonal before any alignment is built. `0` = off.
    ///
    /// A Hough confidence `conf = peak / median(peak over the scanned slopes)` is a RATIO: its numerator
    /// scales with the overlap and its denominator with the read length, so the same
    /// absolute evidence scores lower in longer reads. `z = (obs - e_max)/sqrt(e)` is
    /// a significance and does not have that property. Measured on a D. melanogaster
    /// set (89,162 true overlaps >= 2 kb): `conf` alone rejects **86.8%** of
    /// all true overlaps -- 98% of the 2-6 kb bands and still 76% of 8-12 kb -- while
    /// `null_z`, `extent`, `matched`, `max_similarity` and `max_anchor_density` reject
    /// ZERO, and only ONE true pair fails to become a candidate at all. The pairs
    /// `conf` discards carry a median `z` of 15.8 (8-12 kb) to 18.8 (>20 kb) on
    /// 161-680 anchors. Hence a floor on `z`, not on `conf` (`haricot assemble`: 3).
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
    /// Minimum aligned span in bins (wtdbg2 `min_aln`, 8 = 2048/256).
    pub min_aln_bins: u32,
    /// Max gap-bin fraction of the alignment (wtdbg2 `max_gap`, 0.6).
    pub max_gap: f32,
    /// DEVIATION from wtdbg2 — **null-hypothesis significance gate** (Algorithm S2, line 11;
    /// off in [`Params::r10_4_1`], on in `haricot assemble`). After an overlap's diagonal is
    /// found, require its on-diagonal anchor concentration to clear a hash-frequency-aware
    /// random-seed null: accept only if the Poisson z-score
    /// ([`crate::overlap::dp::diagonal_significance_parts`]) is `>= null_z`. It keeps strong
    /// real diagonals (high anchor count, tight offset) and rejects diffuse / repeat clouds, in
    /// addition to the fixed `max_similarity` / `max_anchor_density` guards.
    pub null_gate: bool,
    /// Minimum diagonal significance z-score to accept an overlap when `null_gate`
    /// is set (`haricot assemble`: 4). Higher = stricter. No effect when `null_gate` is false.
    pub null_z: f32,
    /// The Hough overlapper's slope search (see [`BeamCfg`]).
    pub slope_beam: BeamCfg,
    /// DEVIATION from wtdbg2 — **partner pre-filter** (`--prefilter-nlp L`; `Nominate` of
    /// Algorithm S3). Before
    /// `align_pair`, every gathered pair is scored by its owner with the diagonal statistic
    /// `nlp2n_1` (512-event slope-1 bins, two-bin windows, every key after
    /// the repeat mask, Poisson tail after the window multiplicity; see
    /// [`crate::overlap::prefilter`]) and aligned only when the score reaches `L`. `L` is a
    /// stated chance level (2 = p <= 1e-2 per pair under the uniform null), not a fitted
    /// threshold. `None` (the default) = no pre-filter.
    /// `haricot-detect --production` sets 2.
    pub prefilter_nlp: Option<f64>,
    /// The canonical pair frame (`haricot-detect --production` sets it): every pair is
    /// aligned with the LONGER read as the query (`readorder::query_first`; equal lengths: the
    /// smaller tie rank) and the result is stored in the owner's orientation (transposed exactly
    /// when the owner is
    /// not the query); the pre-filter nominates on the larger of its two frames' scores. A pair's
    /// detection is then a function of the pair, not of which read owns it. Needs the reads'
    /// lengths and tie ranks (`overlap::find_overlap_candidates_frame`). `false` = each pair
    /// aligned in its owner's frame.
    pub symmetric_pairs: bool,
}

/// wtdbg2-style fuzzy de Bruijn graph + layout parameters. Adapted from `wtdbg.c`.
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
    /// product rounds to the absolute floor anyway. `haricot assemble`: 0.20.
    pub min_edge_support: f32,
    /// Record every read's node placements (`graph::threading`) for
    /// the fragment layout ([`crate::layout::frg_contigs`], wtdbg2's second layout level).
    /// `false` (the base value): no placements are recorded. `haricot assemble`: `true`.
    pub frg_layout: bool,
    /// Tip-trimming length, in nodes.
    pub tip_steps: u32,
    /// Bubble-popping search length, in nodes: wtdbg2's popper (`pop_bubble_core_graph`,
    /// `wtdbg-graph.h:2665`, `graph::clean::pop_branched_bubbles`) resolves any region with one
    /// entrance and one exit within this many steps.
    pub bubble_steps: u32,
}

/// Full pipeline configuration.
#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Signal acquisition: sample rate, translocation speed, pA band-pass.
    pub signal: SignalParams,
    /// Event segmentation (RawHash2's dual t-test).
    pub segmentation: SegmentationParams,
    /// Quantization of event values.
    pub quant: QuantParams,
    /// E-mers and minimizer sketching.
    pub sketch: SketchParams,
    /// Seed index, pre-filter and the Hough overlapper.
    pub overlap: OverlapParams,
    /// Fuzzy de Bruijn graph, cleaning and layout.
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
            },
            sketch: SketchParams {
                // roptions.c:7,10 + apply_r10_chemistry k=9 (main.c:140)
                e: 8,
                w: 2,
                k: 9,
                diff: 0,
            },
            overlap: OverlapParams {
                // DEVIATION: wtdbg2's KBM bin is 256 bp (`KBM_BSIZE`, kbm.h:37). Bins
                // here are much finer because the fixture overlaps are short
                // (~457 bp ≈ ~550 events).
                bin_events: 24,
                // wtdbg2 kbm.h:269 (`init_kbmpar`) defaults, in event units, except
                // min_matches: wtdbg2's `min_mat` is 200 bp; 20 events is calibrated
                // against the multi-read null.
                min_matches: 20,
                // Within-cluster overlaps top out near sim≈0.95; repeat-saturation
                // FPs sit at sim≈1.45+. 1.2 splits the gap with margin. See
                // OverlapParams::max_similarity.
                max_similarity: 1.2,
                // Within-cluster density p90≈1347 (max≈4013); cross-cluster repeat
                // FPs ≥2042. 2000 rejected all 5 full-fixture cross FPs, dropping only
                // ~3.4% within (tolerable given overlap redundancy). See
                // OverlapParams::max_anchor_density.
                max_anchor_density: 2000.0,
                // The z floor: off.
                hough_min_z: 0.0,
                max_seed_occ: 1000,
                seed_occ_rel: 0.0,
                min_aln_bins: 8,
                max_gap: 0.6,
                // The null-hypothesis gate: off. Its z = 6 (~6σ over the random-seed null)
                // separated, on a small fixture, the diffuse spurious overlaps (z <= 5.4) from
                // the real diagonals (z >= 27.6).
                null_gate: false,
                null_z: 6.0,
                slope_beam: BeamCfg::default(),
                // The partner pre-filter: off.
                prefilter_nlp: None,
                // The canonical pair frame: off, so every pair is aligned in its owner's frame.
                symmetric_pairs: false,
            },
            graph: GraphParams {
                node_bins: 3,   // K: a node is 3 DP bins (= 72 events at bin_events 24)
                regovl_bins: 1, // wtdbg2 regovl = 1 bin
                min_node_cov: 2,
                node_drop: 0.25,
                max_node_cov: 200,
                // 1 = keep every threaded edge, so a single read's skip edge becomes a
                // graph edge; the relative edge gate (`min_edge_support`) removes weak edges
                // instead. DEVIATION: wtdbg2 auto-sets --edge-min from depth
                // (`wtdbg.c:923-932`: dep<=40 -> 2, >=80 -> 4, else 3).
                min_edge_cov: 1,
                // The relative edge gate: off (`haricot assemble`: 0.20).
                min_edge_support: 0.0,
                frg_layout: false,
                tip_steps: 10,
                bubble_steps: 40,
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r10_preset_is_valid() {
        assert!(Params::r10_4_1().validate().is_ok());
    }

    #[test]
    fn rejects_overflowing_emer() {
        let mut p = Params::r10_4_1();
        p.sketch.e = 30;
        p.quant.q = 3; // 90 bits > 64
        assert!(p.validate().is_err());
    }
}
