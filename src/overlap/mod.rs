//! Read-to-read overlap detection (the wtdbg2 KBM analog).
//!
//! Stages:
//!   1. [`SeedIndex`] — invert seeds by hash, repeat-masked (implemented in [`kbm`]).
//!   2. [`find_overlap_candidates`] — for every read, gather its shared-seed anchors
//!      against every higher-id read, grouped per pair.
//!   3. [`align_pair`] — per pair, [`dp::hough_aln`] (as `haricot assemble` runs it) or the
//!      faithful wtdbg2 banded bin-DP [`dp::bin_dp`], then the repeat guards and the
//!      null gate. Opt-in ([`OverlapParams::prefilter_nlp`]): a pair is first scored by
//!      [`prefilter`] and aligned only when it is nominated.
//!
//! The bin-DP ([`dp`]) is a direct port of wtdbg2's `map_kbm` recurrence
//! (`kbm.h:1638`): seed matches are floored into bin cells and a banded
//! Smith-Waterman-shaped DP integrates the diagonal, accepting on
//! matched/aligned (`min_sim`).

pub mod dp;
pub mod kbm;
mod kernel;
pub mod prefilter;
pub mod read_model;

pub use kbm::{SeedCounts, SeedHit, SeedIndex};
pub use kernel::kernel_name;

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::emer::Seed;
use crate::overlap::dp::BinDp;
use crate::params::{OverlapParams, Strand};

/// A shared-seed match between two reads: `(a_pos on read a, b_pos on read b)`,
/// in event coordinates.
#[derive(Debug, Clone, Copy)]
pub struct Anchor {
    pub a_pos: u32,
    pub b_pos: u32,
}

/// A scored overlap between two reads, carrying the bin alignment so the FBG
/// (`graph::build_graph`) can slice it into K-bins.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Overlap {
    /// Lower read id.
    pub a: u32,
    /// Higher read id.
    pub b: u32,
    /// Relative strand. Always forward: [`anchor_clouds`] pairs same-strand seed
    /// hits only, so reverse-strand overlaps are not produced yet.
    pub strand: Strand,
    /// Total shared seeds (anchors) between the two reads.
    pub n_anchors: u32,
    /// Matched event length on the alignment path (wtdbg2 `mat`).
    pub matched: u32,
    /// Seed-match count on the path (wtdbg2 `cnt`).
    pub cnt: u32,
    /// Aligned span in bins (wtdbg2 `aln`).
    pub aligned_bins: u32,
    /// Bin size (events) this overlap was aligned at — for event↔bin conversion.
    pub bin_events: u32,
    /// Gap bins on the path (wtdbg2 `gap`).
    pub gap: u32,
    /// Slope the DP ran at (1.0 = baseline; ≠1 when slope-corrected).
    pub slope: f64,
    /// Query (read `a`) bin span `(lo, hi)`.
    pub a_bin_span: (u32, u32),
    /// Target (read `b`) bin span `(lo, hi)` in real bin coords.
    pub b_bin_span: (u32, u32),
    /// `(a_bin, b_bin)` cells along the alignment (occupied cells, real coords).
    pub path: Vec<(u32, u32)>,
    /// Diagonal significance at the accepted fit (`dp::diagonal_significance`).
    ///
    /// Persisted so a single permissive overlap pass can be re-thresholded offline at
    /// any operating point, instead of one expensive pass per gate -- the anchor work
    /// is identical for every gate, only the decision differs. It is also the per-edge
    /// evidence weight the graph stage needs: `build_graph` otherwise sees only
    /// coverage COUNTS, so a barely-accepted overlap and an overwhelming one vote
    /// alike.
    ///
    /// `NaN` means "not measured" -- synthesized overlaps (the oracle, tests) carry it,
    /// and every `z <` comparison is false for `NaN`, so those pass any filter.
    pub z: f32,
}

impl Overlap {
    /// Similarity = matched / aligned (wtdbg2 `min_sim`), in event units.
    pub fn similarity(&self) -> f32 {
        let aligned_events = (self.aligned_bins.max(1) * self.bin_events) as f32;
        self.matched as f32 / aligned_events
    }
}

/// Read pair `(a, b)` → the shared-seed anchors between them.
type PairClouds = FxHashMap<(u32, u32), Vec<Anchor>>;

/// Anchors gathered for one read before they are grouped by partner. Above this a read
/// is gathered in several slices of the partner id range, so a read carrying a huge
/// cloud (a long repetitive read against the whole set) cannot make one worker hold
/// gigabytes: 8M anchors is ~160 MB of scratch, and a typical read holds well under 1M.
const MAX_READ_ANCHORS: u64 = 8 << 20;

/// Per-worker scratch for [`gather`], reused across reads.
struct Scratch {
    /// One `(partner, anchor)` per shared seed occurrence, in the read's seed order.
    raw: Vec<(u32, Anchor)>,
    /// The same anchors grouped by partner (stable), and each partner's slice of them.
    anchors: Vec<Anchor>,
    pairs: Vec<(u32, u32, u32)>,
    /// Per-partner counts, then write cursors; zero outside a call.
    cnt: Vec<u32>,
}

impl Scratch {
    fn new(n_reads: u32) -> Self {
        Scratch {
            raw: Vec::new(),
            anchors: Vec::new(),
            pairs: Vec::new(),
            cnt: vec![0; n_reads as usize],
        }
    }
}

/// Gather read `a`'s anchors against every read with a higher id in `[b_lo, b_hi)`,
/// grouped per partner, into `sc.anchors` / `sc.pairs` (partners ascending, `(b, start,
/// end)` slices).
///
/// Enumerates exactly the hit pairs the index holds: for every seed of `a` whose hash
/// survived masking, every same-strand occurrence on a higher-id read — the lower id
/// owns the pair, so each `(i, j)` hit pair is visited once over the whole run. Within a
/// pair the anchors keep the order of `a`'s seeds, so a read sketched in position order
/// (every sketch is) hands [`dp`] clouds sorted by `a_pos`.
fn gather(
    a: u32,
    seeds: &[Seed],
    hits: &[(&[SeedHit], u32)],
    b_lo: u32,
    b_hi: u32,
    sc: &mut Scratch,
) {
    sc.raw.clear();
    for (s, &(entry, start)) in seeds.iter().zip(hits) {
        let (lo, hi) = if b_lo <= a + 1 && b_hi == u32::MAX {
            (start as usize, entry.len())
        } else {
            (
                entry.partition_point(|h| h.read_id() < b_lo.max(a + 1)),
                entry.partition_point(|h| h.read_id() < b_hi),
            )
        };
        for y in &entry[lo..hi] {
            // Reverse-strand support would need revcomp e-mer hashing; not handled.
            if y.strand() == s.strand() {
                sc.raw.push((
                    y.read_id(),
                    Anchor {
                        a_pos: s.pos,
                        b_pos: y.pos,
                    },
                ));
            }
        }
    }
    // Stable counting sort by partner, touching only the partners seen.
    sc.pairs.clear();
    for &(b, _) in &sc.raw {
        let c = &mut sc.cnt[b as usize];
        if *c == 0 {
            sc.pairs.push((b, 0, 0));
        }
        *c += 1;
    }
    sc.pairs.sort_unstable_by_key(|p| p.0);
    let mut off = 0u32;
    for p in sc.pairs.iter_mut() {
        let c = sc.cnt[p.0 as usize];
        *p = (p.0, off, off + c);
        sc.cnt[p.0 as usize] = off;
        off += c;
    }
    sc.anchors.clear();
    sc.anchors
        .resize(sc.raw.len(), Anchor { a_pos: 0, b_pos: 0 });
    for &(b, an) in &sc.raw {
        let c = &mut sc.cnt[b as usize];
        sc.anchors[*c as usize] = an;
        *c += 1;
    }
    for &(b, _, _) in &sc.pairs {
        sc.cnt[b as usize] = 0;
    }
}

/// The index entries read `a`'s seeds hit, each with the offset of its first hit on a
/// higher-id read (entries are sorted by read id), and the anchor count that implies.
/// Masked or unshared hashes contribute an empty entry, so the result stays aligned
/// with `seeds`.
fn read_hits<'i>(
    a: u32,
    seeds: &[Seed],
    index: &'i SeedIndex,
    out: &mut Vec<(&'i [SeedHit], u32)>,
) -> u64 {
    out.clear();
    let mut n = 0u64;
    for s in seeds {
        let entry = index.occurrences(s.hash());
        let start = entry.partition_point(|h| h.read_id() <= a);
        n += (entry.len() - start) as u64;
        out.push((entry, start as u32));
    }
    n
}

/// Run `f` on every `(b, anchors)` cloud of read `a` against a higher-id read, in
/// ascending `b`, slicing the partner range when the read's cloud is too large to
/// gather at once ([`MAX_READ_ANCHORS`]).
fn for_each_cloud<F>(a: u32, seeds: &[Seed], index: &SeedIndex, sc: &mut Scratch, mut f: F)
where
    F: FnMut(&Scratch),
{
    let mut hits: Vec<(&[SeedHit], u32)> = Vec::new();
    let n = read_hits(a, seeds, index, &mut hits);
    let slices = n.div_ceil(MAX_READ_ANCHORS).max(1) as u32;
    if slices == 1 {
        gather(a, seeds, &hits, 0, u32::MAX, sc);
        f(sc);
        return;
    }
    // Partner ids are uncorrelated with anchor load, so equal id ranges are close to
    // equal anchor loads; a slice is still a complete set of pairs.
    let (lo, span) = (a + 1, index.n_reads().saturating_sub(a + 1));
    for k in 0..slices {
        let b_lo = lo + (span as u64 * k as u64 / slices as u64) as u32;
        let b_hi = lo + (span as u64 * (k + 1) as u64 / slices as u64) as u32;
        if b_lo < b_hi {
            gather(a, seeds, &hits, b_lo, b_hi, sc);
            f(sc);
        }
    }
}

/// The per-pair decision: the aligner, wtdbg2-style repeat guards, and the null gate.
///
/// A pure function of the pair's anchors (their order does not matter — every statistic
/// is computed over sorted copies or order-free sums), so the overlap set does not depend
/// on how the pairs were gathered or scheduled.
pub fn align_pair(
    a: u32,
    b: u32,
    anchors: &[Anchor],
    dp: &BinDp,
    params: &OverlapParams,
) -> Option<Overlap> {
    // Cheap candidacy floor: a real alignment needs at least `min_aln_bins` occupied
    // bins, hence at least that many anchors.
    if (anchors.len() as u32) < params.min_aln_bins {
        return None;
    }
    let aln = if params.hough {
        // Global Hough-extent overlapper (detects the diagonal itself). See `dp::hough_aln`.
        dp::hough_aln(anchors, dp, params.hough_min_conf as f64)
    } else {
        dp::bin_dp(anchors, 1.0, dp)
    };
    if !aln.is_overlap() {
        return None;
    }
    // Repeat guards (cross-cluster false-positive filtering). Both are raw-signal analogs
    // of wtdbg2's `kmax` seed-occurrence mask.
    //  (1) saturation: matched event length far exceeds the aligned span;
    //  (2) density: a huge anchor cloud crammed onto a short diagonal.
    // See OverlapParams::{max_similarity, max_anchor_density}.
    let aligned_bins = aln.aligned_bins.max(1);
    let aligned_events = (aligned_bins * params.bin_events) as f32;
    if aln.matched as f32 / aligned_events > params.max_similarity {
        return None;
    }
    if anchors.len() as f32 / aligned_bins as f32 > params.max_anchor_density {
        return None;
    }
    // The null-hypothesis significance (`dp::diagonal_significance`), for EVERY accepted
    // pair so the cache carries it, and the gate on it when `null_gate` is set. The
    // Hough aligner has already computed it when it had a `min_z` floor to apply.
    let parts = aln
        .null
        .or_else(|| dp::diagonal_significance_parts(anchors, aln.slope, params.bin_events as f64));
    let z = parts.map_or(0.0, |p| p.z);
    if params.null_gate && z < params.null_z as f64 {
        return None;
    }
    Some(Overlap {
        a,
        b,
        strand: Strand::Forward,
        n_anchors: anchors.len() as u32,
        matched: aln.matched,
        cnt: aln.cnt,
        aligned_bins: aln.aligned_bins,
        bin_events: params.bin_events,
        gap: aln.gap,
        slope: aln.slope,
        a_bin_span: aln.a_bin_span,
        b_bin_span: aln.b_bin_span,
        path: aln.path,
        z: z as f32,
    })
}

/// Find candidate overlaps: every read's anchors against every higher-id read, aligned
/// pair by pair. `reads[i]` holds the seeds of read `i` (the sketch order), `index` the
/// same seeds inverted; `e` is the e-mer width (`SketchParams::e`), used by the DP for the
/// matched-event length.
///
/// **Streamed per read.** A worker gathers one read's clouds, aligns them, and drops
/// them before the next read, so the live anchor set is a few reads' worth per thread
/// instead of the whole candidate set. Each `(i, j)` hit pair is enumerated exactly once
/// (the lower read id owns the pair), and a pair's alignment is a pure function of its
/// anchors, so the result does not depend on scheduling. (A blocked, sharded, two-pass
/// group-by serialises its fill phase on shard locks.)
///
/// Returns one accepted [`Overlap`] per read pair, sorted by `(a, b)`.
///
/// Without read lengths the opt-in pre-filter cannot run: see
/// [`find_overlap_candidates_with`].
pub fn find_overlap_candidates(
    reads: &[&[Seed]],
    index: &SeedIndex,
    params: &OverlapParams,
    e: u32,
) -> Vec<Overlap> {
    find_overlap_candidates_with(reads, index, params, e, None)
}

/// [`find_overlap_candidates`] with the reads' lengths in retained events (`lens[i]` for read
/// `i`, the coordinate space; `ReadSketch::n_coord_events`), which the opt-in pre-filter
/// ([`OverlapParams::prefilter_nlp`], [`prefilter`]) needs.
///
/// With `prefilter_nlp = None` the lengths are unused and the result is exactly
/// [`find_overlap_candidates`]'s. With `Some(level)`, every gathered pair is scored by its owner
/// (the lower id, the same cloud it would align) and [`align_pair`] runs only on the pairs that
/// reach `level`; the result is the unfiltered candidate set restricted to the nominated pairs,
/// since a pair's alignment depends on its anchors alone.
///
/// # Panics
/// When `prefilter_nlp` is set and `lens` is `None` or shorter than the read set.
pub fn find_overlap_candidates_with(
    reads: &[&[Seed]],
    index: &SeedIndex,
    params: &OverlapParams,
    e: u32,
    lens: Option<&[u32]>,
) -> Vec<Overlap> {
    find_overlap_candidates_frame(reads, index, params, e, lens, None)
}

/// A pair aligned in its canonical frame ([`crate::readorder::query_first`]: the longer read is the
/// query; equal lengths, the smaller tie rank), the result stored in the owner's orientation:
/// `a` (the owner, `a < b`) and `b` with `anchors` as `a` gathered them (`a_pos` on `a`).
/// When `a` is not the query, the cloud is swapped, aligned as `align_pair(b, a, ..)` and the result
/// transposed back ([`crate::readorder::transpose_overlap`]). The cloud is also put in one order,
/// sorted by (query position, target position): the gather lists a pair's anchors in an order that
/// depends on its owner, so without the sort the alignment would be a function of the anchor
/// multiset only as far as every step of it is order-free. The pair's decision and its overlap (up
/// to the transposition) are then those either read's gather produces.
pub fn align_pair_canonical(
    a: u32,
    b: u32,
    anchors: &[Anchor],
    lens: &[u32],
    rank: &[u32],
    dp: &BinDp,
    params: &OverlapParams,
) -> Option<Overlap> {
    let (ai, bi) = (a as usize, b as usize);
    let owner_is_query = crate::readorder::query_first(lens[ai], lens[bi], rank[ai], rank[bi]);
    let mut cloud: Vec<Anchor> = if owner_is_query {
        anchors.to_vec()
    } else {
        anchors
            .iter()
            .map(|x| Anchor {
                a_pos: x.b_pos,
                b_pos: x.a_pos,
            })
            .collect()
    };
    cloud.sort_unstable_by_key(|x| (x.a_pos, x.b_pos));
    if owner_is_query {
        align_pair(a, b, &cloud, dp, params)
    } else {
        align_pair(b, a, &cloud, dp, params).map(|o| crate::readorder::transpose_overlap(&o))
    }
}

/// [`find_overlap_candidates_with`] with the reads' tie ranks (`rank[i]`,
/// [`crate::readorder::rank_by_content`]), which the opt-in canonical pair frame
/// ([`OverlapParams::symmetric_pairs`]) needs with the lengths: every pair is aligned by
/// [`align_pair_canonical`] and the pre-filter scores both frames. With `symmetric_pairs` off the
/// ranks are unused and the result is [`find_overlap_candidates_with`]'s.
///
/// # Panics
/// When `symmetric_pairs` is set and `lens` or `rank` is `None` or shorter than the read set.
pub fn find_overlap_candidates_frame(
    reads: &[&[Seed]],
    index: &SeedIndex,
    params: &OverlapParams,
    e: u32,
    lens: Option<&[u32]>,
    rank: Option<&[u32]>,
) -> Vec<Overlap> {
    use std::sync::atomic::{AtomicU64, Ordering};
    let need = index.n_reads().min(reads.len() as u32) as usize;
    let nominator = params.prefilter_nlp.map(|level| {
        let lens =
            lens.expect("prefilter_nlp needs the reads' lengths (find_overlap_candidates_with)");
        assert!(
            lens.len() >= need,
            "prefilter_nlp: {} lengths for {} reads",
            lens.len(),
            reads.len()
        );
        prefilter::Nominator::fixed(lens, level).symmetric(params.symmetric_pairs)
    });
    let frame: Option<(&[u32], &[u32])> = params.symmetric_pairs.then(|| {
        let (l, r) = (
            lens.expect("symmetric_pairs needs the reads' lengths"),
            rank.expect("symmetric_pairs needs the reads' tie ranks"),
        );
        assert!(
            l.len() >= need && r.len() >= need,
            "symmetric_pairs: {} lengths and {} ranks for {} reads",
            l.len(),
            r.len(),
            reads.len()
        );
        (l, r)
    });
    let (scored, nominated) = (AtomicU64::new(0), AtomicU64::new(0));
    let dp = BinDp::from_params(params, e);
    let n_reads = index.n_reads().min(reads.len() as u32);
    log::debug!(
        "overlap: {n_reads} reads, {} index entries, streamed per read",
        index.len()
    );
    let per_read: Vec<Vec<Overlap>> = (0..n_reads)
        .into_par_iter()
        .map_init(
            || Scratch::new(index.n_reads()),
            |sc, a| {
                let seeds = reads[a as usize];
                debug_assert!(seeds.iter().all(|s| s.read_id() == a));
                let mut out = Vec::new();
                for_each_cloud(a, seeds, index, sc, |sc| {
                    // The pre-filter's counters, per slice of this read's partners rather than
                    // per pair (one shared atomic per pair cost ~2-4% of detection).
                    let nominated_here = AtomicU64::new(0);
                    if nominator.is_some() {
                        scored.fetch_add(sc.pairs.len() as u64, Ordering::Relaxed);
                    }
                    // Pairs in parallel too, so one read with many partners cannot
                    // leave a single worker running at the end of the loop.
                    let found: Vec<Overlap> = sc
                        .pairs
                        .par_iter()
                        .filter_map(|&(b, s, e)| {
                            let anchors = &sc.anchors[s as usize..e as usize];
                            // Opt-in pre-filter: the owner `a` scores the pair on this
                            // same cloud and skips `align_pair` unless it is nominated.
                            if let Some(nom) = &nominator {
                                if !nom.nominated(a, b, anchors) {
                                    return None;
                                }
                                nominated_here.fetch_add(1, Ordering::Relaxed);
                            }
                            match frame {
                                Some((l, rk)) => {
                                    align_pair_canonical(a, b, anchors, l, rk, &dp, params)
                                }
                                None => align_pair(a, b, anchors, &dp, params),
                            }
                        })
                        .collect();
                    nominated.fetch_add(nominated_here.into_inner(), Ordering::Relaxed);
                    out.extend(found);
                });
                out
            },
        )
        .collect();
    if let Some(level) = params.prefilter_nlp {
        // One line on stderr, only when the opt-in pre-filter ran (the default is silent).
        let (n, m) = (scored.into_inner(), nominated.into_inner());
        eprintln!(
            "prefilter: nlp2n_1 >= {level}: {m} of {n} gathered pairs nominated ({:.2}%), the rest not aligned",
            100.0 * m as f64 / n.max(1) as f64
        );
    }
    // Reads come back in order and each read's partners ascend, so this is already
    // sorted; the sort keeps the contract explicit.
    let mut overlaps: Vec<Overlap> = per_read.into_iter().flatten().collect();
    overlaps.sort_by_key(|o| (o.a, o.b));
    overlaps
}

/// The anchor clouds of every pair whose lower read id lies in `a_range`, keyed by pair
/// — exactly what [`find_overlap_candidates`] aligns. For inspection and tests.
pub fn anchor_clouds(
    reads: &[&[Seed]],
    index: &SeedIndex,
    a_range: std::ops::Range<u32>,
) -> PairClouds {
    let mut sc = Scratch::new(index.n_reads());
    let mut out = PairClouds::default();
    for a in a_range.start..a_range.end.min(index.n_reads()).min(reads.len() as u32) {
        for_each_cloud(a, reads[a as usize], index, &mut sc, |sc| {
            for &(b, s, e) in &sc.pairs {
                out.insert((a, b), sc.anchors[s as usize..e as usize].to_vec());
            }
        });
    }
    out
}

/// Build seed lists for many reads and find their overlaps in one call.
pub fn overlaps_from_seeds(reads: &[Vec<Seed>], params: &OverlapParams, e: u32) -> Vec<Overlap> {
    let slices: Vec<&[Seed]> = reads.iter().map(|v| v.as_slice()).collect();
    let index = SeedIndex::build(slices.iter().copied(), params);
    find_overlap_candidates(&slices, &index, params, e)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(hash: u64, read_id: u32, pos: u32, strand: Strand) -> Seed {
        Seed::new(hash, read_id, pos, strand)
    }

    /// End-to-end wiring: two reads sharing a long collinear seed run produce one
    /// accepted overlap with a non-empty bin path.
    #[test]
    fn two_reads_sharing_seeds_overlap() {
        let mut params = crate::params::Params::r10_4_1().overlap;
        // Toy-scale knobs so the synthetic run clears the gate.
        params.bin_events = 4;
        params.min_aln_bins = 4;
        params.min_matches = 2;
        params.min_similarity = 0.0;
        // Disable the repeat guards for this synthetic dense overlap.
        params.max_similarity = f32::INFINITY;
        params.max_anchor_density = f32::INFINITY;
        params.pgap = -2;
        params.pvar = -3;

        // 40 collinear shared seeds (distinct hashes) → slope-1 diagonal.
        let mk = |rid: u32, base: u32| -> Vec<Seed> {
            (0..40)
                .map(|i| Seed::new(100 + i as u64, rid, base + i, Strand::Forward))
                .collect()
        };
        let reads = vec![mk(0, 0), mk(1, 5)];
        let ovls = overlaps_from_seeds(&reads, &params, 2);
        assert_eq!(ovls.len(), 1);
        assert_eq!(ovls[0].a, 0);
        assert_eq!(ovls[0].b, 1);
        assert_eq!(ovls[0].n_anchors, 40);
        assert!(
            ovls[0].aligned_bins >= 4,
            "aligned_bins={}",
            ovls[0].aligned_bins
        );
        assert!(!ovls[0].path.is_empty());
    }

    /// Non-overlapping reads (disjoint hashes) produce no overlap.
    #[test]
    fn disjoint_reads_no_overlap() {
        let mut params = crate::params::Params::r10_4_1().overlap;
        params.bin_events = 4;
        params.min_aln_bins = 4;
        let a: Vec<Seed> = (0..40)
            .map(|i| Seed::new(100 + i, 0, i as u32, Strand::Forward))
            .collect();
        let b: Vec<Seed> = (0..40)
            .map(|i| Seed::new(5000 + i, 1, i as u32, Strand::Forward))
            .collect();
        let ovls = overlaps_from_seeds(&[a, b], &params, 2);
        assert!(ovls.is_empty());
    }

    /// A deterministic pseudo-random multi-read seed set with plenty of repeated
    /// hashes (so within-pair multiplicity `ca * cb > 1` actually occurs) and reads of
    /// unequal length.
    fn synthetic_reads(n_reads: u32) -> Vec<Vec<Seed>> {
        let mut st = 0x1234_5678_9abc_def0u64;
        let mut rnd = move || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            st
        };
        (0..n_reads)
            .map(|rid| {
                let len = 40 + (rid as usize % 7) * 15;
                (0..len)
                    // Small hash alphabet => hashes recur within and across reads.
                    .map(|p| Seed::new(rnd() % 37, rid, p as u32, Strand::Forward))
                    .collect()
            })
            .collect()
    }

    /// Per-read gathering must reproduce the straightforward serial grouping of the same
    /// index: every cross-read same-strand match in the right `(a, b)` pair, self-repeats
    /// and cross-strand hits dropped, none lost or duplicated. Compared
    /// order-insensitively, since [`align_pair`] does not depend on anchor order.
    #[test]
    fn anchor_clouds_match_serial_reference() {
        // Reads 0 and 1 share 2500 collinear seeds; read 2 shares 10 of those hashes on
        // the reverse strand (cross-strand -> dropped), and read 0 carries a hash twice
        // (self-repeat -> dropped).
        let n = 2500u32;
        let mut reads: Vec<Vec<Seed>> = vec![Vec::new(); 3];
        for i in 0..n {
            reads[0].push(seed(1000 + i as u64, 0, i, Strand::Forward));
            reads[1].push(seed(1000 + i as u64, 1, i + 3, Strand::Forward));
        }
        for i in 0..10u32 {
            reads[2].push(seed(1000 + i as u64, 2, i, Strand::Reverse));
        }
        for p in 0..2u32 {
            reads[0].push(seed(9999, 0, 2500 + p, Strand::Forward));
        }
        let mut all = reads.clone();
        all.extend(synthetic_reads(23).into_iter().map(|v| {
            v.into_iter()
                .map(|mut s| {
                    s.set_read_id(s.read_id() + 3);
                    s
                })
                .collect()
        }));

        let params = crate::params::Params::r10_4_1().overlap;
        let slices: Vec<&[Seed]> = all.iter().map(|v| v.as_slice()).collect();
        let index = SeedIndex::build(slices.iter().copied(), &params);
        let got = anchor_clouds(&slices, &index, 0..u32::MAX);

        // Serial reference over the exact same (post-masking) index entries.
        let mut want: PairClouds = FxHashMap::default();
        for (_h, hits) in index.entries() {
            for i in 0..hits.len() {
                for j in (i + 1)..hits.len() {
                    let (x, y) = (hits[i], hits[j]);
                    if x.read_id() == y.read_id() || x.strand() != y.strand() {
                        continue;
                    }
                    let (a, b, a_pos, b_pos) = if x.read_id() < y.read_id() {
                        (x.read_id(), y.read_id(), x.pos, y.pos)
                    } else {
                        (y.read_id(), x.read_id(), y.pos, x.pos)
                    };
                    want.entry((a, b))
                        .or_default()
                        .push(Anchor { a_pos, b_pos });
                }
            }
        }
        let flat = |m: &PairClouds, k: &(u32, u32)| -> Vec<(u32, u32)> {
            let mut v: Vec<(u32, u32)> = m[k].iter().map(|a| (a.a_pos, a.b_pos)).collect();
            v.sort_unstable();
            v
        };
        assert_eq!(got.len(), want.len(), "pair count");
        for k in want.keys() {
            assert_eq!(flat(&got, k), flat(&want, k), "anchors for pair {k:?}");
        }
        assert_eq!(got[&(0, 1)].len(), n as usize);
        assert!(!got.contains_key(&(0, 2)) && !got.contains_key(&(1, 2)));
        // Seeds in position order give clouds sorted by `a_pos`.
        assert!(got[&(0, 1)].is_sorted_by_key(|x| x.a_pos));
    }

    /// Gathering a read in slices of the partner range yields the same pairs and anchors
    /// as gathering it whole: the slicing is a memory bound, not a different enumeration.
    #[test]
    fn sliced_gathering_matches_the_whole_read() {
        let reads = synthetic_reads(23);
        let params = crate::params::Params::r10_4_1().overlap;
        let slices: Vec<&[Seed]> = reads.iter().map(|v| v.as_slice()).collect();
        let index = SeedIndex::build(slices.iter().copied(), &params);
        let mut sc = Scratch::new(index.n_reads());
        let collect = |sc: &Scratch| -> Vec<(u32, Vec<(u32, u32)>)> {
            sc.pairs
                .iter()
                .map(|&(b, s, e)| {
                    let v = sc.anchors[s as usize..e as usize]
                        .iter()
                        .map(|x| (x.a_pos, x.b_pos))
                        .collect();
                    (b, v)
                })
                .collect()
        };
        for a in 0..index.n_reads() {
            let mut hits = Vec::new();
            read_hits(a, slices[a as usize], &index, &mut hits);
            gather(a, slices[a as usize], &hits, 0, u32::MAX, &mut sc);
            let whole = collect(&sc);
            let mut parts = Vec::new();
            let mid = (a + 1 + index.n_reads()) / 2;
            for (lo, hi) in [(a + 1, mid), (mid, index.n_reads())] {
                gather(a, slices[a as usize], &hits, lo, hi, &mut sc);
                parts.extend(collect(&sc));
            }
            assert_eq!(parts, whole, "read {a}");
            assert!(sc.cnt.iter().all(|&c| c == 0), "scratch counts reset");
        }
    }

    /// The detection settings of `haricot assemble`, for the pre-filter tests.
    fn production() -> OverlapParams {
        let mut p = crate::params::Params::r10_4_1().overlap;
        p.bin_events = 192;
        p.hough = true;
        p.hough_min_conf = 0.0;
        p.hough_min_z = 3.0;
        p.null_gate = true;
        p.null_z = 4.0;
        p.min_aln_bins = 4;
        p
    }

    /// Fourteen reads of 2,500-5,500 positions tiling a 12,000-position genome of distinct
    /// hashes (one seed every 2 positions, so overlapping reads share long slope-1
    /// diagonals), each also carrying 40 seeds drawn from 40 shared "repeat" hashes, so
    /// non-overlapping pairs share a few chance anchors. Returns the seeds (in position order)
    /// and each read's length.
    fn tiled_reads() -> (Vec<Vec<Seed>>, Vec<u32>) {
        let mut st = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = move || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            st
        };
        // distinct hashes below 2^32 (a seed's hash is 32-bit): strictly increasing in the position
        let genome: Vec<u64> = (0..12_000u64)
            .map(|g| 1_000_000 + g * 7_919 + rnd() % 7_919)
            .collect();
        let (mut reads, mut lens) = (Vec::new(), Vec::new());
        for rid in 0..14u32 {
            let start = (rid as usize * 797) % 9_000;
            let end = (start + 2_500 + (rid as usize * 331) % 3_000).min(genome.len());
            let mut v: Vec<Seed> = (start..end)
                .step_by(2)
                .map(|g| seed(genome[g], rid, (g - start) as u32, Strand::Forward))
                .collect();
            for _ in 0..40 {
                let p = (rnd() % (end - start) as u64) as u32;
                v.push(seed(rnd() % 40, rid, p, Strand::Forward));
            }
            v.sort_by_key(|s| s.pos);
            reads.push(v);
            lens.push((end - start) as u32);
        }
        (reads, lens)
    }

    /// Every field of an overlap, for exact comparisons.
    type OvlKey = (
        u32,
        u32,
        u32,
        u32,
        u32,
        u32,
        u32,
        u64,
        (u32, u32),
        (u32, u32),
        Vec<(u32, u32)>,
        u32,
    );

    fn keys(v: &[Overlap]) -> Vec<OvlKey> {
        v.iter()
            .map(|o| {
                (
                    o.a,
                    o.b,
                    o.n_anchors,
                    o.matched,
                    o.cnt,
                    o.aligned_bins,
                    o.gap,
                    o.slope.to_bits(),
                    o.a_bin_span,
                    o.b_bin_span,
                    o.path.clone(),
                    o.z.to_bits(),
                )
            })
            .collect()
    }

    /// The serial reference clouds (every same-strand hit pair of the index, owner = the
    /// lower id), independent of [`gather`].
    fn reference_clouds(index: &SeedIndex) -> PairClouds {
        let mut want: PairClouds = FxHashMap::default();
        for (_h, hits) in index.entries() {
            for i in 0..hits.len() {
                for j in (i + 1)..hits.len() {
                    let (x, y) = (hits[i], hits[j]);
                    if x.read_id() == y.read_id() || x.strand() != y.strand() {
                        continue;
                    }
                    let (a, b, a_pos, b_pos) = if x.read_id() < y.read_id() {
                        (x.read_id(), y.read_id(), x.pos, y.pos)
                    } else {
                        (y.read_id(), x.read_id(), y.pos, x.pos)
                    };
                    want.entry((a, b))
                        .or_default()
                        .push(Anchor { a_pos, b_pos });
                }
            }
        }
        want
    }

    /// With the pre-filter off, the detector is unchanged: `find_overlap_candidates_with`
    /// (lengths given, `prefilter_nlp = None`) returns exactly `find_overlap_candidates`'s set,
    /// and so does a pre-filter that nominates every pair (level -inf).
    #[test]
    fn prefilter_off_is_the_plain_candidate_set() {
        let (reads, lens) = tiled_reads();
        let p = production();
        let slices: Vec<&[Seed]> = reads.iter().map(|v| v.as_slice()).collect();
        let index = SeedIndex::build(slices.iter().copied(), &p);
        let base = find_overlap_candidates(&slices, &index, &p, 8);
        assert!(base.len() >= 10, "only {} overlaps", base.len());
        let off = find_overlap_candidates_with(&slices, &index, &p, 8, Some(&lens));
        assert_eq!(keys(&off), keys(&base));
        let mut all = p;
        all.prefilter_nlp = Some(f64::NEG_INFINITY);
        let every = find_overlap_candidates_with(&slices, &index, &all, 8, Some(&lens));
        assert_eq!(keys(&every), keys(&base));
    }

    /// With the pre-filter on, the candidate set equals the per-pair brute force over the
    /// serial reference clouds: a pair is kept when its owner-frame `nlp2n_1` reaches the level
    /// and `align_pair` accepts it.
    #[test]
    fn prefilter_keeps_exactly_the_nominated_pairs() {
        let (reads, lens) = tiled_reads();
        let mut p = production();
        p.prefilter_nlp = Some(2.0);
        let slices: Vec<&[Seed]> = reads.iter().map(|v| v.as_slice()).collect();
        let index = SeedIndex::build(slices.iter().copied(), &p);
        let dp = BinDp::from_params(&p, 8);
        let clouds = reference_clouds(&index);
        let mut want: Vec<Overlap> = Vec::new();
        let (mut nom, mut skipped) = (0, 0);
        let mut cnt = Vec::new();
        for (&(a, b), c) in &clouds {
            let s = prefilter::pair_nlp(c, lens[a as usize], lens[b as usize], &mut cnt);
            if s >= 2.0 {
                nom += 1;
                want.extend(align_pair(a, b, c, &dp, &p));
            } else {
                skipped += 1;
            }
        }
        want.sort_by_key(|o| (o.a, o.b));
        assert!(
            nom >= 10 && skipped >= 10,
            "nominated {nom}, skipped {skipped}"
        );
        let got = find_overlap_candidates_with(&slices, &index, &p, 8, Some(&lens));
        assert_eq!(keys(&got), keys(&want));
    }

    /// The owner (the lower id) scores the pair in its own frame: with a level between a
    /// pair's owner-frame and swapped-frame scores, the pair's presence follows the owner's
    /// score. Each pair appears once, `a < b`.
    #[test]
    fn prefilter_scores_in_the_owner_frame() {
        let (reads, lens) = tiled_reads();
        let p = production();
        let slices: Vec<&[Seed]> = reads.iter().map(|v| v.as_slice()).collect();
        let index = SeedIndex::build(slices.iter().copied(), &p);
        let base = find_overlap_candidates(&slices, &index, &p, 8);
        let clouds = reference_clouds(&index);
        let mut cnt = Vec::new();
        let pick = base.iter().find_map(|o| {
            let c = &clouds[&(o.a, o.b)];
            let sw: Vec<Anchor> = c
                .iter()
                .map(|x| Anchor {
                    a_pos: x.b_pos,
                    b_pos: x.a_pos,
                })
                .collect();
            let so = prefilter::pair_nlp(c, lens[o.a as usize], lens[o.b as usize], &mut cnt);
            let ss = prefilter::pair_nlp(&sw, lens[o.b as usize], lens[o.a as usize], &mut cnt);
            (so != ss).then_some((o.a, o.b, so, ss))
        });
        let (a, b, so, ss) = pick.expect("a pair whose two frames score differently");
        let mut q = p;
        q.prefilter_nlp = Some(0.5 * (so + ss));
        let got = find_overlap_candidates_with(&slices, &index, &q, 8, Some(&lens));
        assert_eq!(got.iter().any(|o| (o.a, o.b) == (a, b)), so > ss);
        assert!(got.iter().all(|o| o.a < o.b));
        let mut ps: Vec<(u32, u32)> = got.iter().map(|o| (o.a, o.b)).collect();
        ps.dedup();
        assert_eq!(ps.len(), got.len());
    }

    /// [`align_pair`] is a pure function of the anchor MULTISET: shuffling a cloud (a
    /// diagonal in noise) changes no field of the result.
    #[test]
    fn align_pair_does_not_depend_on_anchor_order() {
        let mut params = crate::params::Params::r10_4_1().overlap;
        params.hough = true;
        params.hough_min_conf = 0.0;
        params.hough_min_z = 3.0;
        params.null_gate = true;
        params.null_z = 4.0;
        params.min_aln_bins = 4;
        params.bin_events = 16;
        let dp = BinDp::from_params(&params, 8);
        let mut st = 0x9e37_79b9_7f4a_7c15u64;
        let mut rnd = move || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            st
        };
        let mut cloud: Vec<Anchor> = (0..400)
            .map(|i| Anchor {
                a_pos: 50 + (i * 97) / 100,
                b_pos: i,
            })
            .collect();
        cloud.extend((0..600).map(|_| Anchor {
            a_pos: (rnd() % 1200) as u32,
            b_pos: (rnd() % 1200) as u32,
        }));
        let key = |o: &Overlap| {
            (
                o.matched,
                o.cnt,
                o.aligned_bins,
                o.gap,
                o.slope.to_bits(),
                o.a_bin_span,
                o.b_bin_span,
                o.path.clone(),
                o.z.to_bits(),
            )
        };
        let base = align_pair(0, 1, &cloud, &dp, &params).expect("the diagonal is found");
        for _ in 0..5 {
            for i in (1..cloud.len()).rev() {
                cloud.swap(i, (rnd() % (i as u64 + 1)) as usize);
            }
            let got = align_pair(0, 1, &cloud, &dp, &params).expect("still found");
            assert_eq!(key(&got), key(&base));
        }
    }

    /// The canonical pair frame: the same pair gathered by either read (numbering 1: x = 0 owns it;
    /// numbering 2: y = 0 owns it, the cloud swapped) gives the same decision and the same alignment
    /// once both are put in the canonical frame (the longer read y as the query).
    #[test]
    fn canonical_alignment_does_not_depend_on_the_owner() {
        let mut params = crate::params::Params::r10_4_1().overlap;
        params.hough = true;
        params.hough_min_conf = 0.0;
        params.hough_min_z = 3.0;
        params.null_gate = true;
        params.null_z = 4.0;
        params.min_aln_bins = 4;
        params.bin_events = 16;
        let dp = BinDp::from_params(&params, 8);
        let mut st = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = move || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            st
        };
        // x (1,300 positions) overlaps y (2,000) at slope ~1.08 in x = s * y + c
        let mut cloud: Vec<Anchor> = (0..500)
            .map(|i| Anchor {
                a_pos: 100 + (i * 108) / 100,
                b_pos: 300 + i,
            })
            .collect();
        cloud.extend((0..700).map(|_| Anchor {
            a_pos: (rnd() % 1300) as u32,
            b_pos: (rnd() % 2000) as u32,
        }));
        // y's gather lists the anchors in another order (here reversed)
        let swapped: Vec<Anchor> = cloud
            .iter()
            .rev()
            .map(|x| Anchor {
                a_pos: x.b_pos,
                b_pos: x.a_pos,
            })
            .collect();
        // numbering 1: 0 = x, 1 = y; numbering 2: 0 = y, 1 = x (ranks follow the names x < y)
        let (lens1, rank1) = ([1300u32, 2000], [0u32, 1]);
        let (lens2, rank2) = ([2000u32, 1300], [1u32, 0]);
        let o1 = align_pair_canonical(0, 1, &cloud, &lens1, &rank1, &dp, &params);
        let o2 = align_pair_canonical(0, 1, &swapped, &lens2, &rank2, &dp, &params);
        let (o1, o2) = (o1.expect("found from x"), o2.expect("found from y"));
        let c1 = crate::readorder::in_canonical_frame(&o1, &lens1, &rank1);
        let c2 = crate::readorder::in_canonical_frame(&o2, &lens2, &rank2);
        assert_eq!(
            (c1.a, c1.b, c2.a, c2.b),
            (1, 0, 0, 1),
            "y is the query in both"
        );
        let key = |o: &Overlap| {
            (
                o.n_anchors,
                o.matched,
                o.cnt,
                o.aligned_bins,
                o.gap,
                o.a_bin_span,
                o.b_bin_span,
                o.path.clone(),
                o.z.to_bits(),
            )
        };
        assert_eq!(key(&c1), key(&c2));
        assert!((c1.slope - c2.slope).abs() <= 1e-12 * c2.slope.abs());
    }
}
