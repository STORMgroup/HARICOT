//! Seed index — the analog of wtdbg2's KBM (k-bin mapping).
//!
//! wtdbg2 indexes k-mers into 256 bp bins and finds read-to-read matches by
//! looking up shared k-mers (`extern/wtdbg2/kbm.{c,h}`). Here we index e-mer
//! seeds by hash and find reads that share many seeds. This is the core of the
//! costly all-vs-all comparison that the de-Bruijn-style approach aims to
//! avoid doing exhaustively.
//!
//! Scope: this module is only the inverted seed index (with repeat masking) and
//! shared-seed lookup. The banded alignment over the resulting anchors lives in
//! [`crate::overlap::dp`].

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::emer::Seed;
use crate::emer::sketch::{STRAND_BIT, pack_read_strand, strand_of};
use crate::params::{OverlapParams, Strand};

/// One occurrence of a seed on a read, in 8 bytes (it was 12): the position, and the read id
/// with the strand in its top bit (as [`Seed`]).
#[derive(Debug, Clone, Copy)]
pub struct SeedHit {
    pub pos: u32,
    id: u32,
}

impl SeedHit {
    pub fn new(read_id: u32, pos: u32, strand: Strand) -> SeedHit {
        SeedHit {
            pos,
            id: pack_read_strand(read_id, strand),
        }
    }

    #[inline]
    pub fn read_id(&self) -> u32 {
        self.id & !STRAND_BIT
    }

    #[inline]
    pub fn strand(&self) -> Strand {
        strand_of(self.id)
    }
}

/// Inverted index: e-mer hash → occurrences across all reads.
///
/// **Entry invariant: each entry's hits are sorted by `read_id` (ties by `pos`).**
/// Every read's occurrences of a hash therefore form one contiguous run, and the hits
/// on reads with a higher id than a given read are a suffix found by one binary search
/// — which is how [`crate::overlap::find_overlap_candidates`] gives each pair to its
/// lower read.
#[derive(Debug, Default)]
pub struct SeedIndex {
    /// The entries, split by [`shard_of`] the hash into [`SHARDS`] maps built in parallel.
    shards: Vec<FxHashMap<u64, Vec<SeedHit>>>,
    /// Hashes dropped as too frequent (repeat-masked).
    masked: u64,
    /// One past the highest `read_id` seen at build time (before masking).
    n_reads: u32,
}

/// Hash shards of the index: each is filled by one thread, so the build is parallel without
/// a shared map.
const SHARDS: usize = 16;

/// The shard of a hash (its halves folded, so it does not rest on the hash's top bits).
#[inline]
fn shard_of(hash: u64) -> usize {
    ((hash ^ (hash >> 32) ^ (hash >> 16)) as usize) & (SHARDS - 1)
}

/// How often each hash occurs over many seed lists, counted in parallel [`SHARDS`] (counts do not
/// depend on the order hashes are counted in, so this is the serial count).
#[derive(Debug, Default)]
pub struct SeedCounts {
    shards: Vec<FxHashMap<u64, u32>>,
}

impl SeedCounts {
    /// Count the hashes `hash(list, i)`, `i < len(list)`, of every list. Each seed's shard is
    /// tagged once, in parallel; each shard then reads only its own seeds.
    pub fn count<L, N, H>(lists: &[L], len: N, hash: H) -> Self
    where
        L: Sync,
        N: Fn(&L) -> usize + Sync,
        H: Fn(&L, usize) -> u64 + Sync,
    {
        let tags: Vec<Vec<u8>> = lists
            .par_iter()
            .map(|l| (0..len(l)).map(|i| shard_of(hash(l, i)) as u8).collect())
            .collect();
        let shards = (0..SHARDS)
            .into_par_iter()
            .map(|k| {
                let mut map: FxHashMap<u64, u32> = FxHashMap::default();
                for (l, tags) in lists.iter().zip(&tags) {
                    for (i, &t) in tags.iter().enumerate() {
                        if t as usize == k {
                            *map.entry(hash(l, i)).or_insert(0) += 1;
                        }
                    }
                }
                map
            })
            .collect();
        SeedCounts { shards }
    }

    /// Occurrences of `hash` (0 if it was never counted).
    pub fn get(&self, hash: u64) -> u32 {
        self.shards
            .get(shard_of(hash))
            .and_then(|m| m.get(&hash))
            .copied()
            .unwrap_or(0)
    }

    /// Number of distinct hashes counted.
    pub fn len(&self) -> usize {
        self.shards.iter().map(|m| m.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.shards.iter().all(|m| m.is_empty())
    }

    /// Every `(hash, count)`, sorted by hash (a deterministic order, for writing the table out).
    pub fn sorted_entries(&self) -> Vec<(u64, u32)> {
        let mut v: Vec<(u64, u32)> = self
            .shards
            .iter()
            .flat_map(|m| m.iter().map(|(&h, &c)| (h, c)))
            .collect();
        v.sort_unstable();
        v
    }

    /// The counts of `entries` (`(hash, count)`, each hash once): what [`SeedCounts::count`] gave
    /// for the lists those entries were taken from ([`SeedCounts::sorted_entries`]).
    pub fn from_entries(entries: impl IntoIterator<Item = (u64, u32)>) -> Self {
        let mut shards: Vec<FxHashMap<u64, u32>> = vec![FxHashMap::default(); SHARDS];
        for (h, c) in entries {
            shards[shard_of(h)].insert(h, c);
        }
        SeedCounts { shards }
    }
}

/// The repeat-mask cutoff for this input: `max_seed_occ`, or the relative rule
/// ([`relative_seed_cap`]) when [`OverlapParams::seed_occ_rel`] is set.
fn effective_seed_cap(shards: &[FxHashMap<u64, Vec<SeedHit>>], params: &OverlapParams) -> u32 {
    if params.seed_occ_rel > 0.0 {
        return relative_seed_cap(shards, params);
    }
    params.max_seed_occ
}

/// [`OverlapParams::seed_occ_rel`]: `max(max_seed_occ, seed_occ_rel * mean occurrence)`, the
/// mean being seed instances per distinct hash of this input (before masking). A DEVIATION
/// from wtdbg2, whose `ktop` form of the same idea (`-K 1000.05`, `kbm.h:1054-1063`: the
/// occurrence above which the most frequent seeds hold 5% of all instances) was measured
/// and does not transfer to the e-mer alphabet: its top occurrences are low-complexity
/// e-mers in the hundreds of thousands (a cap of 129,495 on A. thaliana chr4 30x, out of
/// memory).
fn relative_seed_cap(shards: &[FxHashMap<u64, Vec<SeedHit>>], params: &OverlapParams) -> u32 {
    let floor = params.max_seed_occ;
    let total: u64 = shards
        .iter()
        .flat_map(|m| m.values())
        .map(|h| h.len() as u64)
        .sum();
    let n_hashes: usize = shards.iter().map(|m| m.len()).sum();
    let mean = total as f64 / n_hashes.max(1) as f64;
    let rel = (params.seed_occ_rel as f64 * mean).round() as u32;
    let cap = floor.max(rel);
    eprintln!(
        "seed index: {} hashes, {} instances, mean occurrence {mean:.1}; {:.2} x mean = {rel} -> repeat-mask cap {cap} (floor {floor})",
        n_hashes, total, params.seed_occ_rel
    );
    cap
}

impl SeedIndex {
    /// Build an index from per-read seed lists, dropping seeds that occur more
    /// than `params.max_seed_occ` times (repeat masking).
    ///
    /// Entries are left sorted by `read_id` (see the struct invariant). Callers
    /// normally feed reads in id order and each read's seeds in position order, so
    /// the check below is a linear scan that never sorts; the sort is the fallback
    /// for a caller that does not.
    pub fn build<'a, I>(reads: I, params: &OverlapParams) -> Self
    where
        I: IntoIterator<Item = &'a [Seed]>,
    {
        let reads: Vec<&[Seed]> = reads.into_iter().collect();
        let n_reads = reads
            .par_iter()
            .flat_map_iter(|seeds| seeds.iter().map(|s| s.read_id() + 1))
            .max()
            .unwrap_or(0);
        // Each shard scans every read in order and keeps its own hashes, so every entry's hits
        // arrive in the order a single serial pass would push them. The shard of each seed is
        // tagged once, in parallel, so the scans read one byte per seed.
        let tags: Vec<Vec<u8>> = reads
            .par_iter()
            .map(|seeds| seeds.iter().map(|s| shard_of(s.hash()) as u8).collect())
            .collect();
        let mut shards: Vec<FxHashMap<u64, Vec<SeedHit>>> = (0..SHARDS)
            .into_par_iter()
            .map(|k| {
                let mut map: FxHashMap<u64, Vec<SeedHit>> = FxHashMap::default();
                for (seeds, tags) in reads.iter().zip(&tags) {
                    for (s, &t) in seeds.iter().zip(tags) {
                        if t as usize == k {
                            map.entry(s.hash()).or_default().push(SeedHit::new(
                                s.read_id(),
                                s.pos,
                                s.strand(),
                            ));
                        }
                    }
                }
                map
            })
            .collect();
        drop(tags);

        let cap = effective_seed_cap(&shards, params);
        let masked: u64 = shards
            .par_iter_mut()
            .map(|map| {
                let before = map.len();
                map.retain(|_, hits| hits.len() as u32 <= cap);
                (before - map.len()) as u64
            })
            .sum();
        // Entries are independent, and even the no-op case is a full linear scan of
        // every hit (823 M on the D. melanogaster both-strand run), so this is worth a
        // `par_iter_mut` rather than a serial walk. Mirrors `FullSeeds::build`.
        shards.par_iter_mut().for_each(|map| {
            for hits in map.values_mut() {
                if !hits.is_sorted_by_key(|h| (h.read_id(), h.pos)) {
                    // Stable, so equal `(read_id, pos)` keys keep insertion order and the
                    // anchor multiset is unchanged.
                    hits.sort_by_key(|h| (h.read_id(), h.pos));
                }
            }
        });

        SeedIndex {
            shards,
            masked,
            n_reads,
        }
    }

    /// One past the highest read id the index was built over.
    pub fn n_reads(&self) -> u32 {
        self.n_reads
    }

    /// Number of distinct (non-masked) seed hashes indexed.
    pub fn len(&self) -> usize {
        self.shards.iter().map(|m| m.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.shards.iter().all(|m| m.is_empty())
    }

    /// Number of hashes dropped by repeat masking.
    pub fn masked(&self) -> u64 {
        self.masked
    }

    /// Occurrences of a given seed hash.
    pub fn occurrences(&self, hash: u64) -> &[SeedHit] {
        self.shards
            .get(shard_of(hash))
            .and_then(|m| m.get(&hash))
            .map_or(&[], |v| v.as_slice())
    }

    /// Iterate over all indexed `(hash, hits)` entries.
    pub fn entries(&self) -> impl Iterator<Item = (&u64, &Vec<SeedHit>)> {
        self.shards.iter().flat_map(|m| m.iter())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Params;

    fn seed(hash: u64, read_id: u32, pos: u32) -> Seed {
        Seed::new(hash, read_id, pos, Strand::Forward)
    }

    /// A count table written out and read back answers every hash as the counted one does.
    #[test]
    fn seed_counts_round_trip_through_their_entries() {
        let lists: Vec<Vec<u64>> = (0..40u64)
            .map(|r| {
                (0..200u64)
                    .map(|i| (r * 7919 + i * i * 31) % 997 + (i % 3) * 1_000_003)
                    .collect()
            })
            .collect();
        let c = SeedCounts::count(&lists, |l| l.len(), |l, i| l[i]);
        let e = c.sorted_entries();
        assert!(e.windows(2).all(|w| w[0].0 < w[1].0));
        let back = SeedCounts::from_entries(e.iter().copied());
        assert_eq!(back.len(), c.len());
        for h in 0..3_100_000u64 {
            assert_eq!(back.get(h), c.get(h), "hash {h}");
        }
        assert_eq!(back.sorted_entries(), e);
    }

    #[test]
    fn index_groups_by_hash() {
        let r0 = vec![seed(10, 0, 0), seed(20, 0, 5)];
        let r1 = vec![seed(10, 1, 3)];
        let params = Params::r10_4_1().overlap;
        let idx = SeedIndex::build([r0.as_slice(), r1.as_slice()], &params);
        assert_eq!(idx.len(), 2);
        assert_eq!(idx.occurrences(10).len(), 2);
        assert_eq!(idx.occurrences(20).len(), 1);
        assert_eq!(idx.occurrences(999).len(), 0);
    }

    #[test]
    fn relative_cap_scales_with_the_mean_occurrence_and_never_lowers_the_floor() {
        // Hashes 1..=4 occur 10 / 5 / 2 / 1 times: 18 instances over 4 hashes, mean 4.5.
        let mut params = Params::r10_4_1().overlap;
        params.max_seed_occ = 3;
        let mut r = Vec::new();
        for (h, n) in [(1u64, 10u32), (2, 5), (3, 2), (4, 1)] {
            r.extend((0..n).map(|i| seed(h, i, 0)));
        }
        // 2 x 4.5 = 9: hash 1 (10) masked, hash 2 (5) kept although the floor 3 would mask it.
        params.seed_occ_rel = 2.0;
        let idx = SeedIndex::build([r.as_slice()], &params);
        assert_eq!(idx.masked(), 1);
        assert_eq!(idx.occurrences(2).len(), 5);
        // 0.5 x 4.5 = 2 is below the floor: the floor (3) applies, masking hashes 1 and 2.
        params.seed_occ_rel = 0.5;
        let idx = SeedIndex::build([r.as_slice()], &params);
        assert_eq!(idx.masked(), 2);
    }

    #[test]
    fn sharded_counts_are_the_serial_counts() {
        let lists: Vec<Vec<u64>> = (0..20u64)
            .map(|i| {
                (0..500u64)
                    .map(|j| (i * 7919 + j * 104_729) % 977 * 0x9E37_79B9)
                    .collect()
            })
            .collect();
        let c = SeedCounts::count(&lists, |l| l.len(), |l, i| l[i]);
        let mut want: FxHashMap<u64, u32> = FxHashMap::default();
        for l in &lists {
            for &h in l {
                *want.entry(h).or_insert(0) += 1;
            }
        }
        assert_eq!(c.len(), want.len());
        assert!(want.iter().all(|(&h, &n)| c.get(h) == n));
        assert_eq!(c.get(1), 0);
    }

    #[test]
    fn repeat_masking_drops_frequent_seeds() {
        let mut params = Params::r10_4_1().overlap;
        params.max_seed_occ = 2;
        // hash 7 occurs 3 times -> masked.
        let r: Vec<Seed> = (0..3).map(|i| seed(7, i, 0)).collect();
        let idx = SeedIndex::build([r.as_slice()], &params);
        assert_eq!(idx.occurrences(7).len(), 0);
        assert_eq!(idx.masked(), 1);
    }
}
