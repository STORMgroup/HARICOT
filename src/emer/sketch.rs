//! E-mer construction, hashing, and (optional) minimizer sketching.
//!
//! An **e-mer** is `e` consecutive quantized events packed into one integer
//! (`q` bits each), analogous to a nucleotide k-mer. Ported from
//! `ri_sketch_reg` / `ri_sketch_min` and `hash64`
//! (`extern/RawHash2/src/rsketch.c`).
//!
//! Intentional deviations from RawHash2 (this is a readable prototype):
//!   - seeds are a plain struct rather than a packed 128-bit `(x, y)` pair;
//!   - the minimizer selection uses a straightforward per-window argmin with
//!     dedup by position, not RawHash2's duplicate-handling rotation, so the two
//!     can differ where window hashes tie. `tests/hash_compat.rs` pins the dense
//!     path (quantize → sig-diff → pack → `hash64`) against the C reference; it
//!     deliberately does not compare minimizer selection.

use crate::emer::quantize::dynamic_quantize;
use crate::params::{Params, Strand};

/// One sampled seed (an e-mer occurrence on a read), packed into 12 bytes: every
/// read's seeds are held for the whole overlap stage, so this is most of detection's memory
/// (472 M seeds on A. thaliana chr1 30x).
///
/// * the hash in 32 bits: every seed hash is `hash64(.., HASH_MASK)` ([`emer_hash`], the
///   strobemers' `select.rs` combination), and [`hash64`] ends with `& mask`, so it is below
///   2^32 ([`HASH_MASK`]); [`Seed::new`] refuses a wider one;
/// * the coordinate `pos` (`u32`);
/// * the read id in 31 bits and the strand in the top bit ([`Seed::new`] refuses a read id of
///   2^31 or more: a pod5 of two billion reads).
///
/// The raw-stream index (`event_pos`, diagnostics and the RawHash2 golden only, never a
/// coordinate) is not stored: it is `retained[pos].event_pos` of the read's retained
/// stream ([`retained_quantized`]), and the e-mers ([`Emer`]) still carry it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seed {
    hash: u32,
    /// **The coordinate**: index of the e-mer's first event in the *retained*
    /// (post-merge) stream. Every downstream coordinate — DP bins, anchors, K-bin
    /// windows, edge offsets — is measured in this space. See
    /// [`retained_quantized`] for why it is not the raw event index.
    pub pos: u32,
    /// The read id (bits 0-30) and the strand (bit 31: reverse).
    id: u32,
}

/// The top bit of a packed read id: the strand ([`Seed`], [`crate::overlap::SeedHit`]).
pub const STRAND_BIT: u32 = 1 << 31;

/// A read id and a strand in one `u32` (the id below 2^31).
pub fn pack_read_strand(read_id: u32, strand: Strand) -> u32 {
    assert!(
        read_id < STRAND_BIT,
        "read id {read_id} does not fit the packed seed (below 2^31)"
    );
    read_id | (strand.bit() << 31)
}

/// The strand of a packed read id.
pub fn strand_of(id: u32) -> Strand {
    if id & STRAND_BIT != 0 {
        Strand::Reverse
    } else {
        Strand::Forward
    }
}

impl Seed {
    /// A seed; the hash must lie below 2^32 (every sketched one does) and the read id below 2^31.
    pub fn new(hash: u64, read_id: u32, pos: u32, strand: Strand) -> Seed {
        assert!(hash <= HASH_MASK, "seed hash {hash:#x} wider than 32 bits");
        Seed {
            hash: hash as u32,
            pos,
            id: pack_read_strand(read_id, strand),
        }
    }

    /// Hash of the packed e-mer value.
    #[inline]
    pub fn hash(&self) -> u64 {
        self.hash as u64
    }

    /// Read identifier.
    #[inline]
    pub fn read_id(&self) -> u32 {
        self.id & !STRAND_BIT
    }

    /// Strand of the read this seed was sampled from.
    #[inline]
    pub fn strand(&self) -> Strand {
        strand_of(self.id)
    }

    /// The same seed on read `read_id` (renumbering).
    pub fn set_read_id(&mut self, read_id: u32) {
        self.id = pack_read_strand(read_id, self.strand());
    }
}

/// RawHash2's 64-bit invertible mixing hash (`hash64`, rsketch.c:7).
pub fn hash64(key: u64, mask: u64) -> u64 {
    let mut key = key;
    key = (!key).wrapping_add(key << 21) & mask; // key = (key << 21) - key - 1
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8) & mask; // key * 265
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4) & mask; // key * 21
    key ^= key >> 28;
    key = key.wrapping_add(key << 31) & mask;
    key
}

/// A dense (per-position) e-mer occurrence before optional minimizing.
#[derive(Clone, Copy, Debug)]
pub struct Emer {
    pub hash: u64,
    /// Coordinate: index of the e-mer's first event in the retained stream.
    pub pos: u32,
    /// The same event's index in the raw segmented stream (diagnostics only).
    pub event_pos: u32,
}

/// A quantized event that survived the sig-diff merge. This is the intermediate
/// stream that both [`dense_emers`] and the content-based selectors in
/// [`crate::emer::select`] (syncmers, strobemers) are built from — and it is the
/// stream that defines the coordinate space (see [`retained_quantized`]).
#[derive(Clone, Copy, Debug)]
pub struct QEvent {
    /// **The coordinate**: this event's index in the retained stream.
    pub pos: u32,
    /// This event's index in the raw segmented stream, before the merge.
    /// Diagnostics and the RawHash2 golden only — **never** a coordinate.
    pub event_pos: u32,
    /// Quantization bucket in `0..n_buckets`.
    pub bucket: u32,
}

/// 32-bit hash mask used throughout sketching (matches RawHash2's `hash64`).
pub const HASH_MASK: u64 = (1u64 << 32) - 1;

/// Quantize the event stream and apply the sig-diff merge, returning the retained
/// quantized events in order. Factored out of [`dense_emers`] so the
/// content-based selectors can share exactly the same front end.
///
/// **This function defines the coordinate space.** `diff >= 0` merges a run of
/// consecutive events that quantize alike into a single retained event, so each
/// `QEvent.pos` is its index in the *retained* stream, not in the raw segmented
/// stream. That distinction is the whole point: segmentation is variable, so the
/// same stretch of DNA yields a different number of raw events on every read (a
/// homopolymer read as 3 events on one read and 7 on another). The merge removes
/// exactly that variability — stamping the pre-merge index would put it straight
/// back into the coordinate, leaving overlapping reads with coordinates that drift
/// apart by the segmentation noise the merge just discarded. The raw index is kept
/// as `event_pos` for diagnostics.
///
/// (RawHash2's `ri_sketch_reg` stamps the raw index, `rsketch.c:214-218`. That is
/// a holdover from when its `diff` was a signal-level sigdiff merge rather than
/// the bucket-level filter it is now, and RawHash2 is a mapper whose chainer
/// tolerates drift. A bin-DP feeding an FBG does not.)
pub fn retained_quantized(events: &[f32], params: &Params) -> Vec<QEvent> {
    let diff = params.sketch.diff;
    let mask_quant = (1u64 << params.quant.q) - 1;
    let mut retained: Vec<QEvent> = Vec::with_capacity(events.len());
    let mut prev_quant: Option<u32> = None;
    for (i, &v) in events.iter().enumerate() {
        let quant = dynamic_quantize(v, &params.quant) & (mask_quant as u32);
        if diff >= 0
            && let Some(pq) = prev_quant
            && (quant as i32 - pq as i32).abs() <= diff
        {
            continue; // collapse near-identical consecutive events
        }
        prev_quant = Some(quant);
        let idx = retained.len() as u32;
        retained.push(QEvent {
            pos: idx,
            event_pos: i as u32,
            bucket: quant,
        });
    }
    retained
}

/// Pack `buckets.len()` quantized events (`q` bits each, oldest first) into one
/// integer and hash it with [`hash64`]. `buckets.len() * q` must be `<= 64`.
///
/// Used for both e-mers and the shorter *s*-mers of the content-based selectors,
/// so a selected syncmer's hash is identical to the dense e-mer's hash at the
/// same position (anchors still join on hash).
pub fn emer_hash(buckets: &[u32], q: u32) -> u64 {
    let mut packed = 0u64;
    for &b in buckets {
        packed = (packed << q) | b as u64;
    }
    hash64(packed, HASH_MASK)
}

/// Quantize the event stream, apply the sig-diff filter, and pack dense e-mers.
///
/// Returns the e-mers in event order. Each e-mer's [`Emer::pos`] is the index of
/// the first (oldest) of its `e` events in the **retained** (post-merge) stream —
/// the coordinate space everything downstream is measured in. Its raw segmented
/// index is kept alongside as [`Emer::event_pos`], which is what RawHash2 stamps
/// (`rsketch.c:214-218`) and what `tests/hash_compat.rs` pins against the golden.
/// See [`retained_quantized`] for why `pos` is not the raw index.
pub fn dense_emers(events: &[f32], params: &Params) -> Vec<Emer> {
    let retained = retained_quantized(events, params);
    dense_emers_from_retained(&retained, params.sketch.e as usize, params.quant.q)
}

/// Pack every length-`e` window of an already-quantized stream into an e-mer.
/// Shared by [`dense_emers`] and [`crate::emer::select`].
pub fn dense_emers_from_retained(retained: &[QEvent], e: usize, q: u32) -> Vec<Emer> {
    if e == 0 || retained.len() < e {
        return Vec::new();
    }
    let buckets: Vec<u32> = retained.iter().map(|r| r.bucket).collect();
    (0..=retained.len() - e)
        .map(|i| Emer {
            hash: emer_hash(&buckets[i..i + e], q),
            pos: retained[i].pos,
            event_pos: retained[i].event_pos,
        })
        .collect()
}

/// Select minimizers: for each window of `w` consecutive e-mers keep the
/// minimum-hash e-mer, de-duplicating by event index.
pub fn minimizers(emers: &[Emer], w: usize) -> Vec<Emer> {
    if w <= 1 {
        return emers.to_vec();
    }
    if emers.len() < w {
        // Fewer e-mers than a window: keep the single global minimum.
        return emers
            .iter()
            .min_by(|a, b| a.hash.cmp(&b.hash).then(a.pos.cmp(&b.pos)))
            .copied()
            .into_iter()
            .collect();
    }

    let mut out = Vec::new();
    let mut last_pos = u32::MAX;
    for start in 0..=emers.len() - w {
        let window = &emers[start..start + w];
        let (off, m) = window
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.hash.cmp(&b.hash).then(a.pos.cmp(&b.pos)))
            .unwrap();
        let _ = off;
        if m.pos != last_pos {
            out.push(*m);
            last_pos = m.pos;
        }
    }
    out
}

/// Sketch a read's events into seeds, using `params.sketch.method` to choose
/// the selection scheme (default [`crate::emer::select::SketchMethod::Minimizer`]).
///
/// All schemes share the same `retained_quantized` front end; see
/// [`crate::emer::select`] for the indel-robust alternatives.
pub fn sketch(events: &[f32], read_id: u32, strand: Strand, params: &Params) -> Vec<Seed> {
    let retained = retained_quantized(events, params);
    sketch_from_retained(&retained, read_id, strand, params)
}

/// Like [`sketch`], but takes an already-quantized stream. A caller that also needs
/// the read's length *in coordinate space* (`retained.len()`, e.g. to size the bin
/// grid) can then quantize once instead of twice.
pub fn sketch_from_retained(
    retained: &[QEvent],
    read_id: u32,
    strand: Strand,
    params: &Params,
) -> Vec<Seed> {
    use crate::emer::select::{SelectCfg, select};
    let cfg = SelectCfg::from_params(params);
    select(retained, params.sketch.method, &cfg)
        .into_iter()
        .map(|m| Seed::new(m.hash, read_id, m.pos, strand))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Params;

    #[test]
    fn hash64_is_deterministic_and_masked() {
        let mask = (1u64 << 32) - 1;
        let h = hash64(123456789, mask);
        assert_eq!(h, hash64(123456789, mask));
        assert!(h <= mask);
    }

    #[test]
    fn distinct_emers_mostly_distinct_hashes() {
        let mask = (1u64 << 32) - 1;
        let a = hash64(0b001_010_011, mask);
        let b = hash64(0b011_010_001, mask);
        assert_ne!(a, b);
    }

    #[test]
    fn dense_sketch_has_expected_count() {
        // No diff collapsing: distinct ramp so consecutive quant values differ.
        let mut p = Params::r10_4_1();
        p.sketch.w = 0; // dense
        p.sketch.diff = -1; // disable sig-diff filter
        let events: Vec<f32> = (0..50).map(|i| (i as f32 - 25.0) * 0.1).collect();
        let seeds = sketch(&events, 0, Strand::Forward, &p);
        // dense e-mers = n_events - e + 1 (when nothing is filtered)
        assert_eq!(seeds.len(), events.len() - p.sketch.e as usize + 1);
        assert!(seeds.iter().all(|s| s.read_id() == 0));
    }

    #[test]
    fn minimizer_sketch_is_subset_and_ordered() {
        let mut p = Params::r10_4_1();
        p.sketch.diff = -1;
        let events: Vec<f32> = (0..200)
            .map(|i| ((i * 7) % 11) as f32 * 0.2 - 1.0)
            .collect();
        let dense_p = Params {
            sketch: crate::params::SketchParams { w: 0, ..p.sketch },
            ..p
        };
        let dense = sketch(&events, 1, Strand::Forward, &dense_p);
        let mins = sketch(&events, 1, Strand::Forward, &p);
        assert!(mins.len() <= dense.len());
        assert!(!mins.is_empty());
        // positions strictly increasing (dedup by pos)
        for pair in mins.windows(2) {
            assert!(pair[0].pos < pair[1].pos);
        }
    }

    /// A seed is 12 bytes, an index hit 8, and every field reads back: hashes up to 2^32 - 1, read
    /// ids up to 2^31 - 1, both strands, renumbering keeps the strand.
    #[test]
    fn packed_seeds_read_back() {
        use crate::overlap::SeedHit;
        assert_eq!(std::mem::size_of::<Seed>(), 12);
        assert_eq!(std::mem::size_of::<SeedHit>(), 8);
        for &(h, r, p) in &[
            (0u64, 0u32, 0u32),
            (HASH_MASK, (1 << 31) - 1, u32::MAX),
            (12345, 7, 99),
        ] {
            for st in [Strand::Forward, Strand::Reverse] {
                let mut s = Seed::new(h, r, p, st);
                assert_eq!((s.hash(), s.read_id(), s.pos, s.strand()), (h, r, p, st));
                s.set_read_id(r / 2);
                assert_eq!(
                    (s.hash(), s.read_id(), s.pos, s.strand()),
                    (h, r / 2, p, st)
                );
                let x = SeedHit::new(r, p, st);
                assert_eq!((x.read_id(), x.pos, x.strand()), (r, p, st));
            }
        }
    }

    #[test]
    #[should_panic(expected = "wider than 32 bits")]
    fn a_wide_hash_is_refused() {
        Seed::new(HASH_MASK + 1, 0, 0, Strand::Forward);
    }

    #[test]
    #[should_panic(expected = "below 2^31")]
    fn a_read_id_of_2_pow_31_is_refused() {
        Seed::new(1, 1 << 31, 0, Strand::Forward);
    }

    #[test]
    fn too_few_events_yield_no_seeds() {
        let p = Params::r10_4_1();
        let events = vec![0.1f32; (p.sketch.e - 1) as usize]; // fewer than e
        assert!(sketch(&events, 0, Strand::Forward, &p).is_empty());
    }

    /// Under the merge (`diff >= 0`, the default) a seed's `pos` is its index
    /// in the **retained** stream, NOT the raw event index, which is kept separately
    /// as `event_pos`.
    ///
    /// This is the only sketch test that pins it: the others set `diff = -1`, which
    /// disables the merge and makes `pos == event_pos` identically. Repeating each
    /// value collapses runs, so the two coordinates must diverge here.
    #[test]
    fn merged_seed_pos_is_the_retained_index_not_the_raw_one() {
        let mut p = Params::r10_4_1();
        p.sketch.diff = 0; // the default: collapse alike-quantizing runs
        // Cycle levels that land in DISTINCT buckets (the fine band is ~0.67 wide,
        // so small steps would collapse everything and leave < e events), repeating
        // each 3x so the merge keeps exactly 1 of every 3.
        const LEVELS: [f32; 5] = [-3.0, -1.0, 0.0, 1.0, 3.0];
        let events: Vec<f32> = (0..40)
            .flat_map(|i| {
                let v = LEVELS[i % LEVELS.len()];
                [v, v, v]
            })
            .collect();

        let retained = retained_quantized(&events, &p);
        assert!(
            retained.len() < events.len(),
            "fixture must actually trigger the merge (got {} of {})",
            retained.len(),
            events.len()
        );
        // pos is dense 0..n over the retained stream; event_pos keeps the raw index.
        for (i, q) in retained.iter().enumerate() {
            assert_eq!(q.pos, i as u32, "pos must be the retained index");
        }
        assert!(
            retained.iter().any(|q| q.pos != q.event_pos),
            "the two coordinates must diverge under the merge — this is the whole \
             point of the refactor"
        );
        assert!(
            retained.windows(2).all(|w| w[0].event_pos < w[1].event_pos),
            "event_pos stays strictly increasing in raw-stream order"
        );

        // And the property propagates to the packed e-mers / seeds.
        let emers = dense_emers(&events, &p);
        assert!(!emers.is_empty());
        assert!(
            emers.iter().any(|e| e.pos != e.event_pos),
            "Emer::pos must carry the retained coordinate, not the raw index"
        );
        for e in &emers {
            assert!(
                (e.pos as usize) < retained.len(),
                "pos indexes the retained stream, so it is bounded by its length"
            );
        }
    }
}
