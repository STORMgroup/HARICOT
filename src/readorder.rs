//! Read order and tie-breaks.
//!
//! **Canonical read numbering** (opt-in, `haricot-detect --canonical-order`): reads numbered by their read id
//! instead of their position in the pod5. Every downstream step keys on the read number: the seed
//! index, the "lower id owns the pair" rule of the per-read pair streaming
//! (`overlap::find_overlap_candidates`), the candidate cache's read table. [`canonicalize`] renumbers
//! the sketches right after sketching: a read's seeds do not depend on its number (the id is only
//! stored, `emer::sketch_from_retained`), so this equals sketching the reads in read-id order. The key
//! is the read id string in byte order (pod5 read ids are UUIDs, so byte order of the canonical
//! lower-case text is the order of the 128-bit values).
//!
//! **The tie rank** (the canonical pair frame and the content-hash ties of `haricot assemble`): where
//! two reads must be put in an order that their lengths or coverages do not decide (the frame of an
//! equal-length pair, the graph build's equal-coverage commits, the contig-end rescue's and the end
//! bridges' score ties), the order is [`rank_by_content`]: by a hash of each read's own signal samples
//! ([`content_hash`]), never its id or its place in the file, so every output is a function of the
//! reads' data alone: invariant to the pod5 order AND to how the reads are named. A read-id rank
//! ([`rank_by_name`]) would make a renaming of the reads change the graph. Two reads with the same
//! content hash (an identical signal, i.e. a duplicated read, or a 64-bit collision) are an input
//! error: the order must be total.

use anyhow::{Result, bail};

use crate::overlap::Overlap;
use crate::pipeline::ReadSketch;

/// Errors on a read id that occurs twice. Every step after the detection keys reads by id (the
/// candidate cache's name table, the feature and score tables, the drop lists, the rescue's accept
/// list), and [`rank_by_name`] gives two equal ids the order of their places in the input, so a
/// repeated id would silently merge two reads in the name-keyed steps and make the canonical
/// tie-breaks depend on the input order. The detection checks it whatever the numbering
/// (`haricot-detect`), not only under [`canonicalize`].
pub fn ensure_distinct_ids<'a>(ids: impl IntoIterator<Item = &'a str>) -> Result<()> {
    let mut seen: rustc_hash::FxHashSet<&str> = rustc_hash::FxHashSet::default();
    for id in ids {
        if !seen.insert(id) {
            bail!("read id {id} occurs twice in the input");
        }
    }
    Ok(())
}

/// A read's signal content hash: 64 bits over its samples (the bit patterns of the calibrated pA
/// values, in order) and their count, and nothing else: not the read id, not the read's place in the
/// file. Dependency-free and cheap (it must not cost detection time): four
/// independent lanes (sample `i` into lane `i % 4`), each folding a sample in by a xor, a multiply
/// and a rotation, then the lanes and the count chained through splitmix64's mixer. Computed once
/// per read while the sketch pass holds its signal (`pipeline::sketch_reads_from` over `f32` reads,
/// `pipeline::sketch_adc_from` over ADC reads: the same pA samples); two different signals share a
/// hash with probability about 2^-64.
pub fn content_hash(signal: &[f32]) -> u64 {
    const K: u64 = 0x9e37_79b9_7f4a_7c15;
    fn fold(h: u64, x: f32) -> u64 {
        (h ^ x.to_bits() as u64).wrapping_mul(K).rotate_left(29)
    }
    fn mix(mut z: u64) -> u64 {
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    let mut h: [u64; 4] = [
        0x243f_6a88_85a3_08d3,
        0x1319_8a2e_0370_7344,
        0xa409_3822_299f_31d0,
        0x082e_fa98_ec4e_6c89,
    ];
    let mut quads = signal.chunks_exact(4);
    for q in &mut quads {
        for (hl, &x) in h.iter_mut().zip(q) {
            *hl = fold(*hl, x);
        }
    }
    for (hl, &x) in h.iter_mut().zip(quads.remainder()) {
        *hl = fold(*hl, x);
    }
    let mut z = mix((signal.len() as u64).wrapping_mul(K));
    for lane in h {
        z = mix(z ^ lane);
    }
    z
}

/// The first two reads (indices into `hashes`) that share a content hash, if any.
pub fn content_duplicate(hashes: &[u64]) -> Option<(usize, usize)> {
    let mut idx: Vec<usize> = (0..hashes.len()).collect();
    idx.sort_unstable_by_key(|&i| (hashes[i], i));
    idx.windows(2)
        .find(|w| hashes[w[0]] == hashes[w[1]])
        .map(|w| (w[0], w[1]))
}

/// Each read's **tie rank**: its place among the reads ordered by [`content_hash`] (`rank[i]` for
/// read `i`, `hashes[i]` its hash). Errors when two reads share a hash: the order would not be total.
pub fn rank_by_content(hashes: &[u64]) -> Result<Vec<u32>> {
    if let Some((x, y)) = content_duplicate(hashes) {
        bail!(
            "reads {x} and {y} have the same signal content hash {:016x}: an identical signal (a \
             duplicated read) or a 64-bit collision; the tie order would not be total",
            hashes[x]
        );
    }
    let mut idx: Vec<usize> = (0..hashes.len()).collect();
    idx.sort_unstable_by_key(|&i| hashes[i]);
    let mut rank = vec![0u32; hashes.len()];
    for (r, &i) in idx.iter().enumerate() {
        rank[i] = r as u32;
    }
    Ok(rank)
}

/// Sort `sketches` by read id and renumber them `0..n` (the sketch's `read_id` and every seed's).
/// Errors on a duplicated read id, which no numbering could make unambiguous.
pub fn canonicalize(sketches: &mut [ReadSketch]) -> Result<()> {
    sketches.sort_by(|x, y| x.name.as_bytes().cmp(y.name.as_bytes()));
    if let Some(w) = sketches.windows(2).find(|w| w[0].name == w[1].name) {
        bail!("read id {} occurs twice in the input", w[0].name);
    }
    for (i, s) in sketches.iter_mut().enumerate() {
        let id = i as u32;
        s.read_id = id;
        for seed in s.seeds.iter_mut() {
            seed.set_read_id(id);
        }
    }
    Ok(())
}

/// Each read's rank in read-id order (`rank[i]` = the place of `names[i]` among the sorted read ids).
/// Not the tie rank: it would make the outputs depend on the read
/// names ([`rank_by_content`] is the tie rank); kept for diagnostics. The identity when the reads
/// are numbered canonically.
pub fn rank_by_name(names: &[String]) -> Vec<u32> {
    let mut idx: Vec<usize> = (0..names.len()).collect();
    idx.sort_by(|&x, &y| names[x].as_bytes().cmp(names[y].as_bytes()));
    let mut rank = vec![0u32; names.len()];
    for (r, &i) in idx.iter().enumerate() {
        rank[i] = r as u32;
    }
    rank
}

/// **The canonical pair frame** (opt-in, `OverlapParams::symmetric_pairs`; on in `haricot assemble`):
/// whether read `x` (length `len_x` in retained events, tie rank `rank_x`, [`rank_by_content`]) is the
/// QUERY of its pair with `y`, i.e. the read whose coordinates are the first axis (`a`) of every
/// computation that needs one geometry: the Hough fit and its gates, the stage-1 / stage-2 features,
/// the fill, the contig-end rescue scan. The longer read is the query; on equal lengths, the smaller
/// tie rank.
/// A rule on the pair's own data, so for any numbering exactly one of `query_first(x, y)` and
/// `query_first(y, x)` holds and the frame of a pair does not depend on which read owns it.
pub fn query_first(len_x: u32, len_y: u32, rank_x: u32, rank_y: u32) -> bool {
    len_x > len_y || (len_x == len_y && rank_x < rank_y)
}

/// The same overlap seen from the other read: ids and spans swapped, every path cell transposed
/// (re-sorted), the slope inverted; the detector's counts (`n_anchors`, `matched`, `cnt`,
/// `aligned_bins`, `gap`) and `z` as they are. Transposing twice restores every field exactly
/// except the slope, which becomes `1 / (1 / s)` (equal to `s` up to its last bits).
pub fn transpose_overlap(o: &Overlap) -> Overlap {
    let mut path: Vec<(u32, u32)> = o.path.iter().map(|&(x, y)| (y, x)).collect();
    path.sort_unstable();
    Overlap {
        a: o.b,
        b: o.a,
        a_bin_span: o.b_bin_span,
        b_bin_span: o.a_bin_span,
        path,
        slope: 1.0 / o.slope,
        ..o.clone()
    }
}

/// `o` (stored with `o.a` < `o.b`, the owner first) in its canonical frame: as it is when `o.a` is
/// the query, transposed otherwise. `lens` / `rank`: every read's length and tie rank.
pub fn in_canonical_frame(o: &Overlap, lens: &[u32], rank: &[u32]) -> Overlap {
    let (a, b) = (o.a as usize, o.b as usize);
    if query_first(lens[a], lens[b], rank[a], rank[b]) {
        o.clone()
    } else {
        transpose_overlap(o)
    }
}

/// Whether `sketches` are already numbered canonically (sorted by read id, numbered `0..n`).
pub fn is_canonical(sketches: &[ReadSketch]) -> bool {
    sketches
        .windows(2)
        .all(|w| w[0].name.as_bytes() < w[1].name.as_bytes())
        && sketches
            .iter()
            .enumerate()
            .all(|(i, s)| s.read_id == i as u32 && s.seeds.iter().all(|x| x.read_id() == i as u32))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::emer::Seed;
    use crate::params::Strand;

    /// A sketch as compared here: id, name, length, seeds as (hash, read id, pos).
    type SketchKey = (u32, String, u32, Vec<(u64, u32, u32)>);

    fn sk(name: &str, id: u32, n: u32) -> ReadSketch {
        ReadSketch {
            read_id: id,
            name: name.to_string(),
            n_samples: 100 + n,
            n_events: 10 + n,
            n_coord_events: 8 + n,
            content: n as u64,
            seeds: (0..n)
                .map(|p| Seed::new((p as u64) * 7 + n as u64, id, p, Strand::Forward))
                .collect(),
        }
    }

    #[test]
    fn any_input_order_gives_the_same_numbering() {
        let names = ["c3", "a1", "b2", "a0", "d4"];
        let base: Vec<ReadSketch> = names
            .iter()
            .enumerate()
            .map(|(i, n)| sk(n, i as u32, i as u32 + 1))
            .collect();
        let mut fwd = base.clone();
        canonicalize(&mut fwd).unwrap();
        assert!(is_canonical(&fwd));
        let mut rev: Vec<ReadSketch> = base.iter().rev().cloned().collect();
        for (i, s) in rev.iter_mut().enumerate() {
            s.read_id = i as u32;
        }
        canonicalize(&mut rev).unwrap();
        let key = |v: &[ReadSketch]| -> Vec<SketchKey> {
            v.iter()
                .map(|s| {
                    (
                        s.read_id,
                        s.name.clone(),
                        s.n_coord_events,
                        s.seeds
                            .iter()
                            .map(|x| (x.hash(), x.read_id(), x.pos))
                            .collect(),
                    )
                })
                .collect()
        };
        assert_eq!(key(&fwd), key(&rev));
        let order: Vec<&str> = fwd.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(order, ["a0", "a1", "b2", "c3", "d4"]);
    }

    #[test]
    fn duplicate_read_ids_are_an_error() {
        let mut v = vec![sk("x", 0, 1), sk("y", 1, 1), sk("x", 2, 1)];
        assert!(canonicalize(&mut v).is_err());
    }

    /// The check without renumbering (the detection runs it under any numbering). Without it a
    /// repeated id ranks by its place in the input, the order dependence the ranks exist to remove.
    #[test]
    fn a_repeated_id_is_an_error_in_any_numbering() {
        assert!(ensure_distinct_ids(["a", "b", "c"]).is_ok());
        assert!(ensure_distinct_ids(std::iter::empty::<&str>()).is_ok());
        let e = ensure_distinct_ids(["b", "a", "b"]).unwrap_err();
        assert!(e.to_string().contains("read id b occurs twice"), "{e}");
        let two = |x: &str, y: &str| vec![x.to_string(), y.to_string()];
        assert_eq!(rank_by_name(&two("x", "x")), vec![0, 1]);
    }

    /// The tie rank is a function of the signal alone: one changed sample or one more sample moves
    /// the hash, the same samples give the same hash whatever the read is called, and a repeated
    /// signal is an error.
    #[test]
    fn the_tie_rank_follows_the_signal_content() {
        let a: Vec<f32> = (0..500).map(|i| 60.0 + (i % 37) as f32 * 0.5).collect();
        let mut b = a.clone();
        b[250] += 0.25;
        let mut c = a.clone();
        c.push(a[0]);
        let (ha, hb, hc) = (content_hash(&a), content_hash(&b), content_hash(&c));
        assert_eq!(ha, content_hash(&a[..]));
        assert!(ha != hb && ha != hc && hb != hc);
        assert_ne!(content_hash(&[]), content_hash(&[0.0]));
        let hs = [hb, ha, hc];
        let rank = rank_by_content(&hs).unwrap();
        let mut order: Vec<usize> = (0..hs.len()).collect();
        order.sort_by_key(|&i| rank[i]);
        assert!(
            order.windows(2).all(|w| hs[w[0]] < hs[w[1]]),
            "ranks follow the hashes"
        );
        let mut r = rank;
        r.sort_unstable();
        assert_eq!(r, [0, 1, 2]);
        assert_eq!(content_duplicate(&[ha, hb, ha]), Some((0, 2)));
        let e = rank_by_content(&[ha, hb, ha]).unwrap_err();
        assert!(e.to_string().contains("same signal content hash"), "{e}");
    }

    #[test]
    fn ranks_follow_the_read_ids() {
        fn own(v: &[&str]) -> Vec<String> {
            v.iter().map(|s| s.to_string()).collect()
        }
        assert_eq!(rank_by_name(&own(&["c", "a", "d", "b"])), vec![2, 0, 3, 1]);
        assert_eq!(rank_by_name(&own(&["a", "b", "c"])), vec![0, 1, 2]);
    }

    #[test]
    fn exactly_one_read_of_a_pair_is_the_query() {
        for (lx, ly, rx, ry) in [(10, 9, 0, 1), (9, 10, 0, 1), (10, 10, 0, 1), (10, 10, 1, 0)] {
            assert_ne!(query_first(lx, ly, rx, ry), query_first(ly, lx, ry, rx));
        }
        assert!(query_first(10, 9, 5, 0), "the longer read is the query");
        assert!(query_first(10, 10, 0, 1), "equal lengths: the smaller rank");
    }

    fn ovl(a: u32, b: u32) -> Overlap {
        Overlap {
            a,
            b,
            strand: Strand::Forward,
            n_anchors: 40,
            matched: 900,
            cnt: 31,
            aligned_bins: 12,
            bin_events: 192,
            gap: 2,
            slope: 0.93,
            a_bin_span: (3, 16),
            b_bin_span: (0, 11),
            path: vec![(3, 0), (4, 1), (6, 2), (7, 4), (9, 6), (12, 8), (16, 11)],
            z: 7.5,
        }
    }

    /// One pair stored by either owner has the same canonical frame: numbering 1 stores it as
    /// (0, 1) with read 0 = "x"; numbering 2 numbers x as 1, so stores the transposed record (0, 1)
    /// with read 0 = "y". Both canonical frames are the same alignment of the same two reads.
    #[test]
    fn the_canonical_frame_does_not_depend_on_the_owner() {
        let o1 = ovl(0, 1); // reads: 0 = x (len 20), 1 = y (len 30): y is the query
        let (lens1, rank1) = (vec![20u32, 30], vec![0u32, 1]);
        let mut o2 = transpose_overlap(&o1); // the same alignment stored by the other owner
        (o2.a, o2.b) = (0, 1); // numbering 2: 0 = y, 1 = x
        let (lens2, rank2) = (vec![30u32, 20], vec![1u32, 0]);
        let c1 = in_canonical_frame(&o1, &lens1, &rank1);
        let c2 = in_canonical_frame(&o2, &lens2, &rank2);
        // c1 is in numbering 1's ids, c2 in numbering 2's: the query is y in both
        assert_eq!((c1.a, c1.b), (1, 0));
        assert_eq!((c2.a, c2.b), (0, 1));
        assert_eq!(c1.path, c2.path);
        assert_eq!(
            (c1.a_bin_span, c1.b_bin_span),
            (c2.a_bin_span, c2.b_bin_span)
        );
        assert!((c1.slope - c2.slope).abs() <= 1e-12 * c2.slope.abs());
        assert_eq!(
            (c1.matched, c1.gap, c1.z.to_bits()),
            (c2.matched, c2.gap, c2.z.to_bits())
        );
        let back = transpose_overlap(&transpose_overlap(&o1));
        assert_eq!(
            (back.a, back.b, back.path.clone()),
            (o1.a, o1.b, o1.path.clone())
        );
    }
}
