//! An evidence-based READ MODEL: clip each read at the junctions its own partners reveal.
//!
//! A port of wtdbg2's default-on read clipping, `clip_read_core` / `clip_read_algo`
//! (`extern/wtdbg2/wtdbg.h:741-884`, turned off by `--no-read-clip`, `wtdbg.c:108,346`). A read
//! that is contiguous in the sample has partners whose alignments run across every position of
//! it; a chimera (two loci) or a fold-back (a read that turns back on the other strand) has a
//! position where the partners of one side all stop while the partners of the other side all
//! start, and nobody runs across. The read is cut there and only its best-supported segment
//! keeps its overlaps. Nothing here looks at a reference: the junctions are the reads' own
//! partner evidence, so a split that only the reference sees (a real join in the sample) is
//! spanned by partners and never cut.
//!
//! Per read, from its overlaps (bin spans on the read and on the partner):
//!
//! * **Spur ends** (`clip_read_core`, `wtdbg.h:857-868`, same-strand branch). The left end of
//!   an overlap is a *spur* when both reads continue past it by more than `margin` bins (the
//!   alignment stops while both molecules go on); likewise the right end. An overlap with a
//!   spur at both ends is an internal match -- a repeat copy in wtdbg2's experience -- and,
//!   as there (`wtdbg.h:869-871`), it is neither coverage nor a spur.
//! * **Containment** (`wtdbg.h:850-855`): a read covered end to end (within `margin`) by one
//!   overlap is contiguous by that partner's evidence and is never clipped.
//! * **Spanning coverage.** `S(j)`, at the boundary between bins `j-1` and `j`, counts the
//!   overlaps (not internal matches) that extend more than `margin` bins on BOTH sides of it
//!   (cover bins `j-1-margin ..= j+margin`). So two sides whose ends overlap by up to
//!   `2 * margin` bins of end noise do not count as spanning each other.
//! * **Junction.** A maximal run of boundaries `[ja, jb]` with `S = 0`, flanked on the left by
//!   at least one overlap whose right end is a spur within `margin` bins of where spanning
//!   stopped (`hi` in `ja-1 ± margin`) and on the right by at least one whose left end is a
//!   spur within `margin` of where it resumes (`lo` in `jb ± margin`). The `± margin` window is
//!   wtdbg2's one-bin spur pile (`wtdbg.h:766-771`, `KBM_BIN_SIZE` on each side). A run of any
//!   width is allowed, so an unaligned stretch at the junction (an adapter, a hairpin) is one
//!   junction, not none.
//! * **Clip** (`wtdbg.h:814-818`). The junctions cut the read into segments; an overlap belongs
//!   to the segment its span lies in (within `margin` of the segment's ends); the read keeps
//!   the segment with the most partner support (the summed spans of its overlaps on the read)
//!   and every overlap of the read outside it is dropped, including any that crosses a
//!   junction. An overlap survives when neither of its reads drops it. All reads are decided
//!   from the same input, so the result does not depend on read order.
//!
//! Deviations from wtdbg2:
//! 1. No `min_dep = 2` coverage clip (`wtdbg.h:747-773,836`): a read with no junction is
//!    never clipped. At 6x a depth-1 read end is usually its only dovetail, the same single
//!    adjacency `min_edge_cov >= 2` destroyed (measured at 6x).
//! 2. The junction test is the exact spanning count `S = 0`, where wtdbg2 compares the sweep
//!    depth at one spur event with the spur pile (`wtdbg.h:801-810`); that depth includes one
//!    side's spurs and depends on the sort order of tied events, and `S` is the quantity it
//!    stands for. Hence also no `avg_dep` support threshold: one spur on each side, with no
//!    partner across, is the evidence.
//! 3. Several junctions per read (wtdbg2 cuts at the single strongest one); the kept segment
//!    is chosen by support in consistent units (wtdbg2 compares `c->pos - clps[0]`, which
//!    mixes bp and bins, `wtdbg.h:814`).
use crate::overlap::Overlap;

/// Parameters of [`read_model`]: only wtdbg2's overhang tolerance.
#[derive(Debug, Clone, Copy)]
pub struct ReadModelParams {
    /// Bins a read must continue past an alignment end for that end to be a spur, and the
    /// end-noise tolerance of spanning and of the spur pile: wtdbg2's `max_overhang`
    /// (256 bp = 1 bin, `wtdbg.c:459,962`), which `clip_read_core` takes as `margin`
    /// (`wtdbg.h:835`).
    pub margin: u32,
}

impl Default for ReadModelParams {
    fn default() -> Self {
        Self { margin: 1 }
    }
}

/// One read the model cut.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadClip {
    pub read: u32,
    /// Read length in whole bins (`counts / bin_events`).
    pub n_bins: u32,
    /// The junctions: runs `[ja, jb]` of bin boundaries with no spanning partner.
    pub junctions: Vec<(u32, u32)>,
    /// Index of the kept segment (0 = leftmost; there are `junctions.len() + 1`).
    pub kept: usize,
    /// The kept segment's bin limits (an overlap is in it when `lo >= .0 && hi <= .1`).
    pub kept_bins: (i64, i64),
    /// Overlaps of this read it dropped.
    pub dropped: usize,
    /// Of those, overlaps in no segment: across a junction, or inside a junction's gap.
    pub crossing: usize,
}

/// What [`read_model`] decided.
#[derive(Debug, Clone, Default)]
pub struct ReadModelOutcome {
    /// Per input overlap: does it survive?
    pub keep: Vec<bool>,
    /// The reads cut, by read id.
    pub clips: Vec<ReadClip>,
    /// Reads with at least one overlap.
    pub reads: usize,
    /// Reads that showed a junction but were protected by a containing overlap.
    pub contained_protected: usize,
}

impl ReadModelOutcome {
    pub fn dropped(&self) -> usize {
        self.keep.iter().filter(|&&k| !k).count()
    }
    /// Reads cut at two or more junctions (three or more segments).
    pub fn multi(&self) -> usize {
        self.clips.iter().filter(|c| c.junctions.len() >= 2).count()
    }
    pub fn junctions(&self) -> usize {
        self.clips.iter().map(|c| c.junctions.len()).sum()
    }
}

/// One overlap as one of its reads sees it.
#[derive(Debug, Clone, Copy)]
struct End {
    idx: usize,
    lo: i64,
    hi: i64,
    lspur: bool,
    rspur: bool,
}

/// The decision for one read: its junctions and kept segment, or `None` when it is whole.
struct Decision {
    junctions: Vec<(i64, i64)>,
    contained: bool,
}

/// Apply the read model to an overlap set. `counts[r]` is read `r`'s event count (the
/// cache's first field); spans are in each overlap's `bin_events`.
pub fn read_model(ovl: &[Overlap], counts: &[u32], p: &ReadModelParams) -> ReadModelOutcome {
    let m = p.margin as i64;
    let n_reads = counts.len();
    // Per-read lists of (overlap index, is_a), as a CSR.
    let mut start = vec![0usize; n_reads + 1];
    for o in ovl {
        start[o.a as usize + 1] += 1;
        start[o.b as usize + 1] += 1;
    }
    for r in 0..n_reads {
        start[r + 1] += start[r];
    }
    let mut fill = start.clone();
    let mut ends = vec![(0usize, false); start[n_reads]];
    for (i, o) in ovl.iter().enumerate() {
        for (r, is_a) in [(o.a, true), (o.b, false)] {
            ends[fill[r as usize]] = (i, is_a);
            fill[r as usize] += 1;
        }
    }

    let mut keep = vec![true; ovl.len()];
    let mut out = ReadModelOutcome::default();
    let mut list: Vec<End> = Vec::new();
    for r in 0..n_reads {
        let mine = &ends[start[r]..start[r + 1]];
        if mine.is_empty() {
            continue;
        }
        out.reads += 1;
        list.clear();
        let mut n_r = 0i64;
        for &(i, is_a) in mine {
            let o = &ovl[i];
            let bin = o.bin_events.max(1) as i64;
            let (span, pspan, partner) = if is_a {
                (o.a_bin_span, o.b_bin_span, o.b)
            } else {
                (o.b_bin_span, o.a_bin_span, o.a)
            };
            n_r = counts[r] as i64 / bin;
            let n_p = counts[partner as usize] as i64 / bin;
            let (lo, hi) = (span.0 as i64, span.1 as i64);
            let (plo, phi) = (pspan.0 as i64, pspan.1 as i64);
            // `wtdbg.h:861-863` (same strand): x = both continue on the left, y = on the right.
            let lspur = lo > m && plo > m;
            let rspur = n_r - 1 - hi > m && n_p - 1 - phi > m;
            list.push(End {
                idx: i,
                lo,
                hi,
                lspur,
                rspur,
            });
        }
        let d = decide(&list, n_r, m);
        if d.junctions.is_empty() {
            continue;
        }
        if d.contained {
            out.contained_protected += 1;
            continue;
        }
        let clip = clip_read(r as u32, n_r, &list, &d.junctions, m);
        for e in &list {
            if !in_segment(e, clip.kept_bins) {
                keep[e.idx] = false;
            }
        }
        out.clips.push(clip);
    }
    out.keep = keep;
    out
}

/// Junctions of one read (`clip_read_algo` as adapted; see the module doc).
fn decide(list: &[End], n_r: i64, m: i64) -> Decision {
    // `wtdbg.h:850-855`: contained in one partner -> never clipped.
    let contained = list.iter().any(|e| e.lo <= m && e.hi + 1 + m >= n_r);
    let ext = list
        .iter()
        .map(|e| e.hi + 1)
        .max()
        .unwrap_or(0)
        .max(n_r)
        .max(1) as usize;
    // Spanning coverage over boundaries 0..=ext, and the spur ends by position.
    let mut diff = vec![0i64; ext + 2];
    let mut rs_end = vec![0i64; ext + 1];
    let mut ls_start = vec![0i64; ext + 1];
    for e in list.iter().filter(|e| !(e.lspur && e.rspur)) {
        let (from, to) = (e.lo + 1 + m, e.hi - m);
        if from <= to {
            diff[from as usize] += 1;
            diff[(to + 1) as usize] -= 1;
        }
        if e.rspur {
            rs_end[e.hi as usize] += 1;
        }
        if e.lspur {
            ls_start[e.lo as usize] += 1;
        }
    }
    let window = |v: &[i64], c: i64| -> i64 {
        let (a, b) = ((c - m).max(0), (c + m).min(v.len() as i64 - 1));
        if a > b {
            0
        } else {
            v[a as usize..=b as usize].iter().sum()
        }
    };
    let mut junctions = Vec::new();
    let (mut s, mut run): (i64, Option<i64>) = (0, None);
    // boundaries 1..=ext-1 are interior; S(0) and S(ext) are not boundaries of the read
    s += diff[0];
    for j in 1..ext as i64 {
        s += diff[j as usize];
        let zero = s == 0;
        match (zero, run) {
            (true, None) => run = Some(j),
            (false, Some(ja)) => {
                let jb = j - 1;
                if window(&rs_end, ja - 1) > 0 && window(&ls_start, jb) > 0 {
                    junctions.push((ja, jb));
                }
                run = None;
            }
            _ => {}
        }
    }
    if let Some(ja) = run {
        let jb = ext as i64 - 1;
        if window(&rs_end, ja - 1) > 0 && window(&ls_start, jb) > 0 {
            junctions.push((ja, jb));
        }
    }
    Decision {
        junctions,
        contained,
    }
}

fn in_segment(e: &End, lim: (i64, i64)) -> bool {
    e.lo >= lim.0 && e.hi <= lim.1
}

/// Cut one read at its junctions and keep its best-supported segment.
fn clip_read(read: u32, n_r: i64, list: &[End], junctions: &[(i64, i64)], m: i64) -> ReadClip {
    // Segment k runs from the previous junction's resumption to the next junction's stop,
    // each widened by `margin` (an end within `margin` of a junction is end noise).
    let k_n = junctions.len() + 1;
    let lim = |k: usize| -> (i64, i64) {
        let lo = if k == 0 {
            i64::MIN
        } else {
            junctions[k - 1].1 - m
        };
        let hi = if k + 1 == k_n {
            i64::MAX
        } else {
            junctions[k].0 - 1 + m
        };
        (lo, hi)
    };
    let mut support = vec![0i64; k_n];
    let mut crossing = 0usize;
    for e in list {
        match (0..k_n).find(|&k| in_segment(e, lim(k))) {
            Some(k) => support[k] += e.hi - e.lo + 1,
            None => crossing += 1,
        }
    }
    // Most support; ties to the leftmost (deterministic).
    let kept = (0..k_n)
        .max_by_key(|&k| (support[k], std::cmp::Reverse(k)))
        .unwrap_or(0);
    let kept_bins = lim(kept);
    let dropped = list.iter().filter(|e| !in_segment(e, kept_bins)).count();
    ReadClip {
        read,
        n_bins: n_r.max(0) as u32,
        junctions: junctions
            .iter()
            .map(|&(a, b)| (a as u32, b as u32))
            .collect(),
        kept,
        kept_bins,
        dropped,
        crossing,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Strand;

    const BIN: u32 = 192;
    /// Partner reads are 80 bins long unless stated.
    const P_LEN: u32 = 80;

    /// An overlap between reads `a < b` with the given inclusive bin spans.
    fn ovl(a: u32, b: u32, sa: (u32, u32), sb: (u32, u32)) -> Overlap {
        let n = sa.1 - sa.0 + 1;
        Overlap {
            a,
            b,
            strand: Strand::Forward,
            n_anchors: n,
            matched: n * BIN,
            cnt: n,
            aligned_bins: n,
            bin_events: BIN,
            gap: 0,
            slope: 1.0,
            a_bin_span: sa,
            b_bin_span: sb,
            path: (0..n).map(|k| (sa.0 + k, sb.0 + k)).collect(),
            z: f32::NAN,
        }
    }

    /// Read 0's bins `lo..=hi` against partner `p`, which STARTS at the alignment (its bin 0)
    /// and runs on past `hi`: the right end is a spur when read 0 continues too.
    fn stops_at(p: u32, lo: u32, hi: u32) -> Overlap {
        ovl(0, p, (lo, hi), (0, hi - lo))
    }

    /// Read 0's bins `lo..=hi` against partner `p`, which ENDS at the alignment (its last
    /// bin) and extends before `lo`: the left end is a spur when read 0 extends too.
    fn starts_at(p: u32, lo: u32, hi: u32) -> Overlap {
        ovl(0, p, (lo, hi), (P_LEN - 1 - (hi - lo), P_LEN - 1))
    }

    fn run(o: &[Overlap], bins: &[u32]) -> ReadModelOutcome {
        let counts: Vec<u32> = bins.iter().map(|&b| b * BIN).collect();
        read_model(o, &counts, &ReadModelParams::default())
    }

    /// Read 0 of 100 bins, partners 1.. of `P_LEN`.
    fn lens(n_partners: usize) -> Vec<u32> {
        std::iter::once(100)
            .chain(std::iter::repeat_n(P_LEN, n_partners))
            .collect()
    }

    /// A CHIMERA (read 0 = locus 1 on bins 0..50, locus 2 on bins 50..100) is cut at bin 50:
    /// four partners on the left stop there while continuing, three on the right start
    /// there, nobody spans. The better-supported left keeps its overlaps, the right's go.
    #[test]
    fn chimera_is_split_and_the_weaker_locus_dropped() {
        let mut o: Vec<Overlap> = (1..=4).map(|p| stops_at(p, 5 + 5 * p, 49)).collect();
        o.extend((5..=7).map(|p| starts_at(p, 50, 115 - 5 * p)));
        let r = run(&o, &lens(7));
        assert_eq!(r.clips.len(), 1, "only the chimera is cut: {:?}", r.clips);
        let c = &r.clips[0];
        assert_eq!(c.read, 0);
        assert_eq!(c.junctions.len(), 1);
        let (ja, jb) = c.junctions[0];
        assert!(
            ja <= 50 && jb >= 50,
            "junction at bin 50: {:?}",
            c.junctions
        );
        assert_eq!(c.kept, 0, "the left locus has more support");
        assert_eq!(r.keep, [true, true, true, true, false, false, false]);
        assert_eq!((c.dropped, c.crossing), (3, 0));
    }

    /// A FOLD-BACK: read 0 (80 bins) is forward on bins 0..60 (five partners, two of them
    /// across the read's start), then turns back on bins 62..80 after a 2-bin hairpin with
    /// no placement. The one tail partner (another fold-back's tail, extending to ITS left)
    /// starts at bin 62; the forward partners that reach the fold stop there while
    /// continuing. The tail overlap is clipped; the forward side is untouched.
    #[test]
    fn fold_back_tail_is_clipped() {
        let o = vec![
            ovl(0, 1, (0, 30), (49, 79)), // dovetail over the read's start
            ovl(0, 2, (0, 45), (34, 79)), // dovetail over the read's start
            stops_at(3, 20, 59),
            stops_at(4, 30, 59),
            stops_at(5, 10, 59),
            ovl(0, 6, (62, 79), (40, 57)), // the tail partner
        ];
        let r = run(&o, &[80; 7]);
        assert_eq!(r.clips.len(), 1, "{:?}", r.clips);
        let c = &r.clips[0];
        assert_eq!((c.read, c.kept), (0, 0));
        assert_eq!(r.keep, [true, true, true, true, true, false]);
    }

    /// A CONTIGUOUS read with a REFERENCE-ONLY split. In the reference's view the read
    /// breaks at bin 50 and partners on each side stop / start there (as the per-segment
    /// oracle reports them); in the sample the read is contiguous, so some partners run
    /// straight across bin 50. Nothing is clipped -- and without those partners the same
    /// read would be.
    #[test]
    fn contiguous_read_with_reference_only_split_is_kept_whole() {
        let mut o: Vec<Overlap> = (1..=3).map(|p| stops_at(p, 5 + 5 * p, 49)).collect();
        o.extend((4..=6).map(|p| starts_at(p, 50, 110 - 5 * p)));
        let split_only = o.clone();
        o.push(ovl(0, 7, (30, 75), (0, 45))); // across the split
        o.push(ovl(0, 8, (40, 99), (5, 64))); // across the split, to the read's end
        let r = run(&o, &lens(8));
        assert!(r.clips.is_empty(), "{:?}", r.clips);
        assert!(r.keep.iter().all(|&k| k));
        assert_eq!(run(&split_only, &lens(6)).clips.len(), 1);
    }

    /// End noise: the two sides may share up to `2 * margin` bins and still be a junction;
    /// sharing more is spanning evidence.
    #[test]
    fn sides_sharing_two_bins_are_a_junction_three_are_not() {
        let mk = |left_hi: u32, right_lo: u32| {
            vec![
                stops_at(1, 10, left_hi),
                stops_at(2, 15, left_hi),
                starts_at(3, right_lo, 90),
                starts_at(4, right_lo, 85),
            ]
        };
        // left ends at 50, right starts at 49: bins 49 and 50 are shared
        assert_eq!(run(&mk(50, 49), &lens(4)).clips.len(), 1);
        // left ends at 51, right starts at 49: three shared bins -> they span each other
        assert!(run(&mk(51, 49), &lens(4)).clips.is_empty());
    }

    /// A low-depth BRIDGE: read 0 is the only read joining partner 1 (whose read ENDS inside
    /// it) and partner 2 (whose read STARTS inside it). Nobody spans the hand-over, but
    /// neither end is a spur, so it is not a junction.
    #[test]
    fn dovetail_handover_is_not_a_junction() {
        let o = vec![ovl(0, 1, (0, 49), (30, 79)), ovl(0, 2, (52, 99), (0, 47))];
        let r = run(&o, &lens(2));
        assert!(r.clips.is_empty(), "{:?}", r.clips);
    }

    /// A read covered end to end by one partner (wtdbg2's containment, `wtdbg.h:850-855`) is
    /// never clipped, whatever its other partners show.
    #[test]
    fn contained_read_is_never_clipped() {
        let mut o = vec![
            stops_at(1, 10, 49),
            stops_at(2, 15, 49),
            starts_at(3, 50, 90),
            starts_at(4, 50, 85),
        ];
        o.push(ovl(0, 5, (0, 99), (40, 139)));
        let mut bins = lens(4);
        bins.push(200);
        let r = run(&o, &bins);
        assert!(r.clips.is_empty());
        assert!(r.keep.iter().all(|&k| k));
    }

    /// An internal match (spur at both ends) across a junction neither hides the junction
    /// nor survives it.
    #[test]
    fn internal_match_across_a_junction_is_dropped() {
        let o = vec![
            stops_at(1, 10, 49),
            stops_at(2, 15, 49),
            stops_at(3, 20, 49),
            starts_at(4, 50, 90),
            ovl(0, 5, (40, 60), (20, 40)), // both reads continue on both sides
        ];
        let r = run(&o, &lens(5));
        assert_eq!(r.clips.len(), 1, "{:?}", r.clips);
        assert_eq!(r.clips[0].kept, 0);
        assert_eq!(r.keep, [true, true, true, false, false]);
        assert_eq!(r.clips[0].crossing, 1);
    }

    /// Two junctions: the middle segment is the best supported and the read keeps it.
    #[test]
    fn two_junctions_keep_the_best_segment() {
        let o = vec![
            stops_at(1, 5, 29),   // segment A: bins 0..30
            starts_at(2, 30, 55), // segment B: bins 30..70
            stops_at(3, 45, 69),
            starts_at(4, 30, 60),
            starts_at(5, 70, 95), // segment C: bins 70..100
        ];
        let r = run(&o, &lens(5));
        assert_eq!(r.clips.len(), 1, "{:?}", r.clips);
        let c = &r.clips[0];
        assert_eq!(c.junctions.len(), 2);
        assert_eq!(c.kept, 1);
        assert_eq!(r.keep, [false, true, true, true, false]);
        assert_eq!(r.multi(), 1);
        assert_eq!(r.junctions(), 2);
    }

    /// The decision is per read and from the input alone: the partners of a chimera keep
    /// their own overlaps with each other, and are not cut by the chimera's spurs.
    #[test]
    fn partners_are_not_clipped() {
        let mut o: Vec<Overlap> = (1..=3).map(|p| stops_at(p, 5 + 5 * p, 49)).collect();
        o.extend((4..=5).map(|p| starts_at(p, 50, 110 - 5 * p)));
        o.push(ovl(1, 2, (10, 79), (5, 74))); // partners 1 and 2 overlap each other
        let r = run(&o, &lens(5));
        assert_eq!(r.clips.len(), 1);
        assert_eq!(r.clips[0].read, 0);
        assert_eq!(r.keep, [true, true, true, false, false, true]);
    }
}
