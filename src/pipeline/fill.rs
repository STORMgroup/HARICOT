//! The fill, a step of `haricot-layout`: each cached overlap's path straightened, using nothing but the
//! path itself.
//!
//! `hough_aln` builds its path by collapsing each query bin to its MODAL target bin (`dp::hough_aln`). At
//! ~2 anchors per bin that mode is one or two anchors deep, so each cell carries the full per-bin noise and
//! the path wanders around the true diagonal even when the fit is right; and a query bin without an anchor
//! is a hole. `graph::build` slices the path into K-bin node windows. The fill rewrites each path from a
//! line fitted to the path itself: no truth, no reference, nothing the detector did not already have.
//!
//! A Theil-Sen fit, and a cell for EVERY query bin of the path's span, interpolated from the line.
//! Against the signal-exact truth, detected paths were about 32% holes (99% of overlaps have one) while
//! 93% of the overlaps lie within one bin of the true line, so the missing cells can be interpolated with
//! no new information. Every overlap is rewritten in the canonical pair frame, along the longer read's
//! bins as the query axis ([`straighten_in_frame`]).

use crate::overlap::Overlap;
use rayon::prelude::*;

/// Theil-Sen: the median of pairwise slopes, then the median intercept. Costs
/// `O(n^2)` in the cell count (paths here are tens of cells), and unlike least squares
/// a single badly-placed cell cannot drag the line.
pub fn theil(path: &[(u32, u32)]) -> Option<(f64, f64)> {
    if path.len() < 3 {
        return None;
    }
    let mut sl: Vec<f64> = Vec::with_capacity(path.len() * path.len() / 2);
    for i in 0..path.len() {
        for j in (i + 1)..path.len() {
            let dx = path[j].0 as f64 - path[i].0 as f64;
            if dx.abs() > 0.5 {
                sl.push((path[j].1 as f64 - path[i].1 as f64) / dx);
            }
        }
    }
    if sl.is_empty() {
        return None;
    }
    // The median by selection: the value a full sort puts at `len / 2`, in linear time.
    let k = sl.len() / 2;
    let m = *sl
        .select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap())
        .1;
    Some((m, median_intercept(path, m)))
}

/// The median of the path cells' intercepts `b - m * a` at slope `m` (the upper median).
pub fn median_intercept(path: &[(u32, u32)], m: f64) -> f64 {
    let mut ic: Vec<f64> = path.iter().map(|&(a, b)| b as f64 - m * a as f64).collect();
    let k = ic.len() / 2;
    *ic.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap())
        .1
}

/// Theil-Sen over points `(x, y)`, the stage-1 and stage-2 features' line fit
/// (`haricot-pair-features`, `haricot-graph-features`): the median of the pairwise slopes over
/// pairs with `|dx| > 1e-9`, then the median intercept `y - s * x` (both the upper median, by a
/// full sort). `None` for fewer than 2 points or no pair with distinct `x`.
pub fn theil_points(pts: &[(f64, f64)]) -> Option<(f64, f64)> {
    if pts.len() < 2 {
        return None;
    }
    let mut sl = Vec::with_capacity(pts.len() * (pts.len() - 1) / 2);
    for i in 0..pts.len() {
        for j in i + 1..pts.len() {
            let dx = pts[j].0 - pts[i].0;
            if dx.abs() > 1e-9 {
                sl.push((pts[j].1 - pts[i].1) / dx);
            }
        }
    }
    if sl.is_empty() {
        return None;
    }
    sl.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let s = sl[sl.len() / 2];
    let mut ic: Vec<f64> = pts.iter().map(|p| p.1 - s * p.0).collect();
    ic.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some((s, ic[ic.len() / 2]))
}

/// What one or more overlaps' rewrites did (summed for the report line).
#[derive(Default, Clone, Copy)]
pub struct Stats {
    pub fitted: u64,
    pub moved: u64,
    pub cells: u64,
    pub shifted: u64,
    pub filled: u64,
}

impl Stats {
    fn add(self, o: Stats) -> Stats {
        Stats {
            fitted: self.fitted + o.fitted,
            moved: self.moved + o.moved,
            cells: self.cells + o.cells,
            shifted: self.shifted + o.shifted,
            filled: self.filled + o.filled,
        }
    }
}

/// Rewrite one overlap's path (and the spans it moves).
pub fn straighten(o: &mut Overlap) -> Stats {
    let mut st = Stats::default();
    let Some((m, c)) = theil(&o.path) else {
        return st;
    };
    st.fitted += 1;
    let mut any = false;
    for cell in &mut o.path {
        st.cells += 1;
        let want = (m * cell.0 as f64 + c).round();
        if !want.is_finite() || want < 0.0 {
            continue;
        }
        let d = want as i64 - cell.1 as i64;
        if d == 0 {
            continue;
        }
        cell.1 = want as u32;
        st.shifted += 1;
        any = true;
    }
    {
        // Every query bin in the span gets a cell, interpolated from the line. The
        // query span is unchanged, so `aligned_bins` is unchanged; only the holes
        // that `hough_aln`'s modal-bin collapse left behind are closed.
        let (alo, ahi) = o
            .path
            .iter()
            .fold((u32::MAX, 0u32), |(l, h), &(a, _)| (l.min(a), h.max(a)));
        let mut dense = Vec::with_capacity((ahi - alo + 1) as usize);
        let mut push = |q: u32| {
            let want = (m * q as f64 + c).round();
            if want.is_finite() && want >= 0.0 {
                dense.push((q, want as u32));
            }
        };
        (alo..=ahi).for_each(&mut push);
        if dense.len() > o.path.len() {
            st.filled += (dense.len() - o.path.len()) as u64;
            o.path = dense;
            any = true;
        }
    }
    if any {
        st.moved += 1;
        // The b-span is the only derived field the rewrite can invalidate; the cell
        // count, the a-span, `matched` and `cnt` are all untouched by construction.
        let (lo, hi) = o
            .path
            .iter()
            .fold((u32::MAX, 0u32), |(l, h), &(_, b)| (l.min(b), h.max(b)));
        o.b_bin_span = (lo, hi);
    }
    st
}

/// [`straighten`] in the overlap's canonical pair frame (`frame` = every read's length and tie rank):
/// when the owner `o.a` is not the query (`readorder::query_first`: the longer read), the overlap is
/// transposed, rewritten, and transposed back, so the fill runs along the longer read's bins whichever
/// read owns the pair. The rewrite does not touch the slope, which is restored exactly (no `1 / (1 / s)`).
pub fn straighten_in_frame(o: &mut Overlap, (lens, rank): (&[u32], &[u32])) -> Stats {
    use crate::readorder::{query_first, transpose_overlap};
    if !query_first(
        lens[o.a as usize],
        lens[o.b as usize],
        rank[o.a as usize],
        rank[o.b as usize],
    ) {
        let slope = o.slope;
        let mut t = transpose_overlap(o);
        let st = straighten(&mut t);
        *o = transpose_overlap(&t);
        o.slope = slope;
        st
    } else {
        straighten(o)
    }
}

/// Every overlap rewritten from its own path alone ([`straighten_in_frame`]), in parallel; the stats
/// summed.
pub fn fill_all(ovl: &mut [Overlap], frame: (&[u32], &[u32])) -> Stats {
    ovl.par_iter_mut()
        .map(|o| straighten_in_frame(o, frame))
        .collect::<Vec<Stats>>()
        .into_iter()
        .fold(Stats::default(), Stats::add)
}

/// The report line of a run over `n` overlaps.
pub fn report(n: usize, st: &Stats) -> String {
    let Stats {
        fitted,
        moved,
        cells,
        shifted,
        filled,
    } = *st;
    format!(
        "{n} overlaps, {fitted} fitted (theil), {moved} moved; {shifted}/{cells} cells shifted \
         ({:.1}%), {filled} cells filled",
        100.0 * shifted as f64 / cells.max(1) as f64,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The same pair stored by either owner is filled along the same (longer) read,
    /// so both rewrites are the same alignment seen from the two owners.
    #[test]
    fn the_fill_does_not_depend_on_the_owner() {
        use crate::readorder::{in_canonical_frame, transpose_overlap};
        let path = vec![
            (2, 0),
            (3, 1),
            (6, 3),
            (7, 4),
            (11, 6),
            (12, 7),
            (16, 9),
            (20, 12),
        ];
        let o1 = Overlap {
            a: 0,
            b: 1,
            strand: crate::params::Strand::Forward,
            n_anchors: 40,
            matched: 900,
            cnt: 31,
            aligned_bins: 13,
            bin_events: 192,
            gap: 5,
            slope: 1.0 / 0.62,
            a_bin_span: (2, 20),
            b_bin_span: (0, 12),
            path,
            z: 7.5,
        };
        // numbering 1: read 0 = x, the longer read (the owner is the query); numbering 2
        // renumbers x as 1 and stores the transposed record
        let mut o2 = transpose_overlap(&o1);
        (o2.a, o2.b) = (0, 1);
        let (lens1, rank1): (&[u32], &[u32]) = (&[40, 30], &[0, 1]);
        let (lens2, rank2): (&[u32], &[u32]) = (&[30, 40], &[1, 0]);
        let (mut r1, mut r2) = (o1.clone(), o2.clone());
        straighten_in_frame(&mut r1, (lens1, rank1));
        straighten_in_frame(&mut r2, (lens2, rank2));
        assert_eq!(
            r2.slope.to_bits(),
            o2.slope.to_bits(),
            "the slope is not rewritten"
        );
        let c1 = in_canonical_frame(&r1, lens1, rank1);
        let c2 = in_canonical_frame(&r2, lens2, rank2);
        assert_eq!(
            c1.path.len(),
            19,
            "every bin of the longer read's span gets a cell"
        );
        assert_eq!(c1.path, c2.path);
        assert_eq!(
            (c1.a_bin_span, c1.b_bin_span),
            (c2.a_bin_span, c2.b_bin_span)
        );
        // in the owner's frame, numbering 2's owner would be filled along the shorter read: another path
        let mut r3 = o2.clone();
        straighten(&mut r3);
        assert_ne!(in_canonical_frame(&r3, lens2, rank2).path, c1.path);
    }
}
