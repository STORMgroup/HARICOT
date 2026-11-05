//! The fill, a step of `haricot-layout`: each cached overlap's path straightened, using nothing but the
//! path itself (`#[path]`-included like `cache.rs`).
//!
//! `hough_aln` builds its path by collapsing each query bin to its MODAL target bin (`dp::hough_aln`). At
//! ~2 anchors per bin that mode is one or two anchors deep, so each cell carries the full per-bin noise and
//! the path wanders around the true diagonal even when the fit is right; and a query bin without an anchor
//! is a hole. `graph::build` slices the path into K-bin node windows. The fill rewrites each path from a
//! line fitted to the path itself: no truth, no reference, nothing the detector did not already have.
//!
//! `haricot assemble` uses [`DRIVER`]: a Theil-Sen fit, and a cell for EVERY query bin of the path's span,
//! interpolated from the line. Against the signal-exact truth, detected paths were about 32% holes (99% of
//! overlaps have one) while 93% of the overlaps lie within one bin of the true line, so the missing cells
//! can be interpolated with no new information. In the canonical pair frame (`--symmetric`,
//! `haricot-layout`'s default) every overlap is rewritten along the longer read's bins as the query axis
//! ([`straighten_in_frame`]).
//!
//! [`Cfg`] also holds the alternatives compared during development: least squares (`lsq`) and Siegel's
//! repeated median (`rm`; breakdown 50% vs Theil-Sen's ~29%) fits; a cap on how far a cell may move
//! (`max_shift`); cutting the path at every run of >= `split_hole` empty query bins and keeping only the
//! piece with the most cells (`hough_aln` has no notion of where an alignment ENDS, so stray in-band
//! anchors past the true end extend the span); and filling only the holes shorter than `fill_max_hole`
//! bins.
#![allow(dead_code)]

use anyhow::Result;
use haricot::overlap::Overlap;
use rayon::prelude::*;
use std::io::Write;

/// Ordinary least squares through the path cells.
pub fn lsq(path: &[(u32, u32)]) -> Option<(f64, f64)> {
    let n = path.len() as f64;
    if path.len() < 3 {
        return None;
    }
    let (sx, sy) = path
        .iter()
        .fold((0.0, 0.0), |(x, y), &(a, b)| (x + a as f64, y + b as f64));
    let (mx, my) = (sx / n, sy / n);
    let (mut num, mut den) = (0.0, 0.0);
    for &(a, b) in path {
        let dx = a as f64 - mx;
        num += dx * (b as f64 - my);
        den += dx * dx;
    }
    (den > 0.0).then(|| (num / den, my - (num / den) * mx))
}

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

pub fn median_intercept(path: &[(u32, u32)], m: f64) -> f64 {
    let mut ic: Vec<f64> = path.iter().map(|&(a, b)| b as f64 - m * a as f64).collect();
    let k = ic.len() / 2;
    *ic.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap())
        .1
}

/// Siegel's repeated median: for each cell the median slope to every other cell, then
/// the median of those; the intercept as Theil-Sen's. Breaks down at 50% off-line cells
/// where Theil-Sen breaks down at ~29%, at the same `O(n^2)` cost.
pub fn repeated_median(path: &[(u32, u32)]) -> Option<(f64, f64)> {
    if path.len() < 3 {
        return None;
    }
    let mut per: Vec<f64> = Vec::with_capacity(path.len());
    let mut sl: Vec<f64> = Vec::with_capacity(path.len());
    for i in 0..path.len() {
        sl.clear();
        for (j, q) in path.iter().enumerate() {
            let dx = q.0 as f64 - path[i].0 as f64;
            if j != i && dx.abs() > 0.5 {
                sl.push((q.1 as f64 - path[i].1 as f64) / dx);
            }
        }
        if !sl.is_empty() {
            sl.sort_by(|a, b| a.partial_cmp(b).unwrap());
            per.push(sl[sl.len() / 2]);
        }
    }
    if per.is_empty() {
        return None;
    }
    per.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = per[per.len() / 2];
    Some((m, median_intercept(path, m)))
}

/// Cells more than one bin off the line `(m, c)`.
pub fn off_line(path: &[(u32, u32)], (m, c): (f64, f64)) -> usize {
    path.iter()
        .filter(|&&(a, b)| (b as f64 - (m * a as f64 + c)).abs() > 1.0)
        .count()
}

/// Split `path` (sorted by query bin) at every run of >= `h` empty query bins and keep
/// the piece with the most cells (ties: the wider query span, then the first). Returns
/// the number of cells dropped.
pub fn keep_main_piece(path: &mut Vec<(u32, u32)>, h: u32) -> usize {
    if path.len() < 2 {
        return 0;
    }
    let mut pieces: Vec<(usize, usize)> = Vec::new(); // [start, end) into path
    let mut start = 0usize;
    for i in 1..path.len() {
        if path[i].0 - path[i - 1].0 > h {
            pieces.push((start, i));
            start = i;
        }
    }
    pieces.push((start, path.len()));
    if pieces.len() == 1 {
        return 0;
    }
    let width = |p: &(usize, usize)| path[p.1 - 1].0 - path[p.0].0;
    let mut best = pieces[0];
    for p in &pieces[1..] {
        let (n, nb) = (p.1 - p.0, best.1 - best.0);
        if n > nb || (n == nb && width(p) > width(&best)) {
            best = *p;
        }
    }
    let dropped = path.len() - (best.1 - best.0);
    *path = path[best.0..best.1].to_vec();
    dropped
}

/// The options of one run.
pub struct Cfg<'a> {
    pub mode: &'a str,
    pub max_shift: i64,
    pub fill: bool,
    pub split_hole: Option<u32>,
    pub fill_max_hole: Option<u32>,
}

/// What one or more overlaps' rewrites did (summed for the report line).
#[derive(Default, Clone, Copy)]
pub struct Stats {
    pub fitted: u64,
    pub moved: u64,
    pub cells: u64,
    pub shifted: u64,
    pub capped: u64,
    pub filled: u64,
    pub trimmed: u64,
    pub trim_cells: u64,
}

impl Stats {
    pub fn add(self, o: Stats) -> Stats {
        Stats {
            fitted: self.fitted + o.fitted,
            moved: self.moved + o.moved,
            cells: self.cells + o.cells,
            shifted: self.shifted + o.shifted,
            capped: self.capped + o.capped,
            filled: self.filled + o.filled,
            trimmed: self.trimmed + o.trimmed,
            trim_cells: self.trim_cells + o.trim_cells,
        }
    }
}

/// Rewrite one overlap's path (and the spans it moves), writing its `--census` row if asked.
pub fn straighten(
    o: &mut Overlap,
    cfg: &Cfg,
    census: Option<&mut std::io::BufWriter<std::fs::File>>,
    names: &[String],
) -> Result<Stats> {
    let Cfg {
        mode,
        max_shift,
        fill,
        split_hole,
        fill_max_hole,
    } = *cfg;
    let mut st = Stats::default();
    let mut census = census;
    if let Some(h) = split_hole {
        o.path.sort_unstable();
        let dropped = keep_main_piece(&mut o.path, h);
        if dropped > 0 {
            st.trimmed += 1;
            st.trim_cells += dropped as u64;
            // The span is what `graph::build` slices into node windows, so it must
            // follow the cells. `n_anchors`/`matched`/`cnt`/`z` still describe the
            // untrimmed detection (the cache carries no anchors to recount).
            let fold = |f: fn(&(u32, u32)) -> u32| {
                o.path
                    .iter()
                    .map(f)
                    .fold((u32::MAX, 0u32), |(l, h), x| (l.min(x), h.max(x)))
            };
            o.a_bin_span = fold(|c| c.0);
            o.b_bin_span = fold(|c| c.1);
            let aligned =
                (o.a_bin_span.1 - o.a_bin_span.0 + 1).min(o.b_bin_span.1 - o.b_bin_span.0 + 1);
            o.aligned_bins = aligned;
            o.gap = aligned.saturating_sub(o.path.len() as u32);
        }
    }
    if let Some(w) = census.as_mut() {
        // --census: the path as detected (after any --split-hole), both robust fits.
        if let (Some(t), Some(r)) = (theil(&o.path), repeated_median(&o.path)) {
            let (alo, ahi) = o
                .path
                .iter()
                .fold((u32::MAX, 0u32), |(l, h), &(a, _)| (l.min(a), h.max(a)));
            let at = |(m, c): (f64, f64), q: u32| m * q as f64 + c;
            let dend = (at(t, alo) - at(r, alo))
                .abs()
                .max((at(t, ahi) - at(r, ahi)).abs());
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.2}",
                names[o.a as usize],
                names[o.b as usize],
                o.path.len(),
                off_line(&o.path, t),
                off_line(&o.path, r),
                t.0,
                r.0,
                dend
            )?;
        }
    }
    let Some((m, c)) = (match mode {
        "lsq" => lsq(&o.path),
        "rm" => repeated_median(&o.path),
        _ => theil(&o.path),
    }) else {
        return Ok(st);
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
        if d.abs() > max_shift {
            st.capped += 1;
            continue;
        }
        cell.1 = want as u32;
        st.shifted += 1;
        any = true;
    }
    if fill {
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
        match fill_max_hole {
            None => (alo..=ahi).for_each(&mut push),
            Some(h) => {
                // Interpolate only the holes shorter than `h`; keep the cells themselves.
                let mut qs: Vec<u32> = o.path.iter().map(|x| x.0).collect();
                qs.sort_unstable();
                qs.dedup();
                for (i, &q) in qs.iter().enumerate() {
                    push(q);
                    if let Some(&nx) = qs.get(i + 1)
                        && nx - q - 1 < h
                    {
                        (q + 1..nx).for_each(&mut push);
                    }
                }
            }
        }
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
    Ok(st)
}

/// [`straighten`] in the overlap's canonical pair frame (`--symmetric`; `frame` =
/// every read's length and tie rank): when the owner `o.a` is not the query
/// (`readorder::query_first`: the longer read), the overlap is transposed, rewritten, and transposed
/// back, so the fill runs along the longer read's bins whichever read owns the pair. The rewrite does
/// not touch the slope, which is restored exactly (no `1 / (1 / s)`). `None` = [`straighten`].
pub fn straighten_in_frame(
    o: &mut Overlap,
    cfg: &Cfg,
    census: Option<&mut std::io::BufWriter<std::fs::File>>,
    names: &[String],
    frame: Option<(&[u32], &[u32])>,
) -> Result<Stats> {
    use haricot::readorder::{query_first, transpose_overlap};
    match frame {
        Some((lens, rank))
            if !query_first(
                lens[o.a as usize],
                lens[o.b as usize],
                rank[o.a as usize],
                rank[o.b as usize],
            ) =>
        {
            let slope = o.slope;
            let mut t = transpose_overlap(o);
            let st = straighten(&mut t, cfg, census, names)?;
            *o = transpose_overlap(&t);
            o.slope = slope;
            Ok(st)
        }
        _ => straighten(o, cfg, census, names),
    }
}

/// The fill of `haricot assemble`: Theil-Sen, no shift cap,
/// every hole of the query span filled, no trimming.
pub const DRIVER: Cfg<'static> = Cfg {
    mode: "theil",
    max_shift: i64::MAX,
    fill: true,
    split_hole: None,
    fill_max_hole: None,
};

/// Every overlap rewritten from its own path alone ([`straighten_in_frame`], no census), in
/// parallel; the stats summed.
pub fn fill_all(
    ovl: &mut [Overlap],
    cfg: &Cfg,
    names: &[String],
    frame: Option<(&[u32], &[u32])>,
) -> Result<Stats> {
    Ok(ovl
        .par_iter_mut()
        .map(|o| straighten_in_frame(o, cfg, None, names, frame))
        .collect::<Result<Vec<Stats>>>()?
        .into_iter()
        .fold(Stats::default(), Stats::add))
}

/// The report line of a run over `n` overlaps.
pub fn report(n: usize, st: &Stats, cfg: &Cfg) -> String {
    let Stats {
        fitted,
        moved,
        cells,
        shifted,
        capped,
        filled,
        trimmed,
        trim_cells,
    } = *st;
    format!(
        "{n} overlaps, {fitted} fitted ({}), {moved} moved; {shifted}/{cells} cells shifted \
         ({:.1}%), {capped} capped at --max-shift, {filled} cells filled; {trimmed} trimmed at \
         holes >= {} ({trim_cells} cells dropped)",
        cfg.mode,
        100.0 * shifted as f64 / cells.max(1) as f64,
        cfg.split_hole.map_or("-".to_string(), |h| h.to_string()),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 30 cells on b = a + 5, and 16 (35%) on a steep overrun tail: past Theil-Sen's ~29%
    /// breakdown, inside the repeated median's 50%.
    fn tailed() -> Vec<(u32, u32)> {
        let mut p: Vec<(u32, u32)> = (0..30).map(|a| (a, a + 5)).collect();
        p.extend((30..46).map(|a| (a, 3 * a - 55)));
        p
    }

    #[test]
    fn repeated_median_holds_the_line_where_theil_sen_tilts() {
        let p = tailed();
        let (mr, cr) = repeated_median(&p).unwrap();
        assert!(
            (mr - 1.0).abs() < 1e-9 && (cr - 5.0).abs() < 1e-9,
            "rm {mr} {cr}"
        );
        assert_eq!(off_line(&p, (mr, cr)), 15, "the tail minus its first cell");
        let (mt, _) = theil(&p).unwrap();
        assert!(
            mt > 1.05,
            "Theil-Sen should be pulled by the tail, got {mt}"
        );
    }

    /// `--symmetric`: the same pair stored by either owner is filled along the same (longer) read,
    /// so both rewrites are the same alignment seen from the two owners.
    #[test]
    fn the_fill_does_not_depend_on_the_owner() {
        use haricot::readorder::{in_canonical_frame, transpose_overlap};
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
            strand: haricot::params::Strand::Forward,
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
        let names = ["x".to_string(), "y".to_string()];
        let cfg = Cfg {
            mode: "theil",
            max_shift: i64::MAX,
            fill: true,
            split_hole: None,
            fill_max_hole: None,
        };
        let (mut r1, mut r2) = (o1.clone(), o2.clone());
        straighten_in_frame(&mut r1, &cfg, None, &names, Some((lens1, rank1))).unwrap();
        straighten_in_frame(&mut r2, &cfg, None, &names, Some((lens2, rank2))).unwrap();
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
        // without the frame, numbering 2's owner is filled along the shorter read: another path
        let mut r3 = o2.clone();
        straighten_in_frame(&mut r3, &cfg, None, &names, None).unwrap();
        assert_ne!(in_canonical_frame(&r3, lens2, rank2).path, c1.path);
    }

    #[test]
    fn both_fits_agree_on_a_clean_line() {
        let p: Vec<(u32, u32)> = (0..20).map(|a| (a, 2 * a + 1)).collect();
        assert_eq!(theil(&p), repeated_median(&p));
        assert_eq!(repeated_median(&p[..2]), None);
    }
}
