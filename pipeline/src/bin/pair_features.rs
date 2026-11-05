//! Per-candidate features for a learned overlap decision, with EXACT labels.
//!
//! Candidates are the detection's permissive cache (`haricot-detect --production`): it holds most
//! true overlaps, so the remaining problem is one decision over it, and `z` alone cannot make it
//! because its scale moves with read length. Four feature families, every one expressed as a ratio
//! or fraction wherever a raw count would carry the training genome's depth or read length:
//!
//! * **geometry** -- `z`, anchors, spans, slope, overhangs, containment, lengths.
//! * **path shape** -- hole fraction and longest hole run, Theil-Sen residuals, anchor density
//!   and its evenness along the path, re-collected from the pair's own seeds.
//! * **read context / repeat** -- what one pair cannot see: the candidate DEPTH over the
//!   overlap against each read's own typical depth (a repeat copy is visited by the other
//!   copies' reads too), the global occurrence of the pair's shared seeds against the
//!   dataset median, and whether third reads place A and B consistently.
//! * **signal-level DTW** -- banded DTW on the CONTINUOUS event levels (not the quantized
//!   alphabet) along the pair's fitted line, window by window, against a decoy that uses the
//!   same windows of A and a far offset in B. A paralog is a good alignment at the e-mer level
//!   but differs at the base level; the continuous levels are the only place that shows.
//!
//! Labels come from the signal-exact oracle (`haricot-truth-oracle`): `cls` is `true` (in the strict
//! oracle, path within 1 bin for >= 80% of cells), `misplaced` (in the oracle, path off),
//! `short` (a real intersection the strict oracle does not hold: UNKNOWN, not false) or
//! `false` (no placed intersection).
//!
//! Usage:
//!   haricot-pair-features READS.pod5 CANDIDATES.bin [--oracle ORACLE.bin --truth truth.tsv] --out F.tsv
//!                [--key NAME] [--sample-true 0.25] [--band 64] [--win 200] [--max-win 48]
//!                [--seed-occ-rel R [--seed-occ-floor 2000]] [--only PAIRS.tsv] [--no-symmetric]
//!                [--no-dtw] [--save-occ OCC.bin] [--occ OCC.bin (with --only)]
//! An unknown flag, a flag without its value and a value that does not parse are errors
//! (`cliargs.rs`).
//!
//! `--save-occ OCC.bin` writes the seed-occurrence table over every read and its two dataset statistics
//! (`occtable.rs`). `--occ OCC.bin`, with `--only`, takes them from such a table (written over the same
//! pod5) and segments only the reads of the listed pairs: the rows are the ones the full pass writes,
//! since a row reads the two reads' levels and seeds, the occurrence counts and the candidate set alone.
//!
//! The canonical pair frame (the default; `--symmetric` is accepted): every value is computed with
//! the pair in its canonical frame (`haricot::readorder::query_first`: the longer read as `a`; equal
//! lengths, the smaller tie rank, the reads' content-hash order), so no feature depends on which read
//! owns the pair; `name_a` / `name_b` stay in the cache's orientation (owner first) and so do the
//! labels (a cache without content hashes is refused); the `a_lo` .. `line_c` columns are the
//! canonical frame's, and so are the signal-DTW windows (along the canonical `a`, into `b`'s band).
//! `--no-symmetric`: every value in the cache's orientation (the owner's frame).
//!
//! `--no-dtw`: the signal-DTW columns are not computed (`dtw_windows` 0, the others NaN), in either
//! frame, for forests trained without them.
//!
//! `--only PAIRS.tsv` (lines `a<TAB>b`, read indices in the cache's order, either order): write rows
//! for those candidates only. Everything else -- the seed occurrences, the read context and the
//! triangles -- is still computed over the whole candidate set, so a row is the one the full run
//! would write (`haricot assemble` scores the contig-end rescue's pool this way, its pairs added to
//! the candidates).
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use haricot::emer;
use haricot::overlap::{Overlap, SeedCounts};
use haricot::params::{Params, Strand};
use haricot::signal::{self, AdcRead};
use rayon::prelude::*;

#[path = "../cache.rs"]
mod cache;
#[path = "../cliargs.rs"]
mod cliargs;
#[path = "../occtable.rs"]
mod occtable;

/// The flags that take a value, and the switches (`cliargs.rs`).
const VALUED: &[&str] = &[
    "--oracle",
    "--truth",
    "--out",
    "--key",
    "--sample-true",
    "--band",
    "--win",
    "--max-win",
    "--seed-occ-rel",
    "--seed-occ-floor",
    "--only",
    "--occ",
    "--save-occ",
];
const SWITCHES: &[&str] = &["--symmetric", "--no-symmetric", "--no-dtw"];

/// One read in the pipeline's coordinate space.
struct Read {
    /// Continuous level of each RETAINED event: the mean of the raw events its run merged.
    level: Vec<f32>,
    /// The sketch's seeds, `(hash, pos)` sorted by hash: 8 bytes each (a seed hash is 32-bit,
    /// `emer::Seed`), so the order is the `u64` hashes' order.
    seeds: Vec<(u32, u32)>,
    /// Every level is finite (see [`banded_dtw_pair`]).
    finite: bool,
}

fn build_read(signal_pa: &[f32], params: &Params) -> Read {
    let ev = signal::signal_to_events(signal_pa, params);
    let ret = emer::retained_quantized(&ev, params);
    let mut level = Vec::with_capacity(ret.len());
    for (k, q) in ret.iter().enumerate() {
        let lo = q.event_pos as usize;
        let hi = ret.get(k + 1).map_or(ev.len(), |n| n.event_pos as usize);
        level.push(ev[lo..hi].iter().sum::<f32>() / (hi - lo) as f32);
    }
    let mut seeds: Vec<(u32, u32)> = emer::sketch_from_retained(&ret, 0, Strand::Forward, params)
        .into_iter()
        .map(|s| (s.hash() as u32, s.pos))
        .collect();
    seeds.sort_unstable();
    let finite = level.iter().all(|x| x.is_finite());
    Read {
        level,
        seeds,
        finite,
    }
}

/// Theil-Sen line `b = s * a + c` through event-space points.
fn theil(pts: &[(f64, f64)]) -> Option<(f64, f64)> {
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

/// Mean per-row cost of an open-begin/open-end banded DTW of `x` against `y`, where row `i`
/// may use columns `centre(i) - w ..= centre(i) + w` (clipped to `y`).
///
/// The same recurrence, in the same order of operations, as the straightforward form kept in the
/// tests (`banded_dtw_ref`): the band's columns outside `y` are cut once per row instead of tested
/// per cell, and the previous row is addressed by the band's shift. It runs twice per window, up to
/// 48 windows per candidate, and is most of the features step.
#[cfg(test)]
fn banded_dtw(x: &[f32], y: &[f32], centre: &dyn Fn(usize) -> i64, w: i64) -> Option<f32> {
    let width = (2 * w + 1) as usize;
    let inf = f32::INFINITY;
    let mut prev = vec![inf; width];
    let mut cur = vec![inf; width];
    let mut prev_lo: Option<i64> = None;
    let mut rows = 0usize;
    let ny = y.len() as i64;
    for (i, &xi) in x.iter().enumerate() {
        let lo = centre(i) - w;
        // columns of the band inside `y`: k in [k0, k1)
        let k0 = (-lo).clamp(0, width as i64) as usize;
        let k1 = (ny - lo).clamp(0, width as i64) as usize;
        cur[..k0].fill(inf);
        let mut any = false;
        let mut left = inf; // cur[k - 1]
        match prev_lo {
            // first row: open begin, every cell of its band is a start
            None => {
                for k in k0..k1.max(k0) {
                    let d = (xi - y[(lo + k as i64) as usize]).abs();
                    cur[k] = d;
                    any |= d.is_finite();
                }
            }
            Some(plo) => {
                // prev index of column lo + k is k + sh
                let sh = lo - plo;
                let at = |kk: i64| -> f32 {
                    if kk < 0 || kk >= width as i64 {
                        inf
                    } else {
                        prev[kk as usize]
                    }
                };
                for k in k0..k1.max(k0) {
                    let d = (xi - y[(lo + k as i64) as usize]).abs();
                    let kk = k as i64 + sh;
                    let best = at(kk).min(at(kk - 1)).min(left);
                    let v = if best.is_finite() { d + best } else { inf };
                    cur[k] = v;
                    left = v;
                    any |= v.is_finite();
                }
            }
        }
        cur[k1.max(k0)..].fill(inf);
        if !any {
            return None;
        }
        std::mem::swap(&mut prev, &mut cur);
        prev_lo = Some(lo);
        rows += 1;
    }
    let end = prev.iter().copied().fold(inf, f32::min);
    (end.is_finite() && rows > 0).then(|| end / rows as f32)
}

/// Padding of `inf` cells on each side of a [`banded_dtw_pair`] row: a row whose band moved by
/// `sh` in `[1 - DTW_PAD, DTW_PAD]` reads its previous row without bounds tests.
const DTW_PAD: usize = 2;

/// One lane of [`banded_dtw_pair`]: a row buffer pair padded with `DTW_PAD` `inf` cells.
struct DtwLane {
    prev: Vec<f32>,
    cur: Vec<f32>,
    prev_lo: Option<i64>,
}

impl DtwLane {
    fn new(width: usize) -> Self {
        let n = width + 2 * DTW_PAD;
        DtwLane {
            prev: vec![f32::INFINITY; n],
            cur: vec![f32::INFINITY; n],
            prev_lo: None,
        }
    }

    /// Cell `k` of the current row (`y` its level, `left` the cell before it): [`banded_dtw`]'s
    /// recurrence, operation for operation. `sh` is the band's move since the previous row,
    /// in `[1 - DTW_PAD, DTW_PAD]`.
    /// With `FINITE` (both levels finite, so `d` is), `d + best` is already infinite exactly
    /// when `best` is, and the test drops out of the cell's dependency chain.
    #[inline(always)]
    fn cell<const FINITE: bool>(
        &mut self,
        xi: f32,
        y: f32,
        sh: i64,
        k: usize,
        left: &mut f32,
        any: &mut bool,
    ) {
        let d = (xi - y).abs();
        let kk = (k as i64 + sh) as usize + DTW_PAD;
        let best = self.prev[kk].min(self.prev[kk - 1]).min(*left);
        let v = if FINITE || best.is_finite() {
            d + best
        } else {
            f32::INFINITY
        };
        self.cur[k + DTW_PAD] = v;
        *left = v;
        *any |= v.is_finite();
    }

    /// One whole row by [`banded_dtw`]'s general form (any band move, the first row).
    fn row_general(&mut self, xi: f32, y: &[f32], lo: i64, width: usize) -> bool {
        let inf = f32::INFINITY;
        let ny = y.len() as i64;
        let k0 = (-lo).clamp(0, width as i64) as usize;
        let k1 = ((ny - lo).clamp(0, width as i64) as usize).max(k0);
        let (prev, cur) = (&self.prev, &mut self.cur);
        cur[DTW_PAD..DTW_PAD + k0].fill(inf);
        let mut any = false;
        let mut left = inf;
        match self.prev_lo {
            None => {
                for k in k0..k1 {
                    let d = (xi - y[(lo + k as i64) as usize]).abs();
                    cur[k + DTW_PAD] = d;
                    any |= d.is_finite();
                }
            }
            Some(plo) => {
                let sh = lo - plo;
                let at = |kk: i64| -> f32 {
                    if kk < 0 || kk >= width as i64 {
                        inf
                    } else {
                        prev[kk as usize + DTW_PAD]
                    }
                };
                for k in k0..k1 {
                    let d = (xi - y[(lo + k as i64) as usize]).abs();
                    let kk = k as i64 + sh;
                    let best = at(kk).min(at(kk - 1)).min(left);
                    let v = if best.is_finite() { d + best } else { inf };
                    cur[k + DTW_PAD] = v;
                    left = v;
                    any |= v.is_finite();
                }
            }
        }
        cur[DTW_PAD + k1..DTW_PAD + width].fill(inf);
        any
    }

    fn end_row(&mut self, lo: i64) {
        std::mem::swap(&mut self.prev, &mut self.cur);
        self.prev_lo = Some(lo);
    }

    fn end(&self, width: usize, rows: usize) -> Option<f32> {
        let end = self.prev[DTW_PAD..DTW_PAD + width]
            .iter()
            .copied()
            .fold(f32::INFINITY, f32::min);
        (end.is_finite() && rows > 0).then(|| end / rows as f32)
    }
}

/// [`banded_dtw`] of `x` against `y` along two centre lines at once (a window's line and its
/// decoy), `Some` iff both are: the same values bit for bit. Each lane's cells follow one
/// dependency chain along the row; running two independent chains in one loop lets the core
/// overlap them. A row whose band stays inside `y` and moved by at most `DTW_PAD` in both lanes
/// takes the joint loop; any other row, each lane's general form. The first lane to find no
/// finite cell in a row ends both, since the caller uses a window only when both align.
fn banded_dtw_pair(
    x: &[f32],
    y: &[f32],
    ca: &dyn Fn(usize) -> i64,
    cb: &dyn Fn(usize) -> i64,
    w: i64,
    finite: bool,
) -> Option<(f32, f32)> {
    if finite {
        dtw_pair::<true>(x, y, ca, cb, w)
    } else {
        dtw_pair::<false>(x, y, ca, cb, w)
    }
}

/// [`banded_dtw_pair`]; `FINITE`: every value of `x` and `y` is finite.
fn dtw_pair<const FINITE: bool>(
    x: &[f32],
    y: &[f32],
    ca: &dyn Fn(usize) -> i64,
    cb: &dyn Fn(usize) -> i64,
    w: i64,
) -> Option<(f32, f32)> {
    let width = (2 * w + 1) as usize;
    let ny = y.len() as i64;
    let (mut la, mut lb) = (DtwLane::new(width), DtwLane::new(width));
    let inside = |lo: i64| lo >= 0 && lo + width as i64 <= ny;
    let near = |lo: i64, plo: Option<i64>| {
        plo.is_some_and(|p| (1 - DTW_PAD as i64..=DTW_PAD as i64).contains(&(lo - p)))
    };
    let mut rows = 0usize;
    for (i, &xi) in x.iter().enumerate() {
        let (loa, lob) = (ca(i) - w, cb(i) - w);
        let (ok_a, ok_b) =
            if inside(loa) && inside(lob) && near(loa, la.prev_lo) && near(lob, lb.prev_lo) {
                let (sha, shb) = (loa - la.prev_lo.unwrap(), lob - lb.prev_lo.unwrap());
                let (ya, yb) = (
                    &y[loa as usize..loa as usize + width],
                    &y[lob as usize..lob as usize + width],
                );
                let (mut left_a, mut left_b) = (f32::INFINITY, f32::INFINITY);
                let (mut any_a, mut any_b) = (false, false);
                for (k, (&ya, &yb)) in ya.iter().zip(yb).enumerate() {
                    la.cell::<FINITE>(xi, ya, sha, k, &mut left_a, &mut any_a);
                    lb.cell::<FINITE>(xi, yb, shb, k, &mut left_b, &mut any_b);
                }
                (any_a, any_b)
            } else {
                (
                    la.row_general(xi, y, loa, width),
                    lb.row_general(xi, y, lob, width),
                )
            };
        if !ok_a || !ok_b {
            return None;
        }
        la.end_row(loa);
        lb.end_row(lob);
        rows += 1;
    }
    la.end(width, rows).zip(lb.end(width, rows))
}

/// The signal-DTW windows of one candidate: per window start `s0` in `starts`, [`banded_dtw_pair`]
/// of A's `win` levels from `s0` (`xa[s0..s0 + win]`) against `y` (B's levels), along the line
/// `b = s a + c` (`line`) and along its decoy (the line moved by `shift`, moved back by `2 shift`
/// past `y`'s end), `band` cells either side. With finite levels on an AVX2 CPU the windows run
/// through the batched kernel [`dtw16`], which gives the same values bit for bit; the windows it
/// does not take, and every window elsewhere, run
/// [`banded_dtw_pair`].
#[allow(clippy::too_many_arguments)]
fn dtw_windows(
    xa: &[f32],
    y: &[f32],
    starts: &[usize],
    win: usize,
    line: (f64, f64),
    shift: i64,
    band: i64,
    finite: bool,
) -> Vec<Option<(f32, f32)>> {
    let (s, c) = line;
    let nb = y.len() as i64;
    let cen = |s0: usize, i: usize| (s * (s0 + i) as f64 + c).round() as i64;
    let dec = |s0: usize, i: usize| {
        let j = cen(s0, i) + shift;
        if j >= nb { j - 2 * shift } else { j }
    };
    let mut out: Vec<Option<(f32, f32)>> = vec![None; starts.len()];
    #[allow(unused_mut)]
    let mut todo: Vec<usize> = (0..starts.len()).collect();
    #[cfg(target_arch = "x86_64")]
    if finite && dtw16::available() {
        todo = dtw16::windows(xa, y, starts, win, band, &cen, &dec, &mut out);
    }
    for w in todo {
        let s0 = starts[w];
        out[w] = banded_dtw_pair(
            &xa[s0..s0 + win],
            y,
            &|i| cen(s0, i),
            &|i| dec(s0, i),
            band,
            finite,
        );
    }
    out
}

/// [`banded_dtw_pair`] for finite levels, batched for AVX2: the line and decoy problems of up to 8
/// windows of one candidate in the 16 lanes of two `__m256` (lanes 0-7 the lines, 8-15 the
/// decoys), every lane running the scalar form's recurrence on its own problem, operation for
/// operation: `d = |x - y|`, `best = min(min(prev[kk], prev[kk - 1]), left)`, `v = d + best`, the
/// last row's minimum over the row count. The values are the scalar form's bit for bit: every cell
/// is `+0` or larger, never `NaN` or `-0`, so `min` returns the same bits whichever operand wins a
/// tie; `inf` stands for every cell outside `y` and for every cell no previous cell reaches, as in
/// the scalar form; the vector add, subtract and sign clear are the scalar ones lane by lane.
///
/// Lanes differ in where their bands sit. Row `i` of a lane reads the previous row's cells
/// `k + sh` and `k + sh - 1` (`sh` its band's move since the previous row) and its `y` strip from
/// `lo - base` (`base` the lane's lowest band start): each distinct move and strip offset of a row
/// is loaded once and blended in by lane. A move in `1 - PAD ..= PAD` reads the padded previous
/// row. A move past the whole band (`sh > width`, `sh <= -width`: the decoy's wrap) reaches no
/// previous cell, so that row and every later one is `inf` (the scalar form returns `None` there).
/// A window with any other move is left to [`banded_dtw_pair`].
#[cfg(target_arch = "x86_64")]
mod dtw16 {
    use std::arch::x86_64::*;
    use std::sync::OnceLock;

    /// `inf` cells either side of a DP row: the band moves the kernel reads directly.
    pub(super) const PAD: usize = 8;

    /// AVX2 on this CPU; decided once.
    pub(super) fn available() -> bool {
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::arch::is_x86_feature_detected!("avx2"))
    }

    /// A move that leaves the band: no cell of the new row has a previous cell.
    fn leaves(sh: i64, width: usize) -> bool {
        sh > width as i64 || sh <= -(width as i64)
    }

    /// The first row whose band leaves the previous one, `lo.len()` if none.
    fn first_leave(lo: &[i64], width: usize) -> usize {
        (1..lo.len())
            .find(|&i| leaves(lo[i] - lo[i - 1], width))
            .unwrap_or(lo.len())
    }

    /// Every move of the band starts `lo` is one the kernel reads.
    fn movable(lo: &[i64], width: usize) -> bool {
        lo.windows(2).all(|p| {
            let sh = p[1] - p[0];
            (1 - PAD as i64..=PAD as i64).contains(&sh) || leaves(sh, width)
        })
    }

    /// One lane's problem: A's levels for the window, the band start per row, and the first row
    /// whose band leaves the previous one (`rows` if none).
    struct Lane<'a> {
        x: &'a [f32],
        lo: &'a [i64],
        dies: usize,
    }

    /// [`super::dtw_windows`]'s batched part: fills `out` for the windows whose bands only move in
    /// ways the kernel reads, returns the others (in order) for the scalar form. `cen` / `dec`:
    /// the line's and the decoy's centre in `y` at row `i` of the window starting at `s0`.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn windows(
        xa: &[f32],
        y: &[f32],
        starts: &[usize],
        win: usize,
        band: i64,
        cen: &dyn Fn(usize, usize) -> i64,
        dec: &dyn Fn(usize, usize) -> i64,
        out: &mut [Option<(f32, f32)>],
    ) -> Vec<usize> {
        if band < 0 || win == 0 {
            return (0..starts.len()).collect();
        }
        let width = (2 * band + 1) as usize;
        let mut rest = Vec::new();
        let mut ok: Vec<(usize, Vec<i64>, Vec<i64>)> = Vec::new();
        for (w, &s0) in starts.iter().enumerate() {
            let lt: Vec<i64> = (0..win).map(|i| cen(s0, i) - band).collect();
            let ld: Vec<i64> = (0..win).map(|i| dec(s0, i) - band).collect();
            if movable(&lt, width) && movable(&ld, width) {
                ok.push((w, lt, ld));
            } else {
                rest.push(w);
            }
        }
        for chunk in ok.chunks(8) {
            let lanes: [Option<Lane>; 16] = std::array::from_fn(|l| {
                chunk.get(l % 8).map(|(w, lt, ld)| {
                    let lo = if l < 8 { lt.as_slice() } else { ld.as_slice() };
                    Lane {
                        x: &xa[starts[*w]..starts[*w] + win],
                        lo,
                        dies: first_leave(lo, width),
                    }
                })
            });
            // SAFETY: `available()` found AVX2 on this CPU (the only caller checks it).
            let ends = unsafe { batch(&lanes, win, width, y) };
            for (q, (w, _, _)) in chunk.iter().enumerate() {
                let (et, ed) = (ends[q], ends[8 + q]);
                out[*w] =
                    (et.is_finite() && ed.is_finite()).then(|| (et / win as f32, ed / win as f32));
            }
        }
        rest
    }

    /// The distinct values of one row's per-lane quantity over the alive lanes of one half, each
    /// with the mask of the lanes holding it (the first one is the blend's start; with no alive
    /// lane, the value 0: an in-range read whose lanes are all discarded).
    struct Sel {
        n: usize,
        e: [(i64, __m256); 8],
    }

    impl Sel {
        /// `vals[q]`: lane `8h + q`'s value, read where `alive` has the lane.
        #[target_feature(enable = "avx2")]
        #[allow(unused_unsafe)]
        fn of(alive: u16, h: usize, vals: &[i64; 8]) -> Sel {
            let mut v = [0i64; 8];
            let mut bits = [[0i32; 8]; 8];
            let mut n = 0usize;
            for (q, &x) in vals.iter().enumerate() {
                if alive & (1 << (8 * h + q)) == 0 {
                    continue;
                }
                let at = match v[..n].iter().position(|&u| u == x) {
                    Some(at) => at,
                    None => {
                        v[n] = x;
                        n += 1;
                        n - 1
                    }
                };
                bits[at][q] = -1;
            }
            // SAFETY: AVX2 is enabled here; each load reads one whole `[i32; 8]`.
            unsafe {
                let zero = _mm256_setzero_ps();
                let mut e = [(0i64, zero); 8];
                for ((slot, &x), b) in e.iter_mut().zip(&v).zip(&bits).take(n) {
                    *slot = (
                        x,
                        _mm256_castsi256_ps(_mm256_loadu_si256(b.as_ptr() as *const __m256i)),
                    );
                }
                Sel { n: n.max(1), e }
            }
        }
    }

    /// The 16 lanes over `rows` rows of `width` cells: each lane's last-row minimum (`inf` where
    /// its problem has no alignment). Unused lanes (`None`) are `inf` from row 0.
    #[target_feature(enable = "avx2")]
    #[allow(unused_unsafe, clippy::needless_range_loop)]
    fn batch(lanes: &[Option<Lane>; 16], rows: usize, width: usize, y: &[f32]) -> [f32; 16] {
        let nb = y.len() as i64;
        // Each used lane's `y` strip starts at its lowest band start over the rows it is alive.
        let mut base = [0i64; 16];
        let mut reach = 0usize;
        for (l, ln) in lanes.iter().enumerate() {
            if let Some(ln) = ln {
                let alive = &ln.lo[..ln.dies];
                let mn = alive.iter().copied().min().unwrap_or(0);
                let mx = alive.iter().copied().max().unwrap_or(0);
                base[l] = mn;
                reach = reach.max((mx - mn) as usize);
            }
        }
        let span = reach + width;
        let rowlen = width + 2 * PAD;
        // SAFETY: AVX2 is enabled here; every pointer read or written is a whole local array.
        unsafe {
            let inf = _mm256_set1_ps(f32::INFINITY);
            let sign = _mm256_set1_ps(-0.0);
            // yi[j * 2 + h], lane q: y[base + j] of lane 8h + q (inf outside y and for unused lanes)
            let mut yi: Vec<__m256> = Vec::with_capacity(span * 2);
            for j in 0..span {
                for h in 0..2 {
                    let v: [f32; 8] = std::array::from_fn(|q| {
                        let l = 8 * h + q;
                        let p = base[l] + j as i64;
                        if lanes[l].is_some() && (0..nb).contains(&p) {
                            y[p as usize]
                        } else {
                            f32::INFINITY
                        }
                    });
                    yi.push(_mm256_loadu_ps(v.as_ptr()));
                }
            }
            let mut prev: Vec<__m256> = vec![inf; rowlen * 2];
            let mut cur: Vec<__m256> = vec![inf; rowlen * 2];
            let mut dead = 0u16;
            for i in 0..rows {
                // lanes whose row is all `inf` from here: a band that leaves, an unused lane
                let mut dies = 0u16;
                for (l, ln) in lanes.iter().enumerate() {
                    let d = match ln {
                        Some(ln) => ln.dies == i,
                        None => i == 0,
                    };
                    dies |= (d as u16) << l;
                }
                let alive = !(dead | dies);
                let mut xv = [inf; 2];
                let mut km = [inf; 2];
                let mut ysel = [Sel::of(0, 0, &[0; 8]), Sel::of(0, 1, &[0; 8])];
                let mut psel = [Sel::of(0, 0, &[0; 8]), Sel::of(0, 1, &[0; 8])];
                for h in 0..2 {
                    let mut xs = [0f32; 8];
                    let mut kb = [0i32; 8];
                    let mut off = [0i64; 8];
                    let mut mv = [0i64; 8];
                    for q in 0..8 {
                        let l = 8 * h + q;
                        kb[q] = -(((dies >> l) & 1) as i32);
                        if let Some(ln) = &lanes[l] {
                            xs[q] = ln.x[i];
                            if alive & (1 << l) != 0 {
                                off[q] = ln.lo[i] - base[l];
                                if i > 0 {
                                    mv[q] = ln.lo[i] - ln.lo[i - 1];
                                }
                            }
                        }
                    }
                    xv[h] = _mm256_loadu_ps(xs.as_ptr());
                    km[h] = _mm256_castsi256_ps(_mm256_loadu_si256(kb.as_ptr() as *const __m256i));
                    ysel[h] = Sel::of(alive, h, &off);
                    psel[h] = Sel::of(alive, h, &mv);
                }
                let dies_h = [(dies & 0xff) != 0, (dies >> 8) != 0];
                if i == 0 {
                    for k in 0..width {
                        for h in 0..2 {
                            let s = &ysel[h];
                            let mut yv = yi[(s.e[0].0 as usize + k) * 2 + h];
                            for &(o, m) in &s.e[1..s.n] {
                                yv = _mm256_blendv_ps(yv, yi[(o as usize + k) * 2 + h], m);
                            }
                            let mut d = _mm256_andnot_ps(sign, _mm256_sub_ps(xv[h], yv));
                            if dies_h[h] {
                                d = _mm256_blendv_ps(d, inf, km[h]);
                            }
                            cur[(PAD + k) * 2 + h] = d;
                        }
                    }
                } else {
                    let mut left = [inf; 2];
                    for k in 0..width {
                        for h in 0..2 {
                            let s = &ysel[h];
                            let mut yv = yi[(s.e[0].0 as usize + k) * 2 + h];
                            for &(o, m) in &s.e[1..s.n] {
                                yv = _mm256_blendv_ps(yv, yi[(o as usize + k) * 2 + h], m);
                            }
                            let d = _mm256_andnot_ps(sign, _mm256_sub_ps(xv[h], yv));
                            let p = &psel[h];
                            let kk = (PAD as i64 + k as i64 + p.e[0].0) as usize;
                            let mut sel = _mm256_min_ps(prev[kk * 2 + h], prev[(kk - 1) * 2 + h]);
                            for &(sh, m) in &p.e[1..p.n] {
                                let kk = (PAD as i64 + k as i64 + sh) as usize;
                                let mm = _mm256_min_ps(prev[kk * 2 + h], prev[(kk - 1) * 2 + h]);
                                sel = _mm256_blendv_ps(sel, mm, m);
                            }
                            if dies_h[h] {
                                sel = _mm256_blendv_ps(sel, inf, km[h]);
                            }
                            let v = _mm256_add_ps(d, _mm256_min_ps(sel, left[h]));
                            cur[(PAD + k) * 2 + h] = v;
                            left[h] = v;
                        }
                    }
                }
                std::mem::swap(&mut prev, &mut cur);
                dead |= dies;
            }
            let mut e = [inf; 2];
            for k in 0..width {
                for (h, eh) in e.iter_mut().enumerate() {
                    *eh = _mm256_min_ps(*eh, prev[(PAD + k) * 2 + h]);
                }
            }
            let mut out = [0f32; 16];
            _mm256_storeu_ps(out.as_mut_ptr(), e[0]);
            _mm256_storeu_ps(out.as_mut_ptr().add(8), e[1]);
            out
        }
    }
}

fn quantile(v: &mut [f32], q: f64) -> f32 {
    if v.is_empty() {
        return f32::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() - 1) as f64 * q).round() as usize]
}

/// `-log10 P(X >= k)` for `X ~ Poisson(lam)`: how surprising a candidate depth `k` is at an
/// expected depth `lam`. A doubled depth is weak evidence of a repeat at 6x (5 -> 10) and
/// overwhelming at 30x (30 -> 60); a plain ratio calls both 2.0, this does not.
fn pois_tail(k: f64, lam: f64) -> f64 {
    let k = k.round().max(0.0) as u64;
    let lam = lam.max(0.5);
    if (k as f64) <= lam {
        return 0.0;
    }
    // Sum the upper tail from k, in log space, until the terms stop mattering.
    let ln_fact = |n: u64| -> f64 { (1..=n).map(|i| (i as f64).ln()).sum() };
    let mut lt = -lam + k as f64 * lam.ln() - ln_fact(k);
    let mut acc = lt;
    let mut i = k;
    loop {
        i += 1;
        lt += lam.ln() - (i as f64).ln();
        let hi = acc.max(lt);
        acc = hi + ((acc - hi).exp() + (lt - hi).exp()).ln();
        if lt < acc - 30.0 || i > k + 10_000 {
            break;
        }
    }
    (-acc / std::f64::consts::LN_10).min(300.0)
}

fn median_u32(v: &mut [u32]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_unstable();
    v[v.len() / 2] as f64
}

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let usage = "usage: haricot-pair-features READS.pod5 CAND.bin [--oracle O.bin --truth T.tsv] --out F.tsv";
    let args = cliargs::Cli::parse(&argv, VALUED, SWITCHES, usage)?;
    args.positional(2, usage)?;
    let (pod5, cand_p) = (&args.pos[0], &args.pos[1]);
    // The oracle and the truth only LABEL the rows (`cls`, `place_*`, `ov_bp`); no forest reads
    // them, so a production run (no reference) leaves them out and every row is `unlabeled`.
    let orc_p = args.opt("--oracle");
    let truth_p = args.opt("--truth");
    anyhow::ensure!(
        orc_p.is_some() == truth_p.is_some(),
        "--oracle and --truth go together"
    );
    let out_p = args.opt("--out").context("--out is required")?;
    let key = args.opt("--key").unwrap_or_else(|| "?".into());
    let sample_true: f64 = args.num("--sample-true", 1.0)?;
    let band: i64 = args.num("--band", 64)?;
    let win: usize = args.num("--win", 200)?;
    let max_win: usize = args.num("--max-win", 48)?;
    let params = Params::r10_4_1();
    params.validate().map_err(anyhow::Error::msg)?;
    // `--seed-occ-rel R [--seed-occ-floor F]`: the detection cap was max(F, R x mean occurrence)
    // (`OverlapParams::seed_occ_rel`, F = the detection's --max-seed-occ, 2000). The capped-anchor
    // features count anchors at occurrence <= the preset's max_seed_occ (1000, half the detection
    // cap), so they scale by the SAME factor, max(1, R x mean / F): unchanged wherever detection
    // kept its floor (every training set), scale-free above it. Off (0): 1000.
    let seed_occ_rel: f64 = args.num("--seed-occ-rel", 0.0)?;
    let seed_occ_floor: f64 = args.num("--seed-occ-floor", 2000.0)?;
    let only: Option<std::collections::HashSet<(u32, u32)>> = match args.opt("--only") {
        Some(p) => {
            let mut set = std::collections::HashSet::new();
            for line in std::fs::read_to_string(&p)
                .with_context(|| format!("reading {p}"))?
                .lines()
            {
                let f: Vec<&str> = line.split('\t').collect();
                if let (Some(x), Some(y)) = (f.first(), f.get(1))
                    && let (Ok(x), Ok(y)) = (x.parse::<u32>(), y.parse::<u32>())
                {
                    set.insert((x.min(y), x.max(y)));
                }
            }
            Some(set)
        }
        None => None,
    };

    let cache::Cache {
        counts,
        names,
        ovl: cand,
        content,
    } = cache::read(cand_p)?;
    // The canonical pair frame (the default; `--no-symmetric` for the owner's): every candidate's
    // features are computed in its canonical pair frame (`readorder::query_first`: the longer read as
    // `a`), so they do not depend on which read owns the pair; the labels and the row's names stay in
    // the cache's orientation.
    let symmetric = args.symmetric()?;
    let no_dtw = args.has("--no-dtw");
    // `--occ OCC.bin` (with `--only`): the occurrence table of a run over the same pod5
    // (`--save-occ`); only the reads of the `--only` pairs are segmented
    let occ_p = args.opt("--occ");
    anyhow::ensure!(
        occ_p.is_none() || only.is_some(),
        "--occ needs --only: without it every read is segmented anyway"
    );
    let cand_c: Vec<Overlap> = if symmetric {
        // the tie rank: the reads' signal content-hash order (a cache without hashes is refused)
        let rank = cache::tie_rank(content.as_deref())?;
        cand.iter()
            .map(|o| haricot::readorder::in_canonical_frame(o, &counts, &rank))
            .collect()
    } else {
        Vec::new()
    };
    // the overlaps the features read: canonical, or the cache's with --no-symmetric
    let fc: &[Overlap] = if symmetric { &cand_c } else { &cand };
    let oracle = match &orc_p {
        Some(p) => {
            let cache::Cache {
                names: onames,
                ovl: oracle,
                ..
            } = cache::read(p)?;
            if onames != names {
                bail!("oracle and candidate caches index reads differently");
            }
            oracle
        }
        None => Vec::new(),
    };
    let bin = cand.first().map(|o| o.bin_events).unwrap_or(192);
    eprintln!(
        "[{key}] {} candidates over {} reads, bin {bin}",
        cand.len(),
        names.len()
    );

    // Wall time per phase, reported on the last line: signal processing of every read, the
    // read context, the per-candidate loop (DTW included; `--max-win 1` is DTW nearly off).
    let t_start = Instant::now();

    // ---- reads, in cache order ---------------------------------------------------------
    let idx: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let mut reads: Vec<Option<Read>> = (0..names.len()).map(|_| None).collect();
    // with --occ, only the reads of the --only pairs (no row reads the others)
    let mut need: Option<std::collections::HashSet<usize>> = None;
    if let (Some(_), Some(set)) = (&occ_p, &only) {
        let ids = set.iter().flat_map(|&(x, y)| [x, y]);
        need = Some(ids.map(|x| x as usize).collect());
    }
    let wanted = |i: usize| need.as_ref().is_none_or(|n| n.contains(&i));
    let flush = |chunk: &mut Vec<AdcRead>, reads: &mut Vec<Option<Read>>| {
        let built: Vec<(usize, Read)> = chunk
            .par_iter()
            .filter_map(|r| {
                let &i = idx.get(r.id.as_str())?;
                wanted(i).then(|| (i, build_read(&r.signal_pa(), &params)))
            })
            .collect();
        for (i, r) in built {
            reads[i] = Some(r);
        }
        chunk.clear();
    };
    signal::for_each_adc_chunk(std::path::Path::new(pod5), None, 2 << 30, |mut chunk| {
        flush(&mut chunk, &mut reads);
        Ok(())
    })?;
    for (i, r) in reads.iter().enumerate() {
        match r {
            Some(r) if r.level.len() == counts[i] as usize => {}
            Some(r) => bail!(
                "{}: {} events vs cache {}",
                names[i],
                r.level.len(),
                counts[i]
            ),
            None if !wanted(i) => {}
            None => bail!("{} missing from the pod5", names[i]),
        }
    }
    let reads: Vec<Read> = reads
        .into_iter()
        .map(|r| {
            r.unwrap_or_else(|| Read {
                level: Vec::new(),
                seeds: Vec::new(),
                finite: true,
            })
        })
        .collect();
    let t_reads = t_start.elapsed().as_secs_f64();

    // Global seed occurrence (uncapped), and its dataset median over seed instances: counted over
    // every read, or (--occ) read from the table a full pass wrote.
    let (occ, occ_med, occ_mean) = match &occ_p {
        Some(p) => {
            let t = occtable::load(p)?;
            (t.counts, t.median, t.mean)
        }
        None => {
            let occ = SeedCounts::count(&reads, |r| r.seeds.len(), |r, i| r.seeds[i].0 as u64);
            let mut inst: Vec<u32> = reads
                .par_iter()
                .flat_map_iter(|r| r.seeds.iter().map(|&(h, _)| occ.get(h as u64)))
                .collect();
            let occ_med = median_u32(&mut inst).max(1.0);
            let occ_mean = inst.len() as f64 / occ.len().max(1) as f64;
            (occ, occ_med, occ_mean)
        }
    };
    if let Some(p) = args.opt("--save-occ") {
        occtable::save(&p, &occ, occ_med, occ_mean)?;
    }
    let occ_cap = (params.overlap.max_seed_occ as f64
        * (seed_occ_rel * occ_mean / seed_occ_floor).max(1.0))
    .round() as u32;
    eprintln!(
        "[{key}] mean occurrence {occ_mean:.1}; capped-anchor features at occurrence <= {occ_cap}"
    );

    // ---- read context from the whole candidate set --------------------------------------
    // Candidate depth per bin of every read, and each read's median non-zero depth.
    let mut depth: Vec<Vec<u32>> = counts
        .iter()
        .map(|&c| vec![0u32; (c / bin + 1) as usize])
        .collect();
    for o in &cand {
        for b in o.a_bin_span.0..=o.a_bin_span.1 {
            if let Some(x) = depth[o.a as usize].get_mut(b as usize) {
                *x += 1;
            }
        }
        for b in o.b_bin_span.0..=o.b_bin_span.1 {
            if let Some(x) = depth[o.b as usize].get_mut(b as usize) {
                *x += 1;
            }
        }
    }
    let depth_med: Vec<f64> = depth
        .iter()
        .map(|d| {
            let mut v: Vec<u32> = d.iter().copied().filter(|&x| x > 0).collect();
            median_u32(&mut v)
        })
        .collect();
    let mut all_med: Vec<u32> = depth_med
        .iter()
        .filter(|&&x| x > 0.0)
        .map(|&x| x as u32)
        .collect();
    let data_depth = median_u32(&mut all_med).max(1.0);

    // Per-pair fitted line in events, for the triangle test and the DTW (in the features' frame).
    let lines: Vec<Option<(f64, f64)>> = fc
        .par_iter()
        .map(|o| {
            let pts: Vec<(f64, f64)> = o
                .path
                .iter()
                .map(|&(a, b)| ((a * bin + bin / 2) as f64, (b * bin + bin / 2) as f64))
                .collect();
            theil(&pts)
        })
        .collect();
    // Adjacency: read -> (partner, candidate index).
    let mut adj: Vec<Vec<(u32, u32)>> = vec![Vec::new(); names.len()];
    for (k, o) in cand.iter().enumerate() {
        adj[o.a as usize].push((o.b, k as u32));
        adj[o.b as usize].push((o.a, k as u32));
    }
    for v in adj.iter_mut() {
        v.sort_unstable();
    }
    // Position of read X's event 0 in read Y's coordinates, from candidate k.
    let origin_in = |k: usize, x: u32| -> Option<f64> {
        let (s, c) = lines[k]?;
        let o = &fc[k];
        if o.a == x {
            Some(c) // b = s*a + c  -> a = 0 maps to c in b
        } else {
            (s.abs() > 1e-6).then(|| -c / s) // b = 0 in a's coordinates
        }
    };

    // ---- labels ----------------------------------------------------------------------------
    let orc: HashMap<(u32, u32), &Overlap> = oracle.iter().map(|o| ((o.a, o.b), o)).collect();
    let mut truth: HashMap<(u32, u32), i64> = HashMap::new();
    let truth_text = match &truth_p {
        Some(p) => std::fs::read_to_string(p)?,
        None => String::new(),
    };
    for line in truth_text.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if let (Some(&x), Some(&y)) = (idx.get(f[2]), idx.get(f[3])) {
            truth.insert((x.min(y) as u32, x.max(y) as u32), f[5].parse()?);
        }
    }

    let t_context = t_start.elapsed().as_secs_f64();
    // ---- per-candidate features --------------------------------------------------------------
    let seed_of = |k: usize| -> u64 {
        let o = &cand[k];
        (o.a as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ (o.b as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F)
    };
    let rows: Vec<Option<String>> = (0..cand.len())
        .into_par_iter()
        .map(|k| {
            let o = &cand[k];
            if let Some(set) = &only
                && !set.contains(&(o.a.min(o.b), o.a.max(o.b)))
            {
                return None;
            }
            // Label first, so the sample can skip work.
            let (cls, place_ok, place_med) = match orc.get(&(o.a, o.b)) {
                Some(x) => {
                    let m: HashMap<u32, u32> = x.path.iter().copied().collect();
                    let mut d: Vec<u32> = o
                        .path
                        .iter()
                        .filter_map(|&(ab, bb)| m.get(&ab).map(|&bo| bb.abs_diff(bo)))
                        .collect();
                    let ok = d.iter().filter(|&&e| e <= 1).count() as f64 / d.len().max(1) as f64;
                    let med = median_u32(&mut d);
                    (
                        if ok >= 0.8 && !d.is_empty() {
                            "true"
                        } else {
                            "misplaced"
                        },
                        ok,
                        med,
                    )
                }
                None if truth.contains_key(&(o.a, o.b)) => ("short", f64::NAN, f64::NAN),
                None if orc_p.is_none() => ("unlabeled", f64::NAN, f64::NAN),
                None => ("false", f64::NAN, f64::NAN),
            };
            if cls == "true" && sample_true < 1.0 {
                let r = (seed_of(k) >> 11) as f64 / (1u64 << 53) as f64;
                if r >= sample_true {
                    return None;
                }
            }
            // From here on the features' frame: the candidate in its canonical frame, or as the cache
            // stores it with --no-symmetric. `o0` (the cache's) names the row and keys the truth.
            let o0 = o;
            let o = &fc[k];
            let (a, b) = (o.a as usize, o.b as usize);
            let ra = &reads[a];
            let rb = &reads[b];
            let na_bins = (counts[a] / bin).max(1) as f64;
            let nb_bins = (counts[b] / bin).max(1) as f64;
            let (a0, a1) = o.a_bin_span;
            let (b0, b1) = o.b_bin_span;
            let span_a = (a1 - a0 + 1) as f64;
            let span_b = (b1 - b0 + 1) as f64;

            // geometry
            let left = a0.min(b0) as f64;
            let right = (na_bins - 1.0 - a1 as f64)
                .min(nb_bins - 1.0 - b1 as f64)
                .max(0.0);
            let oh = left + right;
            let shorter = na_bins.min(nb_bins);
            let a_cont = (a0 <= 1 && a1 as f64 >= na_bins - 2.0) as u8;
            let b_cont = (b0 <= 1 && b1 as f64 >= nb_bins - 2.0) as u8;

            // path shape
            let hole_frac = o.gap as f64 / span_a;
            let mut cells: Vec<u32> = o.path.iter().map(|c| c.0).collect();
            cells.sort_unstable();
            cells.dedup();
            let max_hole = cells.windows(2).map(|w| w[1] - w[0] - 1).max().unwrap_or(0) as f64;
            let (resid_med, resid_max, resid_frac1, slope) = match lines[k] {
                Some((s, c)) => {
                    let mut r: Vec<f32> = o
                        .path
                        .iter()
                        .map(|&(ab, bb)| {
                            let x = (ab * bin + bin / 2) as f64;
                            (((bb * bin + bin / 2) as f64 - (s * x + c)) / bin as f64).abs() as f32
                        })
                        .collect();
                    let frac1 =
                        r.iter().filter(|&&e| e > 1.0).count() as f64 / r.len().max(1) as f64;
                    let mx = r.iter().copied().fold(0.0f32, f32::max);
                    (quantile(&mut r, 0.5) as f64, mx as f64, frac1, s)
                }
                None => (f64::NAN, f64::NAN, f64::NAN, o.slope),
            };

            // seeds shared in the band around the line
            let (mut n_cap, mut n_all, mut n_high) = (0u32, 0u32, 0u32);
            let mut per_bin = vec![0u32; (a1 - a0 + 1) as usize];
            let mut occs: Vec<u32> = Vec::new();
            if let Some((s, c)) = lines[k] {
                let (mut i, mut j) = (0usize, 0usize);
                let (sa, sb) = (&ra.seeds, &rb.seeds);
                while i < sa.len() && j < sb.len() {
                    if sa[i].0 < sb[j].0 {
                        i += 1;
                    } else if sa[i].0 > sb[j].0 {
                        j += 1;
                    } else {
                        let h = sa[i].0;
                        let (ie, je) = (
                            i + sa[i..].iter().take_while(|x| x.0 == h).count(),
                            j + sb[j..].iter().take_while(|x| x.0 == h).count(),
                        );
                        let oc = occ.get(h as u64);
                        for x in &sa[i..ie] {
                            let ab = x.1 / bin;
                            if ab < a0 || ab > a1 {
                                continue;
                            }
                            let exp = s * x.1 as f64 + c;
                            for y in &sb[j..je] {
                                if (y.1 as f64 - exp).abs() <= 1.5 * bin as f64 {
                                    n_all += 1;
                                    occs.push(oc);
                                    if oc > occ_cap {
                                        n_high += 1;
                                    } else {
                                        n_cap += 1;
                                        per_bin[(ab - a0) as usize] += 1;
                                    }
                                }
                            }
                        }
                        i = ie;
                        j = je;
                    }
                }
            }
            let occ_rel_med = median_u32(&mut occs.clone()) / occ_med;
            let occ_rel_q90 = {
                occs.sort_unstable();
                occs.get(((occs.len().max(1) - 1) as f64 * 0.9) as usize)
                    .copied()
                    .unwrap_or(0) as f64
                    / occ_med
            };
            let bins_hit = per_bin.iter().filter(|&&x| x > 0).count() as f64 / per_bin.len() as f64;
            let dens_mean = n_cap as f64 / per_bin.len() as f64;
            let dens_cv = {
                let m = dens_mean;
                let v = per_bin.iter().map(|&x| (x as f64 - m).powi(2)).sum::<f64>()
                    / per_bin.len() as f64;
                if m > 0.0 { v.sqrt() / m } else { f64::NAN }
            };

            // read context: candidate depth over the overlap vs each read's own typical depth
            let mean_depth = |r: usize, lo: u32, hi: u32| -> f64 {
                let d = &depth[r];
                let hi = (hi as usize).min(d.len() - 1);
                let lo = (lo as usize).min(hi);
                d[lo..=hi].iter().map(|&x| x as f64).sum::<f64>() / (hi - lo + 1) as f64
            };
            let dep_a = mean_depth(a, a0, a1);
            let dep_b = mean_depth(b, b0, b1);
            let ctx_rel_a = dep_a / depth_med[a].max(1.0);
            let ctx_rel_b = dep_b / depth_med[b].max(1.0);
            let ctx_abs = dep_a.max(dep_b) / data_depth;
            let ctx_pois = pois_tail(dep_a, depth_med[a].max(1.0))
                .max(pois_tail(dep_b, depth_med[b].max(1.0)));
            let ctx_pois_abs = pois_tail(dep_a.max(dep_b), data_depth);
            let deg_rel = (adj[a].len() as f64 / na_bins).max(adj[b].len() as f64 / nb_bins)
                / (cand.len() as f64 * 2.0
                    / counts.iter().map(|&c| (c / bin).max(1) as f64).sum::<f64>());

            // triangle consistency through common partners
            let (mut n_common, mut n_cons) = (0u32, 0u32);
            if let Some(dab) = origin_in(k, o.a) {
                let (pa, pb) = (&adj[a], &adj[b]);
                let (mut i, mut j) = (0usize, 0usize);
                while i < pa.len() && j < pb.len() {
                    if pa[i].0 < pb[j].0 {
                        i += 1;
                    } else if pa[i].0 > pb[j].0 {
                        j += 1;
                    } else {
                        let cc = pa[i].0;
                        n_common += 1;
                        // With slopes ~1, off(A->B) = off(A->C) + off(C->B), where off(X->Y)
                        // is where X's origin lands in Y's coordinates.
                        let (ka, kb) = (pa[i].1 as usize, pb[j].1 as usize);
                        if let (Some(ac), Some(cb)) = (origin_in(ka, o.a), origin_in(kb, cc))
                            && (ac + cb - dab).abs() <= 3.0 * bin as f64
                        {
                            n_cons += 1;
                        }
                        i += 1;
                        j += 1;
                    }
                }
            }
            let tri_frac = if n_common > 0 {
                n_cons as f64 / n_common as f64
            } else {
                f64::NAN
            };
            let tri_rel = n_cons as f64 / (adj[a].len().min(adj[b].len()).max(1)) as f64;

            // signal-level DTW, window by window, true line vs a far decoy offset in B (not with
            // --no-dtw: no windows, the columns read NaN)
            let (mut w_true, mut w_ratio) = (Vec::<f32>::new(), Vec::<f32>::new());
            let mut dec_all = Vec::<f32>::new();
            if let Some((s, c)) = lines[k].filter(|_| !no_dtw) {
                let (lo, hi) = (
                    (a0 * bin) as usize,
                    (((a1 + 1) * bin) as usize).min(ra.level.len()),
                );
                let nwin = (hi.saturating_sub(lo)) / win;
                let step = (nwin as f64 / max_win as f64).max(1.0);
                let nb = rb.level.len() as i64;
                let shift = (nb / 2).max(4 * band);
                let mut starts = Vec::new();
                let mut t = 0.0f64;
                while (t as usize) < nwin {
                    starts.push(lo + (t as usize) * win);
                    t += step;
                }
                let finite = ra.finite && rb.finite;
                let line = (s, c);
                let r = dtw_windows(
                    &ra.level, &rb.level, &starts, win, line, shift, band, finite,
                );
                for (tc, dc) in r.into_iter().flatten() {
                    w_true.push(tc);
                    w_ratio.push(tc / dc.max(1e-6));
                    dec_all.push(dc);
                }
            }
            let nw = w_true.len();
            let dtw_mean = if nw > 0 {
                w_true.iter().sum::<f32>() / nw as f32
            } else {
                f32::NAN
            };
            let dec_mean = if nw > 0 {
                dec_all.iter().sum::<f32>() / nw as f32
            } else {
                f32::NAN
            };
            let good = w_ratio.iter().filter(|&&r| r < 0.7).count() as f32 / nw.max(1) as f32;
            let r_q10 = quantile(&mut w_ratio.clone(), 0.1);
            let r_q50 = quantile(&mut w_ratio.clone(), 0.5);
            let r_q90 = quantile(&mut w_ratio.clone(), 0.9);
            let r_max = w_ratio.iter().copied().fold(f32::NAN, f32::max);

            let ov_bp = truth.get(&(o0.a, o0.b)).copied().unwrap_or(0);
            let (line_s, line_c) = lines[k].unwrap_or((f64::NAN, f64::NAN));
            Some(format!(
                "{key}\t{}\t{}\t{cls}\t{place_ok:.3}\t{place_med}\t{ov_bp}\t\
                 {:.3}\t{}\t{}\t{}\t{:.4}\t{na_bins}\t{nb_bins}\t{span_a}\t{span_b}\t{:.4}\t\
                 {left}\t{right}\t{:.4}\t{a_cont}\t{b_cont}\t\
                 {hole_frac:.4}\t{max_hole}\t{resid_med:.3}\t{resid_max:.3}\t{resid_frac1:.4}\t\
                 {:.4}\t{n_cap}\t{n_all}\t{:.4}\t{:.4}\t{dens_cv:.4}\t{bins_hit:.4}\t\
                 {occ_rel_med:.3}\t{occ_rel_q90:.3}\t\
                 {ctx_rel_a:.3}\t{ctx_rel_b:.3}\t{ctx_abs:.3}\t{deg_rel:.3}\t{ctx_pois:.3}\t\
                 {ctx_pois_abs:.3}\t\
                 {n_common}\t{n_cons}\t{tri_frac:.4}\t{tri_rel:.4}\t\
                 {nw}\t{dtw_mean:.4}\t{dec_mean:.4}\t{:.4}\t{good:.4}\t{r_q10:.4}\t{r_q50:.4}\t\
                 {r_q90:.4}\t{r_max:.4}\t{a0}\t{a1}\t{b0}\t{b1}\t{line_s:.5}\t{line_c:.1}",
                names[o0.a as usize],
                names[o0.b as usize],
                o.z,
                o.n_anchors,
                o.aligned_bins,
                o.gap,
                slope,
                span_a / shorter,
                oh / span_a,
                (slope.ln()).abs(),
                dens_mean,
                n_high as f64 / n_all.max(1) as f64,
                if dec_mean > 0.0 {
                    dtw_mean / dec_mean
                } else {
                    f32::NAN
                },
            ))
        })
        .collect();

    let t_pairs = t_start.elapsed().as_secs_f64();
    let mut w = std::io::BufWriter::new(std::fs::File::create(&out_p)?);
    writeln!(
        w,
        "key\tname_a\tname_b\tcls\tplace_ok\tplace_med\tov_bp\t\
         z\tn_anchors\taligned_bins\tgap\tslope\ta_bins\tb_bins\tspan_a\tspan_b\tov_frac_short\t\
         oh_left\toh_right\toh_frac\ta_contained\tb_contained\t\
         hole_frac\tmax_hole\tresid_med\tresid_max\tresid_frac1\t\
         abs_log_slope\tanc_cap\tanc_all\tanc_per_bin\thigh_occ_frac\tdens_cv\tbins_hit\t\
         occ_rel_med\tocc_rel_q90\t\
         ctx_rel_a\tctx_rel_b\tctx_abs\tdeg_rel\tctx_pois\tctx_pois_abs\t\
         tri_common\ttri_cons\ttri_frac\ttri_rel\t\
         dtw_windows\tdtw_mean\tdtw_decoy\tdtw_ratio\tdtw_good_frac\tdtw_r_q10\tdtw_r_q50\t\
         dtw_r_q90\tdtw_r_max\ta_lo\ta_hi\tb_lo\tb_hi\tline_s\tline_c"
    )?;
    let mut n = 0usize;
    let mut by: HashMap<String, usize> = HashMap::new();
    for r in rows.into_iter().flatten() {
        let cls = r.split('\t').nth(3).unwrap().to_string();
        *by.entry(cls).or_insert(0) += 1;
        writeln!(w, "{r}")?;
        n += 1;
    }
    eprintln!(
        "[{key}] wrote {n} rows -> {out_p}; {:?}; occ median {occ_med}, depth median {data_depth}",
        by
    );
    let t_end = t_start.elapsed().as_secs_f64();
    eprintln!(
        "[{key}] PHASE reads {t_reads:.1} s, context {:.1} s, per-candidate {:.1} s, write {:.1} s",
        t_context - t_reads,
        t_pairs - t_context,
        t_end - t_pairs
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cli(l: &str) -> Result<cliargs::Cli> {
        let v: Vec<String> = l.split_whitespace().map(String::from).collect();
        cliargs::Cli::parse(&v, VALUED, SWITCHES, "usage")
    }

    /// The calls in the tree parse, with their frames: `haricot assemble`'s stage-1 and rescue
    /// features, `training/train.sh`'s, one in the owner's frame, a replay's.
    #[test]
    fn the_calls_in_the_tree_parse() {
        for (l, sym) in [
            (
                "P C --key k --out f.tsv --seed-occ-rel 17 --save-occ occ.bin",
                true,
            ),
            (
                "P C --key k --only pool.tsv --occ occ.bin --out f.tsv",
                true,
            ),
            ("P C --oracle o.bin --truth t.tsv --key k --out f.tsv", true),
            (
                "P C --oracle o.bin --truth t.tsv --key k --out f.tsv --no-symmetric",
                false,
            ),
            (
                "P C --key k --only pool.tsv --symmetric --occ occ.bin --out f.tsv",
                true,
            ),
        ] {
            let a = cli(l).unwrap();
            a.positional(2, "").unwrap();
            assert_eq!(a.symmetric().unwrap(), sym, "{l}");
            assert!(a.num("--seed-occ-rel", 0.0).is_ok() && a.num("--band", 64i64).is_ok());
        }
        assert!(
            cli("P C --out f --seed-occ-rel x17")
                .unwrap()
                .num("--seed-occ-rel", 0.0)
                .is_err()
        );
        assert!(cli("P C --out f --symetric").is_err());
        assert!(cli("P C D --out f").unwrap().positional(2, "").is_err());
    }

    /// The straightforward form of [`banded_dtw`], the reference its output must equal.
    fn banded_dtw_ref(x: &[f32], y: &[f32], centre: &dyn Fn(usize) -> i64, w: i64) -> Option<f32> {
        let width = (2 * w + 1) as usize;
        let inf = f32::INFINITY;
        let mut prev = vec![inf; width];
        let mut cur = vec![inf; width];
        let mut prev_lo = i64::MIN;
        let mut rows = 0usize;
        for (i, &xi) in x.iter().enumerate() {
            let lo = centre(i) - w;
            let mut any = false;
            for k in 0..width {
                let j = lo + k as i64;
                if j < 0 || j >= y.len() as i64 {
                    cur[k] = inf;
                    continue;
                }
                let d = (xi - y[j as usize]).abs();
                // (i-1, j), (i-1, j-1) from the previous row; (i, j-1) from this row.
                let from_prev = |jj: i64| -> f32 {
                    let kk = jj - prev_lo;
                    if prev_lo == i64::MIN || kk < 0 || kk >= width as i64 {
                        inf
                    } else {
                        prev[kk as usize]
                    }
                };
                let best =
                    from_prev(j)
                        .min(from_prev(j - 1))
                        .min(if k > 0 { cur[k - 1] } else { inf });
                // Open begin: the first row starts anywhere in its band.
                cur[k] = if i == 0 || !best.is_finite() {
                    if i == 0 { d } else { inf }
                } else {
                    d + best
                };
                any |= cur[k].is_finite();
            }
            if !any {
                return None;
            }
            std::mem::swap(&mut prev, &mut cur);
            prev_lo = lo;
            rows += 1;
        }
        let end = prev.iter().copied().fold(inf, f32::min);
        (end.is_finite() && rows > 0).then(|| end / rows as f32)
    }

    /// The two-lane DTW is the two single ones, `Some` iff both are, bit for bit: over bands
    /// that run off both ends of `y`, drift at slopes 0.7-1.4, jump, and (the decoy) wrap.
    #[test]
    fn banded_dtw_pair_matches_two_singles() {
        let mut h: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut rnd = move || {
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            (h >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut n = 0;
        for trial in 0..600 {
            let (nx, ny) = (1 + (rnd() * 250.0) as usize, 1 + (rnd() * 1500.0) as usize);
            let x: Vec<f32> = (0..nx).map(|_| (rnd() * 4.0 - 2.0) as f32).collect();
            let y: Vec<f32> = (0..ny).map(|_| (rnd() * 4.0 - 2.0) as f32).collect();
            let (s, c) = (0.7 + rnd() * 0.7, rnd() * ny as f64 * 1.4 - 0.2 * ny as f64);
            let w = 1 + (rnd() * 70.0) as i64;
            let jump = (rnd() * nx as f64) as usize;
            let cen =
                move |i: usize| (s * i as f64 + c).round() as i64 + if i > jump { 17 } else { 0 };
            let shift = (ny as i64 / 2).max(4 * w);
            let dec = move |i: usize| {
                let j = cen(i) + shift;
                if j >= ny as i64 { j - 2 * shift } else { j }
            };
            let pair = banded_dtw_pair(&x, &y, &cen, &dec, w, true);
            assert_eq!(
                pair.map(|(a, b)| (a.to_bits(), b.to_bits())),
                banded_dtw_pair(&x, &y, &cen, &dec, w, false)
                    .map(|(a, b)| (a.to_bits(), b.to_bits())),
                "trial {trial}: finite and general forms"
            );
            let singles = banded_dtw(&x, &y, &cen, w).zip(banded_dtw(&x, &y, &dec, w));
            assert_eq!(
                pair.map(|(a, b)| (a.to_bits(), b.to_bits())),
                singles.map(|(a, b)| (a.to_bits(), b.to_bits())),
                "trial {trial}"
            );
            n += pair.is_some() as usize;
        }
        assert!(n > 100, "only {n} trials produced an alignment pair");
    }

    /// The batched AVX2 kernel is the pair form, window by window, bit for bit: continuous and
    /// tied (quantised) levels, slopes 0.2-9.4 (moves of 0 to 10 cells per row; past `PAD` the
    /// window goes to the scalar form), bands from 1 to 201 cells running off both ends of `y`,
    /// decoys that wrap inside a window, 1 to 20 windows per candidate (partial batches), and the
    /// dispatcher over all of them.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn batched_dtw_is_the_pair_form_bit_for_bit() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            eprintln!("no AVX2 on this CPU: the batched DTW kernel is not exercised");
            return;
        }
        let mut h: u64 = 0x5851_f42d_4c95_7f2d;
        let mut rnd = move || {
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            (h >> 11) as f64 / (1u64 << 53) as f64
        };
        let bits = |r: Option<(f32, f32)>| r.map(|(a, b)| (a.to_bits(), b.to_bits()));
        let (mut batched, mut some, mut scalar) = (0usize, 0usize, 0usize);
        for trial in 0..400 {
            let tied = trial % 3 == 0;
            let level = |r: f64| {
                let v = (r * 4.0 - 2.0) as f32;
                if tied { (v * 4.0).round() / 4.0 } else { v }
            };
            let nxa = 1 + (rnd() * 3000.0) as usize;
            let ny = 1 + (rnd() * 3000.0) as usize;
            let xa: Vec<f32> = (0..nxa).map(|_| level(rnd())).collect();
            let y: Vec<f32> = (0..ny).map(|_| level(rnd())).collect();
            let s = match trial % 5 {
                0 => 1.0,
                1 => 0.2 + rnd() * 0.6,
                2 => 1.4 + rnd() * 8.0,
                _ => 0.7 + rnd() * 0.7,
            };
            let c = (rnd() * 1.4 - 0.2) * ny as f64 - s * nxa as f64 * rnd();
            let band = [0i64, 1, 5, 64, 64, 64, 100][trial % 7];
            let win = [1usize, 7, 200, 200, 37][(trial / 5) % 5].min(nxa);
            let nw = 1 + (rnd() * 20.0) as usize;
            let starts: Vec<usize> = (0..nw)
                .map(|_| (rnd() * (nxa - win + 1) as f64) as usize)
                .collect();
            let shift = (ny as i64 / 2).max(4 * band);
            let nb = ny as i64;
            let cen = |s0: usize, i: usize| (s * (s0 + i) as f64 + c).round() as i64;
            let dec = |s0: usize, i: usize| {
                let j = cen(s0, i) + shift;
                if j >= nb { j - 2 * shift } else { j }
            };
            let want: Vec<Option<(f32, f32)>> = starts
                .iter()
                .map(|&s0| {
                    let x = &xa[s0..s0 + win];
                    banded_dtw_pair(x, &y, &|i| cen(s0, i), &|i| dec(s0, i), band, true)
                })
                .collect();
            let mut got = vec![None; starts.len()];
            let rest = dtw16::windows(&xa, &y, &starts, win, band, &cen, &dec, &mut got);
            for (w, want) in want.iter().enumerate() {
                if rest.contains(&w) {
                    scalar += 1;
                    continue;
                }
                assert_eq!(
                    bits(got[w]),
                    bits(*want),
                    "trial {trial} window {w} (batched)"
                );
                batched += 1;
                some += want.is_some() as usize;
            }
            let all = dtw_windows(&xa, &y, &starts, win, (s, c), shift, band, true);
            for (w, want) in want.iter().enumerate() {
                assert_eq!(
                    bits(all[w]),
                    bits(*want),
                    "trial {trial} window {w} (dispatch)"
                );
            }
        }
        assert!(
            batched > 500 && some > 100 && scalar > 10,
            "coverage: {batched} batched windows ({some} aligned), {scalar} left to the scalar form"
        );
    }

    /// Bit-identical to the reference on random levels, over bands that run off both ends of
    /// `y`, drift, jump, and (the decoy) wrap.
    #[test]
    fn banded_dtw_matches_the_reference() {
        let mut h: u64 = 0x2545_f491_4f6c_dd1d;
        let mut rnd = move || {
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            (h >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut n = 0;
        for trial in 0..400 {
            let (nx, ny) = (1 + (rnd() * 250.0) as usize, 1 + (rnd() * 900.0) as usize);
            let x: Vec<f32> = (0..nx).map(|_| (rnd() * 4.0 - 2.0) as f32).collect();
            let y: Vec<f32> = (0..ny).map(|_| (rnd() * 4.0 - 2.0) as f32).collect();
            let (s, c) = (0.7 + rnd() * 0.7, rnd() * ny as f64 * 1.4 - 0.2 * ny as f64);
            let w = 1 + (rnd() * 70.0) as i64;
            let jump = (rnd() * nx as f64) as usize;
            let cen =
                move |i: usize| (s * i as f64 + c).round() as i64 + if i > jump { 17 } else { 0 };
            let (a, b) = (banded_dtw(&x, &y, &cen, w), banded_dtw_ref(&x, &y, &cen, w));
            assert_eq!(a.map(f32::to_bits), b.map(f32::to_bits), "trial {trial}");
            n += a.is_some() as usize;
        }
        assert!(n > 100, "only {n} trials produced an alignment");
    }
}
