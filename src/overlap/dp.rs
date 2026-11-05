//! Banded bin-DP — a faithful port of wtdbg2's KBM alignment.
//!
//! This is the **clean wtdbg2 baseline**: the recurrence, penalties, band, and
//! acceptance gate are ported directly from `extern/wtdbg2/kbm.h`, with no
//! raw-signal adaptation baked in. wtdbg2 finds seed matches at base-pair
//! resolution, floors each match's endpoints to a 256 bp **bin** to address a DP
//! cell `(⌊q/256⌋, ⌊t/256⌋)`, and runs a Smith-Waterman-shaped DP over the bin
//! grid. The cell score is `kmat + kcnt` (matched length + seed count); a banded
//! DP integrates the weak per-cell signal into a diagonal.
//!
//! Mapping into the e-mer domain (the only domain change — *not* an algorithm
//! change): the bin is `bin_events` events instead of 256 bp; an anchor is a
//! shared e-mer match `(a_event, b_event)`; a cell's `kcnt` is the number of
//! matches in it and `kmat` is their matched **event** length (`min(kcnt·e,
//! bin_events)`, each e-mer covering `e` events, saturated at the bin).
//!
//! The DP assumes a slope-1 main diagonal exactly as wtdbg2: the `max_bvar` band
//! and `max_bgap` gap bound budget *cumulative drift*, not a different slope. `bin_dp`
//! takes the slope as an argument and rescales the target coordinate by it before
//! binning; `1.0` is the pure wtdbg2 baseline. Off-slope diagonals are what the
//! Hough overlapper [`hough_aln`] finds (`haricot assemble` runs it).
//!
//! Citations: recurrence `_dp_cal_spare_row_kbm` (`kbm.h:1638`), backtrace +
//! acceptance gate `_backtrace_map_kbm` (`kbm.h:1765`, gate `kbm.h:1836`),
//! floor-to-bin (`kbm.h:2098`).

use std::cell::RefCell;
use std::collections::HashMap;

use rustc_hash::FxHashSet;

use super::Anchor;
use super::kernel::{self, densest_sorted};
use crate::params::{BeamCfg, OverlapParams};

/// Bin-DP configuration, in **bins** (and the e-mer width `e`). Faithful wtdbg2
/// defaults live in [`OverlapParams`]; [`from_params`](BinDp::from_params) builds
/// this from them.
#[derive(Debug, Clone, Copy)]
pub struct BinDp {
    /// Bin size in events (wtdbg2 `KBM_BSIZE` = 256 bp). The PAYLOAD resolution:
    /// it floors anchors into the bin `path` the FBG slicer consumes.
    pub bin_events: f64,
    /// Events per e-mer — used for the matched-event length `kmat`.
    pub e: u32,
    /// Per-gap-bin score `pgap` (wtdbg2 −7).
    pub pgap: i64,
    /// Per-bin diagonal-deviation penalty `pvar` (wtdbg2 −21).
    pub pvar: i64,
    /// Max consecutive gap bins on the diagonal `max_bgap` (wtdbg2 4).
    pub max_bgap: i64,
    /// DP band half-width `max_bvar` (wtdbg2 4).
    pub max_bvar: i64,
    /// Minimum aligned bins `min_aln` (wtdbg2 8 = 2048/256).
    pub min_aln_bins: i64,
    /// Minimum matched length `min_mat` (wtdbg2 200 bp; here in events).
    pub min_mat: i64,
    /// Minimum matched/aligned ratio `min_sim` (wtdbg2 0.05).
    pub min_sim: f64,
    /// Max gap-bin fraction of the alignment `max_gap` (wtdbg2 0.6).
    pub max_gap: f64,
    /// Max span-length variation `aln_var` (wtdbg2 0.25).
    pub aln_var: f64,
    /// The Hough slope search (see [`BeamCfg`]); the bin-DP path ignores it.
    pub slope_beam: BeamCfg,
    /// Significance floor `z` at the Hough early exit, from
    /// [`crate::params::OverlapParams::hough_min_z`]. `0` = off.
    pub min_z: f64,
    /// Compute the Hough confidence exactly even for pairs that provably fail the
    /// `min_conf` gate.
    ///
    /// `conf = peak / median(background)` with the median floored at 1, so
    /// `conf <= peak`. A pair whose winning in-window peak is already below `min_conf`
    /// is therefore rejected whatever the background says, and the 17-slope background
    /// scan — which exists only to compute that median — is dead work. `hough_aln`
    /// skips it for those pairs and reports `conf = peak`, an upper bound.
    ///
    /// With the gate OFF (`min_conf <= 0`, the permissive `z`-primary gate) no pair can
    /// fail it, so the background is skipped for EVERY pair, accepted ones included, and
    /// each reports `conf = peak`. On E. coli 6x that scan ran for all 586,986 pairs and
    /// decided none of them.
    ///
    /// This is invisible to the pipeline: `conf` gates and is never stored on an
    /// [`crate::overlap::Overlap`]. It *is* visible to any caller that reads or gates on
    /// the returned `conf` itself — set this there to keep its numbers exact.
    pub exact_conf: bool,
}

impl BinDp {
    /// Build from the overlap params + the e-mer width.
    ///
    /// The penalties are kept at wtdbg2's literal `pgap=−7`/`pvar=−21` rather than
    /// scaled to `bin_events`: relaxing them to help the slope-1 baseline integrate
    /// the weak per-e-mer-bin signal introduced cross-cluster false positives on the
    /// `overlap_clusters` fixture, and did not recover the debug pair (a slope-1 band
    /// cannot span its ≈1.28 diagonal, which is a geometry limit rather than a
    /// scoring one). The lever that did recover it is the Hough overlapper's slope
    /// search — see [`hough_aln`].
    pub fn from_params(p: &OverlapParams, e: u32) -> Self {
        BinDp {
            bin_events: p.bin_events as f64,
            e,
            pgap: p.pgap as i64,
            pvar: p.pvar as i64,
            max_bgap: p.max_bgap as i64,
            max_bvar: p.max_bvar as i64,
            min_aln_bins: p.min_aln_bins as i64,
            min_mat: p.min_matches as i64,
            min_sim: p.min_similarity as f64,
            max_gap: p.max_gap as f64,
            aln_var: p.aln_var as f64,
            slope_beam: p.slope_beam,
            min_z: p.hough_min_z as f64,
            // Production skips the background scan for provably-rejected pairs; the
            // diagnostics turn it back on. See the field docs.
            exact_conf: false,
        }
    }
}

/// Result of the bin-DP for one read pair. Carries the bin `path` so the FBG
/// (Phase 2) can slice the alignment into K-bins. An *empty* result
/// (`aligned_bins == 0`, empty `path`) means no overlap passed the gate.
/// Why a pair did **not** become an overlap: which gate ended it.
///
/// [`BinAln`] already carries `conf` on rejected results so a caller can tell
/// "measured, and it failed" from "never measured". This carries the same
/// information for the *decision*: which floor actually ended the pair. Without it
/// a caller has to infer the reason from which fields happen to be zero, which is
/// ambiguous.
///
/// Only [`hough_aln`] sets this; [`bin_dp`] leaves the default.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Reject {
    /// No diagonal was fitted, and no floor was reached: the default, and what
    /// `bin_dp` always reports.
    #[default]
    NotMeasured,
    /// Fewer than [`SLOPE_MIN_ANCHORS`] anchors in the cloud, so no line could be fitted.
    TooFewAnchors,
    /// `conf < min_conf` — the Hough significance gate.
    Conf,
    /// `z < hough_min_z` — the Hough *significance* gate (see
    /// [`crate::params::OverlapParams::hough_min_z`]).
    Signif,
    /// The level-0 peak fell under the first filter ([`BeamCfg::screen_peak`] or
    /// [`BeamCfg::screen_z`]).
    Screen,
    /// Fewer than [`SLOPE_MIN_ANCHORS`] anchors landed within `±tol` of the fitted line.
    OnLine,
    /// `aligned_bins < min_aln_bins` — wtdbg2's extent floor (`kbm.h` `min_aln`).
    Extent,
    /// `matched < min_mat` — wtdbg2's matched-length floor.
    Matched,
    /// `gap > max_gap · aligned_bins` — wtdbg2's gap-fraction floor (`kbm.h:1836`).
    ///
    /// Only [`hough_aln`] reports it: [`bin_dp`] ANDs the same condition into its own
    /// acceptance and has no separate reason for it.
    Gappy,
    /// The detector accepted it. Guards downstream in `align_clouds`
    /// (`max_similarity`, `max_anchor_density`, the null gate) may still reject.
    Accepted,
}

impl Reject {
    /// Stable lowercase tag, for diagnostics.
    pub fn tag(self) -> &'static str {
        match self {
            Reject::NotMeasured => "not_measured",
            Reject::TooFewAnchors => "few_anchors",
            Reject::Conf => "conf",
            Reject::Signif => "hough_z",
            Reject::Screen => "screen",
            Reject::OnLine => "online",
            Reject::Extent => "extent",
            Reject::Matched => "matched",
            Reject::Gappy => "gappy",
            Reject::Accepted => "accepted",
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct BinAln {
    /// Matched event length on the path (wtdbg2 `mat`).
    pub matched: u32,
    /// Seed-match count on the path (wtdbg2 `cnt`).
    pub cnt: u32,
    /// Aligned span in bins, `min(qspan, tspan)+1` (wtdbg2 `aln`).
    pub aligned_bins: u32,
    /// Gap bins on the path (wtdbg2 `gap`).
    pub gap: u32,
    /// `matched / (aligned_bins · bin_events)` (wtdbg2 `min_sim` numerator/denominator).
    pub similarity: f32,
    /// The slope the DP ran at (1.0 for the baseline).
    pub slope: f64,
    /// Offset of the fitted diagonal, `a_pos − slope·b_pos` ([`hough_aln`] only;
    /// 0.0 from [`bin_dp`], which is handed its slope and never fits a line).
    pub intercept: f64,
    /// Hough confidence (winning peak ÷ median peak over the scanned slopes) that
    /// [`hough_aln`] gated on. **Carried on rejected results too**, so a caller can
    /// tell "measured, and it failed" from "never measured" — a defaulted 0.0 is
    /// indistinguishable from a genuine zero. 0.0 from [`bin_dp`], which computes no
    /// confidence.
    pub conf: f64,
    /// Query (read A) bin span `(lo, hi)`.
    pub a_bin_span: (u32, u32),
    /// Target (read B) bin span `(lo, hi)` in **real** (un-rescaled) bin coords.
    pub b_bin_span: (u32, u32),
    /// `(a_bin, b_bin)` cells along the alignment, in real coords (occupied cells
    /// only). The per-cell correspondence the FBG K-bin slicer consumes.
    pub path: Vec<(u32, u32)>,
    /// Which floor ended this pair. Diagnostic only — nothing in the pipeline branches
    /// on it. See [`Reject`].
    pub reason: Reject,
    /// The null terms at the fitted slope, when [`hough_aln`] already computed them for
    /// its `min_z` floor. `None` otherwise (and always from [`bin_dp`]).
    pub null: Option<NullParts>,
}

impl BinAln {
    /// Whether this is a real (accepted) overlap.
    pub fn is_overlap(&self) -> bool {
        self.aligned_bins > 0 && !self.path.is_empty()
    }
}

/// One DP cell's running state (wtdbg2 `kbm_cell_t`).
#[derive(Debug, Clone, Copy)]
struct Cell {
    score: i64,
    mat: i64,
    cnt: i64,
    gap: i64,
    var: i64,
    live: bool,
}

const DEAD: Cell = Cell {
    score: 0,
    mat: 0,
    cnt: 0,
    gap: 0,
    var: 0,
    live: false,
};

/// Run the banded bin-DP over a read pair's anchors and return the best
/// alignment that passes wtdbg2's acceptance gate (or an empty [`BinAln`]).
///
/// `slope` rescales the target coordinate (1.0 = the wtdbg2 baseline, which is what
/// the pipeline passes).
pub fn bin_dp(anchors: &[Anchor], slope: f64, dp: &BinDp) -> BinAln {
    let bin = dp.bin_events.max(1.0);
    let bin_sz = bin as i64;

    // Floor each anchor to a cell. Row coord `trow` is slope-corrected; we also
    // remember a representative *real* b-bin per cell for the output path.
    let mut occ: HashMap<(u32, u32), (u32, u32)> = HashMap::new(); // (qcol,trow) -> (cnt, real_b)
    let (mut max_q, mut max_t) = (0u32, 0u32);
    for a in anchors {
        let qcol = (a.a_pos as f64 / bin) as u32;
        let trow = (a.b_pos as f64 * slope / bin) as u32;
        let real_b = (a.b_pos as f64 / bin) as u32;
        let cell = occ.entry((qcol, trow)).or_insert((0, real_b));
        cell.0 += 1;
        cell.1 = cell.1.min(real_b);
        max_q = max_q.max(qcol);
        max_t = max_t.max(trow);
    }
    let empty = || BinAln {
        slope,
        ..Default::default()
    };
    if occ.is_empty() {
        return empty();
    }
    let nq = max_q as usize + 1;
    let nt = max_t as usize + 1;

    // Rolling two-row DP (wtdbg2 keeps two rows; `kbm.h:1649`). `back` records
    // the chosen move per cell for backtrace: -1 dead, 0 path-start, 1 D, 2 H, 3 V.
    let mut prev = vec![DEAD; nq];
    let mut cur = vec![DEAD; nq];
    let mut back = vec![-1i8; nq * nt];
    let mut best: (i64, usize, usize) = (i64::MIN, 0, 0);

    for t in 0..nt {
        for c in cur.iter_mut() {
            *c = DEAD;
        }
        for q in 0..nq {
            let occ_cell = occ.get(&(q as u32, t as u32)).copied();
            let is_gap = occ_cell.is_none();
            let (kcnt, kmat) = match occ_cell {
                Some((cnt, _)) => {
                    let cnt = cnt as i64;
                    (cnt, (cnt * dp.e as i64).min(bin_sz))
                }
                None => (0, 0),
            };
            let score_cell = if is_gap { dp.pgap } else { kmat + kcnt };

            let d_pred = if q > 0 { prev[q - 1] } else { DEAD };
            let h_pred = if q > 0 { cur[q - 1] } else { DEAD };
            let v_pred = prev[q];

            // Diagonal (`bt=0`): resets var to 0, no var penalty (`kbm.h:1692`).
            let cand_d = Cell {
                score: d_pred.score + score_cell,
                mat: d_pred.mat + kmat,
                cnt: d_pred.cnt + kcnt,
                gap: d_pred.gap,
                var: 0,
                live: true,
            };
            // Horizontal (`bt=1`, query insertion): pvar charged as |var| grows,
            // refunded as it shrinks; var++ (`kbm.h:1678`).
            let h_pen = if h_pred.var < 0 { -dp.pvar } else { dp.pvar };
            let cand_h = Cell {
                score: h_pred.score + score_cell + h_pen,
                mat: h_pred.mat + kmat,
                cnt: h_pred.cnt + kcnt,
                gap: h_pred.gap,
                var: h_pred.var + 1,
                live: true,
            };
            // Vertical (`bt=2`, target insertion): symmetric; var-- (`kbm.h:1702`).
            let v_pen = if v_pred.var > 0 { -dp.pvar } else { dp.pvar };
            let cand_v = Cell {
                score: v_pred.score + score_cell + v_pen,
                mat: v_pred.mat + kmat,
                cnt: v_pred.cnt + kcnt,
                gap: v_pred.gap,
                var: v_pred.var - 1,
                live: true,
            };

            // Merge order (`kbm.h:1699`): start H; D wins ties (>=); V if strictly >.
            let (mut chosen, mut mv, mut pred_live) = (cand_h, 2i8, h_pred.live);
            if cand_d.score >= chosen.score {
                chosen = cand_d;
                mv = 1;
                pred_live = d_pred.live;
            }
            if cand_v.score > chosen.score {
                chosen = cand_v;
                mv = 3;
                pred_live = v_pred.live;
            }
            chosen.gap = if is_gap { chosen.gap + 1 } else { 0 };

            // Keep iff score>0, gap band, var band (`kbm.h:1734`).
            if chosen.score > 0 && chosen.gap <= dp.max_bgap && chosen.var.abs() <= dp.max_bvar {
                cur[q] = chosen;
                back[t * nq + q] = if pred_live { mv } else { 0 };
                // A path end is a non-gap cell (`kbm.h:1738` records mat>=min_mat;
                // we track the best-scoring end and gate after backtrace).
                if !is_gap && chosen.score > best.0 {
                    best = (chosen.score, q, t);
                }
            } else {
                cur[q] = DEAD;
                back[t * nq + q] = -1;
            }
        }
        std::mem::swap(&mut prev, &mut cur);
    }

    if best.0 == i64::MIN {
        return empty();
    }

    // Backtrace from the best end to its path start (`_backtrace_map_kbm`).
    let (mut q, mut t) = (best.1, best.2);
    let mut cells_on_path: Vec<(u32, u32, bool)> = Vec::new(); // (qcol, trow, occupied)
    loop {
        let occupied = occ.contains_key(&(q as u32, t as u32));
        cells_on_path.push((q as u32, t as u32, occupied));
        match back[t * nq + q] {
            1 if q > 0 && t > 0 => {
                q -= 1;
                t -= 1;
            }
            2 if q > 0 => q -= 1,
            3 if t > 0 => t -= 1,
            _ => break, // 0 (start), -1, or would underflow
        }
    }
    cells_on_path.reverse();

    // Recompute matched/cnt/gap and spans along the path. `trow` (corrected) feeds
    // the aligned-bins / length-consistency gate; `real_b` feeds the output.
    let (mut matched, mut cnt, mut gap) = (0i64, 0i64, 0u32);
    let (mut qlo, mut qhi, mut tlo, mut thi) = (u32::MAX, 0u32, u32::MAX, 0u32);
    let (mut blo, mut bhi) = (u32::MAX, 0u32);
    let mut path = Vec::new();
    for &(qc, tr, occupied) in &cells_on_path {
        qlo = qlo.min(qc);
        qhi = qhi.max(qc);
        tlo = tlo.min(tr);
        thi = thi.max(tr);
        if occupied {
            let (cc, real_b) = occ[&(qc, tr)];
            matched += (cc as i64 * dp.e as i64).min(bin_sz);
            cnt += cc as i64;
            blo = blo.min(real_b);
            bhi = bhi.max(real_b);
            path.push((qc, real_b));
        } else {
            gap += 1;
        }
    }
    let qspan = (qhi - qlo) as i64;
    let tspan = (thi - tlo) as i64;
    let aligned_bins = qspan.min(tspan) + 1;
    let aligned_events = aligned_bins as f64 * bin;
    let similarity = matched as f64 / aligned_events.max(1.0);

    // wtdbg2 acceptance gate (`kbm.h:1836`): all ANDed.
    let accept = aligned_bins >= dp.min_aln_bins
        && matched >= dp.min_mat
        && matched as f64 >= aligned_events * dp.min_sim
        && gap as f64 <= aligned_bins as f64 * dp.max_gap
        && (qspan - tspan).abs() as f64 <= (dp.aln_var * aligned_bins as f64).max(1.0);
    if !accept || path.is_empty() {
        return empty();
    }

    BinAln {
        matched: matched as u32,
        cnt: cnt as u32,
        aligned_bins: aligned_bins as u32,
        gap,
        similarity: similarity as f32,
        slope,
        // The banded DP is *handed* its slope and never fits a line, so it has
        // neither an intercept nor a Hough confidence to report.
        intercept: 0.0,
        conf: 0.0,
        a_bin_span: (qlo, qhi),
        b_bin_span: (blo, bhi),
        path,
        // The banded DP does not fit a line, so it has no `Reject` axis to report.
        reason: Reject::NotMeasured,
        null: None,
    }
}

/// Minimum anchors to attempt a slope estimate.
const SLOPE_MIN_ANCHORS: usize = 8;
/// A fitted overlap diagonal and its Hough confidence against the wrong-slope background.
#[derive(Debug, Clone, Copy)]
pub struct DiagFit {
    pub slope: f64,
    pub intercept: f64,
    /// `peak / median(peaks)` — the [`OverlapParams::hough_min_conf`] statistic.
    pub conf: f64,
    /// Stopped by the first filter ([`BeamCfg::screen_peak`]) after level 0: `slope` is
    /// the best coarse slope and `conf` its peak, and nothing was refined.
    pub screened: bool,
}

/// Find the densest diagonal `a_pos ≈ slope·b_pos + intercept` in an anchor cloud: the
/// slope at which the anchors are most **concentrated** onto one diagonal (the beam of
/// [`beam_search_slope`], scoring the densest offset window), then a least-squares refit of
/// slope *and* intercept over just that diagonal's anchors.
///
/// Concentration, not chain length, is the criterion on purpose: in a noisy cloud the
/// longest non-decreasing chain can follow a long, shallow ramp through off-diagonal
/// noise rather than the dense true diagonal (measured on the debug pair, a chain-length
/// estimator reported slope ≈ 2 / 0.3 where the concentration scan finds ≈ 1.28; see
/// [`hough_aln`]). `tol_events` is the diagonal half-width.
///
/// Returns `None` only when there are too few anchors. `conf_floor` is the caller's
/// `min_conf` gate, used only to skip work that cannot change the verdict (see
/// [`BinDp::exact_conf`]); pass `0.0` for the exact confidence. The slope is clamped to
/// `[0.2, 5.0]`.
pub fn estimate_diagonal_fit(
    anchors: &[Anchor],
    tol_events: u32,
    beam: &BeamCfg,
    conf_floor: f64,
) -> Option<DiagFit> {
    if anchors.len() < SLOPE_MIN_ANCHORS {
        return None;
    }
    Some(fit_cloud(
        &Cloud::new(anchors),
        tol_events,
        beam,
        conf_floor,
    ))
}

/// [`estimate_diagonal_fit`] of a cloud of at least [`SLOPE_MIN_ANCHORS`] anchors.
fn fit_cloud(cloud: &Cloud, tol_events: u32, beam: &BeamCfg, conf_floor: f64) -> DiagFit {
    let anchors = cloud.anchors;
    let tol = tol_events.max(1) as f64;
    // Locate the densest diagonal `best = (slope, peak, offset)` and score its
    // confidence (winning peak ÷ median peak) against the wrong-slope background. The
    // background is taken from a coarse scan over the FULL slope range, so the
    // `min_conf` gate is not inflated by a window-only (mostly near-diagonal) median.
    let (best, bg_peaks, screened, finalists) = beam_search_slope(cloud, tol, beam, conf_floor);
    if screened {
        return DiagFit {
            slope: best.0,
            intercept: best.2,
            conf: best.1 as f64,
            screened: true,
        };
    }
    let conf_of = |peak: usize| match &bg_peaks {
        Some(peaks) => peak_conf(peak, peaks.clone()),
        // Background skipped because the peak alone already fails the gate. `conf <=
        // peak`, so reporting the bound keeps the caller's `conf < min_conf` test giving
        // the same answer.
        None => peak as f64,
    };
    // `BeamCfg::select_z`: refit every finalist basin as the winner would be, and keep the
    // one the gate scores highest (ties keep the finalist order, i.e. the higher count).
    if finalists.len() > 1 {
        let mut pick: Option<(f64, f64, f64, usize)> = None; // (z, slope, intercept, peak)
        for f in &finalists {
            let (slope, intercept) = refit_line(anchors, f.0, f.2, tol);
            let z = significance_parts(cloud, slope, tol).z;
            if pick.is_none_or(|p| z > p.0) {
                pick = Some((z, slope, intercept, f.1));
            }
        }
        let (_, slope, intercept, peak) = pick.expect("finalists is non-empty");
        return DiagFit {
            slope,
            intercept,
            conf: conf_of(peak),
            screened: false,
        };
    }
    let (slope, intercept) = refit_line(anchors, best.0, best.2, tol);
    DiagFit {
        slope,
        intercept,
        conf: conf_of(best.1),
        screened: false,
    }
}

/// Refine slope + intercept with a least-squares fit over only the anchors on the diagonal
/// `a = m0·b + c0` (offset within ±tol of it), so neither is pulled by off-diagonal noise.
/// Too few anchors on it (or a degenerate fit) keep the line itself. The slope is clamped to
/// `[0.2, 5.0]`.
fn refit_line(anchors: &[Anchor], m0: f64, c0: f64, tol: f64) -> (f64, f64) {
    let (mut n, mut sb, mut sa, mut sbb, mut sab) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for an in anchors {
        let (x, y) = (an.b_pos as f64, an.a_pos as f64);
        if (y - m0 * x - c0).abs() <= tol {
            n += 1.0;
            sb += x;
            sa += y;
            sbb += x * x;
            sab += x * y;
        }
    }
    let den = n * sbb - sb * sb;
    if n < SLOPE_MIN_ANCHORS as f64 || den.abs() < 1e-9 {
        (m0.clamp(0.2, 5.0), c0)
    } else {
        let slope = ((n * sab - sb * sa) / den).clamp(0.2, 5.0);
        (slope, (sa - slope * sb) / n)
    }
}

/// DEVIATION from wtdbg2: a global, gap-agnostic Hough-extent overlapper, the
/// alternative to [`bin_dp`] selected by [`OverlapParams::hough`].
///
/// Where [`bin_dp`] integrates the diagonal with a *banded* Smith-Waterman that can
/// stop at the first gap budget violation, this finds the overlap diagonal with the
/// 1-D Hough scan ([`estimate_diagonal_fit`]) and takes all anchors within `±tol`
/// of that line as the alignment. Two consequences feed the FBG:
///   * it recovers the whole diagonal (no early stop at gaps) → denser, fuller bin
///     `path` → fewer FBG fragmentation breaks;
///   * it rejects diffuse repeat / cross-cluster anchor clouds on Hough
///     significance, which a banded DP can still walk a thin chain through.
///
/// It emits the bin `path` the FBG slicer consumes. The on-line anchors are
/// floored to `bin_events` cells exactly as `bin_dp` floors matches; a cell's
/// matched length is the wtdbg2 `min(count·e, bin)` (`kbm.h:1638`).
///
/// `tol`, the Hough line half-width in events, is `bin_events`: detection and the FBG
/// grid share one resolution. Accepts iff the Hough confidence ≥ `min_conf`, the null
/// `z` ≥ `dp.min_z`, and wtdbg2's extent / matched / gap floors pass; otherwise returns
/// an empty [`BinAln`] stamped with the [`Reject`] reason.
pub fn hough_aln(anchors: &[Anchor], dp: &BinDp, min_conf: f64) -> BinAln {
    let bin = dp.bin_events.max(1.0);
    let tol = bin;
    let empty = BinAln::default();

    // 1. Hough detection + significance gate (the specificity lever). The floor lets
    //    the search skip the confidence background where it cannot
    //    change the verdict; with the gate OFF (`min_conf <= 0`) that is every pair,
    //    since `conf` then decides nothing and is not stored (see `BinDp::exact_conf`).
    let conf_floor = if dp.exact_conf {
        0.0
    } else if min_conf <= 0.0 {
        f64::INFINITY
    } else {
        min_conf
    };
    if anchors.len() < SLOPE_MIN_ANCHORS {
        return BinAln {
            reason: Reject::TooFewAnchors,
            ..empty
        };
    }
    let cloud = Cloud::new(anchors);
    let fit = fit_cloud(&cloud, tol as u32, &dp.slope_beam, conf_floor);
    let (slope, intercept, conf) = (fit.slope, fit.intercept, fit.conf);
    if fit.screened {
        return BinAln {
            slope,
            intercept,
            conf,
            reason: Reject::Screen,
            ..empty
        };
    }
    if conf < min_conf {
        // Rejected, but report WHAT was measured: the fitted diagonal and the
        // confidence it failed on. `is_overlap()` is still false (no bins, no path),
        // so the pipeline is unaffected — this only stops the rejected result from
        // being indistinguishable from a defaulted one.
        return BinAln {
            slope,
            intercept,
            conf,
            reason: Reject::Conf,
            ..empty
        };
    }

    // 1b. Significance gate, on the SAME fit and at the SAME early exit as `conf`.
    //     `conf` is a ratio and is length-biased; `z` is a significance. On a
    //     D. melanogaster set `conf` rejected 86.8% of true overlaps at a median `z`
    //     of 15.8-18.8 -- see `OverlapParams::hough_min_z`.
    //     The terms are handed back on the accepted result, so the caller's null gate
    //     (the same statistic at the same slope) does not compute them a second time.
    let mut null = None;
    if dp.min_z > 0.0 {
        let parts = Some(significance_parts(&cloud, slope, tol));
        if parts.map_or(0.0, |p| p.z) < dp.min_z {
            return BinAln {
                slope,
                intercept,
                conf,
                reason: Reject::Signif,
                ..empty
            };
        }
        null = parts;
    }

    // 2. On-line anchors (within ±tol of the line). Tally, per query bin, the target
    //    bins they fall in. The ±tol band is ~3 bins wide, so a query bin can see
    //    several target bins; emitting every one of them gave the FBG a fat,
    //    phase-shifted tiling that tangled into spurious branches. So the band is
    //    kept only for *detection*, and each query bin is collapsed to its single
    //    best-supported target bin — a thin monotone path that still spans the whole
    //    diagonal without the fat-band tangling.
    let mut per_a: HashMap<u32, HashMap<u32, u32>> = HashMap::new();
    let mut online = 0usize;
    for an in anchors {
        if (an.a_pos as f64 - slope * an.b_pos as f64 - intercept).abs() <= tol {
            online += 1;
            let abn = (an.a_pos as f64 / bin) as u32;
            let bbn = (an.b_pos as f64 / bin) as u32;
            *per_a.entry(abn).or_default().entry(bbn).or_default() += 1;
        }
    }
    if online < SLOPE_MIN_ANCHORS || per_a.is_empty() {
        return BinAln {
            slope,
            intercept,
            conf,
            reason: Reject::OnLine,
            ..empty
        };
    }

    // 3. Collapse each query bin to its modal target bin (ties → smaller b for
    //    determinism), forming the thin path. Matched per query bin is the wtdbg2
    //    `min(count·e, bin)` over all on-line anchors in that bin (`kmat`, kbm.h:1638).
    let bin_sz = bin as i64;
    let (mut matched, mut cnt) = (0i64, 0i64);
    let (mut alo, mut ahi, mut blo, mut bhi) = (u32::MAX, 0u32, u32::MAX, 0u32);
    let mut path: Vec<(u32, u32)> = Vec::with_capacity(per_a.len());
    for (&abn, bs) in &per_a {
        let (&best_b, _) = bs
            .iter()
            .max_by(|x, y| x.1.cmp(y.1).then(y.0.cmp(x.0)))
            .unwrap();
        let c: i64 = bs.values().map(|&v| v as i64).sum();
        matched += (c * dp.e as i64).min(bin_sz);
        cnt += c;
        path.push((abn, best_b));
        alo = alo.min(abn);
        ahi = ahi.max(abn);
        blo = blo.min(best_b);
        bhi = bhi.max(best_b);
    }
    path.sort_unstable();
    let a_bins = ahi - alo + 1;
    let b_bins = bhi - blo + 1;
    let aligned_bins = a_bins.min(b_bins);
    let gap = aligned_bins.saturating_sub(path.len() as u32);

    // 4. wtdbg2's floors. `max_gap` was omitted here on the assumption that the Hough gate
    // did that work; it does not — 52.9% of the false overlaps at `z >= 12` exceed it
    // against 1.85% of the true. Applying wtdbg2's own value: refNG50 87,617 -> 100,161,
    // neutral on both bacteria. `min_sim` and
    // `aln_var` stay unapplied — they are weaker and would be new operating points.
    let too_gappy = gap as f64 > dp.max_gap * aligned_bins.max(1) as f64;

    if (aligned_bins as i64) < dp.min_aln_bins || matched < dp.min_mat || too_gappy {
        return BinAln {
            slope,
            intercept,
            conf,
            // Report WHICH floor: they fail on different pairs (a long thin diagonal
            // clears extent but not matched; a short dense one the reverse; a repeat
            // smear clears both and fails the gap).
            reason: if (aligned_bins as i64) < dp.min_aln_bins {
                Reject::Extent
            } else if matched < dp.min_mat {
                Reject::Matched
            } else {
                Reject::Gappy
            },
            ..empty
        };
    }
    let aligned_events = aligned_bins as f64 * bin;
    let similarity = (matched as f64 / aligned_events.max(1.0)) as f32;

    BinAln {
        matched: matched as u32,
        cnt: cnt as u32,
        aligned_bins,
        gap,
        similarity,
        slope,
        intercept,
        conf,
        a_bin_span: (alo, ahi),
        b_bin_span: (blo, bhi),
        path,
        reason: Reject::Accepted,
        null,
    }
}

/// Densest `win`-wide window over a SORTED sequence: `(count, lo, hi)` of the
/// window's extreme members — the peak [`diagonal_significance_parts`] scores, found the
/// way [`diagonal_peak`] finds it.
fn densest_window(sorted: &[f64], win: f64) -> (usize, f64, f64) {
    let (mut i, mut best, mut lo, mut hi) = (0usize, 0usize, sorted[0], sorted[0]);
    for j in 0..sorted.len() {
        while sorted[j] - sorted[i] > win {
            i += 1;
        }
        if j - i + 1 > best {
            best = j - i + 1;
            lo = sorted[i];
            hi = sorted[j];
        }
    }
    (best, lo, hi)
}

/// NOT a wtdbg2 port — a novel opt-in overlap-precision gate. The **null-hypothesis
/// significance** of the best diagonal in an anchor cloud, as a Poisson z-score.
///
/// # The test
/// A real overlap piles its shared-seed anchors onto one diagonal `a = m·b + c₀`;
/// a diffuse / repeat / cross-locus collision scatters them. We ask: is the
/// observed on-diagonal concentration more than a random-seed null would produce?
///
/// * **Observed.** `obs` = anchors in the densest `±tol` offset window at slope `m`
///   (the same window [`diagonal_peak`] / [`hough_aln`] call "on-line"; `c₀` is its
///   center).
/// * **Expected under the null.** The design baseline is the analytic Poisson form
///   `E = N · band_area / box_area` — with `N` anchors uniform in the
///   `a_range × b_range` bounding box, a `±tol` band around one diagonal covers a
///   `2·tol / a_range` fraction of every column, so `E = N · 2·tol / a_range`.
///   That **uniform-box** assumption is exactly what low e-mer entropy breaks:
///   a handful of frequent hashes (entropy ≈2.53/3 bits here) dump a
///   Cartesian block of anchors onto a few repeated a- and b-positions, so the
///   cloud is far from uniform and the box null under-predicts the on-band count,
///   inflating significance.
///
/// # Hash-frequency-aware null (what this actually computes)
/// We replace the uniform box with the **observed marginals**: draw an a-position
/// and a b-position *independently* from the anchors' own position histograms —
/// random pairing that preserves both marginals but destroys the joint diagonal.
/// The expected on-band count is then
///
/// ```text
///   E = (1/N) · #{ (i,j) : | a_i − m·b_j − c₀ | ≤ tol }
/// ```
///
/// (each of the `N²` random pairings has weight `1/N`; the `N` true pairings are a
/// negligible `1/N` of them). Because a repeated position enters the histogram with
/// its full multiplicity, a repeat block *raises* `E` — random pairing keeps
/// reconstructing near-diagonal offsets — so its spurious concentration is
/// discounted. This is the same intent as weighting anchors by `1/occ`, achieved
/// from the marginals alone (no hash plumbing needed), and it **reduces to the
/// box null** `E ≈ N·2·tol/a_range` when the marginals are spread. Computed in
/// `O(N log N)` by binary-searching the sorted `{m·b_j}` for each `a_i`.
///
/// # Score — extreme-value corrected
/// `obs` is the count in the **densest** window, i.e. a MAXIMUM over the
/// `M_eff ≈ c_range / (2·tol)` window positions scanned. Comparing that maximum to
/// the single-window mean `E` is a selection bias: even pure scatter beats `E` by an
/// extreme-value margin, and the margin grows as `√(2·ln M_eff)` — so an uncorrected
/// `z = (obs−E)/√E` reads several "sigma" on noise and, worse, drifts with read
/// length / seed density (more anchors, longer reads → larger `M_eff` → larger null
/// `z`). A fixed threshold is then not comparable across pairs (measured: pure-null
/// `z` ran 1.5–4 and up to 9 before this correction; see `null_gate_calibration_probe`).
///
/// We therefore compare `obs` to the **expected null maximum**, not the mean:
///
/// ```text
///   E_max = the value where P(Poisson(E) ≥ E_max) = 1 / M_eff      (typical max of
///                                                                   M_eff windows)
///   z     = (obs − E_max) / √max(E, 1)
/// ```
///
/// The exact Poisson upper tail (not a Gaussian approximation) keeps it honest at the
/// small-`E`, sparse-cloud end, where a Gaussian floor misbehaves. Under the null `obs ≈
/// E_max`, so `z ≈ 0` regardless of `N` or read length; a genuine diagonal keeps
/// `obs ≫ E_max` and stays at `z ≫ 0`. `M_eff = 1` (a cloud tighter than one window)
/// degrades gracefully to `E_max = E`. Deterministic (no RNG).
///
/// Returns `0.0` (never significant) for fewer than [`SLOPE_MIN_ANCHORS`] anchors.
/// This augments — never replaces — wtdbg2's `min_sim`/extent gates and the
/// `max_similarity`/`max_anchor_density` repeat guards; it is the raw-signal analog
/// of masking repeats by seed occurrence (`kmax`, `kbm.h`), phrased as a per-pair
/// significance test rather than a global mask.
pub fn diagonal_significance(anchors: &[Anchor], slope: f64, tol: f64) -> f64 {
    diagonal_significance_parts(anchors, slope, tol).map_or(0.0, |p| p.z)
}

/// The individual terms behind [`diagonal_significance`]'s score.
///
/// The score collapses four quantities into one number, which makes the gate
/// impossible to calibrate after the fact: a given `z` can come from a small `obs`
/// against a tiny `e`, or a large one against a large `e`, and those want different
/// thresholds. Emitting the terms lets an offline study fit the null directly
/// (exact Poisson tail p-values, length-adaptive or FDR-controlled thresholds)
/// instead of re-deriving them from `z` alone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NullParts {
    /// Anchors in the densest `±tol` offset window — the observed statistic.
    pub obs: u32,
    /// Single-window null mean `E` under marginal-preserving random pairing.
    pub e: f64,
    /// Expected null *maximum* over the `m_eff` scanned windows.
    pub e_max: f64,
    /// Number of window positions the maximum was taken over.
    pub m_eff: f64,
    /// The score itself: `(obs − e_max) / √max(e, 1)`.
    pub z: f64,
}

/// [`diagonal_significance`] with its intermediate terms exposed. Returns `None`
/// for fewer than [`SLOPE_MIN_ANCHORS`] anchors (where the score is defined as 0).
pub fn diagonal_significance_parts(anchors: &[Anchor], slope: f64, tol: f64) -> Option<NullParts> {
    if anchors.len() < SLOPE_MIN_ANCHORS {
        return None;
    }
    Some(significance_parts(&Cloud::new(anchors), slope, tol))
}

/// [`diagonal_significance_parts`] of a cloud of at least [`SLOPE_MIN_ANCHORS`] anchors.
fn significance_parts(cloud: &Cloud, slope: f64, tol: f64) -> NullParts {
    let (anchors, n) = (cloud.anchors, cloud.anchors.len());
    let tol = tol.max(1.0);
    if slope.is_nan() || slope < 0.0 {
        return significance_parts_sorted(anchors, slope, tol);
    }

    // Observed: densest ±tol offset window over c_i = a_i − slope·b_i (the sliding
    // window of `diagonal_peak`), and its center c0.
    let (obs, c0) = peak_at(cloud, slope, tol);
    let (mut c_min, mut c_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for c in cloud.offsets(slope) {
        c_min = c_min.min(c);
        c_max = c_max.max(c);
    }

    // Null expectation under marginal-preserving random pairing:
    //   E = (1/N) · #{ (i,j) : |a_i − slope·b_j − c0| ≤ tol }.
    // For each a_i, count the {slope·b_j} inside [a_i−c0−tol, a_i−c0+tol]. At `slope >= 0`
    // both `slope·b` and the interval ends are monotone in the integer positions (rounding
    // is monotone), so sorting the positions orders them: one merge-like sweep over the
    // sorted `a` counts what a binary search per anchor would. The anchors arrive in
    // `a` order from the gather, so usually only `b` is sorted.
    let mut bv: Vec<u32> = anchors.iter().map(|x| x.b_pos).collect();
    sort_u32(&mut bv);
    let (mut l, mut r, mut m_pairs) = (0usize, 0usize, 0u64);
    let mut sweep = |a: u32| {
        let center = a as f64 - c0;
        let lo = center - tol;
        let hi = center + tol;
        while l < n && slope * (bv[l] as f64) < lo {
            l += 1;
        }
        while r < n && slope * (bv[r] as f64) <= hi {
            r += 1;
        }
        m_pairs += (r - l) as u64;
    };
    if anchors.is_sorted_by_key(|x| x.a_pos) {
        anchors.iter().for_each(|x| sweep(x.a_pos));
    } else {
        let mut av: Vec<u32> = anchors.iter().map(|x| x.a_pos).collect();
        sort_u32(&mut av);
        av.into_iter().for_each(&mut sweep);
    }
    null_parts(obs, m_pairs, n, c_max - c_min, 2.0 * tol)
}

/// Sort positions: LSD radix sort on 11-bit digits for large inputs (skipping digits
/// every value shares), `sort_unstable` below. Either gives the one sorted sequence.
fn sort_u32(v: &mut Vec<u32>) {
    const D: u32 = 11;
    const NB: usize = 1 << D;
    if v.len() < 512 {
        v.sort_unstable();
        return;
    }
    let (min, max) = v
        .iter()
        .fold((u32::MAX, 0), |(lo, hi), &x| (lo.min(x), hi.max(x)));
    let span = max - min;
    let mut tmp = vec![0u32; v.len()];
    let mut shift = 0u32;
    while shift < 32 && (span >> shift) > 0 {
        let mut cnt = [0u32; NB];
        for &x in v.iter() {
            cnt[(((x - min) >> shift) as usize) & (NB - 1)] += 1;
        }
        let mut sum = 0u32;
        for c in cnt.iter_mut() {
            let t = *c;
            *c = sum;
            sum += t;
        }
        for &x in v.iter() {
            let d = (((x - min) >> shift) as usize) & (NB - 1);
            tmp[cnt[d] as usize] = x;
            cnt[d] += 1;
        }
        std::mem::swap(v, &mut tmp);
        shift += D;
    }
}

/// [`diagonal_significance_parts`] by sorting the offsets and binary-searching the
/// products: the definition, and the path for a negative slope.
fn significance_parts_sorted(anchors: &[Anchor], slope: f64, tol: f64) -> NullParts {
    let n = anchors.len();
    let win = 2.0 * tol;
    // Unstable sorts are exact here for the reason given at `peak_sorted`.
    let (mut keys, mut cs, mut sb) = (Vec::new(), Vec::new(), Vec::new());
    sort_total_into(
        anchors
            .iter()
            .map(|a| a.a_pos as f64 - slope * a.b_pos as f64),
        &mut keys,
        &mut cs,
    );
    let (obs, c_lo, c_hi) = densest_window(&cs, win);
    let c0 = (c_lo + c_hi) / 2.0;
    // For each a_i count sorted {slope·b_j} inside [a_i−c0−tol, a_i−c0+tol].
    sort_total_into(
        anchors.iter().map(|a| slope * a.b_pos as f64),
        &mut keys,
        &mut sb,
    );
    let mut m_pairs: u64 = 0;
    for a in anchors {
        let center = a.a_pos as f64 - c0;
        let lo = center - tol;
        let hi = center + tol;
        let l = sb.partition_point(|&v| v < lo);
        let r = sb.partition_point(|&v| v <= hi);
        m_pairs += (r - l) as u64;
    }
    null_parts(obs, m_pairs, n, cs[cs.len() - 1] - cs[0], win)
}

/// The score from its counts: `obs` in the densest window, `m_pairs` random pairings in
/// it, over `n` anchors whose offsets span `c_range`.
fn null_parts(obs: usize, m_pairs: u64, n: usize, c_range: f64, win: f64) -> NullParts {
    let e = m_pairs as f64 / n as f64;
    // Extreme-value correction: `obs` is the max over M_eff scanned windows, so
    // compare it to the expected null MAXIMUM, not the single-window mean `E`.
    let m_eff = (c_range / win).max(1.0);
    let e_max = expected_null_max(e, m_eff);
    NullParts {
        obs: obs as u32,
        e,
        e_max,
        m_eff,
        z: (obs as f64 - e_max) / e.max(1.0).sqrt(),
    }
}

/// `ln k!` — summed exactly for small `k`, Stirling with two correction terms above.
/// Accurate to well under 1e-9 everywhere, which is far finer than a tail reported as
/// `-log10`.
fn ln_factorial(k: u32) -> f64 {
    const SMALL: u32 = 16;
    if k <= SMALL {
        return (2..=k).map(|i| (i as f64).ln()).sum();
    }
    let k = k as f64;
    k * k.ln() - k + 0.5 * (std::f64::consts::TAU * k).ln() + 1.0 / (12.0 * k)
        - 1.0 / (360.0 * k * k * k)
}

/// `−log10` of the probability that the densest of `m_eff` windows reaches `obs`
/// anchors by chance — the exact form of what [`NullParts::z`] approximates.
///
/// **Measured, not adopted as a replacement for `z`**: within a dataset it keeps the same pairs at recall ≤0.70, and across datasets the
/// cut for a fixed precision spreads 9.7× against `z`'s 6.8× — `−log10 p` is near-linear in
/// `obs` with a dataset-dependent slope. Used by the partner pre-filter and the contig-end
/// rescue; not on the alignment's hot path.
///
/// `z = (obs − e_max)/√max(e, 1)` is a Gaussian standardisation of a Poisson maximum,
/// and the `max(e, 1)` clamp means that wherever `e < 1` — few anchors, i.e. a short
/// alignment — it stops being a standardised deviate and becomes a raw count. Measured
/// on the dmel arm, precision at a fixed `z` then swings with alignment extent and does
/// so *backwards*: at `z ∈ [4, 6)` it is 0.567 below 10 bins and 0.081 at 25-40, because
/// a long sparse repeat smear accumulates extent without significance. A tail
/// probability carries its own multiple-testing correction and has no such dependence,
/// which is also what makes it comparable across datasets.
///
/// Both tails are taken in log space: the interesting pairs sit where `p` underflows.
pub fn tail_nlp(obs: u32, e: f64, m_eff: f64) -> f64 {
    if obs == 0 {
        return 0.0;
    }
    let e = e.max(1e-12);
    let m = m_eff.max(1.0);

    // ln P(X ≥ obs), X ~ Poisson(e). Sum upward from `obs` as ratios to the first term,
    // so the series stays near unity and only the logarithm is ever taken of a small
    // number. Converges geometrically once `e < obs`, which is the regime of interest.
    let k = obs as f64;
    let (mut sum, mut term, mut i) = (1.0f64, 1.0f64, 0.0f64);
    while term > 1e-17 * sum && i < 10_000.0 {
        i += 1.0;
        term *= e / (k + i);
        sum += term;
    }
    let ln_p = -e + k * e.ln() - ln_factorial(obs) + sum.ln();
    if ln_p >= 0.0 {
        return 0.0;
    }

    // P(any of m windows) = 1 − (1 − p)^m. `ln_1p`/`exp_m1` keep this exact for small
    // `p`; when `p` itself underflows to zero the union bound `m·p` is exact to every
    // digit that survives, so take it in logs.
    let p = ln_p.exp();
    if p == 0.0 {
        return -(ln_p + m.ln()) / std::f64::consts::LN_10;
    }
    let any = -(m * (-p).ln_1p()).exp_m1();
    if any <= 0.0 {
        return -(ln_p + m.ln()) / std::f64::consts::LN_10;
    }
    -any.log10()
}

/// An upper bound on [`tail_nlp`]`(obs, e, m_eff)` without its series, for deciding that a score
/// falls short of a level (the partner pre-filter, `prefilter::reaches`).
///
/// [`tail_nlp`] computes `ln p = ln_pmf + ln(sum)` with `ln_pmf = -e + k ln e - ln k!` (the same
/// expression, evaluated here operation for operation) and `sum >= 1`, so its `ln p >= ln_pmf`
/// exactly in floating point; and `1 - (1 - p)^m >= max(p, m p / (1 + m p))` for `m >= 1` (its
/// union-bound branches are larger still). So `-log10 max(p_pmf, m p_pmf / (1 + m p_pmf))` bounds
/// the score from above, up to the last-place rounding of `exp`, `ln_1p`, `exp_m1` and `log10`
/// (relative ~1e-15; callers keep a margin far above it). `0` where [`tail_nlp`] returns `0`
/// (`obs == 0`, or `ln_pmf >= 0`); `inf` where the point mass underflows (no bound).
pub fn tail_nlp_upper(obs: u32, e: f64, m_eff: f64) -> f64 {
    if obs == 0 {
        return 0.0;
    }
    let e = e.max(1e-12);
    let m = m_eff.max(1.0);
    let k = obs as f64;
    let ln_pmf = -e + k * e.ln() - ln_factorial(obs);
    if ln_pmf >= 0.0 {
        return 0.0;
    }
    let p = ln_pmf.exp();
    if p == 0.0 {
        return f64::INFINITY;
    }
    let x = m * p;
    -(p.max(x / (1.0 + x))).log10()
}

/// `e` above which the forward walk from `P(X = 0)` is abandoned for the
/// mode-anchored recurrence. `exp(-e)` goes denormal near `e = 709` and is exactly
/// zero past `e ≈ 745`; the cutoff sits well below that so the walk stays on the
/// values it has always been validated against, and only the regime it cannot
/// represent at all takes the other branch.
const NULL_MAX_DIRECT_E: f64 = 500.0;

/// Typical maximum window count under the null: the value `k` at which the Poisson
/// upper tail `P(X ≥ k | mean = e)` crosses `1 / m_eff` — i.e. the count you expect
/// the densest of `m_eff` independent windows to reach by chance. Uses the exact
/// Poisson tail (correct at small `e`, where a Gaussian extreme-value form fails) and
/// linearly interpolates the crossing so the result is continuous in `e`. Returns `e`
/// itself when `m_eff ≤ 1` (no maximization to correct for).
fn expected_null_max(e: f64, m_eff: f64) -> f64 {
    if m_eff <= 1.0 || e <= 0.0 {
        return e;
    }
    if e > NULL_MAX_DIRECT_E {
        return expected_null_max_dense(e, m_eff);
    }
    let target = 1.0 / m_eff;
    // Walk k = 0,1,2,… tracking pmf(k) and the upper tail P(X ≥ k).
    let mut pmf = (-e).exp(); // P(X = 0)
    let mut upper = 1.0; // P(X ≥ 0)
    let mut k = 0usize;
    loop {
        let upper_next = (upper - pmf).max(0.0); // P(X ≥ k+1)
        if upper_next < target {
            // Crossing lies in (k, k+1]; interpolate on the tail probability.
            let frac = (upper - target) / (upper - upper_next).max(1e-300);
            return k as f64 + frac.clamp(0.0, 1.0);
        }
        k += 1;
        upper = upper_next;
        pmf *= e / k as f64;
        if k > 1_000_000 {
            return k as f64;
        }
    }
}

/// [`expected_null_max`] for large `e`, where starting the walk at `P(X = 0)` is
/// impossible: `exp(-e)` underflows to zero, `upper − pmf` never decreases, and the
/// walk runs to its iteration cap and returns ~10⁶ — an `E_max` so large that `z`
/// goes hugely negative and the null gate rejects the pair outright. Dense repeat
/// clouds are exactly where `e` gets large, so the failure was silent and
/// selective rather than loud.
///
/// Same quantity, computed from the **mode** outward instead of from zero. Working
/// with `s(k) = pmf(k) / pmf(mode)` — built by the Poisson recurrence
/// `s(k+1) = s(k)·e/(k+1)` upward and `s(k−1) = s(k)·k/e` downward — keeps every
/// term near unity, so nothing underflows regardless of `e`. Tail probabilities are
/// ratios `Σ_{j≥k} s(j) / Σ_j s(j)`, in which the unknown `pmf(mode)` cancels.
fn expected_null_max_dense(e: f64, m_eff: f64) -> f64 {
    let target = 1.0 / m_eff;
    let mode = e.floor();
    // Terms below this fraction of the mode cannot move a tail ratio we compare
    // against `target ≥ 1/m_eff`; stopping there bounds both loops to O(√e) steps.
    const NEGLIGIBLE: f64 = 1e-25;

    // Total mass, swept both ways from the mode. `k_hi` is where the upper sweep
    // stopped — the top of the support we need to consider.
    let mut total = 1.0;
    let mut s = 1.0;
    let mut k = mode;
    while s > NEGLIGIBLE {
        k += 1.0;
        s *= e / k;
        total += s;
    }
    let k_hi = k;
    s = 1.0;
    k = mode;
    while k > 0.0 && s > NEGLIGIBLE {
        s *= k / e;
        k -= 1.0;
        total += s;
    }

    // Sweep down from `k_hi` accumulating the upper tail until it reaches `target`.
    // `tail` here is P(X ≥ k)·total; the first k that reaches it brackets the
    // crossing in [k, k+1), interpolated exactly as the direct walk does.
    let mut tail = 0.0;
    s = 1.0;
    k = k_hi;
    // Rebuild s(k_hi) by the same upward recurrence so both sweeps agree bit-for-bit.
    let mut j = mode;
    while j < k_hi {
        j += 1.0;
        s *= e / j;
    }
    loop {
        let upper_next = tail / total; // P(X ≥ k+1)
        tail += s;
        let upper = tail / total; // P(X ≥ k)
        if upper >= target {
            let frac = (upper - target) / (upper - upper_next).max(1e-300);
            return k + frac.clamp(0.0, 1.0);
        }
        if k <= 0.0 {
            return 0.0;
        }
        s *= k / e;
        k -= 1.0;
    }
}

/// Hough confidence: the winning `peak` over the median peak across the scanned
/// slopes. A real diagonal towers over the wrong-slope background; noise does not.
fn peak_conf(peak: usize, mut peaks: Vec<usize>) -> f64 {
    peaks.sort_unstable();
    let median = peaks[peaks.len() / 2].max(1) as f64;
    peak as f64 / median
}

/// Full slope range for the confidence background (the "typical wrong-slope peak").
/// Sampled coarsely: the median is a robust statistic of a smooth, low background, so a
/// sparse sample matches the dense one at a fraction of the sorts. NOT windowed — a
/// window-only background would be mostly near-diagonal, inflating the median and
/// wrongly rejecting good overlaps.
const BG_LO: f64 = 0.2;
const BG_HI: f64 = 5.0;
const BG_STEP: f64 = 0.3;
/// Finest slope step; the beam refines down to this.
const FINE_STEP: f64 = 0.02;

/// Evaluate [`diagonal_peak`] at every slope in `[lo, hi]` stepping by `step`, one
/// `(slope, peak, offset)` candidate per sample. (`step` is accumulated `m += step` so
/// every caller gets the same grid values, bit for bit.)
fn scan_candidates(cloud: &Cloud, tol: f64, lo: f64, hi: f64, step: f64) -> Vec<(f64, usize, f64)> {
    let mut ms = Vec::new();
    let mut m = lo;
    while m <= hi + 1e-9 {
        ms.push(m);
        m += step;
    }
    peaks_at(cloud, &ms, tol)
}

/// Geometric slope-step schedule from `coarse` down to `fine` over `levels` levels
/// (coarsest first; the finest is pinned exactly to `fine`). Collapses to a single
/// `fine` level when `levels <= 1` or `coarse <= fine`.
fn step_schedule(coarse: f64, fine: f64, levels: u32) -> Vec<f64> {
    let levels = levels.max(1) as usize;
    if levels == 1 || coarse <= fine {
        return vec![coarse.min(fine)];
    }
    let ratio = (fine / coarse).powf(1.0 / (levels - 1) as f64);
    let mut steps: Vec<f64> = (0..levels).map(|i| coarse * ratio.powi(i as i32)).collect();
    *steps.last_mut().unwrap() = fine; // pin the finest level exactly to Full's step
    steps
}

/// The top-`k` **distinct** basins among candidate `(slope, peak, offset)` triples:
/// taken by descending peak (ties → smaller slope, for determinism), skipping any whose
/// slope is within `min_sep` of an already-chosen basin — so they are separate diagonals
/// rather than several samples of one apex.
fn top_k_distinct(
    mut cands: Vec<(f64, usize, f64)>,
    k: usize,
    min_sep: f64,
) -> Vec<(f64, usize, f64)> {
    cands.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.total_cmp(&b.0)));
    let mut kept: Vec<(f64, usize, f64)> = Vec::with_capacity(k);
    for c in cands {
        if kept.len() >= k {
            break;
        }
        if kept.iter().all(|&(s, _, _)| (s - c.0).abs() >= min_sep) {
            kept.push(c);
        }
    }
    kept
}

/// The first filter's null for [`BeamCfg::screen_z`]: the mean count of one `2·tol` offset
/// window if the cloud's anchors were spread uniformly over its offset range, which at
/// slope ~1 spans `a_ext + b_ext` events. One pass over the anchors, no sort.
fn screen_lambda(anchors: &[Anchor], tol: f64) -> f64 {
    let (mut amin, mut amax, mut bmin, mut bmax) = (u32::MAX, 0u32, u32::MAX, 0u32);
    for x in anchors {
        amin = amin.min(x.a_pos);
        amax = amax.max(x.a_pos);
        bmin = bmin.min(x.b_pos);
        bmax = bmax.max(x.b_pos);
    }
    let range = (amax.saturating_sub(amin) as u64 + bmax.saturating_sub(bmin) as u64).max(1);
    anchors.len() as f64 * 2.0 * tol / range as f64
}

/// A level-0 peak standardised against [`screen_lambda`]: `(peak - lam) / sqrt(max(lam, 1))`.
pub fn screen_z(peak: usize, lam: f64) -> f64 {
    (peak as f64 - lam) / lam.max(1.0).sqrt()
}

/// Stage 0b's statistic ([`BeamCfg::screen_null`]): the final gate's `z` taken at a level-0
/// slope. `+inf` below [`SLOPE_MIN_ANCHORS`], where no line is fitted and the gate never
/// runs, so the screen never decides such a pair.
fn level0_null_z(anchors: &[Anchor], slope: f64, tol: f64) -> f64 {
    diagonal_significance_parts(anchors, slope, tol).map_or(f64::INFINITY, |p| p.z)
}

/// Prior-windowed coarse→fine beam over slope. The peak is sought only inside
/// `[win_lo, win_hi]` — the bin-DP/FBG cannot use far-off-slope diagonals, and excluding
/// them removes the competing off-window peaks that make a single-winner coarse pass
/// lock onto the wrong basin. Each level keeps the top-`k` *distinct* basins and refines
/// each at the next finer step: the coarse level only has to **nominate** the true basin
/// into the beam (easy, and it degrades gracefully as the step coarsens), while the fine
/// level **ranks** the shortlist by true apex.
///
/// Returns the winning `(slope, peak, offset)`, the full-range background peaks (for
/// the confidence median), whether the first filter ([`BeamCfg::screen_peak`])
/// stopped the search after level 0 -- then the winner is the best level-0 sample -- and,
/// with [`BeamCfg::select_z`], the last level's top-`k` distinct basins (highest count
/// first; empty otherwise), among which [`fit_cloud`] picks by the gate's `z`.
#[allow(clippy::type_complexity)]
fn beam_search_slope(
    cloud: &Cloud,
    tol: f64,
    beam: &BeamCfg,
    conf_floor: f64,
) -> (
    (f64, usize, f64),
    Option<Vec<usize>>,
    bool,
    Vec<(f64, usize, f64)>,
) {
    let anchors = cloud.anchors;
    let lo = (beam.win_lo as f64).max(BG_LO);
    let hi = (beam.win_hi as f64).min(BG_HI).max(lo);
    let coarse = (beam.coarse_step as f64).max(FINE_STEP);
    let steps = step_schedule(coarse, FINE_STEP, beam.levels);
    let k = (beam.beam_k as usize).max(1);
    let sep = steps[0]; // "distinct basin" = at least one coarse step apart

    // First filter, part 1: a window holds at most every anchor, so a cloud too small to
    // reach the floor even if all of it sat in one window is screened before level 0.
    let screen = beam.screen_peak as usize;
    let screen_lam = (beam.screen_z > 0.0).then(|| screen_lambda(anchors, tol));
    if (screen > 0 && anchors.len() < screen)
        || screen_lam.is_some_and(|lam| screen_z(anchors.len(), lam) < beam.screen_z as f64)
    {
        return ((lo, anchors.len(), 0.0), None, true, Vec::new());
    }

    // Level 0: coarse scan of the window; keep the top-k distinct basins. A small cloud is
    // swept here too, its order at the lowest slope kept for the last level, which starts there.
    let level0 = if sweeps(cloud, tol) {
        SWEEP.with_borrow_mut(|sw| {
            let mut out = Vec::new();
            let mut m = lo;
            sw.start(cloud, m);
            while m <= hi + 1e-9 {
                let (peak, c0) = sw.peak(m, tol);
                if out.is_empty() {
                    sw.keep(m);
                }
                out.push((m, peak, c0));
                m += steps[0];
            }
            out
        })
    } else {
        scan_candidates(cloud, tol, lo, hi, steps[0])
    };
    // First filter, part 2: the densest level-0 window. Refinement and the null scan are
    // what a random pair costs, and this stops them before they start.
    let screen_null = beam.screen_null as f64;
    if screen > 0 || screen_lam.is_some() || screen_null > 0.0 {
        let best0 = level0
            .iter()
            .copied()
            .reduce(|a, b| if b.1 > a.1 { b } else { a })
            .unwrap_or((lo, 0, 0.0));
        if best0.1 < screen
            || screen_lam.is_some_and(|lam| screen_z(best0.1, lam) < beam.screen_z as f64)
            || (screen_null > 0.0 && level0_null_z(anchors, best0.0, tol) < screen_null)
        {
            return (best0, None, true, Vec::new());
        }
    }
    let mut basins = top_k_distinct(level0, k, sep);

    // Refine: each finer level rescans ±(previous step) around each basin, dedups
    // overlapping samples, and re-selects the top-k distinct basins.
    for w in 1..steps.len() {
        let (prev, step) = (steps[w - 1], steps[w]);
        // Of the last level only the winner is used: the first sample in
        // `top_k_distinct`'s order (highest peak, ties -> smaller slope). So there a
        // slope is only measured if its peak can reach the best so far -- unless
        // `select_z` ranks the last level's distinct basins too (`keep_all`).
        let last = w + 1 == steps.len();
        let keep_all = !last || beam.select_z;
        let mut refined: Vec<(f64, usize, f64)> = Vec::new();
        let mut top: Option<(f64, usize, f64)> = None;
        let mut seen: FxHashSet<i64> = FxHashSet::default();
        // The level's slopes: each basin's window, deduplicated on a half-step grid so
        // overlapping refine windows don't re-measure (the first window to reach a key keeps it).
        let mut level: Vec<f64> = Vec::new();
        for &(bc, _, _) in &basins {
            let rlo = (bc - prev).max(lo);
            let rhi = (bc + prev).min(hi);
            let mut m = rlo;
            while m <= rhi + 1e-9 {
                let key = (m / (step * 0.5)).round() as i64;
                if seen.insert(key) {
                    level.push(m);
                }
                m += step;
            }
        }
        // A small cloud's last level is swept in ascending slope (`SlopeSweep`); a larger one's
        // slopes are counted together, and of the last level only the winner is measured
        // (`best_peak`). Which sample wins does not depend on the order the samples are measured
        // in.
        let sweep = last && sweeps(cloud, tol);
        if !sweep {
            if keep_all {
                refined = peaks_at(cloud, &level, tol);
            } else {
                top = best_peak(cloud, &level, tol);
            }
        } else if !level.is_empty() {
            let mut swept = level;
            swept.sort_unstable_by(f64::total_cmp);
            SWEEP.with_borrow_mut(|sw| {
                if !sw.resume(swept[0]) {
                    sw.start(cloud, swept[0]);
                }
                for &m in &swept {
                    let (peak, c0) = sw.peak(m, tol);
                    if keep_all {
                        refined.push((m, peak, c0));
                    } else if top.is_none_or(|t| peak > t.1) {
                        // Ascending slopes: a tie keeps the earlier, smaller one.
                        top = Some((m, peak, c0));
                    }
                }
            });
        }
        basins = if keep_all {
            top_k_distinct(refined, k, sep)
        } else {
            top.into_iter().collect()
        };
    }
    let finalists = if beam.select_z {
        basins.clone()
    } else {
        Vec::new()
    };

    // Winner = highest apex among the surviving basins (ties → smaller slope).
    let best = basins
        .into_iter()
        .reduce(|a, b| {
            if b.1 > a.1 || (b.1 == a.1 && b.0 < a.0) {
                b
            } else {
                a
            }
        })
        .unwrap_or((lo, 0, 0.0));

    // The background exists only to supply the median in `peak_conf`, and
    // `conf = peak / median` with `median >= 1` means `conf <= peak`. So a winner whose
    // peak is already under the gate cannot clear it no matter what the background
    // says, and its 17 extra sorts — a quarter of this pair's slope evaluations — are
    // dead work. Skipping them cannot change which pairs are accepted; see
    // `BinDp::exact_conf` for the one place the difference is observable.
    if (best.1 as f64) < conf_floor {
        return (best, None, false, finalists);
    }

    // Background over the full slope range (for conf); coarse but robust.
    let bg_peaks: Vec<usize> = scan_candidates(cloud, tol, BG_LO, BG_HI, BG_STEP)
        .into_iter()
        .map(|(_, p, _)| p)
        .collect();
    (best, Some(bg_peaks), false, finalists)
}

/// Candidate first-filter statistics for one cloud, over the level-0 slope grid of
/// [`beam_search_slope`]: what the screen tests check the live screen against.
#[cfg(test)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ScreenStats {
    /// Span of the cloud along read `a`, in events.
    pub a_ext: u32,
    /// Span of the cloud along read `b`, in events.
    pub b_ext: u32,
    /// Densest `2·tol` offset window, max over the coarse slopes -- what level 0 sees.
    pub pk0: u32,
    /// Densest `2·tol + (step/2)·b_ext` window, max over the coarse slopes. An exact
    /// UPPER BOUND on the `2·tol` peak at every slope within half a coarse step of the
    /// grid (which covers the whole beam window): moving the slope by `d` shifts an
    /// anchor's offset by `d·b`, so a window at the finer slope spreads by at most
    /// `|d|·b_ext` at the grid slope.
    pub pkw: u32,
    /// Sort-free bound on `pkw`: the densest pair of ADJACENT bins of an offset histogram
    /// with bins of that widened width (a window of a bin's width touches at most two).
    pub hb: u32,
    /// The gate's null `z` at the level-0 slope holding `pk0` (the first such slope, as the
    /// beam picks it) -- what [`BeamCfg::screen_null`] decides on.
    pub zn0: f64,
}

/// [`ScreenStats`] for `anchors` at detection half-width `tol` over `beam`'s coarse grid.
#[cfg(test)]
pub fn screen_stats(anchors: &[Anchor], tol: f64, beam: &BeamCfg) -> ScreenStats {
    if anchors.is_empty() {
        return ScreenStats::default();
    }
    let (mut amin, mut amax, mut bmin, mut bmax) = (u32::MAX, 0u32, u32::MAX, 0u32);
    for x in anchors {
        amin = amin.min(x.a_pos);
        amax = amax.max(x.a_pos);
        bmin = bmin.min(x.b_pos);
        bmax = bmax.max(x.b_pos);
    }
    let b_ext = bmax - bmin;
    // The level-0 grid, exactly as `beam_search_slope` builds it.
    let lo = (beam.win_lo as f64).max(BG_LO);
    let hi = (beam.win_hi as f64).min(BG_HI).max(lo);
    let step = (beam.coarse_step as f64).max(FINE_STEP);
    let w0 = 2.0 * tol.max(1.0);
    let ww = w0 + step / 2.0 * b_ext as f64;
    let (mut cs, mut hist) = (Vec::with_capacity(anchors.len()), Vec::<u32>::new());
    let (mut pk0, mut pkw, mut hb) = (0u32, 0u32, 0u32);
    let mut m0 = lo;
    let mut m = lo;
    while m <= hi + 1e-9 {
        cs.clear();
        cs.extend(anchors.iter().map(|a| a.a_pos as f64 - m * a.b_pos as f64));
        let cmin = cs.iter().copied().fold(f64::INFINITY, f64::min);
        let cmax = cs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        hist.clear();
        hist.resize(((cmax - cmin) / ww) as usize + 2, 0);
        for &c in &cs {
            hist[((c - cmin) / ww) as usize] += 1;
        }
        hb = hb.max(hist.windows(2).map(|w| w[0] + w[1]).max().unwrap_or(0));
        cs.sort_unstable_by(f64::total_cmp);
        let p0 = densest_window(&cs, w0).0 as u32;
        if p0 > pk0 {
            (pk0, m0) = (p0, m);
        }
        pkw = pkw.max(densest_window(&cs, ww).0 as u32);
        m += step;
    }
    ScreenStats {
        a_ext: amax - amin,
        b_ext,
        pk0,
        pkw,
        hb,
        zn0: level0_null_z(anchors, m0, tol),
    }
}

/// `f64` -> `i64` preserving [`f64::total_cmp`]'s order (the transform `total_cmp` itself
/// applies to both operands), and its inverse.
#[inline(always)]
fn total_key(x: f64) -> i64 {
    let b = x.to_bits() as i64;
    b ^ ((((b >> 63) as u64) >> 1) as i64)
}

#[inline(always)]
fn from_total_key(k: i64) -> f64 {
    f64::from_bits((k ^ ((((k >> 63) as u64) >> 1) as i64)) as u64)
}

/// Fill `out` with `vals` in the exact sequence `sort_unstable_by(f64::total_cmp)` gives.
///
/// The keys are transformed once and sorted as plain integers, instead of transforming
/// both operands inside every comparison: 1.5x faster on the real anchor clouds, with
/// the same output bit for bit (equal keys are identical values, so the instability of
/// the sort cannot show). The Hough slope scan spends most of the detection stage here.
#[inline]
fn sort_total_into(vals: impl Iterator<Item = f64>, keys: &mut Vec<i64>, out: &mut Vec<f64>) {
    keys.clear();
    keys.extend(vals.map(total_key));
    keys.sort_unstable();
    out.clear();
    out.extend(keys.iter().map(|&k| from_total_key(k)));
}

/// A pair's anchors, their extents, and their positions as `f64` (exact), converted once
/// for every slope the pair is scanned at. At a slope `m >= 0` every offset `a − m·b` lies
/// in `[amin − m·bmax, amax − m·bmin]`, and because `f64` rounding is monotone so does every
/// computed offset: the slope scan knows the offsets' range before it computes one.
#[derive(Debug, Clone)]
pub struct Cloud<'a> {
    pub anchors: &'a [Anchor],
    amin: u32,
    amax: u32,
    bmin: u32,
    bmax: u32,
    af: Vec<f64>,
    bf: Vec<f64>,
}

impl<'a> Cloud<'a> {
    pub fn new(anchors: &'a [Anchor]) -> Self {
        let (mut amin, mut amax, mut bmin, mut bmax) = (u32::MAX, 0, u32::MAX, 0);
        for x in anchors {
            amin = amin.min(x.a_pos);
            amax = amax.max(x.a_pos);
            bmin = bmin.min(x.b_pos);
            bmax = bmax.max(x.b_pos);
        }
        Cloud {
            anchors,
            amin,
            amax,
            bmin,
            bmax,
            af: anchors.iter().map(|x| x.a_pos as f64).collect(),
            bf: anchors.iter().map(|x| x.b_pos as f64).collect(),
        }
    }

    /// The offsets `a − m·b`, computed exactly as everywhere else.
    fn offsets(&self, m: f64) -> impl Iterator<Item = f64> + '_ {
        self.af.iter().zip(&self.bf).map(move |(&a, &b)| a - m * b)
    }
}

/// Offset cells per window in [`diagonal_peak`]'s histogram: finer cells prune more
/// windows but cost more per pass. Measured on real clouds (`hough_bench`), 16 is best
/// for full-length reads and 12 for clouds of a few hundred anchors.
const PEAK_CELLS_BIG: usize = 16;
const PEAK_CELLS_SMALL: usize = 12;
/// Clouds above this size use [`PEAK_CELLS_BIG`].
const PEAK_BIG: usize = 1024;
/// Clouds up to this size are sorted whole; the histogram does not pay below it.
const PEAK_SORT_MAX: usize = 96;
/// Cells per anchor at most: a sparse cloud's window gets fewer cells (to [`PEAK_MIN_CELLS`];
/// below that it is sorted whole), since a pass over the cells costs as much as one over the
/// anchors and the cells of a sparse cloud are mostly empty.
const PEAK_MAX_CELLS_PER_ANCHOR: usize = 2;
const PEAK_MIN_CELLS: usize = 2;

/// [`diagonal_peak`]'s buffers, reused: it runs ~40 times per aligned pair and is the
/// whole run's dominant cost, so a fresh allocation per call is not free.
#[derive(Default)]
struct PeakBuf {
    /// Per anchor, in anchor order: its offset and its cell.
    vals: Vec<f64>,
    cells: Vec<u32>,
    /// One slope's counts, then their prefix sums (`cum[t]` = anchors in cells `< t`); `cum2`
    /// the second count (see the count loop); `cur` the scatter's write cursors.
    cum: Vec<u32>,
    cum2: Vec<u32>,
    cur: Vec<u32>,
    /// The cells a window reaching the threshold can end in, and the cells such windows
    /// reach back to.
    cand: Vec<u32>,
    need: Vec<bool>,
    /// The offsets grouped by cell, each cell sorted when a window first needs it.
    by_cell: Vec<f64>,
    /// The whole-sort path's keys and sorted offsets.
    keys: Vec<i64>,
    cs: Vec<f64>,
    /// [`peaks_at`] / [`best_peak`]: one plan per slope, every slope's prefix sums back to
    /// back, and the order the slopes are finished in.
    plans: Vec<Plan>,
    hist: Vec<u32>,
    order: Vec<usize>,
}

thread_local! {
    static PEAK_BUF: RefCell<PeakBuf> = RefCell::new(PeakBuf::default());
}

/// [`diagonal_peak`] by sorting every offset: the definition the histogram path
/// reproduces, and the path for small clouds.
fn peak_sorted(
    offsets: impl Iterator<Item = f64>,
    win: f64,
    keys: &mut Vec<i64>,
    cs: &mut Vec<f64>,
) -> (usize, f64) {
    // Unstable is exact here: only the sorted VALUES are read, never which anchor an
    // equal value came from.
    sort_total_into(offsets, keys, cs);
    densest_sorted(cs, win)
}

/// One slope's histogram in [`diagonal_peak`]: the cell map (`q` cells per window), `nc` cells,
/// or `sorted` when the slope is answered by sorting every offset. Every offset lies in
/// `[lo, hi]` (see [`Cloud`]), so its cell lies in `[0, top]` with `top = (hi - lo) * inv < nc`;
/// and `q <= 2 n win / (hi - lo)` ([`cells_per_window`]) gives `top <= 2 n`, so a cell fits an
/// `i32`, as `kernel` converts it. Once counted ([`count_plans`]) the slope's prefix sums sit at
/// `hist[base..=base + nc]`, with `lb` a lower bound on its peak and `ub` an upper bound on
/// every window.
#[derive(Debug, Clone, Copy)]
struct Plan {
    cells: kernel::Cells,
    q: usize,
    nc: usize,
    sorted: bool,
    base: usize,
    lb: u32,
    ub: u32,
}

impl Plan {
    fn new(cloud: &Cloud, m: f64, tol: f64) -> Plan {
        let (n, win) = (cloud.anchors.len(), 2.0 * tol);
        let lo = cloud.amin as f64 - m * cloud.bmax as f64;
        let hi = cloud.amax as f64 - m * cloud.bmin as f64;
        let q = cells_per_window(n, hi - lo, win);
        let inv = q as f64 / win;
        let top = (hi - lo) * inv;
        let sorted = n <= PEAK_SORT_MAX || m.is_nan() || m < 0.0 || q < PEAK_MIN_CELLS;
        Plan {
            cells: kernel::Cells { m, lo, inv },
            q,
            nc: if sorted { 0 } else { top as usize + 1 },
            sorted,
            base: 0,
            lb: 0,
            ub: u32::MAX,
        }
    }
}

/// [`finish`]'s buffers.
struct Work<'a> {
    vals: &'a mut Vec<f64>,
    cells: &'a mut Vec<u32>,
    cand: &'a mut Vec<u32>,
    need: &'a mut Vec<bool>,
    cur: &'a mut Vec<u32>,
    by_cell: &'a mut Vec<f64>,
}

/// For slope `m`, the densest diagonal's anchor count and the center offset of that
/// diagonal: sort the offsets `c = a − m·b`, slide a `2·tol`-wide window, and take the
/// first window (in offset order) that reaches the largest count. `None` iff that count
/// is under `floor`: a caller that only needs the slopes reaching `floor` lets the scan
/// skip the others. With `floor = 0` it always answers.
///
/// Every answer is the full sort's, bit for bit, without sorting every offset. The
/// offsets are counted into `q` cells per window width (the cell is a monotone function
/// of the offset, so cells are ordered):
///   * a window touches at most `q + 2` consecutive cells, so their sum bounds every
///     window that ends in the last of them (`ub`);
///   * `q − 1` consecutive cells always fit in one window, so the largest such sum is a
///     lower bound on the peak (`lb`).
///
/// Only a cell whose `ub` reaches `lb`, `floor` and the running best can hold the first
/// maximum. Those cells and the `q + 1` before each are sorted (on real clouds, 7% of
/// the offsets) and the window slid over them; sorted cells concatenate to the sorted
/// offsets, so each window's count and ends are the full sort's. A slope whose largest
/// `ub` is under `floor` stops after the count. The bounds leave a whole cell of slack
/// for `f64` rounding, which is ~1e-11 events here against a cell of at least 24.
///
/// The anchor loops run in `kernel` (an AVX2 copy where the CPU has it, the same arithmetic);
/// [`peaks_at`] and [`best_peak`] count many slopes in one pass.
fn diagonal_peak(cloud: &Cloud, m: f64, tol: f64, floor: usize) -> Option<(usize, f64)> {
    let win = 2.0 * tol;
    PEAK_BUF.with_borrow_mut(|buf| {
        let n = cloud.anchors.len();
        let p = Plan::new(cloud, m, tol);
        if p.sorted {
            let r = peak_sorted(cloud.offsets(m), win, &mut buf.keys, &mut buf.cs);
            return (r.0 >= floor).then_some(r);
        }
        let PeakBuf {
            vals,
            cells,
            cum,
            cum2,
            cur,
            cand,
            need,
            by_cell,
            ..
        } = buf;

        // The offsets and their cells, then the count per cell; `cum` becomes the cells'
        // prefix sums.
        vals.resize(n, 0.0);
        cells.resize(n, 0);
        kernel::cells(&cloud.af, &cloud.bf, p.cells, vals, cells);
        // Counted into two arrays, alternate anchors to each: a diagonal's anchors arrive in
        // order and fall in one cell, and one array would make each increment wait on the last.
        cum.clear();
        cum.resize(p.nc + 1, 0);
        cum2.clear();
        cum2.resize(p.nc + 1, 0);
        let mut pairs = cells.chunks_exact(2);
        for t in &mut pairs {
            cum[t[0] as usize + 1] += 1;
            cum2[t[1] as usize + 1] += 1;
        }
        for &t in pairs.remainder() {
            cum[t as usize + 1] += 1;
        }
        for (c, &c2) in cum.iter_mut().zip(cum2.iter()) {
            *c += c2;
        }
        let (lb, ub) = kernel::bounds(cum, p.q);
        if (ub as usize) < floor {
            return None;
        }
        let w = Work {
            vals,
            cells,
            cand,
            need,
            cur,
            by_cell,
        };
        finish(cloud, &p, floor, win, cum, lb, false, w)
    })
}

/// [`diagonal_peak`] after the counts, for a slope whose upper bound reaches `floor`: `cum` is its
/// prefix sums and `lb` its lower bound (`kernel::bounds`). `w.vals` / `w.cells` hold its offsets
/// and cells, or are filled here when `fresh` (the batched counts do not store them).
#[allow(clippy::too_many_arguments)]
fn finish(
    cloud: &Cloud,
    p: &Plan,
    floor: usize,
    win: f64,
    cum: &[u32],
    lb: u32,
    fresh: bool,
    w: Work,
) -> Option<(usize, f64)> {
    let n = cloud.anchors.len();
    let (q, nc) = (p.q, p.nc);
    let Work {
        vals,
        cells,
        cand,
        need,
        cur,
        by_cell,
    } = w;
    if fresh {
        vals.resize(n, 0.0);
        cells.resize(n, 0);
        kernel::cells(&cloud.af, &cloud.bf, p.cells, vals, cells);
    }
    let back = q + 1; // a window ending in cell k starts in a cell >= k - back
    let need0 = (lb as usize).max(floor);
    // The cells a window reaching `need0` can end in, listed without a branch per cell (a
    // sparse cloud's cells are empty or not at random, and that test was most of the
    // kernel's mispredictions).
    cand.clear();
    cand.resize(nc, 0);
    let mut n_cand = 0usize;
    for k in 0..nc {
        let ub = (cum[k + 1] - cum[k.saturating_sub(back)]) as usize;
        cand[n_cand] = k as u32;
        n_cand += ((cum[k + 1] > cum[k]) & (ub >= need0)) as usize;
    }
    let cand = &cand[..n_cand];

    // Group the offsets by cell. Above a floor few windows qualify, so in a large
    // cloud only the cells they reach are filled; the others are never read.
    cur.clear();
    cur.extend_from_slice(&cum[..nc]);
    if by_cell.len() < n {
        by_cell.resize(n, 0.0);
    }
    if floor > 0 && n > PEAK_BIG {
        need.clear();
        need.resize(nc, false);
        let mut marked_to = 0usize;
        for &k in cand {
            let k = k as usize;
            need[k.saturating_sub(back).max(marked_to)..=k].fill(true);
            marked_to = k + 1;
        }
        for (&v, &t) in vals.iter().zip(cells.iter()) {
            if need[t as usize] {
                let slot = &mut cur[t as usize];
                by_cell[*slot as usize] = v;
                *slot += 1;
            }
        }
    } else {
        for (&v, &t) in vals.iter().zip(cells.iter()) {
            let slot = &mut cur[t as usize];
            by_cell[*slot as usize] = v;
            *slot += 1;
        }
    }

    // Ascending cells, the windows ending in each, as `densest_sorted` scans them: the window
    // ending at `j` holds more than `best` offsets iff `by_cell[j] - by_cell[j - best] <= win`,
    // and it starts no lower than cell `k - back`, whose position bounds `j - best`. Only a
    // count above the running best is recorded, so the first window reaching the maximum is the
    // one kept. Every position read lies in a sorted cell: the cells from `k - back` to `k` are
    // sorted before `k`.
    let (mut best, mut best_c) = (0usize, 0.0f64);
    let mut sorted_to = 0usize;
    for &k in cand {
        let k = k as usize;
        let (s0, s1) = (cum[k] as usize, cum[k + 1] as usize);
        let k0 = k.saturating_sub(back);
        let lo = cum[k0] as usize;
        if s1 - lo < best + 1 {
            continue;
        }
        for t in k0.max(sorted_to)..=k {
            by_cell[cum[t] as usize..cum[t + 1] as usize].sort_unstable_by(f64::total_cmp);
        }
        sorted_to = k + 1;
        for j in s0..s1 {
            while j >= lo + best && by_cell[j] - by_cell[j - best] <= win {
                best += 1;
                best_c = (by_cell[j + 1 - best] + by_cell[j]) / 2.0;
            }
        }
    }
    (best >= floor).then_some((best, best_c))
}

/// [`diagonal_peak`] where every slope is needed.
fn peak_at(cloud: &Cloud, m: f64, tol: f64) -> (usize, f64) {
    diagonal_peak(cloud, m, tol, 0).expect("a zero floor always answers")
}

/// Clouds up to this many anchors have their slopes counted together ([`count_plans`]); larger
/// ones one slope at a time, as their histograms (~2 n cells per slope, up to ~60 slopes) would
/// outgrow the caches. Above it lie 0.02-0.08% of the pair loop's CPU (PAIRHIST, round 3).
const BATCH_MAX: usize = 1 << 14;

/// Plan every slope of `ms` into `buf.plans`, count the histogram ones four slopes per pass over
/// the anchors (`kernel::count`), each into its own run of `buf.hist`, and bound each one
/// (`kernel::bounds`): the counts [`diagonal_peak`] takes slope by slope, and the same bounds.
fn count_plans(cloud: &Cloud, ms: &[f64], tol: f64, buf: &mut PeakBuf) {
    let PeakBuf { plans, hist, .. } = buf;
    plans.clear();
    let mut total = 0usize;
    for &m in ms {
        let mut p = Plan::new(cloud, m, tol);
        if !p.sorted {
            p.base = total;
            total += p.nc + 1;
        }
        plans.push(p);
    }
    hist.clear();
    hist.resize(total, 0);
    let none = kernel::Cells {
        m: 0.0,
        lo: 0.0,
        inv: 0.0,
    };
    let (mut cs, mut base, mut g) = ([none; 4], [0u32; 4], 0usize);
    for p in plans.iter().filter(|p| !p.sorted) {
        cs[g] = p.cells;
        base[g] = p.base as u32;
        g += 1;
        if g == 4 {
            kernel::count(&cloud.af, &cloud.bf, &cs, &base, hist);
            g = 0;
        }
    }
    if g > 0 {
        kernel::count(&cloud.af, &cloud.bf, &cs[..g], &base[..g], hist);
    }
    for p in plans.iter_mut().filter(|p| !p.sorted) {
        (p.lb, p.ub) = kernel::bounds(&mut hist[p.base..=p.base + p.nc], p.q);
    }
}

/// Planned slope `k`'s answer at `floor`, exactly [`diagonal_peak`]'s.
fn finish_plan(
    cloud: &Cloud,
    k: usize,
    floor: usize,
    tol: f64,
    buf: &mut PeakBuf,
) -> Option<(usize, f64)> {
    let (p, win) = (buf.plans[k], 2.0 * tol);
    if p.sorted {
        let r = peak_sorted(cloud.offsets(p.cells.m), win, &mut buf.keys, &mut buf.cs);
        return (r.0 >= floor).then_some(r);
    }
    if (p.ub as usize) < floor {
        return None;
    }
    let PeakBuf {
        vals,
        cells,
        cand,
        need,
        cur,
        by_cell,
        hist,
        ..
    } = buf;
    let w = Work {
        vals,
        cells,
        cand,
        need,
        cur,
        by_cell,
    };
    finish(
        cloud,
        &p,
        floor,
        win,
        &hist[p.base..=p.base + p.nc],
        p.lb,
        true,
        w,
    )
}

/// [`peak_at`] at every slope of `ms`, in order: `(slope, peak, center)`.
fn peaks_at(cloud: &Cloud, ms: &[f64], tol: f64) -> Vec<(f64, usize, f64)> {
    if cloud.anchors.len() > BATCH_MAX {
        return ms
            .iter()
            .map(|&m| {
                let (peak, c0) = peak_at(cloud, m, tol);
                (m, peak, c0)
            })
            .collect();
    }
    PEAK_BUF.with_borrow_mut(|buf| {
        count_plans(cloud, ms, tol, buf);
        (0..ms.len())
            .map(|k| {
                let (peak, c0) =
                    finish_plan(cloud, k, 0, tol, buf).expect("a zero floor always answers");
                (ms[k], peak, c0)
            })
            .collect()
    })
}

/// The slope of `ms` with the highest peak (ties: the smaller slope), with its peak and center:
/// what measuring the slopes one by one with [`diagonal_peak`] at the best count so far finds.
/// That answer does not depend on the order the slopes are measured in -- a slope reaching the
/// running best is answered exactly, any other cannot win -- so here they are counted together
/// and finished in descending order of their upper bounds: the best count rises early, and every
/// slope whose bound is under it is never finished.
fn best_peak(cloud: &Cloud, ms: &[f64], tol: f64) -> Option<(f64, usize, f64)> {
    let better = |t: (f64, usize, f64), m: f64, peak: usize| {
        peak > t.1 || (peak == t.1 && m.total_cmp(&t.0).is_lt())
    };
    if cloud.anchors.len() > BATCH_MAX {
        let mut top: Option<(f64, usize, f64)> = None;
        for &m in ms {
            if let Some((peak, c0)) = diagonal_peak(cloud, m, tol, top.map_or(0, |t| t.1))
                && top.is_none_or(|t| better(t, m, peak))
            {
                top = Some((m, peak, c0));
            }
        }
        return top;
    }
    PEAK_BUF.with_borrow_mut(|buf| {
        count_plans(cloud, ms, tol, buf);
        let mut order = std::mem::take(&mut buf.order);
        order.clear();
        order.extend(0..ms.len());
        let plans = &buf.plans;
        order.sort_unstable_by_key(|&k| {
            (
                std::cmp::Reverse(plans[k].ub),
                std::cmp::Reverse(plans[k].lb),
                k,
            )
        });
        let mut top: Option<(f64, usize, f64)> = None;
        for &k in &order {
            let floor = top.map_or(0, |t| t.1);
            if (buf.plans[k].ub as usize) < floor {
                break; // the bounds descend: no later slope reaches the best count
            }
            let m = ms[k];
            if let Some((peak, c0)) = finish_plan(cloud, k, floor, tol, buf)
                && top.is_none_or(|t| better(t, m, peak))
            {
                top = Some((m, peak, c0));
            }
        }
        buf.order = order;
        top
    })
}

/// Tests only: measure every cloud with the histograms instead of [`SlopeSweep`]. The answers are
/// the same (both are the whole sort's).
fn sweep_off() -> bool {
    #[cfg(test)]
    if let Some(v) = NO_SWEEP.with(|c| c.get()) {
        return v;
    }
    false
}

#[cfg(test)]
thread_local! {
    static NO_SWEEP: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Clouds up to this size are refined by [`SlopeSweep`].
const SWEEP_MAX: usize = 256;
const SWEEP_SPARSE_MAX: usize = 512;
/// Whether the beam sweeps `cloud` ([`SlopeSweep`]) instead of measuring each slope with
/// [`diagonal_peak`]: up to [`SWEEP_MAX`] anchors, and up to [`SWEEP_SPARSE_MAX`] when the cloud is
/// sparse (its window would get fewer than the full histogram cells). Beyond that the histogram's
/// floor skips more than a sweep saves, and a dense cloud's diagonals reorder between slopes
/// (measured on real clouds of six datasets: the dense 256-511-anchor clouds of T. denticola 30x
/// sweep 8% slower than they histogram, the sparse ones of A. thaliana as fast).
fn sweeps(cloud: &Cloud, tol: f64) -> bool {
    if sweep_off() {
        return false;
    }
    let n = cloud.anchors.len();
    if n <= SWEEP_MAX {
        return true;
    }
    let span = (cloud.amax - cloud.amin) as f64 + (cloud.bmax - cloud.bmin) as f64;
    n <= SWEEP_SPARSE_MAX && cells_per_window(n, span, 2.0 * tol) < PEAK_CELLS_SMALL
}

/// [`peak_at`] along ascending slopes. Between slopes 0.02 apart a small cloud's offsets
/// change order rarely (0.01-0.45 swaps per anchor on the A. thaliana and dmel clouds below
/// 512 anchors, `hough_bench inversions`), so the previous slope's order, re-evaluated, is
/// nearly sorted and an insertion pass orders it in `O(n + swaps)` where a fresh sort costs
/// `O(n log n)`. The positions are kept permuted into that order, so each slope's offsets are
/// computed in one sequential pass. Any sort gives the one sorted sequence of values, so each
/// answer is the full sort's, bit for bit.
#[derive(Default)]
struct SlopeSweep {
    /// Positions (as `f64`) and offsets, in ascending offset at the last slope.
    pa: Vec<f64>,
    pb: Vec<f64>,
    c: Vec<f64>,
    keys: Vec<i64>,
    idx: Vec<u32>,
    /// A slope and the positions in their order there, to resume from ([`Self::keep`]).
    kept: Option<f64>,
    kept_pa: Vec<f64>,
    kept_pb: Vec<f64>,
}

impl SlopeSweep {
    /// Remember the current order, measured at slope `m`.
    fn keep(&mut self, m: f64) {
        self.kept = Some(m);
        self.kept_pa.clone_from(&self.pa);
        self.kept_pb.clone_from(&self.pb);
    }

    /// Return to the order [`Self::keep`] remembered, if it was at slope `m` (exactly).
    fn resume(&mut self, m: f64) -> bool {
        if self.kept.is_none_or(|k| k.to_bits() != m.to_bits()) {
            return false;
        }
        self.pa.clone_from(&self.kept_pa);
        self.pb.clone_from(&self.kept_pb);
        true
    }

    /// Sort the cloud afresh at slope `m`.
    fn start(&mut self, cloud: &Cloud, m: f64) {
        let n = cloud.anchors.len();
        self.kept = None;
        self.keys.clear();
        self.keys.extend(cloud.offsets(m).map(total_key));
        self.idx.clear();
        self.idx.extend(0..n as u32);
        let keys = &self.keys;
        self.idx.sort_unstable_by_key(|&i| keys[i as usize]);
        self.pa.clear();
        self.pa
            .extend(self.idx.iter().map(|&i| cloud.af[i as usize]));
        self.pb.clear();
        self.pb
            .extend(self.idx.iter().map(|&i| cloud.bf[i as usize]));
        self.c.clear();
        self.c.resize(n, 0.0);
    }

    /// The densest window at slope `m` (at or above the last slope, for speed only):
    /// `kernel::sweep`.
    fn peak(&mut self, m: f64, tol: f64) -> (usize, f64) {
        kernel::sweep(&mut self.pa, &mut self.pb, &mut self.c, m, 2.0 * tol)
    }
}

thread_local! {
    static SWEEP: RefCell<SlopeSweep> = RefCell::new(SlopeSweep::default());
}

/// Cells per window for `n` offsets spread over `span`: [`PEAK_CELLS_BIG`] or
/// [`PEAK_CELLS_SMALL`], fewer where the cells would outnumber the anchors.
fn cells_per_window(n: usize, span: f64, win: f64) -> usize {
    let q = if n > PEAK_BIG {
        PEAK_CELLS_BIG
    } else {
        PEAK_CELLS_SMALL
    };
    let fit = (PEAK_MAX_CELLS_PER_ANCHOR * n) as f64 * win / span.max(1.0);
    if fit >= q as f64 { q } else { fit as usize }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Gate attribution has to know WHICH floor ended a pair, so `hough_aln` stamps
    /// [`Reject`] on every exit, rather than leaving it to be inferred from which fields
    /// happen to be zero.
    #[test]
    fn hough_stamps_the_reject_reason_on_every_exit() {
        assert_eq!(BinAln::default().reason, Reject::NotMeasured);

        // Below SLOPE_MIN_ANCHORS: no line can be fitted at all.
        let few: Vec<Anchor> = (0..4).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        assert_eq!(hough_aln(&few, &dp(), 3.0).reason, Reject::TooFewAnchors);

        // A clean diagonal is accepted.
        let clean: Vec<Anchor> = (0..120).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        let ok = hough_aln(&clean, &dp(), 3.0);
        assert!(ok.is_overlap());
        assert_eq!(ok.reason, Reject::Accepted);

        // The SAME cloud under an unreachable min_conf fails on the conf gate -- and
        // still reports the confidence it was measured at, which is the whole reason
        // `conf` is carried on rejected results.
        let gated = hough_aln(&clean, &dp(), 1e9);
        assert!(!gated.is_overlap());
        assert_eq!(gated.reason, Reject::Conf);
        assert!(gated.conf > 0.0, "conf must survive a conf rejection");

        // Tags are what the dump writes; keep them stable.
        assert_eq!(Reject::Conf.tag(), "conf");
        assert_eq!(Reject::Accepted.tag(), "accepted");
        assert_eq!(Reject::Gappy.tag(), "gappy");
    }

    /// wtdbg2's `max_gap` (`kbm.h:1836`) applies on the Hough path too. It did not,
    /// which let through the sparse smears a repeat produces: of the false overlaps
    /// surviving to `z >= 12` on the dmel arm, 52.9% exceeded it against 1.85% of the
    /// true ones.
    ///
    /// The cloud below is a perfect diagonal whose anchors sit 16 events apart. With
    /// `bin_events = 4` that occupies every fourth bin, so the path covers a quarter of
    /// the extent and `gap/aligned_bins ≈ 0.75`. It clears the extent and matched floors
    /// and sits exactly on its fitted line, so nothing but the gap condition can reject
    /// it.
    #[test]
    fn hough_rejects_a_gappy_smear() {
        let sparse: Vec<Anchor> = (0..60)
            .map(|k| Anchor {
                a_pos: k * 16,
                b_pos: k * 16,
            })
            .collect();

        let mut p = dp();
        p.max_gap = 0.6;
        let tight = hough_aln(&sparse, &p, 3.0);
        assert!(
            !tight.is_overlap(),
            "a 75%-gap smear must not be an overlap"
        );
        assert_eq!(tight.reason, Reject::Gappy);

        // With the condition relaxed past the cloud's gap fraction the SAME cloud is
        // accepted, so the rejection is the gap floor and not one of the others.
        p.max_gap = 0.95;
        let loose = hough_aln(&sparse, &p, 3.0);
        assert!(loose.is_overlap(), "reason was {:?}", loose.reason);
        assert!(
            loose.gap as f64 > 0.6 * loose.aligned_bins as f64,
            "gap {} of {} bins should exceed the 0.6 the tight arm rejected",
            loose.gap,
            loose.aligned_bins
        );

        // A dense diagonal of the same extent is unaffected either way.
        let clean: Vec<Anchor> = (0..240).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        p.max_gap = 0.6;
        assert!(hough_aln(&clean, &p, 3.0).is_overlap());
    }

    fn dp() -> BinDp {
        // Baseline-ish knobs sized for the synthetic tests: bin=4 events, e=2, and
        // lightened penalties so a clean diagonal survives at the toy scale (the
        // real preset penalties are exercised on the fixture, not here).
        BinDp {
            bin_events: 4.0,
            e: 2,
            pgap: -2,
            pvar: -3,
            max_bgap: 4,
            max_bvar: 4,
            min_aln_bins: 8,
            min_mat: 4,
            min_sim: 0.05,
            max_gap: 0.6,
            aln_var: 0.25,
            slope_beam: BeamCfg::default(),
            // Off: these toy clouds are sized for the geometry tests, not for the
            // significance gate.
            min_z: 0.0,
            exact_conf: false,
        }
    }

    /// A clean slope-1 collinear anchor run is recovered end to end.
    #[test]
    fn slope1_collinear_recovered() {
        let anchors: Vec<Anchor> = (0..80).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        let aln = bin_dp(&anchors, 1.0, &dp());
        assert!(aln.is_overlap(), "clean diagonal should be recovered");
        // ~20 bins of 4 events; near-full coverage.
        assert!(aln.aligned_bins >= 8, "aligned_bins={}", aln.aligned_bins);
        assert!(aln.similarity > 0.3, "similarity={}", aln.similarity);
    }

    /// A slope-2 diagonal is under-recovered by the slope-1 baseline (it walks off
    /// the band), but recovered when `bin_dp` is handed the diagonal's slope — pinning
    /// both the wtdbg2 behavior and the rescaling.
    #[test]
    fn slope2_needs_rescaling() {
        // a advances twice as fast as b: a = 2*i, b = i. Diagonal slope (a vs b) = 2.
        let anchors: Vec<Anchor> = (0..80)
            .map(|i| Anchor {
                a_pos: 2 * i,
                b_pos: i,
            })
            .collect();

        let baseline = bin_dp(&anchors, 1.0, &dp());
        assert!(
            !baseline.is_overlap(),
            "slope-1 baseline should under-recover a slope-2 diagonal (aligned_bins={}, sim={})",
            baseline.aligned_bins,
            baseline.similarity
        );

        let corrected = bin_dp(&anchors, 2.0, &dp());
        assert!(
            corrected.is_overlap(),
            "rescaling by the true slope should recover the overlap"
        );
        assert!(
            corrected.aligned_bins >= 8,
            "aligned_bins={}",
            corrected.aligned_bins
        );
    }

    /// An anti-collinear (reversed) anchor set is a null — nothing should pass.
    #[test]
    fn anticollinear_finds_nothing() {
        let anchors: Vec<Anchor> = (0..80)
            .map(|i| Anchor {
                a_pos: i,
                b_pos: 79 - i,
            })
            .collect();
        let aln = bin_dp(&anchors, 1.0, &dp());
        assert!(!aln.is_overlap(), "reversed diagonal must not be accepted");
    }

    #[test]
    fn empty_anchors_is_empty() {
        let aln = bin_dp(&[], 1.0, &dp());
        assert!(!aln.is_overlap());
        assert_eq!(aln.aligned_bins, 0);
    }

    /// Deterministic LCG so the noise tests are reproducible.
    fn lcg(s: &mut u64) -> u32 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*s >> 33) as u32
    }

    /// `diagonal_peak` (reused buffer, unstable sort) is bit-identical to the plain
    /// form (fresh vector, stable sort) on tie-heavy clouds of every size, at slopes on
    /// and off the grid -- including across calls of different lengths on one thread,
    /// which is what the reused buffer has to survive.
    #[test]
    fn diagonal_peak_matches_the_stable_sort_form() {
        fn reference(anchors: &[Anchor], m: f64, tol: f64) -> (usize, f64) {
            let mut cs: Vec<f64> = anchors
                .iter()
                .map(|a| a.a_pos as f64 - m * a.b_pos as f64)
                .collect();
            cs.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let (mut i, mut peak, mut peak_c) = (0usize, 0usize, cs[0]);
            for j in 0..cs.len() {
                while cs[j] - cs[i] > 2.0 * tol {
                    i += 1;
                }
                if j - i + 1 > peak {
                    peak = j - i + 1;
                    peak_c = (cs[i] + cs[j]) / 2.0;
                }
            }
            (peak, peak_c)
        }
        let mut s = 7u64;
        for &n in &[5usize, 17, 64, 700, 5000, 33] {
            // Positions on a small range so equal offsets are common.
            let anchors: Vec<Anchor> = (0..n)
                .map(|_| Anchor {
                    a_pos: lcg(&mut s) % 400,
                    b_pos: lcg(&mut s) % 400,
                })
                .collect();
            for &m in &[0.2, 0.5, 1.0, 1.02, 1.3, 2.0, 5.0] {
                for &tol in &[1.0, 4.0, 96.0] {
                    let (p, c) = peak_at(&Cloud::new(&anchors), m, tol);
                    let (rp, rc) = reference(&anchors, m, tol);
                    assert_eq!(p, rp, "n={n} m={m} tol={tol}");
                    assert_eq!(c.to_bits(), rc.to_bits(), "n={n} m={m} tol={tol}");
                }
            }
        }
    }

    /// Clouds of every shape the scan meets, for the exactness tests of the histogram
    /// paths: tie-heavy, spread with planted diagonals, and Cartesian repeat blocks.
    fn test_clouds() -> Vec<Vec<Anchor>> {
        let mut s = 11u64;
        let mut out = Vec::new();
        for &n in &[97usize, 150, 400, 1000, 3000, 12000] {
            // Ties: a small box.
            out.push(
                (0..n)
                    .map(|_| Anchor {
                        a_pos: lcg(&mut s) % 300,
                        b_pos: lcg(&mut s) % 300,
                    })
                    .collect(),
            );
            // Spread noise over read-length boxes, with 0-2 diagonals planted in it.
            for (ext, diags) in [(5_000u32, 0usize), (40_000, 1), (60_000, 2)] {
                let mut v: Vec<Anchor> = (0..n)
                    .map(|_| Anchor {
                        a_pos: lcg(&mut s) % ext,
                        b_pos: lcg(&mut s) % ext,
                    })
                    .collect();
                for d in 0..diags {
                    let slope = [0.93, 1.21][d];
                    for i in 0..n / 3 {
                        let b = (lcg(&mut s) % (ext / 2)) + i as u32 % 7;
                        let a = (slope * b as f64) as u32 + 300 * d as u32 + lcg(&mut s) % 150;
                        v.push(Anchor { a_pos: a, b_pos: b });
                    }
                }
                out.push(v);
            }
            // A repeat block: few a positions against many b positions.
            let mut v = Vec::new();
            for i in 0..n {
                v.push(Anchor {
                    a_pos: 1000 + (i as u32 % 5) * 37,
                    b_pos: lcg(&mut s) % 20_000,
                });
            }
            out.push(v);
        }
        out
    }

    /// The histogram path of `diagonal_peak` and its `floor` are the whole sort's answer,
    /// bit for bit: `Some` exactly when the peak reaches the floor, with the same count and
    /// the same center.
    #[test]
    fn diagonal_peak_cells_and_floor_match_the_sort() {
        let (mut keys, mut cs) = (Vec::new(), Vec::new());
        for anchors in test_clouds() {
            let cloud = Cloud::new(&anchors);
            for &m in &[0.2, 0.7, 0.9, 0.98, 1.0, 1.1, 1.22, 1.3, 1.4, 2.3, 5.0] {
                for &tol in &[4.0, 96.0, 192.0] {
                    let offsets = anchors.iter().map(|a| a.a_pos as f64 - m * a.b_pos as f64);
                    let (rp, rc) = peak_sorted(offsets, 2.0 * tol, &mut keys, &mut cs);
                    let (p, c) = peak_at(&cloud, m, tol);
                    let n = anchors.len();
                    assert_eq!(p, rp, "n={n} m={m} tol={tol}");
                    assert_eq!(c.to_bits(), rc.to_bits(), "n={n} m={m} tol={tol}");
                    for floor in [rp.saturating_sub(3), rp - 1, rp, rp + 1, rp + 40] {
                        let r = diagonal_peak(&cloud, m, tol, floor);
                        if floor <= rp {
                            let (p, c) = r.expect("peak reaches the floor");
                            assert_eq!(p, rp, "n={n} m={m} tol={tol} floor={floor}");
                            assert_eq!(c.to_bits(), rc.to_bits(), "n={n} m={m} floor={floor}");
                        } else {
                            assert!(r.is_none(), "n={n} m={m} tol={tol} floor={floor}");
                        }
                    }
                }
            }
        }
    }

    /// `densest_sorted` is the slide it replaces, count and center bit for bit, on sorted
    /// offsets with ties, runs of equal values and every window width.
    #[test]
    fn densest_sorted_matches_the_slide() {
        fn slide(cs: &[f64], win: f64) -> (usize, f64) {
            let (mut i, mut peak, mut peak_c) = (0usize, 0usize, cs[0]);
            for j in 0..cs.len() {
                while cs[j] - cs[i] > win {
                    i += 1;
                }
                if j - i + 1 > peak {
                    peak = j - i + 1;
                    peak_c = (cs[i] + cs[j]) / 2.0;
                }
            }
            (peak, peak_c)
        }
        let mut s = 3u64;
        for &n in &[1usize, 2, 3, 7, 40, 300, 5000] {
            for &spread in &[3u32, 50, 1000, 100_000] {
                let mut cs: Vec<f64> = (0..n)
                    .map(|_| (lcg(&mut s) % spread) as f64 * 0.37 - 900.0)
                    .collect();
                cs.sort_by(f64::total_cmp);
                for &win in &[0.0, 0.37, 1.0, 8.0, 384.0] {
                    let (a, b) = (densest_sorted(&cs, win), slide(&cs, win));
                    assert_eq!(a.0, b.0, "n={n} spread={spread} win={win}");
                    assert_eq!(
                        a.1.to_bits(),
                        b.1.to_bits(),
                        "n={n} spread={spread} win={win}"
                    );
                }
            }
        }
    }

    /// The radix sort gives `sort_unstable`'s sequence at every size and spread.
    #[test]
    fn sort_u32_matches_sort_unstable() {
        let mut s = 5u64;
        for &(n, span) in &[
            (10usize, 7u32),
            (600, 3),
            (600, 5000),
            (5000, 1 << 20),
            (4000, u32::MAX),
        ] {
            let v: Vec<u32> = (0..n).map(|_| 17 + lcg(&mut s) % span.max(1)).collect();
            let mut want = v.clone();
            want.sort_unstable();
            let mut got = v;
            sort_u32(&mut got);
            assert_eq!(got, want, "n={n} span={span}");
        }
    }

    /// The slope sweep answers every slope as the whole sort does, bit for bit, along
    /// ascending slopes with gaps and after a fresh start.
    #[test]
    fn slope_sweep_matches_the_sort() {
        let (mut keys, mut cs) = (Vec::new(), Vec::new());
        let mut sw = SlopeSweep::default();
        for anchors in test_clouds() {
            let cloud = Cloud::new(&anchors);
            for (start, step) in [(0.7, 0.02), (0.9, 0.02), (1.1, 0.013), (0.2, 0.3)] {
                let mut m: f64 = start;
                sw.start(&cloud, m);
                for it in 0..24 {
                    if it % 5 != 3 {
                        let offsets = anchors.iter().map(|a| a.a_pos as f64 - m * a.b_pos as f64);
                        let (rp, rc) = peak_sorted(offsets, 2.0 * 96.0, &mut keys, &mut cs);
                        let (p, c) = sw.peak(m, 96.0);
                        assert_eq!(p, rp, "n={} m={m}", anchors.len());
                        assert_eq!(c.to_bits(), rc.to_bits(), "n={} m={m}", anchors.len());
                    }
                    m += step;
                }
            }
        }
    }

    /// The merge-sweep significance equals the sort-and-search definition in every term.
    #[test]
    fn significance_parts_match_the_sorted_form() {
        for anchors in test_clouds() {
            for &m in &[0.2, 0.71, 1.0, 1.137, 1.4, 5.0] {
                for &tol in &[1.0, 96.0, 192.0] {
                    let fast = diagonal_significance_parts(&anchors, m, tol).unwrap();
                    let slow = significance_parts_sorted(&anchors, m, tol.max(1.0));
                    let n = anchors.len();
                    assert_eq!(fast.obs, slow.obs, "n={n} m={m} tol={tol}");
                    for (x, y) in [
                        (fast.e, slow.e),
                        (fast.e_max, slow.e_max),
                        (fast.m_eff, slow.m_eff),
                        (fast.z, slow.z),
                    ] {
                        assert_eq!(x.to_bits(), y.to_bits(), "n={n} m={m} tol={tol}");
                    }
                }
            }
        }
    }

    /// The live first filter decides exactly what the screen statistic predicts: a pair is
    /// `Reject::Screen` iff its level-0 peak `ScreenStats::pk0` is under the floor, and a
    /// clean diagonal above it is still accepted.
    #[test]
    fn hough_screen_rejects_exactly_below_the_level0_peak() {
        let mut d = dp();
        d.slope_beam.screen_peak = 15;
        let tol = d.bin_events.max(1.0);
        let mut s = 21u64;
        let (mut screened, mut kept) = (0, 0);
        for trial in 0..60u32 {
            let n_diag = (trial % 30) as usize;
            let mut anchors: Vec<Anchor> = (0..n_diag)
                .map(|i| {
                    let b = 200 + i as u32 * 90;
                    Anchor {
                        // On the level-0 grid (1.1): at this toy `tol` (4 events) an
                        // off-grid diagonal smears across coarse windows.
                        a_pos: 700 + (1.1 * b as f64) as u32 + lcg(&mut s) % 3,
                        b_pos: b,
                    }
                })
                .collect();
            anchors.extend((0..(10 + trial as usize * 7)).map(|_| Anchor {
                a_pos: lcg(&mut s) % 4000,
                b_pos: lcg(&mut s) % 4000,
            }));
            let pk0 = screen_stats(&anchors, tol, &d.slope_beam).pk0;
            let aln = hough_aln(&anchors, &d, 0.0);
            assert_eq!(
                aln.reason == Reject::Screen,
                pk0 < 15,
                "trial {trial}: pk0 {pk0}, {aln:?}"
            );
            if aln.reason == Reject::Screen {
                screened += 1;
            } else {
                kept += 1;
            }
        }
        assert!(screened > 0 && kept > 0, "{screened} screened, {kept} kept");
        // Fewer anchors than the floor: screened without any scan.
        let few: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                a_pos: i * 50,
                b_pos: i * 50,
            })
            .collect();
        assert_eq!(hough_aln(&few, &d, 0.0).reason, Reject::Screen);
        // A clean diagonal well above the floor is untouched by it.
        let clean: Vec<Anchor> = (0..120).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        let mut off = d;
        off.slope_beam.screen_peak = 0;
        let (a, b) = (hough_aln(&clean, &d, 0.0), hough_aln(&clean, &off, 0.0));
        assert!(a.is_overlap());
        assert_eq!(
            (a.slope, a.aligned_bins, a.path.clone()),
            (b.slope, b.aligned_bins, b.path)
        );
    }

    /// The density-normalised filter decides exactly what the screen formula predicts from
    /// `ScreenStats` (`lam = n * 2*tol / (a_ext + b_ext)`), across clouds
    /// from bare scatter to a clean diagonal.
    #[test]
    fn hough_screen_z_matches_the_screen_formula() {
        let mut d = dp();
        d.slope_beam.screen_z = 5.0;
        let tol = d.bin_events.max(1.0);
        let mut s = 33u64;
        let (mut screened, mut kept) = (0, 0);
        for trial in 0..80u32 {
            let n_diag = (trial % 40) as usize;
            let mut anchors: Vec<Anchor> = (0..n_diag)
                .map(|i| {
                    let b = 100 + i as u32 * 60;
                    Anchor {
                        a_pos: 500 + (0.9 * b as f64) as u32 + lcg(&mut s) % 3,
                        b_pos: b,
                    }
                })
                .collect();
            // At least SLOPE_MIN_ANCHORS, below which no line is fitted and the screen
            // is never reached.
            anchors.extend(
                (0..(SLOPE_MIN_ANCHORS + trial as usize * 3)).map(|_| Anchor {
                    a_pos: lcg(&mut s) % 3000,
                    b_pos: lcg(&mut s) % 3000,
                }),
            );
            let st = screen_stats(&anchors, tol, &d.slope_beam);
            let lam = anchors.len() as f64 * 2.0 * tol / ((st.a_ext + st.b_ext).max(1) as f64);
            let want = (st.pk0 as f64 - lam) / lam.max(1.0).sqrt() < 5.0;
            let got = hough_aln(&anchors, &d, 0.0).reason == Reject::Screen;
            assert_eq!(got, want, "trial {trial}: pk0 {} lam {lam:.2}", st.pk0);
            if got {
                screened += 1;
            } else {
                kept += 1;
            }
        }
        assert!(screened > 0 && kept > 0, "{screened} screened, {kept} kept");
    }

    /// The stage-0b screen rejects exactly when the screen statistic `ScreenStats::zn0` (the
    /// gate's null `z` at the level-0 slope) is under the threshold.
    #[test]
    fn hough_screen_null_matches_the_screen_statistic() {
        let mut d = dp();
        d.slope_beam.screen_null = 2.0;
        let tol = d.bin_events.max(1.0);
        let mut s = 71u64;
        let (mut screened, mut kept) = (0, 0);
        for trial in 0..80u32 {
            let n_diag = (trial % 40) as usize;
            let mut anchors: Vec<Anchor> = (0..n_diag)
                .map(|i| {
                    let b = 100 + i as u32 * 60;
                    Anchor {
                        a_pos: 500 + (0.9 * b as f64) as u32 + lcg(&mut s) % 3,
                        b_pos: b,
                    }
                })
                .collect();
            anchors.extend(
                (0..(SLOPE_MIN_ANCHORS + trial as usize * 3)).map(|_| Anchor {
                    a_pos: lcg(&mut s) % 3000,
                    b_pos: lcg(&mut s) % 3000,
                }),
            );
            let st = screen_stats(&anchors, tol, &d.slope_beam);
            let want = st.zn0 < 2.0;
            let got = hough_aln(&anchors, &d, 0.0).reason == Reject::Screen;
            assert_eq!(got, want, "trial {trial}: zn0 {:.3}", st.zn0);
            if got {
                screened += 1;
            } else {
                kept += 1;
            }
        }
        assert!(screened > 0 && kept > 0, "{screened} screened, {kept} kept");
    }

    /// `ScreenStats::pkw` really bounds the peak at every fine slope the beam can visit,
    /// and `hb` bounds `pkw`, on clouds with an embedded diagonal at slopes on and between
    /// the coarse grid points.
    #[test]
    fn screen_stats_bound_the_fine_peak() {
        let beam = BeamCfg::default();
        let tol = 192.0;
        let mut s = 5u64;
        for &(slope, n_diag, n_bg) in &[
            (1.0, 40, 60),
            (0.83, 25, 400),
            (1.37, 60, 30),
            (0.7, 5, 2000),
        ] {
            let mut anchors: Vec<Anchor> = (0..n_diag)
                .map(|i| {
                    let b = 300 + i * 70;
                    Anchor {
                        a_pos: (500.0 + slope * b as f64) as u32 + lcg(&mut s) % 90,
                        b_pos: b,
                    }
                })
                .collect();
            anchors.extend((0..n_bg).map(|_| Anchor {
                a_pos: lcg(&mut s) % 9000,
                b_pos: lcg(&mut s) % 9000,
            }));
            let st = screen_stats(&anchors, tol, &beam);
            assert!(st.pk0 <= st.pkw && st.pkw <= st.hb, "{st:?}");
            let mut m = beam.win_lo as f64;
            while m <= beam.win_hi as f64 + 1e-9 {
                let (p, _) = peak_at(&Cloud::new(&anchors), m, tol);
                assert!(p as u32 <= st.pkw, "slope {m}: peak {p} > bound {}", st.pkw);
                m += 0.01;
            }
        }
    }

    /// The concentration estimator recovers an off-1 slope (0.8) even when it is
    /// half-buried in off-diagonal noise.
    #[test]
    fn fit_recovers_noisy_off1_slope() {
        let mut anchors: Vec<Anchor> = (0..300)
            .map(|i| Anchor {
                a_pos: (i as f64 * 0.8) as u32 + 5,
                b_pos: i,
            })
            .collect();
        let mut s = 999u64;
        anchors.extend((0..300).map(|_| Anchor {
            a_pos: lcg(&mut s) % 300,
            b_pos: lcg(&mut s) % 300,
        }));
        let slope = estimate_diagonal_fit(&anchors, 8, &BeamCfg::default(), 0.0)
            .unwrap()
            .slope;
        assert!(
            (slope - 0.8).abs() < 0.1,
            "should recover ~0.8 through noise, got {slope}"
        );
    }

    /// The confidence: pure scatter has no dominant diagonal, so its peak does not stand
    /// above the wrong-slope background (the preset `hough_min_conf` is 6).
    #[test]
    fn fit_confidence_rejects_pure_noise() {
        let mut s = 12345u64;
        let anchors: Vec<Anchor> = (0..500)
            .map(|_| Anchor {
                a_pos: lcg(&mut s) % 1000,
                b_pos: lcg(&mut s) % 1000,
            })
            .collect();
        let conf = estimate_diagonal_fit(&anchors, 8, &BeamCfg::default(), 0.0)
            .unwrap()
            .conf;
        assert!(
            conf < 6.0,
            "pure noise must not yield a confident diagonal ({conf})"
        );
    }

    /// The fit recovers both the off-1 slope AND the intercept through noise
    /// (a = 0.8·b + 5); pure noise does not clear the confidence gate.
    #[test]
    fn fit_recovers_slope_and_intercept() {
        let mut anchors: Vec<Anchor> = (0..300)
            .map(|i| Anchor {
                a_pos: (i as f64 * 0.8) as u32 + 5,
                b_pos: i,
            })
            .collect();
        let mut s = 999u64;
        anchors.extend((0..300).map(|_| Anchor {
            a_pos: lcg(&mut s) % 300,
            b_pos: lcg(&mut s) % 300,
        }));
        let fit = estimate_diagonal_fit(&anchors, 8, &BeamCfg::default(), 0.0).unwrap();
        let (slope, intercept) = (fit.slope, fit.intercept);
        assert!((slope - 0.8).abs() < 0.1, "slope ~0.8, got {slope}");
        assert!(
            (intercept - 5.0).abs() < 15.0,
            "intercept ~5, got {intercept}"
        );

        let mut s = 12345u64;
        let noise: Vec<Anchor> = (0..500)
            .map(|_| Anchor {
                a_pos: lcg(&mut s) % 1000,
                b_pos: lcg(&mut s) % 1000,
            })
            .collect();
        let conf = estimate_diagonal_fit(&noise, 8, &BeamCfg::default(), 0.0)
            .unwrap()
            .conf;
        assert!(
            conf < 6.0,
            "pure noise must yield no confident diagonal ({conf})"
        );
    }

    /// The Hough overlapper recovers the FULL extent of a diagonal even when the
    /// banded DP would stall at a long gap. A slope-1 run with a 20-event hole in
    /// the middle still yields one alignment spanning the whole extent, with a
    /// non-empty bin path.
    #[test]
    fn hough_recovers_diagonal_across_a_gap() {
        let mut anchors: Vec<Anchor> = (0..40).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        anchors.extend((60..120).map(|i| Anchor { a_pos: i, b_pos: i }));
        let aln = hough_aln(&anchors, &dp(), 3.0);
        assert!(aln.is_overlap(), "Hough should recover the gapped diagonal");
        // Extent spans the whole 0..120 events (~30 bins of 4), across the hole.
        assert!(aln.aligned_bins >= 25, "aligned_bins={}", aln.aligned_bins);
        assert!(!aln.path.is_empty());
        assert_eq!(aln.a_bin_span.0, 0);
    }

    /// The Hough significance gate rejects a diffuse (non-overlapping) anchor cloud:
    /// pure scatter has no dominant diagonal, so confidence stays low and nothing is
    /// emitted.
    #[test]
    fn hough_rejects_diffuse_cloud() {
        let mut s = 12345u64;
        let noise: Vec<Anchor> = (0..500)
            .map(|_| Anchor {
                a_pos: lcg(&mut s) % 1000,
                b_pos: lcg(&mut s) % 1000,
            })
            .collect();
        let aln = hough_aln(&noise, &dp(), 6.0);
        assert!(
            !aln.is_overlap(),
            "diffuse cloud must not pass the Hough gate"
        );
    }

    /// An off-1 (event-rate-mismatched) diagonal is recovered without any external
    /// slope correction: the Hough scan finds the slope itself and lays the on-line
    /// anchors onto the path.
    #[test]
    fn hough_recovers_off1_slope() {
        let mut anchors: Vec<Anchor> = (0..120)
            .map(|i| Anchor {
                a_pos: (i as f64 * 1.3) as u32,
                b_pos: i,
            })
            .collect();
        // Sprinkle off-diagonal noise the gate must ignore.
        let mut s = 7u64;
        anchors.extend((0..60).map(|_| Anchor {
            a_pos: lcg(&mut s) % 160,
            b_pos: lcg(&mut s) % 120,
        }));
        let aln = hough_aln(&anchors, &dp(), 3.0);
        assert!(
            aln.is_overlap(),
            "Hough should recover the slope-1.3 diagonal"
        );
        assert!(
            (aln.slope - 1.3).abs() < 0.15,
            "slope ~1.3, got {}",
            aln.slope
        );
    }

    /// The beam keeps a confidence gate, so it must reject a diffuse cloud.
    #[test]
    fn beam_rejects_diffuse_cloud() {
        let mut s = 12345u64;
        let noise: Vec<Anchor> = (0..500)
            .map(|_| Anchor {
                a_pos: lcg(&mut s) % 1000,
                b_pos: lcg(&mut s) % 1000,
            })
            .collect();
        let aln = hough_aln(&noise, &dp(), 6.0);
        assert!(
            !aln.is_overlap(),
            "coarse-to-fine must still reject a diffuse cloud"
        );
    }

    /// The prior window is enforced: when the densest diagonal is OUTSIDE `[win_lo,
    /// win_hi]` but a usable one sits inside, the beam picks the in-window diagonal —
    /// the one the downstream graph can actually use.
    #[test]
    fn beam_enforces_prior_window() {
        // Usable slope-1 diagonal (60 anchors) …
        let mut anchors: Vec<Anchor> = (0..60).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        // … plus a denser but unusable slope-3 diagonal (120 anchors), position-separated.
        anchors.extend((200..320).map(|b| Anchor {
            a_pos: 3 * b,
            b_pos: b,
        }));
        let beam = BeamCfg::default(); // window [0.7, 1.4]
        let s_beam = estimate_diagonal_fit(&anchors, 4, &beam, 0.0)
            .unwrap()
            .slope;
        assert!(
            (s_beam - 1.0).abs() < 0.1,
            "beam stays in the prior window on the usable diagonal, got {s_beam}"
        );
    }

    /// Top-k recovers a basin a single-winner coarse pass loses. The true diagonal
    /// (slope 1.0, 200 anchors) sits BETWEEN coarse grid points, so it is under-measured
    /// there (~80), while a weaker competitor (slope 1.3, 90 anchors) sits ON a grid
    /// point and out-ranks it at coarse resolution. `beam_k=1` refines only the
    /// competitor and returns 1.3; `beam_k>=3` also refines the true basin and recovers
    /// its real apex (200 at slope 1.0).
    #[test]
    fn beam_topk_recovers_true_basin() {
        let mut anchors: Vec<Anchor> = (0..200).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        anchors.extend((300..390).map(|b| Anchor {
            a_pos: (13 * b) / 10, // ≈ slope 1.3
            b_pos: b,
        }));
        let base = BeamCfg::default(); // coarse grid 0.7, 0.9, 1.1, 1.3
        let k1 = BeamCfg { beam_k: 1, ..base };
        let k3 = BeamCfg { beam_k: 3, ..base };
        let s1 = estimate_diagonal_fit(&anchors, 4, &k1, 0.0).unwrap().slope;
        let s3 = estimate_diagonal_fit(&anchors, 4, &k3, 0.0).unwrap().slope;
        assert!(
            (s1 - 1.3).abs() < 0.1,
            "k=1 loses the true basin to the on-grid competitor, got {s1}"
        );
        assert!(
            (s3 - 1.0).abs() < 0.1,
            "k=3 recovers the true dominant diagonal, got {s3}"
        );
    }

    /// A clean diagonal is highly significant; pure scatter is not — the core
    /// contract of the null gate (keep strong diagonals, reject diffuse clouds).
    #[test]
    fn significance_separates_diagonal_from_noise() {
        let diag: Vec<Anchor> = (0..120).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        let z_diag = diagonal_significance(&diag, 1.0, 4.0);
        assert!(
            z_diag > 8.0,
            "clean diagonal should be significant, z={z_diag}"
        );

        let mut s = 12345u64;
        let noise: Vec<Anchor> = (0..500)
            .map(|_| Anchor {
                a_pos: lcg(&mut s) % 1000,
                b_pos: lcg(&mut s) % 1000,
            })
            .collect();
        let z_noise = diagonal_significance(&noise, 1.0, 8.0);
        assert!(
            z_noise < 3.0,
            "diffuse cloud must be insignificant, z={z_noise}"
        );
    }

    /// Hash-frequency awareness: a low-entropy REPEAT GRID (a few a/b positions
    /// producing a Cartesian block of anchors) must NOT look significant, even
    /// though many anchors share the c=0 offset. A uniform-box null would be fooled;
    /// the marginal null discounts it because random pairing reconstructs the grid.
    #[test]
    fn significance_discounts_repeat_grid() {
        // 5 a-positions × 5 b-positions on a 40-event pitch: a 25-anchor grid whose
        // diagonal (a==b) holds 5, but that is exactly what random pairing predicts.
        let mut grid: Vec<Anchor> = Vec::new();
        for ai in 0..5u32 {
            for bi in 0..5u32 {
                grid.push(Anchor {
                    a_pos: ai * 40,
                    b_pos: bi * 40,
                });
            }
        }
        let z_grid = diagonal_significance(&grid, 1.0, 8.0);
        assert!(
            z_grid < 3.0,
            "a repeat grid must not pass as a diagonal, z={z_grid}"
        );

        // A TRUE diagonal with the same number of anchors and the same extent IS
        // significant — showing the gate rejects the grid on structure, not size.
        let diag: Vec<Anchor> = (0..25u32)
            .map(|i| Anchor {
                a_pos: i * 8,
                b_pos: i * 8,
            })
            .collect();
        let z_diag = diagonal_significance(&diag, 1.0, 8.0);
        assert!(z_diag > 3.0, "a true diagonal should pass, z={z_diag}");
    }

    /// The mode-anchored branch must reproduce the direct walk wherever both are
    /// usable. They are compared below the `NULL_MAX_DIRECT_E` cutoff, across the
    /// `m_eff` range real clouds produce, so the branch taken for dense clouds is
    /// pinned to the one that has always been in production.
    #[test]
    fn dense_null_max_matches_direct_walk() {
        for &e in &[1.0, 5.0, 20.0, 75.0, 150.0, 300.0, 499.0] {
            for &m_eff in &[1.5, 4.0, 30.0, 500.0, 20_000.0] {
                let direct = expected_null_max(e, m_eff);
                let dense = expected_null_max_dense(e, m_eff);
                assert!(
                    (direct - dense).abs() < 1e-6,
                    "e={e} m_eff={m_eff}: direct={direct} dense={dense}"
                );
            }
        }
    }

    /// Regression: past `e ≈ 745`, `exp(-e)` is exactly 0 in f64, so the direct walk's
    /// upper tail never decreased and it ran to its 10^6 iteration cap. `E_max` then
    /// came back as ~10^6 against an `obs` of a few thousand, driving `z` to ~-33000
    /// and silently rejecting the pair. Dense repeat clouds are precisely where `e`
    /// gets that large, so the failure selected against the clouds the gate exists
    /// to judge.
    #[test]
    fn null_max_survives_large_e() {
        for &e in &[600.0, 900.0, 5_000.0, 250_000.0] {
            let em = expected_null_max(e, 100.0);
            // The Poisson is near-Gaussian here, so the 1/m_eff point sits a few
            // sd above the mean — never below it, and never orders of magnitude off.
            let sd = e.sqrt();
            assert!(
                em > e && em < e + 8.0 * sd,
                "e={e}: E_max={em} outside [{e}, {}]",
                e + 8.0 * sd
            );
        }
    }

    /// Validate the analytic marginal null against a seeded permutation null: the
    /// analytic `E` must track the mean on-band count when b-positions are shuffled.
    #[test]
    fn analytic_null_matches_permutation() {
        // A noisy off-1 diagonal plus scatter — a realistic cloud.
        let mut anchors: Vec<Anchor> = (0..200)
            .map(|i| Anchor {
                a_pos: (i as f64 * 0.9) as u32 + 3,
                b_pos: i,
            })
            .collect();
        let mut s = 4242u64;
        anchors.extend((0..200).map(|_| Anchor {
            a_pos: lcg(&mut s) % 200,
            b_pos: lcg(&mut s) % 200,
        }));
        let slope = 0.9;
        let tol = 6.0;
        let win = 2.0 * tol;

        // Analytic E, recovered from the z formula: z = (obs - E)/sqrt(max(E,1)).
        // Recompute obs (densest window) and back out E.
        let mut cs: Vec<f64> = anchors
            .iter()
            .map(|a| a.a_pos as f64 - slope * a.b_pos as f64)
            .collect();
        cs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let (mut i, mut obs, mut clo, mut chi) = (0usize, 0usize, cs[0], cs[0]);
        for j in 0..cs.len() {
            while cs[j] - cs[i] > win {
                i += 1;
            }
            if j - i + 1 > obs {
                obs = j - i + 1;
                clo = cs[i];
                chi = cs[j];
            }
        }
        let c0 = (clo + chi) / 2.0;
        let z = diagonal_significance(&anchors, slope, tol);
        let _ = obs; // obs is used to locate c0; E is reconstructed directly below.
        // Direct analytic E (same formula the function uses):
        let sb: Vec<f64> = {
            let mut v: Vec<f64> = anchors.iter().map(|a| slope * a.b_pos as f64).collect();
            v.sort_by(|x, y| x.partial_cmp(y).unwrap());
            v
        };
        let mut m_pairs = 0u64;
        for a in &anchors {
            let center = a.a_pos as f64 - c0;
            let l = sb.partition_point(|&v| v < center - tol);
            let r = sb.partition_point(|&v| v <= center + tol);
            m_pairs += (r - l) as u64;
        }
        let e_direct = m_pairs as f64 / anchors.len() as f64;

        // Empirical permutation null: shuffle b-positions, count best-window density.
        let n = anchors.len();
        let mut perm_sum = 0.0f64;
        let k = 200;
        let mut rng = 987654321u64;
        let a_vals: Vec<f64> = anchors.iter().map(|a| a.a_pos as f64).collect();
        let b_vals: Vec<f64> = anchors.iter().map(|a| slope * a.b_pos as f64).collect();
        for _ in 0..k {
            // Fisher-Yates shuffle of an index permutation of b_vals.
            let mut perm: Vec<usize> = (0..n).collect();
            for idx in (1..n).rev() {
                let j = (lcg(&mut rng) as usize) % (idx + 1);
                perm.swap(idx, j);
            }
            // On-band count at the SAME window [c0-tol, c0+tol].
            let mut cnt = 0usize;
            for (ii, &av) in a_vals.iter().enumerate() {
                let off = av - b_vals[perm[ii]];
                if (off - c0).abs() <= tol {
                    cnt += 1;
                }
            }
            perm_sum += cnt as f64;
        }
        let e_perm = perm_sum / k as f64;

        // Analytic E and permutation-mean E must agree to a few percent.
        assert!(
            (e_direct - e_perm).abs() <= 0.15 * e_perm.max(1.0) + 2.0,
            "analytic E={e_direct:.2} vs permutation E={e_perm:.2} (z={z:.2})"
        );
        // And the diagonal must still read as significant.
        assert!(
            z > 5.0,
            "noisy real diagonal should stay significant, z={z}"
        );
    }

    /// `tail_nlp` against the same probability summed in 80-digit decimal arithmetic.
    /// The cases straddle the branches: `e < 1` (where `z`'s clamp bites), `e ≈ obs`
    /// (no significance at all), and a tail far past f64 underflow.
    #[test]
    fn tail_nlp_matches_exact_arithmetic() {
        for &(obs, e, m_eff, want) in &[
            (10u32, 1.0f64, 100.0f64, 4.953018f64),
            (5, 2.0, 50.0, 0.030074),
            (30, 3.5, 1000.0, 14.569837),
            (3, 0.5, 10.0, 0.869956),
            (100, 20.0, 5000.0, 32.758344),
            (8, 8.0, 2.0, 0.099728),
        ] {
            let got = tail_nlp(obs, e, m_eff);
            assert!(
                (got - want).abs() < 1e-5,
                "tail_nlp({obs}, {e}, {m_eff}) = {got}, want {want}"
            );
        }
    }

    /// `tail_nlp_upper` is never below `tail_nlp` (beyond rounding), from `e` far below `obs`
    /// through `e` far above it, one window to thousands, up to underflowing tails.
    #[test]
    fn tail_nlp_upper_bounds_tail_nlp() {
        let es: Vec<f64> = (0..48).map(|i| 1e-3 * 1.3f64.powi(i)).collect();
        let obs: Vec<u32> = (0..=60)
            .chain([75, 90, 120, 200, 350, 600, 1000, 2500])
            .collect();
        let mut n = 0;
        for &o in &obs {
            for &e in &es {
                for &m in &[1.0f64, 2.0, 7.5, 80.0, 1000.0, 50_000.0] {
                    let (t, u) = (tail_nlp(o, e, m), tail_nlp_upper(o, e, m));
                    assert!(
                        t <= u + 1e-9 * u.max(1.0),
                        "tail_nlp({o}, {e}, {m}) = {t} > upper {u}"
                    );
                    n += (u.is_finite() && u > 0.0) as usize;
                }
            }
        }
        assert!(n > 5_000, "only {n} finite positive bounds");
        assert_eq!(tail_nlp_upper(0, 3.0, 10.0), 0.0);
        assert!(tail_nlp_upper(3000, 1.0, 10.0).is_infinite());
    }

    /// The property the statistic exists for: unlike `z`, it does not drift with the
    /// number of windows scanned. Doubling `m_eff` costs a fixed log10(2), never more,
    /// and a longer alignment at proportionally more anchors keeps the same tail.
    #[test]
    fn tail_nlp_is_a_probability_not_a_count() {
        // A pure multiple-testing change moves it by exactly the union-bound amount.
        let a = tail_nlp(40, 2.0, 1000.0);
        let b = tail_nlp(40, 2.0, 2000.0);
        assert!((a - b - (2.0f64).log10()).abs() < 1e-6, "{a} vs {b}");
        // Monotone in the observation, and zero when the observation is unremarkable.
        assert!(tail_nlp(41, 2.0, 1000.0) > a);
        assert_eq!(tail_nlp(0, 2.0, 1000.0), 0.0);
        assert!(tail_nlp(1, 50.0, 10.0) < 1e-9);
        // Underflowing tails stay finite and keep growing.
        let deep = tail_nlp(400, 1.0, 1000.0);
        assert!(deep.is_finite() && deep > 300.0, "{deep}");
    }

    /// DIAGNOSTIC PROBE (calibration): how does `diagonal_significance` behave on
    /// PURE-NULL clouds (random scatter, no diagonal) as anchor count N and event
    /// range R vary? A well-calibrated significance gives z≈0 for every null cloud
    /// regardless of N/R. Any growth of z with N or R is selection-bias inflation
    /// (`obs` is a MAX over sliding windows; `E` is the single-window mean at the
    /// observed peak c0), which makes a fixed z threshold non-comparable across read
    /// pairs of different length / seed density. Run:
    ///   cargo test -p haricot null_gate_calibration_probe -- --nocapture
    #[test]
    fn null_gate_calibration_probe() {
        let tol = 192.0; // production value: bin_events
        let slope = 1.0;
        let seeds = 400usize;
        println!(
            "\n  {:>5} {:>7} | {:>10} {:>8} {:>8} {:>8} | {:>7} {:>7}",
            "N", "R", "z_mean", "z_p50", "z_p95", "z_max", "E_mean", "obs_mn"
        );
        for &r in &[2000u32, 10000, 40000, 100000] {
            for &n in &[16usize, 64, 256, 1024] {
                if n as u32 > r {
                    continue;
                }
                let mut zs = Vec::with_capacity(seeds);
                let mut s = 0x1234_5678u64 ^ ((n as u64) << 20) ^ ((r as u64) << 3);
                for _ in 0..seeds {
                    let anchors: Vec<Anchor> = (0..n)
                        .map(|_| Anchor {
                            a_pos: lcg(&mut s) % r,
                            b_pos: lcg(&mut s) % r,
                        })
                        .collect();
                    zs.push(diagonal_significance(&anchors, slope, tol));
                }
                zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let mean = zs.iter().sum::<f64>() / zs.len() as f64;
                let pct = |f: f64| zs[((zs.len() as f64 * f) as usize).min(zs.len() - 1)];
                println!(
                    "  {:>5} {:>7} | {:>10.2} {:>8.2} {:>8.2} {:>8.2}",
                    n,
                    r,
                    mean,
                    pct(0.5),
                    pct(0.95),
                    *zs.last().unwrap()
                );
                // Calibration guard: after the extreme-value correction the pure-null
                // score must sit near 0 for EVERY (N, R) — no drift with read length
                // or seed density. (Pre-correction this p95 ran to 5.3 and grew with R.)
                assert!(
                    pct(0.95) < 3.0 && mean.abs() < 1.5,
                    "null z not calibrated at N={n}, R={r}: mean={mean:.2}, p95={:.2}",
                    pct(0.95)
                );
            }
        }
        // Contrast: a genuine planted diagonal (half on-line, half scatter) must read
        // as strongly significant at every scale.
        for &r in &[5000u32, 20000, 80000] {
            let mut s = 7u64;
            let step = r / 200;
            let mut anchors: Vec<Anchor> = (0..200)
                .map(|i| Anchor {
                    a_pos: i * step,
                    b_pos: i * step,
                })
                .collect();
            anchors.extend((0..200).map(|_| Anchor {
                a_pos: lcg(&mut s) % r,
                b_pos: lcg(&mut s) % r,
            }));
            let z = diagonal_significance(&anchors, 1.0, tol);
            println!("  planted diagonal R={r:>6} (200 on + 200 scatter): z = {z:.1}");
            // A real diagonal must stay far above the null band at every scale.
            assert!(z > 20.0, "planted diagonal under-scored at R={r}: z={z:.1}");
        }
    }

    /// Path cells are unique and sorted (the FBG slicer keys on them).
    #[test]
    fn hough_path_is_sorted_and_unique() {
        let anchors: Vec<Anchor> = (0..120).map(|i| Anchor { a_pos: i, b_pos: i }).collect();
        let aln = hough_aln(&anchors, &dp(), 3.0);
        assert!(aln.is_overlap());
        for w in aln.path.windows(2) {
            assert!(w[0] < w[1], "path must be strictly increasing & unique");
        }
    }

    /// The confidence shortcut (skip the background scan when `peak < min_conf`) is
    /// sound **only because** `peak_conf` floors the median at 1, making `conf <= peak`.
    /// A median below 1 would put `conf` *above* the peak and the shortcut would begin
    /// discarding pairs that should pass. That floor is load-bearing for correctness,
    /// not merely a divide-by-zero guard — pin it here, including the all-zero
    /// background that is the only way to reach it.
    #[test]
    fn conf_never_exceeds_the_peak() {
        let peak = 6usize;
        for peaks in [
            vec![0usize; 17],          // background never fires -> median floored to 1
            vec![0, 0, 0, 0, 0, 1, 2], // median is 0 before the floor
            vec![1usize; 17],
            vec![5, 7, 9, 11, 13, 2, 3],
        ] {
            let conf = peak_conf(peak, peaks.clone());
            assert!(
                conf <= peak as f64,
                "conf {conf} exceeded peak {peak} for background {peaks:?} — the \
                 min_conf shortcut in `beam_search_slope` is unsound without conf <= peak"
            );
        }
    }

    /// With the gate OFF (`min_conf <= 0`) every pair skips the confidence background.
    /// That must change no accept decision and no field of an accepted result except
    /// `conf` itself, which becomes the peak bound (never below the exact ratio).
    #[test]
    fn gate_off_skips_background_and_changes_only_conf() {
        let mut s = 0x9E37_79B9_7F4A_7C15u64;
        let mut rnd = |m: u32| -> u32 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as u32 % m.max(1)
        };
        let (mut accepted, mut differed) = (0, 0);
        for on_diag in [0usize, 5, 12, 40, 160] {
            for scatter in [0usize, 50, 200] {
                let mut anchors: Vec<Anchor> = (0..on_diag)
                    .map(|i| Anchor {
                        a_pos: (i as u32) * 3,
                        b_pos: (i as u32) * 3,
                    })
                    .collect();
                anchors.extend((0..scatter).map(|_| Anchor {
                    a_pos: rnd(1000),
                    b_pos: rnd(1000),
                }));
                let fast = dp();
                let mut exact = fast;
                exact.exact_conf = true;
                let a = hough_aln(&anchors, &fast, 0.0);
                let b = hough_aln(&anchors, &exact, 0.0);
                assert_eq!(a.is_overlap(), b.is_overlap());
                assert_eq!(a.reason, b.reason);
                assert_eq!(a.path, b.path);
                assert_eq!(a.a_bin_span, b.a_bin_span);
                assert_eq!(a.b_bin_span, b.b_bin_span);
                assert_eq!(a.matched, b.matched);
                assert_eq!(a.cnt, b.cnt);
                assert_eq!(a.aligned_bins, b.aligned_bins);
                assert_eq!(a.gap, b.gap);
                assert_eq!(a.slope.to_bits(), b.slope.to_bits());
                assert_eq!(a.intercept.to_bits(), b.intercept.to_bits());
                assert!(a.conf >= b.conf, "bound {} below exact {}", a.conf, b.conf);
                accepted += a.is_overlap() as usize;
                differed += (a.conf.to_bits() != b.conf.to_bits()) as usize;
            }
        }
        assert!(accepted > 0, "nothing accepted — test would be vacuous");
        assert!(
            differed > 0,
            "background never skipped — test would be vacuous"
        );
    }

    /// `BinDp::exact_conf` must be a pure performance switch: skipping the confidence
    /// background for pairs whose peak already fails the gate cannot change which pairs
    /// `hough_aln` accepts, nor any field it reports about an accepted one.
    ///
    /// Sweeps planted diagonals of varying density mixed with scatter so both branches
    /// are exercised — asserted at the end, since a run where the shortcut never fired
    /// would pass vacuously.
    #[test]
    fn conf_shortcut_changes_no_accept_decision() {
        let mut s = 0x243F_6A88_85A3_08D3u64;
        let mut rnd = |m: u32| -> u32 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as u32 % m.max(1)
        };
        let (mut skipped, mut accepted) = (0, 0);
        for on_diag in [0usize, 2, 5, 8, 12, 20, 40, 80, 160] {
            for scatter in [0usize, 10, 50, 200] {
                let mut anchors: Vec<Anchor> = (0..on_diag)
                    .map(|i| Anchor {
                        a_pos: (i as u32) * 3,
                        b_pos: (i as u32) * 3,
                    })
                    .collect();
                anchors.extend((0..scatter).map(|_| Anchor {
                    a_pos: rnd(1000),
                    b_pos: rnd(1000),
                }));

                let fast = dp(); // exact_conf = false (production)
                let mut exact = fast;
                exact.exact_conf = true;
                let min_conf = 6.0;

                let a = hough_aln(&anchors, &fast, min_conf);
                let b = hough_aln(&anchors, &exact, min_conf);
                assert_eq!(
                    a.is_overlap(),
                    b.is_overlap(),
                    "accept decision diverged at on_diag={on_diag} scatter={scatter}"
                );
                if a.is_overlap() {
                    accepted += 1;
                    // Accepted pairs always evaluate the background, so every field —
                    // `conf` included — must match bit for bit.
                    assert_eq!(a.matched, b.matched);
                    assert_eq!(a.cnt, b.cnt);
                    assert_eq!(a.aligned_bins, b.aligned_bins);
                    assert_eq!(a.gap, b.gap);
                    assert_eq!(a.path, b.path);
                    assert_eq!(a.a_bin_span, b.a_bin_span);
                    assert_eq!(a.b_bin_span, b.b_bin_span);
                    assert_eq!(a.slope.to_bits(), b.slope.to_bits());
                    assert_eq!(a.conf.to_bits(), b.conf.to_bits());
                } else if a.conf.to_bits() != b.conf.to_bits() {
                    // The one permitted difference: a rejected pair may report the
                    // upper bound `peak` in place of the exact ratio. It must still be
                    // an upper bound, and must still fail the gate.
                    skipped += 1;
                    assert!(
                        a.conf >= b.conf,
                        "shortcut conf {} below exact {}",
                        a.conf,
                        b.conf
                    );
                    assert!(
                        a.conf < min_conf,
                        "shortcut conf {} cleared the gate",
                        a.conf
                    );
                }
            }
        }
        assert!(skipped > 0, "shortcut never fired — test would be vacuous");
        assert!(
            accepted > 0,
            "nothing was ever accepted — test would be vacuous"
        );
    }

    /// The detection settings of `haricot assemble`, with `select_z` as given.
    fn production_dp(select_z: bool) -> BinDp {
        let mut op = crate::params::Params::r10_4_1().overlap;
        op.hough = true;
        op.hough_min_conf = 0.0;
        op.hough_min_z = 3.0;
        op.null_gate = true;
        op.null_z = 4.0;
        op.min_aln_bins = 4;
        op.max_seed_occ = 2000;
        op.bin_events = 192;
        op.slope_beam.select_z = select_z;
        BinDp::from_params(&op, 8)
    }

    fn xorshift(mut s: u64) -> impl FnMut() -> u64 {
        move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        }
    }

    fn mirror(cloud: &[Anchor]) -> Vec<Anchor> {
        cloud
            .iter()
            .map(|x| Anchor {
                a_pos: x.b_pos,
                b_pos: x.a_pos,
            })
            .collect()
    }

    /// A true overlap `a = b + 5000`
    /// plus a dense repeat block (a in [0, 200) x b in [16000, 20000)). The block's raw window
    /// count wins at the window edge, where `z` is ~0; with `select_z` the slope is chosen by
    /// `z`, so the diagonal (z > 20) is found in BOTH read orders, at slope ~1.
    #[test]
    fn select_z_finds_a_diagonal_an_edge_basin_outcounts() {
        let mut r = xorshift(0x2545_f491_4f6c_dd1d);
        let mut cloud: Vec<Anchor> = (0..15_000u32)
            .step_by(40)
            .map(|b| Anchor {
                a_pos: b + 5000,
                b_pos: b,
            })
            .collect();
        for _ in 0..3000 {
            cloud.push(Anchor {
                a_pos: (r() % 200) as u32,
                b_pos: 16_000 + (r() % 4000) as u32,
            });
        }
        assert!(diagonal_significance(&cloud, 1.0, 192.0) > 20.0);
        let d = production_dp(true);
        for (name, c) in [("as given", cloud.clone()), ("mirrored", mirror(&cloud))] {
            let x = hough_aln(&c, &d, 0.0);
            assert!(
                x.is_overlap() && (x.slope - 1.0).abs() < 0.05,
                "{name}: {:?} at slope {:.4}",
                x.reason,
                x.slope
            );
        }
    }

    /// With one clear diagonal (and uniform noise) the z pick is the count winner: `select_z`
    /// changes nothing, the fitted line and the verdict are identical.
    #[test]
    fn select_z_keeps_a_single_diagonal_fit() {
        let mut r = xorshift(0x9e37_79b9_7f4a_7c15);
        for trial in 0..20u64 {
            let slope = 0.85 + 0.015 * trial as f64;
            let mut cloud: Vec<Anchor> = (0..8000u32)
                .step_by(30)
                .map(|b| Anchor {
                    a_pos: (slope * b as f64) as u32 + 3000,
                    b_pos: b,
                })
                .collect();
            for _ in 0..200 {
                cloud.push(Anchor {
                    a_pos: (r() % 12_000) as u32,
                    b_pos: (r() % 8000) as u32,
                });
            }
            let (x, y) = (
                hough_aln(&cloud, &production_dp(false), 0.0),
                hough_aln(&cloud, &production_dp(true), 0.0),
            );
            assert!(x.is_overlap(), "trial {trial}: {:?}", x.reason);
            assert_eq!(x.slope.to_bits(), y.slope.to_bits(), "trial {trial}");
            assert_eq!(
                x.intercept.to_bits(),
                y.intercept.to_bits(),
                "trial {trial}"
            );
            assert_eq!(x.path, y.path, "trial {trial}");
        }
    }

    /// Every field of an alignment, floats by their bits.
    fn aln_key(a: &BinAln) -> String {
        let null = a.null.map(|p| {
            (
                p.obs,
                p.e.to_bits(),
                p.e_max.to_bits(),
                p.m_eff.to_bits(),
                p.z.to_bits(),
            )
        });
        format!(
            "{:?} {:?} {:?} {:?} {:?}",
            (a.matched, a.cnt, a.aligned_bins, a.gap, a.reason),
            (a.a_bin_span, a.b_bin_span, a.similarity.to_bits()),
            (a.slope.to_bits(), a.intercept.to_bits(), a.conf.to_bits()),
            null,
            a.path
        )
    }

    /// `peaks_at` is `peak_at` at every slope, and `best_peak` is what measuring the slopes one
    /// by one in their order at the best count so far finds (the beam's last level before the
    /// batched counts) -- with every slope list shape the beam builds, under both kernels.
    #[test]
    fn best_peak_and_peaks_at_are_the_slope_by_slope_answers() {
        let (mut keys, mut cs) = (Vec::new(), Vec::new());
        let mut grid = Vec::new();
        let mut m = 0.7f64;
        while m <= 1.4 + 1e-9 {
            grid.push(m);
            m += 0.02;
        }
        let lists = [
            grid,
            vec![1.3, 0.75, 1.0000000000000002, 1.0, 0.9],
            vec![0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 5.0],
            vec![1.1],
        ];
        for anchors in test_clouds() {
            let cloud = Cloud::new(&anchors);
            let n = anchors.len();
            for ms in &lists {
                for avx in [false, true] {
                    kernel::force(Some(avx));
                    for &tol in &[96.0, 192.0] {
                        let all = peaks_at(&cloud, ms, tol);
                        assert_eq!(all.len(), ms.len());
                        for (&(m, p, c), &want_m) in all.iter().zip(ms) {
                            let offsets =
                                anchors.iter().map(|a| a.a_pos as f64 - m * a.b_pos as f64);
                            let (rp, rc) = peak_sorted(offsets, 2.0 * tol, &mut keys, &mut cs);
                            assert_eq!(m.to_bits(), want_m.to_bits());
                            assert_eq!(p, rp, "n={n} m={m} tol={tol} avx={avx}");
                            assert_eq!(c.to_bits(), rc.to_bits(), "n={n} m={m} avx={avx}");
                        }
                        let mut seq: Option<(f64, usize, f64)> = None;
                        for &m in ms {
                            if let Some((p, c)) =
                                diagonal_peak(&cloud, m, tol, seq.map_or(0, |t| t.1))
                                && seq.is_none_or(|t| {
                                    p > t.1 || (p == t.1 && m.total_cmp(&t.0).is_lt())
                                })
                            {
                                seq = Some((m, p, c));
                            }
                        }
                        let bits = |t: Option<(f64, usize, f64)>| {
                            t.map(|(m, p, c)| (m.to_bits(), p, c.to_bits()))
                        };
                        assert_eq!(
                            bits(best_peak(&cloud, ms, tol)),
                            bits(seq),
                            "n={n} tol={tol} avx={avx} slopes={}",
                            ms.len()
                        );
                    }
                }
            }
        }
        kernel::force(None);
    }

    /// `hough_aln` answers bit for bit alike under the baseline and the AVX2 kernels, and with
    /// small clouds swept or counted, on clouds of every size the detector
    /// meets, with and without a diagonal, with and without `select_z`.
    #[test]
    fn hough_aln_is_the_same_under_every_kernel_and_path() {
        let (d, dz) = (production_dp(false), production_dp(true));
        let mut clouds = test_clouds();
        let mut s = 31u64;
        for &n in &[20usize, 90, 200, 300, 450, 700, 1500, 3000] {
            for diag in [false, true] {
                let ext = 15_000u32;
                let mut v: Vec<Anchor> = (0..n)
                    .map(|_| Anchor {
                        a_pos: lcg(&mut s) % ext,
                        b_pos: lcg(&mut s) % ext,
                    })
                    .collect();
                if diag {
                    for _ in 0..n / 4 {
                        let b = lcg(&mut s) % 9_000;
                        let a = (1.08 * b as f64) as u32 + 2_000 + lcg(&mut s) % 60;
                        v.push(Anchor { a_pos: a, b_pos: b });
                    }
                }
                // The gather hands each cloud over in the order of read `a`'s seeds.
                v.sort_by_key(|x| x.a_pos);
                clouds.push(v);
            }
        }
        for (i, anchors) in clouds.iter().enumerate() {
            let mut first: Option<(String, String)> = None;
            for avx in [false, true] {
                for no_sweep in [false, true] {
                    kernel::force(Some(avx));
                    NO_SWEEP.with(|c| c.set(Some(no_sweep)));
                    let got = (
                        aln_key(&hough_aln(anchors, &d, 0.0)),
                        aln_key(&hough_aln(anchors, &dz, 0.0)),
                    );
                    match &first {
                        None => first = Some(got),
                        Some(f) => assert_eq!(
                            f,
                            &got,
                            "cloud {i} (n={}) avx={avx} no_sweep={no_sweep}",
                            anchors.len()
                        ),
                    }
                }
            }
        }
        kernel::force(None);
        NO_SWEEP.with(|c| c.set(None));
    }
}
