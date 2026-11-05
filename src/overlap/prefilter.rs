//! Partner nomination before alignment: the partner pre-filter (`haricot-detect --production` runs
//! it at level 2).
//!
//! **Why.** [`super::align_pair`] (the Hough slope search) takes ~96% of detection without it, and
//! rejects 99.7-99.9% of the pairs it is given. Every pair the per-read gather
//! produces is scored here first, with a statistic that costs a pass over the pair's anchors and
//! a window scan; only the pairs that reach a stated chance level go on to `align_pair`, which
//! then runs unchanged on the pair's full cloud.
//!
//! **The statistic** is `nlp2n_1`, the significance of the densest slope-1 diagonal window. Its
//! constants:
//! * the pair's anchors are every key that survived the repeat mask (the detector's own cloud;
//!   no sampling, no weights);
//! * each anchor is binned by its slope-1 offset `k = a_pos - b_pos + L_b` in bins of [`DBIN`] =
//!   512 events (`L_x` = read x's length in retained events, the coordinate space);
//! * a window is two adjacent bins (the last bin alone at the end); its expected count is the
//!   pair's anchor total times the window's share of the `L_a x L_b` rectangle (exact trapezoid
//!   cell counts, [`cells_in`]);
//! * the window with the largest `z = (obs - e) / sqrt(max(e, 1))` is kept (first one on ties),
//!   and the score is `-log10` of its Poisson tail with the number of windows as multiplicity
//!   ([`crate::overlap::dp::tail_nlp`]).
//!
//! An overlap with slope in the beam window [0.7, 1.4] whose offset drifts by at most 512 events
//! keeps its anchors in one window; the chance distribution of this score was measured nearly
//! identical across four regimes (p99 2.8-3.3).
//!
//! **Ownership and the rule.** A pair is scored once, by its owner (the LOWER read id, as
//! [`super::find_overlap_candidates`] owns pairs), in the owner frame: anchors `(a_pos, b_pos)`,
//! `L_a`, `L_b`. It is nominated when the score reaches EITHER read's level. With the fixed level
//! of `--prefilter-nlp L` both reads' levels are `L`, so the rule is the owner's single test;
//! per-read levels ([`Nominator::per_read`]) exist for a per-read null and for the tests.
//!
//! **The level** is a stated chance level, not a fitted one: `L = 2` nominates a pair whose
//! densest diagonal window has p <= 1e-2 under the uniform null after the window multiplicity.
//!
//! **Exactness.** Off (`OverlapParams::prefilter_nlp = None`) the detector is unchanged. On, the
//! candidate set is exactly the off set restricted to the nominated pairs, because `align_pair` is
//! a pure function of a pair's anchors and the cloud is the same gather either way.

use super::Anchor;
use crate::overlap::dp;

/// Diagonal bin width (events) of the nomination statistic; a window is two adjacent bins.
pub const DBIN: i64 = 512;

/// Cells `(x, y)` of `[0, l1) x [0, l2)` with `x - y + l2` in `[1, kk]`: the cumulative sum of
/// the trapezoid `c(k) = min(k, a, l1 + l2 - k)`, `a = min(l1, l2)`, for k in `[1, l1 + l2 - 1]`
/// (exact integer arithmetic).
fn cum_cells(kk: i64, l1: i64, l2: i64) -> f64 {
    let (a, b) = (l1.min(l2), l1.max(l2));
    if kk <= 0 || a <= 0 {
        return 0.0;
    }
    let t = l1 + l2;
    let k = kk.min(t - 1);
    let v = if k <= a {
        k * (k + 1) / 2
    } else if k <= b {
        a * (a + 1) / 2 + (k - a) * a
    } else {
        let r = t - k;
        a * (a + 1) / 2 + (b - a) * a + (a - 1) * a / 2 - (r - 1) * r / 2
    };
    v as f64
}

/// Cells with `x - y + l2` in `[k0, k1)`.
fn cells_in(k0: i64, k1: i64, l1: i64, l2: i64) -> f64 {
    cum_cells(k1 - 1, l1, l2) - cum_cells(k0 - 1, l1, l2)
}

/// The densest two-bin window (by z) of per-bin counts `cnt` of `n` anchors over the `l1 x l2`
/// rectangle, scored as `-log10` of its Poisson tail over the number of windows. `NaN` with no
/// anchor.
fn diag_nlp(cnt: &[u32], n: u32, l1: i64, l2: i64) -> f64 {
    diag_best(cnt, n, l1, l2).map_or(f64::NAN, |(obs, e, m)| dp::tail_nlp(obs, e, m))
}

/// [`diag_nlp`]'s window before its tail: `(obs, e, number of windows)` of the densest window by z
/// (first one on ties), `None` with no anchor.
fn diag_best(cnt: &[u32], n: u32, l1: i64, l2: i64) -> Option<(u32, f64, f64)> {
    let nb = cnt.len();
    if n == 0 || nb == 0 {
        return None;
    }
    let area = (l1 * l2).max(1) as f64;
    let nw = nb.saturating_sub(1).max(1);
    let (mut bz, mut bobs, mut be) = (f64::NEG_INFINITY, 0u32, 0.0f64);
    for j in 0..nw {
        let obs = if j + 1 < nb {
            cnt[j] + cnt[j + 1]
        } else {
            cnt[j]
        };
        let k0 = j as i64 * DBIN;
        let frac = cells_in(k0, k0 + 2 * DBIN, l1, l2) / area;
        let e = n as f64 * frac;
        let z = (obs as f64 - e) / e.max(1.0).sqrt();
        if z > bz {
            bz = z;
            bobs = obs;
            be = e;
        }
    }
    Some((bobs, be, nw as f64))
}

/// Margin (in `-log10 p`) between [`dp::tail_nlp_upper`] and a level below which a window is
/// decided without its tail: far above every rounding difference between the two (~1e-13), far
/// below any level's resolution.
const UPPER_MARGIN: f64 = 1e-6;

/// Whether a window's score `dp::tail_nlp(obs, e, m)` reaches `level` (`false` without a window,
/// as `NaN >= level`). The tail's series runs only where its upper bound does not already fall
/// short of the level: the same answer as computing the score, and most chance pairs never pay
/// for the series (`tail_nlp` alone took 10% of detection).
fn reaches(best: Option<(u32, f64, f64)>, level: f64) -> bool {
    match best {
        None => false,
        Some((obs, e, m)) => {
            if dp::tail_nlp_upper(obs, e, m) < level - UPPER_MARGIN {
                return false;
            }
            dp::tail_nlp(obs, e, m) >= level
        }
    }
}

/// Bin `anchors` by their slope-1 offset in the owner frame (`a_pos - b_pos + L_b`) into `cnt`;
/// returns the frame's `(L_1, L_2)`.
fn bin_owner(anchors: &[Anchor], la: u32, lb: u32, cnt: &mut Vec<u32>) -> (i64, i64) {
    let (l1, l2) = (la as i64, lb as i64);
    let nb = ((l1 + l2 - 1).max(0) / DBIN + 1) as usize;
    cnt.clear();
    cnt.resize(nb, 0);
    for x in anchors {
        let k = x.a_pos as i64 - x.b_pos as i64 + l2;
        let j = (k.max(0) / DBIN).min(nb as i64 - 1) as usize;
        cnt[j] += 1;
    }
    (l1, l2)
}

/// [`bin_owner`] in the other read's frame (`b_pos - a_pos + L_a`, the lengths exchanged).
fn bin_transposed(anchors: &[Anchor], la: u32, lb: u32, cnt: &mut Vec<u32>) -> (i64, i64) {
    let (l1, l2) = (lb as i64, la as i64);
    let nb = ((l1 + l2 - 1).max(0) / DBIN + 1) as usize;
    cnt.clear();
    cnt.resize(nb, 0);
    for x in anchors {
        let k = x.b_pos as i64 - x.a_pos as i64 + l2;
        let j = (k.max(0) / DBIN).min(nb as i64 - 1) as usize;
        cnt[j] += 1;
    }
    (l1, l2)
}

/// `nlp2n_1` of one pair's cloud in the owner frame: `anchors` as gathered by the owner `a`
/// (`a_pos` on `a`, `b_pos` on `b`), `la` / `lb` the two reads' lengths in retained events.
/// `cnt` is scratch.
pub fn pair_nlp(anchors: &[Anchor], la: u32, lb: u32, cnt: &mut Vec<u32>) -> f64 {
    let (l1, l2) = bin_owner(anchors, la, lb, cnt);
    diag_nlp(cnt, anchors.len() as u32, l1, l2)
}

/// [`pair_nlp`] of the same cloud in the OTHER read's frame: every anchor's `a_pos` / `b_pos` and the
/// two lengths exchanged, without copying the cloud. Equal to `pair_nlp` of the swapped anchors with
/// `lb, la`, operation for operation.
pub fn pair_nlp_transposed(anchors: &[Anchor], la: u32, lb: u32, cnt: &mut Vec<u32>) -> f64 {
    let (l1, l2) = bin_transposed(anchors, la, lb, cnt);
    diag_nlp(cnt, anchors.len() as u32, l1, l2)
}

/// The nomination rule over one input: every read's length, and its level.
pub struct Nominator<'a> {
    /// Read lengths in retained events, by read id.
    lens: &'a [u32],
    /// The level every read uses, unless `per_read` is set.
    level: f64,
    /// Per-read levels, by read id.
    per_read: Option<&'a [f64]>,
    /// Score both frames and keep the larger (`OverlapParams::symmetric_pairs`): the pair's score
    /// is then the same whichever read owns it.
    symmetric: bool,
}

impl<'a> Nominator<'a> {
    /// Every read at the same level (`--prefilter-nlp L`).
    pub fn fixed(lens: &'a [u32], level: f64) -> Self {
        Nominator {
            lens,
            level,
            per_read: None,
            symmetric: false,
        }
    }

    /// One level per read: a pair is nominated when its score reaches either read's level.
    pub fn per_read(lens: &'a [u32], levels: &'a [f64]) -> Self {
        Nominator {
            lens,
            level: f64::NAN,
            per_read: Some(levels),
            symmetric: false,
        }
    }

    /// With `symmetric`, a pair's score is the larger of its two frames' (owner and partner):
    /// the nomination does not depend on which read owns the pair.
    pub fn symmetric(mut self, symmetric: bool) -> Self {
        self.symmetric = symmetric;
        self
    }

    fn level_of(&self, r: u32) -> f64 {
        self.per_read.map_or(self.level, |v| v[r as usize])
    }

    /// Whether the owner `a` (the lower id) aligns its pair with `b`, given the pair's cloud
    /// as `a` gathered it. A `NaN` score (no anchor) is never nominated.
    ///
    /// The rule is `s >= level(a) || s >= level(b)` with `s` the owner frame's score, or with
    /// `--symmetric` the larger of the two frames' (`max` keeps the other score when one is `NaN`).
    /// It is decided as `s1 >= L || s2 >= L` with `L = min(level(a), level(b))` (the same truth
    /// value, `NaN` included), the second frame only when the first does not already reach `L`,
    /// and each frame's tail only when its upper bound does not already fall short ([`reaches`]).
    pub fn nominated(&self, a: u32, b: u32, anchors: &[Anchor]) -> bool {
        let mut cnt = Vec::new();
        let (la, lb) = (self.lens[a as usize], self.lens[b as usize]);
        let level = self.level_of(a).min(self.level_of(b));
        let n = anchors.len() as u32;
        let (l1, l2) = bin_owner(anchors, la, lb, &mut cnt);
        if reaches(diag_best(&cnt, n, l1, l2), level) {
            return true;
        }
        if !self.symmetric {
            return false;
        }
        let (l1, l2) = bin_transposed(anchors, la, lb, &mut cnt);
        reaches(diag_best(&cnt, n, l1, l2), level)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trapezoid_cells_match_brute_force() {
        for &(l1, l2) in &[(7i64, 4i64), (4, 7), (5, 5), (1, 9), (9, 1), (13, 6)] {
            let t = l1 + l2;
            for k0 in 0..=t + 1 {
                for k1 in k0..=t + 2 {
                    let mut brute = 0i64;
                    for x in 0..l1 {
                        for y in 0..l2 {
                            if (k0..k1).contains(&(x - y + l2)) {
                                brute += 1;
                            }
                        }
                    }
                    assert_eq!(cells_in(k0, k1, l1, l2), brute as f64);
                }
            }
            assert_eq!(cum_cells(t + 5, l1, l2), (l1 * l2) as f64);
        }
    }

    fn mix(mut x: u64) -> u64 {
        x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = x;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// The statistic's behaviour: a diagonal over uniform chance scores far above the
    /// chance cloud alone, which stays near the uniform null.
    #[test]
    fn a_diagonal_scores_above_a_uniform_cloud() {
        let (la, lb) = (20_000u32, 20_000u32);
        let mut x = 12345u64;
        let mut c: Vec<Anchor> = (0..400)
            .map(|_| {
                x = mix(x);
                let a_pos = (x % la as u64) as u32;
                x = mix(x);
                Anchor {
                    a_pos,
                    b_pos: (x % lb as u64) as u32,
                }
            })
            .collect();
        let mut cnt = Vec::new();
        let u = pair_nlp(&c, la, lb, &mut cnt);
        c.extend((0..150u32).map(|i| Anchor {
            a_pos: 5_000 + 20 * i,
            b_pos: 3_000 + 20 * i,
        }));
        let d = pair_nlp(&c, la, lb, &mut cnt);
        assert!(u < 4.0, "uniform cloud {u}");
        assert!(d > u + 5.0, "diagonal {d} vs uniform {u}");
        assert!(pair_nlp(&[], la, lb, &mut cnt).is_nan());
    }

    /// Either read's level nominates: the rule is an OR over the pair's two reads, evaluated on
    /// the one owner-frame score.
    #[test]
    fn either_reads_level_nominates() {
        let lens = [8_000u32, 9_000, 7_000];
        let c: Vec<Anchor> = (0..200u32)
            .map(|i| Anchor {
                a_pos: 1_000 + 20 * i,
                b_pos: 30 * i,
            })
            .collect();
        let s = pair_nlp(&c, lens[0], lens[1], &mut Vec::new());
        assert!(s.is_finite() && s > 1.0);
        let lv = |a: f64, b: f64| [a, b, f64::INFINITY];
        for (levels, want) in [
            (lv(s - 1.0, s + 1.0), true),
            (lv(s + 1.0, s - 1.0), true),
            (lv(s, s + 1.0), true), // reaching a level counts
            (lv(s + 1.0, s + 1.0), false),
        ] {
            let n = Nominator::per_read(&lens, &levels);
            assert_eq!(n.nominated(0, 1, &c), want, "levels {levels:?}");
        }
        assert!(Nominator::fixed(&lens, s).nominated(0, 1, &c));
        assert!(!Nominator::fixed(&lens, s + 1e-9).nominated(0, 1, &c));
        assert!(!Nominator::fixed(&lens, 0.0).nominated(0, 1, &[]));
    }

    /// The nomination (second frame only when the first falls short, each tail only past its upper
    /// bound) is the stated rule on the two frames' full scores, `max(s1, s2) >= min(level_a,
    /// level_b)` (owner frame only without `symmetric`), over chance clouds, weak and strong
    /// diagonals, empty clouds, both frame modes and levels around the scores themselves.
    #[test]
    fn nomination_is_the_rule_on_the_full_scores() {
        let mut x = 4242u64;
        let mut next = |m: u64| {
            x = mix(x);
            x % m
        };
        let (mut yes, mut no) = (0usize, 0usize);
        for trial in 0..3000u64 {
            let la = 600 + next(30_000) as u32;
            let lb = 600 + next(30_000) as u32;
            let n_rand = next(600) as usize;
            let mut c: Vec<Anchor> = (0..n_rand)
                .map(|_| Anchor {
                    a_pos: next(la as u64) as u32,
                    b_pos: next(lb as u64) as u32,
                })
                .collect();
            let n_diag = [0u64, 3, 8, 15, 30, 120][(trial % 6) as usize];
            let (a0, b0) = (next(la as u64 / 2) as u32, next(lb as u64 / 2) as u32);
            c.extend((0..n_diag as u32).map(|i| Anchor {
                a_pos: (a0 + 37 * i).min(la - 1),
                b_pos: (b0 + 35 * i + next(300) as u32).min(lb - 1),
            }));
            let lens = [la, lb];
            let mut cnt = Vec::new();
            let s1 = pair_nlp(&c, la, lb, &mut cnt);
            let s2 = pair_nlp_transposed(&c, la, lb, &mut cnt);
            let mut levels = vec![0.5, 1.0, 2.0, 3.0, 6.0];
            for s in [s1, s2] {
                if s.is_finite() {
                    levels.extend([s - 1e-7, s, s + 1e-7, s - 2e-6, s + 2e-6]);
                }
            }
            for &level in &levels {
                for sym in [false, true] {
                    let want = if sym {
                        s1.max(s2) >= level
                    } else {
                        s1 >= level
                    };
                    let got = Nominator::fixed(&lens, level)
                        .symmetric(sym)
                        .nominated(0, 1, &c);
                    assert_eq!(
                        got, want,
                        "trial {trial} level {level} sym {sym}: {s1} {s2}"
                    );
                    if got { yes += 1 } else { no += 1 }
                }
                let per = [level, level + 1.0];
                let want = s1.max(s2) >= level;
                let got = Nominator::per_read(&lens, &per)
                    .symmetric(true)
                    .nominated(0, 1, &c);
                assert_eq!(got, want, "trial {trial} per-read levels {per:?}");
            }
        }
        assert!(yes > 2_000 && no > 2_000, "{yes} nominated, {no} not");
    }

    /// The symmetric nomination does not depend on which read owns the pair: the cloud gathered by
    /// either read (anchors swapped, ids and lengths swapped) gets the same score, the larger frame.
    #[test]
    fn symmetric_nomination_is_the_same_from_either_read() {
        let (la, lb) = (11_000u32, 7_300u32);
        let mut x = 99u64;
        let mut c: Vec<Anchor> = (0..300)
            .map(|_| {
                x = mix(x);
                let a_pos = (x % la as u64) as u32;
                x = mix(x);
                Anchor {
                    a_pos,
                    b_pos: (x % lb as u64) as u32,
                }
            })
            .collect();
        c.extend((0..90u32).map(|i| Anchor {
            a_pos: 6_000 + 23 * i,
            b_pos: 40 + 21 * i,
        }));
        let sw: Vec<Anchor> = c
            .iter()
            .map(|a| Anchor {
                a_pos: a.b_pos,
                b_pos: a.a_pos,
            })
            .collect();
        let mut cnt = Vec::new();
        let s_ab = pair_nlp(&c, la, lb, &mut cnt);
        let s_ba = pair_nlp(&sw, lb, la, &mut cnt);
        assert_eq!(pair_nlp_transposed(&c, la, lb, &mut cnt), s_ba);
        let (lens_x, lens_y) = ([la, lb], [lb, la]);
        let (lo, hi) = (s_ab.min(s_ba), s_ab.max(s_ba));
        for level in [lo - 0.5, (lo + hi) / 2.0, hi, hi + 0.5] {
            let nx = Nominator::fixed(&lens_x, level).symmetric(true);
            let ny = Nominator::fixed(&lens_y, level).symmetric(true);
            assert_eq!(
                nx.nominated(0, 1, &c),
                ny.nominated(0, 1, &sw),
                "level {level}"
            );
            assert_eq!(nx.nominated(0, 1, &c), hi >= level, "level {level}");
        }
    }
}
