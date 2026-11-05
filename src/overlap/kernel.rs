//! The slope scan's anchor-level loops, compiled twice: for the baseline x86-64 target (SSE2) and
//! with AVX2, the copy chosen once per process at run time ([`avx2`]). CPUs without AVX2
//! (e.g. Xeon E5-26xx v1/v2) run the baseline copy, so the build sets no `target-cpu`.
//!
//! **Exact by construction.** Both copies run the same source with the same IEEE-754 operations
//! in the same order: an offset is `a - m * b` (a multiply, then a subtract, never fused: FMA is not
//! enabled, and Rust never contracts on its own), a cell `((c - lo) * inv)` truncated. Only the
//! instruction width differs. Integer counting and window maxima vectorise freely. A cell is
//! converted through `i32`, which equals the `u32` conversion on every cell a plan admits (they lie
//! in `[0, top]` with `top <= 2 n`, see [`super::dp`]'s `Cloud`); the AVX2 count uses the packed
//! truncating conversion, which agrees with both there.
//!
//! The loops (the slope scan spends ~95% of detection in them):
//! * [`count`]: every anchor's cell at up to four slopes in one pass, each into its own histogram;
//! * [`bounds`]: one slope's prefix sums and its two window bounds;
//! * [`cells`]: one slope's offsets and cells, stored for the windows that need sorting;
//! * [`sweep`]: [`super::dp`]'s `SlopeSweep` step (offsets, insertion pass, densest window).

use std::sync::OnceLock;

/// One slope's cell map, as `diagonal_peak` draws it: the offset `c = a - m * b` falls in cell
/// `((c - lo) * inv)`, truncated.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Cells {
    pub m: f64,
    pub lo: f64,
    pub inv: f64,
}

/// Whether the AVX2 copies run: the CPU has AVX2. Decided once.
pub(crate) fn avx2() -> bool {
    #[cfg(test)]
    if let Some(f) = FORCE.with(|f| f.get()) {
        return f;
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(has_avx2)
}

fn has_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Which copy runs, for the logs: `avx2`, `scalar (forced)` or `scalar (no avx2)`.
pub fn kernel_name() -> &'static str {
    if avx2() {
        "avx2"
    } else if has_avx2() {
        "scalar (forced)"
    } else {
        "scalar (no avx2)"
    }
}

#[cfg(test)]
thread_local! {
    static FORCE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Tests: run this thread's kernels as AVX2 (`Some(true)`, only where the CPU has it), as the
/// baseline copy (`Some(false)`), or as detected (`None`). Returns whether AVX2 now runs.
#[cfg(test)]
pub(crate) fn force(f: Option<bool>) -> bool {
    FORCE.with(|c| c.set(f.map(|x| x && has_avx2())));
    avx2()
}

/// Count every anchor `(af[i], bf[i])` into its cell at each slope of `cs` (one to four): cell `t`
/// of slope `l` is counted at `hist[base[l] + t + 1]`. Every cell must lie in `[0, 2^31)` and every
/// index in `hist`.
pub(crate) fn count(af: &[f64], bf: &[f64], cs: &[Cells], base: &[u32], hist: &mut [u32]) {
    #[cfg(target_arch = "x86_64")]
    if avx2() {
        // SAFETY: `avx2()` found AVX2 on this CPU.
        return unsafe { x86::count(af, bf, cs, base, hist) };
    }
    count_any(af, bf, cs, base, hist)
}

#[inline(always)]
fn count_any(af: &[f64], bf: &[f64], cs: &[Cells], base: &[u32], hist: &mut [u32]) {
    match cs.len() {
        1 => count_lanes::<1>(af, bf, cs, base, hist),
        2 => count_lanes::<2>(af, bf, cs, base, hist),
        3 => count_lanes::<3>(af, bf, cs, base, hist),
        4 => count_lanes::<4>(af, bf, cs, base, hist),
        l => panic!("count: {l} slopes (1 to 4)"),
    }
}

#[inline(always)]
fn count_lanes<const L: usize>(
    af: &[f64],
    bf: &[f64],
    cs: &[Cells],
    base: &[u32],
    hist: &mut [u32],
) {
    let c: [Cells; L] = std::array::from_fn(|l| cs[l]);
    let base: [u32; L] = std::array::from_fn(|l| base[l] + 1);
    for (&a, &b) in af.iter().zip(bf) {
        let t: [u32; L] = std::array::from_fn(|l| {
            let x = a - c[l].m * b;
            base[l] + ((x - c[l].lo) * c[l].inv) as i32 as u32
        });
        for &t in &t {
            hist[t as usize] += 1;
        }
    }
}

/// Prefix sums of one slope's counts in place (`cum[k]` = anchors in cells `< k` afterwards; the
/// counts sit at `cum[t + 1]`), and the largest sums of `q - 1` and of `q + 2` consecutive cells:
/// `diagonal_peak`'s lower bound on the peak and upper bound on every window. A run cut short by
/// an end of `cum` is inside a full one (`cum` is non-decreasing), so only full runs are scanned.
pub(crate) fn bounds(cum: &mut [u32], q: usize) -> (u32, u32) {
    #[cfg(target_arch = "x86_64")]
    if avx2() {
        // SAFETY: `avx2()` found AVX2 on this CPU.
        return unsafe { x86::bounds(cum, q) };
    }
    bounds_any(cum, q)
}

#[inline(always)]
fn bounds_any(cum: &mut [u32], q: usize) -> (u32, u32) {
    let mut run = 0u32;
    for c in cum.iter_mut() {
        run += *c;
        *c = run;
    }
    let cum = &*cum;
    let nc = cum.len() - 1;
    let run_max = |len: usize| -> u32 {
        if len > nc {
            cum[nc]
        } else {
            cum[len..]
                .iter()
                .zip(cum.iter())
                .map(|(h, l)| h - l)
                .max()
                .unwrap_or(0)
        }
    };
    (run_max(q - 1), run_max(q + 2))
}

/// One slope's offsets `a - m * b` into `vals` and their cells into `cells` (`af.len()` of each).
pub(crate) fn cells(af: &[f64], bf: &[f64], c: Cells, vals: &mut [f64], cells: &mut [u32]) {
    #[cfg(target_arch = "x86_64")]
    if avx2() {
        // SAFETY: `avx2()` found AVX2 on this CPU.
        return unsafe { x86::cells(af, bf, c, vals, cells) };
    }
    cells_any(af, bf, c, vals, cells)
}

#[inline(always)]
fn cells_any(af: &[f64], bf: &[f64], c: Cells, vals: &mut [f64], cells: &mut [u32]) {
    for (((v, t), &a), &b) in vals.iter_mut().zip(cells.iter_mut()).zip(af).zip(bf) {
        let x = a - c.m * b;
        *v = x;
        *t = ((x - c.lo) * c.inv) as i32 as u32;
    }
}

/// `SlopeSweep`'s step at slope `m`: the offsets of the positions `pa`, `pb` (kept in ascending
/// offset at the previous slope) into `c`, an insertion pass that sorts them (moving the positions
/// with them), and the densest `win`-wide window ([`densest_sorted`]).
///
/// Offsets are never NaN and never -0.0 (`a - m*b` rounds an exact zero to +0.0), so the numeric
/// order is `total_cmp`'s and a plain comparison sorts them; any sort gives the one sorted
/// sequence of values, which is all the window reads.
pub(crate) fn sweep(
    pa: &mut [f64],
    pb: &mut [f64],
    c: &mut [f64],
    m: f64,
    win: f64,
) -> (usize, f64) {
    #[cfg(target_arch = "x86_64")]
    if avx2() {
        // SAFETY: `avx2()` found AVX2 on this CPU.
        return unsafe { x86::sweep(pa, pb, c, m, win) };
    }
    sweep_any(pa, pb, c, m, win)
}

#[inline(always)]
fn sweep_any(pa: &mut [f64], pb: &mut [f64], c: &mut [f64], m: f64, win: f64) -> (usize, f64) {
    for ((v, &a), &b) in c.iter_mut().zip(pa.iter()).zip(pb.iter()) {
        *v = a - m * b;
    }
    for k in 1..c.len() {
        if c[k - 1] > c[k] {
            let (v, a, b) = (c[k], pa[k], pb[k]);
            let mut t = k;
            while t > 0 && c[t - 1] > v {
                c[t] = c[t - 1];
                pa[t] = pa[t - 1];
                pb[t] = pb[t - 1];
                t -= 1;
            }
            c[t] = v;
            pa[t] = a;
            pb[t] = b;
        }
    }
    densest_sorted(c, win)
}

/// The densest `win`-wide window over the sorted `c`: `(count, (c[i] + c[j]) / 2)` for the first
/// `j`, in order, whose window -- the `i <= j` with `c[j] - c[i] <= win` -- holds the largest
/// count, as the slide
///
/// ```text
///   for j: { while c[j] - c[i] > win { i += 1 }; if j - i + 1 > peak { record } }
/// ```
///
/// finds it. The test `c[j] - c[i] <= win` is monotone in `i` (rounding is monotone), so the
/// window ending at `j` holds more than `k` offsets iff `c[j] - c[j - k] <= win`: one test per
/// `j` against the best count so far, raised while it passes. The first `j` raising the count
/// to its final value is the slide's first maximum, and its window starts `count - 1` back.
/// Each `j` costs one comparison whose branch is almost always the same way, where the slide's
/// inner loop exits after a coin-flip number of steps on a sparse cloud.
#[inline(always)]
pub(crate) fn densest_sorted(c: &[f64], win: f64) -> (usize, f64) {
    let (mut k, mut at) = (1usize, 0usize);
    for j in 1..c.len() {
        while j >= k && c[j] - c[j - k] <= win {
            k += 1;
            at = j;
        }
    }
    (k, (c[at + 1 - k] + c[at]) / 2.0)
}

/// The AVX2 copies: each compiles the shared body with AVX2 enabled (and nothing else: no FMA).
#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::Cells;
    use std::arch::x86_64::*;

    // The intrinsics take no pointers; recent toolchains make them safe inside an AVX2 function,
    // older ones do not, hence the block and the allowance.
    #[target_feature(enable = "avx2")]
    #[allow(unused_unsafe)]
    pub(super) fn count(af: &[f64], bf: &[f64], cs: &[Cells], base: &[u32], hist: &mut [u32]) {
        if cs.len() != 4 {
            return super::count_any(af, bf, cs, base, hist);
        }
        // Four slopes in the four lanes. The same operations as `count_lanes`, lane by lane:
        // `a - m * b` (multiply, then subtract), `(c - lo) * inv`, truncated to i32.
        // SAFETY: AVX2 (and so AVX and SSE4.1) is enabled here; no pointer is read or written.
        unsafe {
            let m = _mm256_setr_pd(cs[0].m, cs[1].m, cs[2].m, cs[3].m);
            let lo = _mm256_setr_pd(cs[0].lo, cs[1].lo, cs[2].lo, cs[3].lo);
            let inv = _mm256_setr_pd(cs[0].inv, cs[1].inv, cs[2].inv, cs[3].inv);
            let base = _mm_setr_epi32(
                (base[0] + 1) as i32,
                (base[1] + 1) as i32,
                (base[2] + 1) as i32,
                (base[3] + 1) as i32,
            );
            for (&a, &b) in af.iter().zip(bf) {
                let c = _mm256_sub_pd(_mm256_set1_pd(a), _mm256_mul_pd(m, _mm256_set1_pd(b)));
                let x = _mm256_mul_pd(_mm256_sub_pd(c, lo), inv);
                let t = _mm_add_epi32(_mm256_cvttpd_epi32(x), base);
                hist[_mm_extract_epi32::<0>(t) as u32 as usize] += 1;
                hist[_mm_extract_epi32::<1>(t) as u32 as usize] += 1;
                hist[_mm_extract_epi32::<2>(t) as u32 as usize] += 1;
                hist[_mm_extract_epi32::<3>(t) as u32 as usize] += 1;
            }
        }
    }

    #[target_feature(enable = "avx2")]
    pub(super) fn bounds(cum: &mut [u32], q: usize) -> (u32, u32) {
        super::bounds_any(cum, q)
    }

    #[target_feature(enable = "avx2")]
    pub(super) fn cells(af: &[f64], bf: &[f64], c: Cells, vals: &mut [f64], cells: &mut [u32]) {
        super::cells_any(af, bf, c, vals, cells)
    }

    #[target_feature(enable = "avx2")]
    pub(super) fn sweep(
        pa: &mut [f64],
        pb: &mut [f64],
        c: &mut [f64],
        m: f64,
        win: f64,
    ) -> (usize, f64) {
        super::sweep_any(pa, pb, c, m, win)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(s: &mut u64) -> u32 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*s >> 33) as u32
    }

    /// A cloud's positions and the cell maps of `diagonal_peak` at four slopes (its offsets lie in
    /// `[lo, hi]`, so its cells in `[0, top]`).
    fn setup(n: usize, ext: u32, s: &mut u64) -> (Vec<f64>, Vec<f64>, Vec<Cells>, Vec<usize>) {
        let af: Vec<f64> = (0..n).map(|_| (lcg(s) % ext) as f64).collect();
        let bf: Vec<f64> = (0..n).map(|_| (lcg(s) % ext) as f64).collect();
        let max = |v: &[f64]| v.iter().copied().fold(0.0f64, f64::max);
        let min = |v: &[f64]| v.iter().copied().fold(f64::INFINITY, f64::min);
        let (mut cs, mut nc) = (Vec::new(), Vec::new());
        for m in [0.7, 0.91, 1.0300000000000002, 1.4] {
            let lo = min(&af) - m * max(&bf);
            let hi = max(&af) - m * min(&bf);
            let inv = 12.0 / 384.0;
            cs.push(Cells { m, lo, inv });
            nc.push(((hi - lo) * inv) as usize + 1);
        }
        (af, bf, cs, nc)
    }

    /// Run `f` under the baseline copy, then under AVX2 when the CPU has it, and return both.
    fn both<T>(mut f: impl FnMut() -> T) -> (T, T) {
        force(Some(false));
        let s = f();
        force(Some(true));
        let v = f();
        force(None);
        (s, v)
    }

    /// The four-lane count, the narrower ones and the baseline copy give one histogram, equal to
    /// counting each slope alone.
    #[test]
    fn count_copies_agree_with_one_slope_at_a_time() {
        let mut s = 9u64;
        for &(n, ext) in &[(1usize, 50u32), (37, 400), (1000, 30_000), (5000, 2_000)] {
            let (af, bf, cs, nc) = setup(n, ext, &mut s);
            let base: Vec<u32> = nc
                .iter()
                .scan(0u32, |b, &k| {
                    let at = *b;
                    *b += k as u32 + 1;
                    Some(at)
                })
                .collect();
            let total = base[3] as usize + nc[3] + 1;
            let mut want = vec![0u32; total];
            for (c, b) in cs.chunks(1).zip(base.chunks(1)) {
                count_lanes::<1>(&af, &bf, c, b, &mut want);
            }
            for lanes in [1usize, 2, 3, 4] {
                let (sc, vx) = both(|| {
                    let mut h = vec![0u32; total];
                    for (c, b) in cs.chunks(lanes).zip(base.chunks(lanes)) {
                        count(&af, &bf, c, b, &mut h);
                    }
                    h
                });
                assert_eq!(sc, want, "n={n} lanes={lanes} baseline");
                assert_eq!(vx, want, "n={n} lanes={lanes} avx2");
            }
        }
    }

    /// Bounds, cells and the sweep step agree bit for bit between the two copies.
    #[test]
    fn bounds_cells_and_sweep_agree_between_copies() {
        let mut s = 4u64;
        for &(n, ext) in &[(8usize, 20u32), (300, 30_000), (2000, 500)] {
            let (af, bf, cs, nc) = setup(n, ext, &mut s);
            for (l, c) in cs.iter().enumerate() {
                let (sc, vx) = both(|| {
                    let (mut v, mut t) = (vec![0.0; n], vec![0u32; n]);
                    cells(&af, &bf, *c, &mut v, &mut t);
                    let mut h = vec![0u32; nc[l] + 1];
                    for &t in &t {
                        h[t as usize + 1] += 1;
                    }
                    let b = [2usize, 5, 12, 16].map(|q| bounds(&mut h.clone(), q));
                    (v.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), t, b)
                });
                assert_eq!(sc, vx, "n={n} m={}", c.m);
            }
            let (sc, vx) = both(|| {
                let (mut pa, mut pb, mut c) = (af.clone(), bf.clone(), vec![0.0; n]);
                let mut m = 0.7;
                let mut out = Vec::new();
                for _ in 0..30 {
                    let (k, at) = sweep(&mut pa, &mut pb, &mut c, m, 384.0);
                    out.push((k, at.to_bits()));
                    m += 0.02;
                }
                out
            });
            assert_eq!(sc, vx, "n={n} sweep");
        }
    }
}
