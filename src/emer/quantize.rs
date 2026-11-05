//! Quantization of normalized event values into a small bucket alphabet.
//!
//! Three schemes, selected by [`QuantParams::method`]:
//!
//!   * [`QuantMethod::Piecewise`] — the faithful port of `dynamic_quantize`
//!     (`extern/RawHash2/src/rsketch.c:18`): the bucket budget is split into a
//!     coarse `low` region, a high-resolution `fine` region around 0, and a
//!     coarse `high` region across the normalized range `[-3, 3]`. For the
//!     R10.4.1 preset (`n_buckets=8`, `fine_range=0.75`) this is 1 low + 6 fine +
//!     1 high bucket. Measured on the debug pair, it leaves the two end buckets
//!     nearly empty (entropy 2.52 / 3.0 bits).
//!
//!   * [`QuantMethod::NormalQuantile`] — *equal-frequency* bins of a standard
//!     normal: cut points `Φ⁻¹(i / n_buckets)`. Since events are z-normalized,
//!     these carry ~equal probability mass, lifting entropy toward
//!     `log2(n_buckets)`. Data-free, so the codebook is shared across reads by
//!     construction — see [`QuantMethod`] on why that matters for matching.
//!
//!   * [`QuantMethod::Codebook`] — bucket from explicit ascending cut points in
//!     [`QuantParams::boundaries`] (e.g. quantiles of a pooled reference event
//!     distribution, frozen once). Populating it with *per-read* boundaries gives
//!     the adaptive-quantizer control.

use crate::params::{QuantMethod, QuantParams};

/// Quantize one normalized event value to a bucket in `0..n_buckets`.
///
/// Dispatches on [`QuantParams::method`]. Every scheme is non-decreasing in
/// `signal` and returns a bucket `< n_buckets` (given `n_buckets >= 1`, which
/// [`crate::params::Params::validate`] enforces; `n_buckets == 0` returns 0).
pub fn dynamic_quantize(signal: f32, q: &QuantParams) -> u32 {
    if q.n_buckets == 0 {
        return 0;
    }
    match q.method {
        QuantMethod::Piecewise => piecewise_quantize(signal, q),
        QuantMethod::NormalQuantile => normal_quantile_quantize(signal, q.n_buckets),
        QuantMethod::Codebook => codebook_quantize(signal, q.boundaries.as_slice(), q.n_buckets),
    }
}

/// RawHash2's fixed piecewise-linear map (faithful port of `dynamic_quantize`,
/// `extern/RawHash2/src/rsketch.c:18`).
fn piecewise_quantize(signal: f32, q: &QuantParams) -> u32 {
    let n_buckets = q.n_buckets;

    // How many buckets fall in the fine region.
    let nb_f = n_buckets as f32 * q.fine_range + 0.5;
    let mut n_fine = nb_f as i32; // truncation toward zero, as in C
    if n_fine < 1 {
        n_fine = 1;
    }
    if n_fine as u32 > n_buckets {
        n_fine = n_buckets as i32;
    }
    let n_fine = n_fine as u32;
    let n_low = (n_buckets - n_fine) / 2;
    let n_high = n_buckets - n_fine - n_low;

    const LO_IN: f32 = -3.0;
    const HI_IN: f32 = 3.0;

    if signal < q.fine_min {
        if n_low == 0 {
            return 0;
        }
        let denom = q.fine_min - LO_IN;
        let t = if denom > 0.0 {
            (signal - LO_IN) / denom
        } else {
            0.0
        };
        bucket(t, n_low)
    } else if signal <= q.fine_max {
        let denom = q.fine_max - q.fine_min;
        let t = if denom > 0.0 {
            (signal - q.fine_min) / denom
        } else {
            0.0
        };
        n_low + bucket(t, n_fine)
    } else {
        if n_high == 0 {
            return n_buckets - 1;
        }
        let denom = HI_IN - q.fine_max;
        let t = if denom > 0.0 {
            (signal - q.fine_max) / denom
        } else {
            1.0
        };
        n_low + n_fine + bucket(t, n_high)
    }
}

/// Equal-frequency quantizer against a standard normal: bucket `b` holds
/// `Φ⁻¹(b/n) <= signal < Φ⁻¹((b+1)/n)`. Computed directly from the closed-form
/// normal quantile so no codebook table is needed and the map is identical on
/// every read.
fn normal_quantile_quantize(signal: f32, n_buckets: u32) -> u32 {
    // Compare against the n-1 interior cut points q_i = Φ⁻¹(i/n).
    let n = n_buckets;
    let mut b = 0u32;
    for i in 1..n {
        let edge = normal_quantile(i as f64 / n as f64) as f32;
        if signal >= edge {
            b = i;
        } else {
            break;
        }
    }
    b.min(n - 1)
}

/// Bucket from explicit ascending cut points: `bounds[b-1] <= signal < bounds[b]`.
/// Returns the number of boundaries below `signal`, clamped to `n_buckets-1`.
fn codebook_quantize(signal: f32, bounds: &[f32], n_buckets: u32) -> u32 {
    // `partition_point` = count of edges <= signal = the bucket index. Half-open
    // on the low side so the lowest bucket catches everything below bounds[0].
    let b = bounds.partition_point(|&e| signal >= e) as u32;
    b.min(n_buckets.saturating_sub(1))
}

/// Standard-normal quantile (probit) `Φ⁻¹(p)` for `p ∈ (0,1)`.
///
/// Acklam's rational approximation (max abs error ~1.15e-9 in the central
/// region, ample for an 8–32-bucket cut table). Out of domain, `p <= 0` returns
/// `-∞` and `p >= 1` returns `+∞`; the callers here only pass interior
/// `i / n_buckets`, so they never see those.
pub fn normal_quantile(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    // Coefficients (Acklam 2003), digits as published.
    #[allow(clippy::excessive_precision)]
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    // Break-points for the central vs tail regions.
    const P_LOW: f64 = 0.02425;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p < P_LOW {
        let qv = (-2.0 * p.ln()).sqrt();
        (((((C[0] * qv + C[1]) * qv + C[2]) * qv + C[3]) * qv + C[4]) * qv + C[5])
            / ((((D[0] * qv + D[1]) * qv + D[2]) * qv + D[3]) * qv + 1.0)
    } else if p <= P_HIGH {
        let qv = p - 0.5;
        let r = qv * qv;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * qv
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let qv = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * qv + C[1]) * qv + C[2]) * qv + C[3]) * qv + C[4]) * qv + C[5])
            / ((((D[0] * qv + D[1]) * qv + D[2]) * qv + D[3]) * qv + 1.0)
    }
}

/// The `n_buckets-1` interior standard-normal cut points (equal-probability bins).
/// Exposed so callers/harnesses can inspect or reuse the `NormalQuantile` table.
pub fn normal_quantile_edges(n_buckets: u32) -> Vec<f32> {
    (1..n_buckets)
        .map(|i| normal_quantile(i as f64 / n_buckets as f64) as f32)
        .collect()
}

/// Map a fractional position `t` to a bucket index in `0..count`.
fn bucket(t: f32, count: u32) -> u32 {
    let t = t.clamp(0.0, 1.0);
    let b = (t * count as f32) as u32;
    b.min(count - 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Boundaries, Params};

    fn q() -> QuantParams {
        Params::r10_4_1().quant
    }

    #[test]
    fn output_always_in_range() {
        let q = q();
        for i in -100..=100 {
            let v = i as f32 / 10.0; // -10.0 ..= 10.0
            let b = dynamic_quantize(v, &q);
            assert!(b < q.n_buckets, "bucket {b} >= {}", q.n_buckets);
        }
    }

    #[test]
    fn monotonic_nondecreasing() {
        let q = q();
        let mut prev = 0;
        for i in -60..=60 {
            let v = i as f32 / 20.0; // -3.0 ..= 3.0
            let b = dynamic_quantize(v, &q);
            assert!(b >= prev, "quantization must be non-decreasing in signal");
            prev = b;
        }
    }

    #[test]
    fn extremes_hit_end_buckets() {
        let q = q();
        assert_eq!(dynamic_quantize(-100.0, &q), 0);
        assert_eq!(dynamic_quantize(100.0, &q), q.n_buckets - 1);
    }

    #[test]
    fn r10_region_split_is_1_6_1() {
        // n_buckets=8, fine_range=0.75 -> 1 low, 6 fine, 1 high.
        let q = q();
        // The single low bucket is index 0; the single high bucket is index 7.
        assert_eq!(dynamic_quantize(-2.5, &q), 0); // below fine_min -> low region
        assert_eq!(dynamic_quantize(2.5, &q), 7); // above fine_max -> high region
    }

    // ---- equal-frequency (normal-quantile) scheme ---------------------------

    /// `Φ⁻¹` round-trips known quantiles: the median is 0, the 0.975 quantile is
    /// ~1.96, and the function is odd about p=0.5.
    #[test]
    fn probit_known_values() {
        assert!((normal_quantile(0.5)).abs() < 1e-6);
        assert!((normal_quantile(0.975) - 1.959964).abs() < 1e-3);
        assert!((normal_quantile(0.84134) - 1.0).abs() < 1e-3);
        assert!((normal_quantile(0.3) + normal_quantile(0.7)).abs() < 1e-6); // odd
    }

    #[test]
    fn normal_quantile_buckets_are_in_range_and_monotone() {
        let mut q = q();
        q.method = QuantMethod::NormalQuantile;
        let mut prev = 0;
        for i in -80..=80 {
            let v = i as f32 / 20.0; // -4..=4
            let b = dynamic_quantize(v, &q);
            assert!(b < q.n_buckets);
            assert!(b >= prev);
            prev = b;
        }
        assert_eq!(dynamic_quantize(-100.0, &q), 0);
        assert_eq!(dynamic_quantize(100.0, &q), q.n_buckets - 1);
    }

    /// For 8 buckets, value 0 must land just below the median cut so it sits in
    /// bucket 3 or 4 (the two central bins) — the whole point of equal-frequency.
    #[test]
    fn normal_quantile_centers_mass() {
        let mut q = q();
        q.method = QuantMethod::NormalQuantile;
        let b0 = dynamic_quantize(-0.0001, &q);
        assert!(
            b0 == 3,
            "value ~0 should be the top of the lower-central bin, got {b0}"
        );
        assert_eq!(dynamic_quantize(0.0, &q), 4);
    }

    /// A standard-normal sample, bucketed by the normal-quantile map, fills every
    /// bucket to within ~1/3 of uniform (contrast the piecewise map, whose end
    /// buckets stay nearly empty).
    #[test]
    fn normal_quantile_is_near_equal_frequency_on_gaussian_sample() {
        let mut q = q();
        q.method = QuantMethod::NormalQuantile;
        // Cheap N(0,1)-ish sample via Box–Muller on a xorshift stream.
        let mut x: u64 = 0xDEAD_BEEF_1234_5678;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 53) as f64
        };
        let n = q.n_buckets as usize;
        let mut counts = vec![0u64; n];
        let mut total = 0u64;
        for _ in 0..50_000 {
            let u1: f64 = next().max(1e-12);
            let u2: f64 = next();
            let z = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
            counts[dynamic_quantize(z as f32, &q) as usize] += 1;
            total += 1;
        }
        let expect = total as f64 / n as f64;
        for (i, &c) in counts.iter().enumerate() {
            let ratio = c as f64 / expect;
            assert!(
                (0.66..1.5).contains(&ratio),
                "bucket {i} ratio {ratio:.2} far from uniform (count {c}, expect {expect:.0})"
            );
        }
    }

    // ---- explicit codebook scheme -------------------------------------------

    #[test]
    fn codebook_buckets_by_thresholds() {
        let mut q = q();
        q.method = QuantMethod::Codebook;
        // 7 ascending edges -> 8 buckets at -3,-2,-1,0,1,2,3.
        q.boundaries = Boundaries::from_edges(&[-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]);
        assert_eq!(dynamic_quantize(-100.0, &q), 0);
        assert_eq!(dynamic_quantize(-2.5, &q), 1);
        assert_eq!(dynamic_quantize(-0.5, &q), 3);
        assert_eq!(dynamic_quantize(0.0, &q), 4); // >= 0.0 edge
        assert_eq!(dynamic_quantize(100.0, &q), 7);
    }
}
