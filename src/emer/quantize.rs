//! Quantization of normalized event values into a small bucket alphabet.
//!
//! The faithful port of `dynamic_quantize` (`extern/RawHash2/src/rsketch.c:18`): the bucket
//! budget is split into a coarse `low` region, a high-resolution `fine` region around 0, and a
//! coarse `high` region across the normalized range `[-3, 3]`. For the R10.4.1 preset
//! (`n_buckets=8`, `fine_range=0.75`) this is 1 low + 6 fine + 1 high bucket. Measured on the
//! debug pair, it leaves the two end buckets nearly empty (entropy 2.52 / 3.0 bits).

use crate::params::QuantParams;

/// Quantize one normalized event value to a bucket in `0..n_buckets`.
///
/// Non-decreasing in `signal`; returns a bucket `< n_buckets` (given `n_buckets >= 1`, which
/// [`crate::params::Params::validate`] enforces; `n_buckets == 0` returns 0).
pub fn dynamic_quantize(signal: f32, q: &QuantParams) -> u32 {
    if q.n_buckets == 0 {
        return 0;
    }
    piecewise_quantize(signal, q)
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

/// Map a fractional position `t` to a bucket index in `0..count`.
fn bucket(t: f32, count: u32) -> u32 {
    let t = t.clamp(0.0, 1.0);
    let b = (t * count as f32) as u32;
    b.min(count - 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Params;

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
}
