//! Z-score normalization with ±3σ outlier removal.
//!
//! Faithful port of `normalize_signal` (`extern/RawHash2/src/revent.c:233`).
//! For a single read the running sums start at zero, so this reduces to a
//! per-read z-score. Note RawHash2 *drops* (does not clip) samples whose
//! normalized value falls outside (-3, 3); the returned vector is therefore
//! generally shorter than the input.

/// Normalize `signal` to zero mean / unit variance and drop |z| >= 3 outliers.
///
/// Returns the filtered, normalized samples. Returns an empty vector if the
/// input is empty or has (near-)zero variance.
pub fn normalize(signal: &[f32]) -> Vec<f32> {
    let Some((mean, std_dev)) = zstats(signal) else {
        return Vec::new();
    };
    let mut out = Vec::with_capacity(signal.len());
    for &x in signal {
        if let Some(z) = keep_z(x, mean, std_dev) {
            out.push(z);
        }
    }
    out
}

/// [`normalize`], plus the index in `signal` of every sample it kept, so a position in
/// the normalized stream can be traced back to the input.
pub fn normalize_kept(signal: &[f32]) -> (Vec<f32>, Vec<u32>) {
    let Some((mean, std_dev)) = zstats(signal) else {
        return (Vec::new(), Vec::new());
    };
    let mut out = Vec::with_capacity(signal.len());
    let mut kept = Vec::with_capacity(signal.len());
    for (i, &x) in signal.iter().enumerate() {
        if let Some(z) = keep_z(x, mean, std_dev) {
            out.push(z);
            kept.push(i as u32);
        }
    }
    (out, kept)
}

/// Mean and population standard deviation, or `None` for an empty or flat signal.
fn zstats(signal: &[f32]) -> Option<(f64, f64)> {
    let n = signal.len();
    if n == 0 {
        return None;
    }
    let mut sum = 0.0f64;
    let mut sum_sq = 0.0f64;
    for &x in signal {
        let x = x as f64;
        sum += x;
        sum_sq += x * x;
    }
    let mean = sum / n as f64;
    let var = sum_sq / n as f64 - mean * mean;
    let std_dev = var.sqrt();
    if !std_dev.is_finite() || std_dev == 0.0 {
        return None;
    }
    Some((mean, std_dev))
}

/// The z-score of `x`, or `None` when RawHash2 would drop it: it keeps strictly inside
/// (-3, 3).
#[inline]
fn keep_z(x: f32, mean: f64, std_dev: f64) -> Option<f32> {
    let z = (x as f64 - mean) / std_dev;
    (z < 3.0 && z > -3.0).then_some(z as f32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input() {
        assert!(normalize(&[]).is_empty());
    }

    #[test]
    fn constant_signal_has_zero_variance() {
        assert!(normalize(&[5.0; 16]).is_empty());
    }

    #[test]
    fn normalized_output_is_roughly_zero_mean_unit_var() {
        // A spread-out signal with no extreme outliers: nothing is dropped.
        let sig: Vec<f32> = (0..100).map(|i| (i as f32 - 50.0) * 0.1).collect();
        let norm = normalize(&sig);
        assert_eq!(norm.len(), sig.len()); // linear ramp, |z| < 3 everywhere

        let mean: f32 = norm.iter().sum::<f32>() / norm.len() as f32;
        assert!(mean.abs() < 1e-3, "mean was {mean}");

        let var: f32 = norm.iter().map(|&x| (x - mean).powi(2)).sum::<f32>() / norm.len() as f32;
        assert!((var - 1.0).abs() < 1e-2, "var was {var}");
    }

    #[test]
    fn drops_extreme_outliers() {
        let mut sig = vec![0.0f32; 200];
        sig[0] = 1000.0; // a single huge spike -> |z| >= 3 -> dropped
        let norm = normalize(&sig);
        assert_eq!(norm.len(), sig.len() - 1);
    }

    #[test]
    fn kept_indices_trace_back_to_the_input() {
        let mut sig: Vec<f32> = (0..300).map(|i| ((i * 37) % 11) as f32).collect();
        sig[7] = 1000.0;
        sig[150] = -1000.0;
        let (norm, kept) = normalize_kept(&sig);
        assert_eq!(norm, normalize(&sig));
        assert_eq!(kept.len(), norm.len());
        assert!(!kept.contains(&7) && !kept.contains(&150));
        assert!(kept.windows(2).all(|w| w[0] < w[1]));
    }
}
