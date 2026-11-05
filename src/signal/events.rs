//! Event segmentation: the dual t-test peak detector.
//!
//! Faithful port of the segmentation in `extern/RawHash2/src/revent.c`
//! (`comp_prefix_prefixsq`, `comp_tstat`, `gen_peaks`,
//! `calculate_mean_of_filtered_segment`, `gen_events`). RawHash2 uses this
//! algorithm for R10.4.1, so the logic mirrors the C closely.
//!
//! Input is the **normalized** signal (see [`crate::signal::normalize`]); the
//! output is one mean current value per detected event ("segment").
//!
//! Two deliberate robustness deviations from the C:
//!   - prefix sums are accumulated in `f64` (`comp_prefix_prefixsq` uses `f32`);
//!   - zero-/negative-length segments are skipped rather than indexed (the C
//!     assumes strictly increasing peak positions).

use crate::params::SegmentationParams;

/// Prefix sums of the signal and of its squares, each of length `n + 1`.
fn prefix_sums(sig: &[f32]) -> (Vec<f64>, Vec<f64>) {
    let n = sig.len();
    let mut sum = vec![0.0f64; n + 1];
    let mut sum_sq = vec![0.0f64; n + 1];
    for i in 0..n {
        let x = sig[i] as f64;
        sum[i + 1] = sum[i] + x;
        sum_sq[i + 1] = sum_sq[i] + x * x;
    }
    (sum, sum_sq)
}

/// Welch-style t-statistic comparing two adjacent windows of length `w_len`.
/// Returns an array of length `n + 1` (matching the C layout); boundary
/// positions are left at zero.
fn comp_tstat(prefix: &[f64], prefix_sq: &[f64], n: usize, w_len: usize) -> Vec<f32> {
    let eta = f32::MIN_POSITIVE as f64;
    let mut tstat = vec![0.0f32; n + 1];
    if n < 2 * w_len || w_len < 2 {
        return tstat;
    }
    let w = w_len as f64;
    // Positions i in w_len..=n-w_len, each from prefix[i - w_len], prefix[i], prefix[i + w_len]:
    // three shifted slices, so the loop has no bounds tests and vectorises (the same arithmetic,
    // in the same order, per position).
    let m = n - 2 * w_len + 1;
    let (lo, mid, hi) = (
        &prefix[..m],
        &prefix[w_len..w_len + m],
        &prefix[2 * w_len..2 * w_len + m],
    );
    let (lo2, mid2, hi2) = (
        &prefix_sq[..m],
        &prefix_sq[w_len..w_len + m],
        &prefix_sq[2 * w_len..2 * w_len + m],
    );
    let out = &mut tstat[w_len..w_len + m];
    for k in 0..m {
        // window 1 = [i - w_len, i); prefix[0] == 0 so the subtraction is safe at i == w_len.
        let sum1 = mid[k] - lo[k];
        let sumsq1 = mid2[k] - lo2[k];
        // window 2 = [i, i + w_len)
        let sum2 = hi[k] - mid[k];
        let sumsq2 = hi2[k] - mid2[k];
        let mean1 = sum1 / w;
        let mean2 = sum2 / w;
        let combined_var = (sumsq1 / w - mean1 * mean1 + sumsq2 / w - mean2 * mean2) / w;
        let combined_var = combined_var.max(eta);
        out[k] = ((mean2 - mean1).abs() / combined_var.sqrt()) as f32;
    }
    // fudge boundaries: zero the last w_len positions of the valid range
    for v in tstat.iter_mut().take(n + 1).skip(n - w_len + 1) {
        *v = 0.0;
    }
    tstat
}

const DEF_PEAK_POS: i64 = -1;

/// Mutable + configured state of one peak detector.
struct Detector {
    threshold: f32,
    window_length: i64,
    masked_to: i64,
    peak_pos: i64,
    peak_value: f32,
    valid_peak: bool,
}

impl Detector {
    fn new(threshold: f32, window_length: u32) -> Self {
        Detector {
            threshold,
            window_length: window_length as i64,
            masked_to: 0,
            peak_pos: DEF_PEAK_POS,
            peak_value: f32::MAX,
            valid_peak: false,
        }
    }
}

/// Run the coupled short/long detectors over their t-statistic tracks and
/// return peak positions (indices into the normalized signal).
fn gen_peaks(
    tstats: &[Vec<f32>],
    det: &mut [Detector],
    n: usize,
    peak_height: f32,
    w0: i64,
) -> Vec<u32> {
    let mut peaks = Vec::new();
    let n_det = det.len();
    for i in 0..n as i64 {
        for k in 0..n_det {
            if det[k].masked_to >= i {
                continue;
            }
            let current_value = tstats[k][i as usize];

            if det[k].peak_pos == DEF_PEAK_POS {
                // CASE 1: have not recorded a maximum yet.
                if current_value < det[k].peak_value {
                    det[k].peak_value = current_value; // deeper trough
                } else if current_value - det[k].peak_value > peak_height {
                    det[k].peak_value = current_value; // qualifying rise
                    det[k].peak_pos = i;
                }
            } else {
                // CASE 2: tracking an existing peak.
                if current_value > det[k].peak_value {
                    det[k].peak_value = current_value;
                    det[k].peak_pos = i;
                }
                // Mask all higher-index detectors past this peak.
                if det[k].peak_value > det[k].threshold {
                    let mask = det[k].peak_pos + w0;
                    for d in det.iter_mut().skip(k + 1) {
                        d.masked_to = mask;
                        d.peak_pos = DEF_PEAK_POS;
                        d.peak_value = f32::MAX;
                        d.valid_peak = false;
                    }
                }
                if det[k].peak_value - current_value > peak_height
                    && det[k].peak_value > det[k].threshold
                {
                    det[k].valid_peak = true;
                }
                if det[k].valid_peak && (i - det[k].peak_pos) > det[k].window_length / 2 {
                    peaks.push(det[k].peak_pos as u32);
                    det[k].peak_pos = DEF_PEAK_POS;
                    det[k].peak_value = current_value;
                    det[k].valid_peak = false;
                }
            }
        }
    }
    peaks
}

/// IQR-filtered mean of a segment (`calculate_mean_of_filtered_segment`).
fn filtered_segment_mean(segment: &[f32], s: &mut Vec<f32>) -> f32 {
    if segment.is_empty() {
        return 0.0;
    }
    s.clear();
    s.extend_from_slice(segment);
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q1 = s[s.len() / 4];
    let q3 = s[3 * s.len() / 4];
    let iqr = q3 - q1;
    let lower = q1 - iqr;
    let upper = q3 + iqr;

    let mut sum = 0.0f32;
    let mut count = 0u32;
    for &x in s.iter() {
        if x >= lower && x <= upper {
            sum += x;
            count += 1;
        }
    }
    if count > 0 { sum / count as f32 } else { 0.0 }
}

/// One detected event: its half-open sample span `[start, end)` in the stream it was
/// segmented from, and its IQR-filtered mean.
///
/// The span is what ties an event to the raw signal, and through a basecaller's move
/// table to a read base. [`detect_events`] keeps only the means.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventSpan {
    pub start: u32,
    pub end: u32,
    pub mean: f32,
}

/// Turn peak positions into per-segment spans and mean current values.
fn gen_events(sig: &[f32], peaks: &[u32], min_seg: u32, max_seg: u32) -> Vec<EventSpan> {
    let s_len = sig.len() as i64;
    let mut events = Vec::new();
    let mut scratch = Vec::new();
    let mut start_idx: i64 = 0;
    for &p in peaks {
        let p = p as i64;
        if !(p > 0 && p < s_len) {
            continue;
        }
        let segment_length = p - start_idx;
        // Robustness vs. the C: guard non-positive lengths (C assumes sorted peaks).
        if segment_length > 0 && segment_length >= min_seg as i64 && segment_length < max_seg as i64
        {
            events.push(EventSpan {
                start: start_idx as u32,
                end: p as u32,
                mean: filtered_segment_mean(&sig[start_idx as usize..p as usize], &mut scratch),
            });
        }
        start_idx = p;
    }
    events
}

/// Detect events in a normalized signal, returning one mean value per event.
pub fn detect_events(norm: &[f32], seg: &SegmentationParams) -> Vec<f32> {
    detect_event_spans(norm, seg)
        .into_iter()
        .map(|e| e.mean)
        .collect()
}

/// [`detect_events`] with each event's sample span in `norm` kept alongside its mean.
pub fn detect_event_spans(norm: &[f32], seg: &SegmentationParams) -> Vec<EventSpan> {
    let n = norm.len();
    if n == 0 {
        return Vec::new();
    }
    let (prefix, prefix_sq) = prefix_sums(norm);
    let tstats = [
        comp_tstat(&prefix, &prefix_sq, n, seg.window_length1 as usize),
        comp_tstat(&prefix, &prefix_sq, n, seg.window_length2 as usize),
    ];
    let mut detectors = [
        Detector::new(seg.threshold1, seg.window_length1),
        Detector::new(seg.threshold2, seg.window_length2),
    ];
    let peaks = gen_peaks(
        &tstats,
        &mut detectors,
        n,
        seg.peak_height,
        seg.window_length1 as i64,
    );
    gen_events(norm, &peaks, seg.min_segment_length, seg.max_segment_length)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Params;

    /// The slice form of `comp_tstat` is the indexed loop it replaced, bit for bit.
    #[test]
    fn comp_tstat_matches_the_indexed_loop() {
        fn reference(prefix: &[f64], prefix_sq: &[f64], n: usize, w_len: usize) -> Vec<f32> {
            let eta = f32::MIN_POSITIVE as f64;
            let mut tstat = vec![0.0f32; n + 1];
            if n < 2 * w_len || w_len < 2 {
                return tstat;
            }
            let w = w_len as f64;
            for i in w_len..=(n - w_len) {
                let sum1 = prefix[i] - prefix[i - w_len];
                let sumsq1 = prefix_sq[i] - prefix_sq[i - w_len];
                let sum2 = prefix[i + w_len] - prefix[i];
                let sumsq2 = prefix_sq[i + w_len] - prefix_sq[i];
                let mean1 = sum1 / w;
                let mean2 = sum2 / w;
                let combined_var = (sumsq1 / w - mean1 * mean1 + sumsq2 / w - mean2 * mean2) / w;
                let combined_var = combined_var.max(eta);
                tstat[i] = ((mean2 - mean1).abs() / combined_var.sqrt()) as f32;
            }
            for v in tstat.iter_mut().take(n + 1).skip(n - w_len + 1) {
                *v = 0.0;
            }
            tstat
        }
        let mut h: u64 = 0x1234_5678_9abc_def1;
        for n in [0usize, 1, 5, 6, 11, 12, 13, 100, 1000, 4097] {
            let sig: Vec<f32> = (0..n)
                .map(|_| {
                    h ^= h << 13;
                    h ^= h >> 7;
                    h ^= h << 17;
                    ((h >> 40) as f32 / (1u64 << 24) as f32 - 0.5) * 6.0
                })
                .collect();
            let (p, q) = prefix_sums(&sig);
            for w in [1usize, 2, 3, 6, 9] {
                let (a, b) = (comp_tstat(&p, &q, n, w), reference(&p, &q, n, w));
                assert_eq!(
                    a.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    b.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    "n={n} w={w}"
                );
            }
        }
    }

    #[test]
    fn empty_signal_yields_no_events() {
        let seg = Params::r10_4_1().segmentation;
        assert!(detect_events(&[], &seg).is_empty());
    }

    #[test]
    fn flat_signal_yields_no_events() {
        let seg = Params::r10_4_1().segmentation;
        // A constant track has no t-stat peaks.
        let flat = vec![0.0f32; 1000];
        assert!(detect_events(&flat, &seg).is_empty());
    }

    #[test]
    fn step_signal_is_segmented() {
        let seg = Params::r10_4_1().segmentation;
        // Several clean level changes -> the detector should find segments and
        // the per-segment means should track the underlying levels.
        let mut sig = Vec::new();
        for &level in &[-1.5f32, 1.5, -0.5, 1.0, -1.0] {
            for _ in 0..60 {
                sig.push(level);
            }
        }
        let events = detect_events(&sig, &seg);
        assert!(
            !events.is_empty(),
            "expected the step signal to produce events"
        );
        // Every event mean must lie within the signal's range.
        for &e in &events {
            assert!((-2.0..=2.0).contains(&e), "event mean {e} out of range");
        }
    }

    #[test]
    fn filtered_mean_ignores_outliers() {
        let seg = [1.0f32, 1.0, 1.0, 1.0, 1000.0];
        let m = filtered_segment_mean(&seg, &mut Vec::new());
        assert!((m - 1.0).abs() < 1e-6, "got {m}");
    }
}
