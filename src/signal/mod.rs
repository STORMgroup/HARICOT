//! Raw-signal processing: I/O, normalization, and event segmentation.
//!
//! Pipeline for a single read (all stages follow RawHash2's R10.4.1 approach):
//! ```text
//!   pA signal ──pa_bandpass──▶ z-score normalize (drop |z|≥3) ──▶ dual t-test segmentation ──▶ event means
//! ```
//! The output (one mean current value per event) is the input to the e-mer
//! sketching stage in [`crate::emer`].

pub mod events;
pub mod io;
pub mod normalize;

pub use events::detect_events;
pub use io::{
    AdcRead, RawRead, for_each_adc_chunk, for_each_raw_chunk, for_each_raw_read, read_fast5,
    read_pod5, read_pod5_limit, read_signals, read_signals_limit,
};
pub use normalize::normalize;

use crate::params::Params;

/// Keep only samples within the inclusive pA band (RawHash2 uses 30..=200 pA).
pub fn pa_bandpass(signal: &[f32], lo: f32, hi: f32) -> Vec<f32> {
    signal
        .iter()
        .copied()
        .filter(|&x| x >= lo && x <= hi)
        .collect()
}

/// Full signal → events pipeline for one read.
///
/// Optional pA band-pass → z-score normalization (dropping |z| ≥ 3 outliers) →
/// dual t-test event segmentation. Returns one mean current value per event.
pub fn signal_to_events(signal_pa: &[f32], params: &Params) -> Vec<f32> {
    let filtered = match params.signal.pa_filter {
        Some((lo, hi)) => pa_bandpass(signal_pa, lo, hi),
        None => signal_pa.to_vec(),
    };
    let norm = normalize::normalize(&filtered);
    let evts = events::detect_events(&norm, &params.segmentation);
    if params.signal.event_merge > 0.0 {
        merge_events(&evts, params.signal.event_merge)
    } else {
        evts
    }
}

/// [`signal_to_events`] with each event's span expressed in **raw signal samples**, i.e.
/// indices into `signal_pa`.
///
/// Both filters ahead of segmentation drop samples (the pA band-pass, then the |z| >= 3
/// cut in [`normalize::normalize`]), so a position in the segmented stream is not a
/// position in the read's signal. This carries the kept indices through both, which is
/// what lets an event be tied to a basecaller move table and from there to a read base.
/// `end` is one past the event's last kept sample; dropped samples inside an event are
/// spanned, not counted.
///
/// Means are identical to [`signal_to_events`]'s. `event_merge > 0` averages events into
/// runs and is not supported here.
pub fn signal_to_event_spans(signal_pa: &[f32], params: &Params) -> Vec<events::EventSpan> {
    assert!(
        params.signal.event_merge <= 0.0,
        "signal_to_event_spans does not track event_merge runs"
    );
    let (filtered, band_kept): (Vec<f32>, Vec<u32>) = match params.signal.pa_filter {
        Some((lo, hi)) => signal_pa
            .iter()
            .enumerate()
            .filter(|&(_, &x)| x >= lo && x <= hi)
            .map(|(i, &x)| (x, i as u32))
            .unzip(),
        None => (signal_pa.to_vec(), (0..signal_pa.len() as u32).collect()),
    };
    let (norm, norm_kept) = normalize::normalize_kept(&filtered);
    let raw = |i: u32| band_kept[norm_kept[i as usize] as usize];
    events::detect_event_spans(&norm, &params.segmentation)
        .into_iter()
        .map(|e| events::EventSpan {
            start: raw(e.start),
            end: raw(e.end - 1) + 1,
            mean: e.mean,
        })
        .collect()
}

/// Rawsamble/RawHash2-style pre-quantization event merge: collapse consecutive events
/// whose z-normalized means stay within `threshold` (z-units) of the running mean into
/// a single (averaged) event. `threshold <= 0` is a no-op. See
/// [`SignalParams::event_merge`](crate::params::SignalParams::event_merge).
pub fn merge_events(events: &[f32], threshold: f32) -> Vec<f32> {
    if threshold <= 0.0 || events.is_empty() {
        return events.to_vec();
    }
    let mut out = Vec::with_capacity(events.len());
    let mut sum = events[0] as f64;
    let mut cnt = 1u32;
    for &v in &events[1..] {
        let mean = sum / cnt as f64;
        if (v as f64 - mean).abs() <= threshold as f64 {
            sum += v as f64;
            cnt += 1;
        } else {
            out.push((sum / cnt as f64) as f32);
            sum = v as f64;
            cnt = 1;
        }
    }
    out.push((sum / cnt as f64) as f32);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bandpass_keeps_only_in_range() {
        let sig = [10.0, 30.0, 100.0, 200.0, 250.0];
        let out = pa_bandpass(&sig, 30.0, 200.0);
        assert_eq!(out, vec![30.0, 100.0, 200.0]);
    }

    #[test]
    fn merge_collapses_similar_runs() {
        // 0.0,0.1,0.2 stay within 0.3 of the running mean → one event; 2.0 breaks it.
        let m = merge_events(&[0.0, 0.1, 0.2, 2.0, 2.1], 0.3);
        assert_eq!(m.len(), 2, "expected two runs, got {m:?}");
        assert!((m[0] - 0.1).abs() < 1e-5, "mean of 0,0.1,0.2 = {}", m[0]);
        assert!((m[1] - 2.05).abs() < 1e-5, "mean of 2.0,2.1 = {}", m[1]);
    }

    #[test]
    fn event_spans_match_events_and_map_to_raw_samples() {
        let params = Params::r10_4_1();
        // Clean level steps inside the pA band, with out-of-band samples and a spike
        // sprinkled in so both filters actually drop something.
        let mut sig = Vec::new();
        for (k, &level) in [60.0f32, 110.0, 80.0, 140.0, 70.0, 120.0]
            .iter()
            .enumerate()
        {
            for j in 0..80 {
                sig.push(if j == 13 && k % 2 == 0 { 250.0 } else { level });
            }
        }
        sig.insert(0, 10.0);
        let means = signal_to_events(&sig, &params);
        let spans = signal_to_event_spans(&sig, &params);
        assert!(!spans.is_empty());
        assert_eq!(spans.iter().map(|e| e.mean).collect::<Vec<_>>(), means);
        for w in spans.windows(2) {
            assert!(w[0].start < w[0].end && w[0].end <= w[1].start);
        }
        for e in &spans {
            // Every event lies on in-band samples at both ends.
            assert!((30.0..=200.0).contains(&sig[e.start as usize]));
            assert!((30.0..=200.0).contains(&sig[e.end as usize - 1]));
        }
    }

    #[test]
    fn merge_off_is_identity() {
        assert_eq!(merge_events(&[0.0, 1.0, 2.0], 0.0), vec![0.0, 1.0, 2.0]);
    }
}
