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
pub use io::{AdcRead, RawRead, for_each_adc_chunk, for_each_raw_chunk, for_each_raw_read};
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
    events::detect_events(&norm, &params.segmentation)
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
/// Means are identical to [`signal_to_events`]'s.
pub fn signal_to_event_spans(signal_pa: &[f32], params: &Params) -> Vec<events::EventSpan> {
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
}
