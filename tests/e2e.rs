//! End-to-end tests on `tests/fixtures/debug_reads.pod5`: the real signal → events → seeds → overlap
//! path on two high-quality R10.4.1 reads. The file is not distributed (`tests/fixtures/README.md`);
//! without it these tests are reported as ignored (`build.rs`).

use std::path::Path;

use haricot::params::Params;
use haricot::{layout, pipeline, signal};

const DEBUG_POD5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/debug_reads.pod5"
);

fn fixture() -> Vec<signal::RawRead> {
    signal::read_pod5(Path::new(DEBUG_POD5)).expect("reading tests/fixtures/debug_reads.pod5")
}

#[test]
#[cfg_attr(
    not(haricot_debug_fixture),
    ignore = "needs tests/fixtures/debug_reads.pod5 (not distributed; tests/fixtures/README.md)"
)]
fn signal_to_seeds_on_debug_reads() {
    let reads = fixture();
    assert_eq!(reads.len(), 2, "the debug set contains two reads");

    let params = Params::r10_4_1();
    let sketches = pipeline::sketch_reads(&reads, &params);

    for s in &sketches {
        // Each read should segment into many events and yield seeds.
        assert!(
            s.n_events > 1_000,
            "read {} produced only {} events",
            s.name,
            s.n_events
        );
        assert!(!s.seeds.is_empty(), "read {} produced no seeds", s.name);
        // Minimizer sketch never yields more seeds than events.
        assert!(s.seeds.len() as u32 <= s.n_events);
    }
}

#[test]
#[cfg_attr(
    not(haricot_debug_fixture),
    ignore = "needs tests/fixtures/debug_reads.pod5 (not distributed; tests/fixtures/README.md)"
)]
fn bin_dp_detects_and_hough_recovers_debug_overlap() {
    // The two debug reads are a known same-strand containment. Its overlap diagonal
    // has slope ≈ 0.78 (a_pos vs b_pos; equivalently ≈1.28 in long/short terms) — the
    // event-rate mismatch of this one pair (it does not generalize: the median slope at
    // scale is ≈1.0). A "slope ≈ 2" figure was previously cited for this pair; it came
    // from a patience-LIS estimator that maximized chain length, and the concentration
    // estimator measures ≈0.78.
    let reads = fixture();

    // (1) Default = the wtdbg2 baseline: slope-1 DP, faithful pgap=−7/pvar=−21.
    // It detects only a short slope-1 fragment of the overlap: the off-slope (~0.78)
    // diagonal drifts out of the ±max_bvar band after a handful of bins.
    let base = Params::r10_4_1();
    let mins = pipeline::overlaps(&pipeline::sketch_reads(&reads, &base), &base);
    assert_eq!(
        mins.len(),
        1,
        "baseline detects the single known overlap (a short fragment)"
    );
    assert_eq!(mins[0].slope, 1.0, "the baseline runs at slope 1");
    assert!(mins[0].aligned_bins >= base.overlap.min_aln_bins);

    // (2) The Hough overlapper finds the diagonal's slope itself: with coarse
    // (wtdbg2-scale) 128-event bins and a dense sketch it recovers the bulk of the
    // containment, at the measured ≈0.78.
    let mut rec = Params::r10_4_1();
    rec.overlap.hough = true;
    rec.overlap.bin_events = 128;
    rec.sketch.method = haricot::emer::select::SketchMethod::Dense;
    let recov = pipeline::overlaps(&pipeline::sketch_reads(&reads, &rec), &rec);
    assert_eq!(recov.len(), 1, "the overlap is recovered");
    assert!(
        recov[0].aligned_bins >= 50,
        "the Hough fit should recover a large span, got {} bins",
        recov[0].aligned_bins
    );
    assert!(
        (recov[0].slope - 0.78).abs() < 0.12,
        "recovered slope should be ~0.78, got {}",
        recov[0].slope
    );
}

#[test]
#[cfg_attr(
    not(haricot_debug_fixture),
    ignore = "needs tests/fixtures/debug_reads.pod5 (not distributed; tests/fixtures/README.md)"
)]
fn assemble_emits_nonempty_gfa() {
    let reads = fixture();
    let params = Params::r10_4_1();
    let graph = pipeline::assemble(&reads, &params);
    let gfa = layout::to_gfa_string(&graph);
    assert!(gfa.starts_with("H\tVN:Z:1.0"));
    // The debug pair's overlap produces a non-empty FBG (a short K-bin chain).
    assert!(
        graph.n_nodes() >= 1,
        "build_graph should emit nodes from the overlap"
    );
    assert!(gfa.contains("\nS\t"), "GFA should contain segment lines");
}

/// The ADC chunks are the raw chunks: the same reads, in the same order, in the same chunks, and each
/// read's pA signal bit for bit; and the streamed sketch (over ADC chunks) is the sketch of the
/// whole read set.
#[test]
#[cfg_attr(
    not(haricot_debug_fixture),
    ignore = "needs tests/fixtures/debug_reads.pod5 (not distributed; tests/fixtures/README.md)"
)]
fn adc_chunks_are_the_raw_chunks() {
    let reads = fixture();
    let path = Path::new(DEBUG_POD5);
    for chunk_bytes in [1usize, 1 << 20, 1 << 30] {
        let mut raw: Vec<Vec<(String, Vec<u32>)>> = Vec::new();
        signal::for_each_raw_chunk(path, None, chunk_bytes, |c| {
            let bits = |r: &signal::RawRead| -> Vec<u32> {
                r.signal_pa.iter().map(|x| x.to_bits()).collect()
            };
            raw.push(c.iter().map(|r| (r.id.clone(), bits(r))).collect());
            Ok(())
        })
        .unwrap();
        let mut adc: Vec<Vec<(String, Vec<u32>)>> = Vec::new();
        signal::for_each_adc_chunk(path, None, chunk_bytes, |c| {
            let bits = |r: &signal::AdcRead| -> Vec<u32> {
                r.signal_pa().iter().map(|x| x.to_bits()).collect()
            };
            adc.push(c.iter().map(|r| (r.id.clone(), bits(r))).collect());
            Ok(())
        })
        .unwrap();
        assert_eq!(adc, raw, "chunk_bytes {chunk_bytes}");
        assert_eq!(raw.iter().map(|c| c.len()).sum::<usize>(), reads.len());
    }
    let params = Params::r10_4_1();
    let whole = pipeline::sketch_reads(&reads, &params);
    let streamed = pipeline::sketch_pod5_streaming(path, None, 1, &params).unwrap();
    assert_eq!(streamed.len(), whole.len());
    for ((s, w), r) in streamed.iter().zip(&whole).zip(&reads) {
        assert_eq!(
            (s.read_id, &s.name, s.n_samples),
            (w.read_id, &w.name, w.n_samples)
        );
        assert_eq!(
            (s.n_events, s.n_coord_events),
            (w.n_events, w.n_coord_events)
        );
        // the tie rank's key: the same pA samples hashed whichever path read them
        let h = haricot::readorder::content_hash(&r.signal_pa);
        assert_eq!((s.content, w.content), (h, h), "{}", r.id);
        let key = |x: &haricot::emer::Seed| (x.hash(), x.read_id(), x.pos, x.strand());
        assert_eq!(
            s.seeds.iter().map(key).collect::<Vec<_>>(),
            w.seeds.iter().map(key).collect::<Vec<_>>()
        );
    }
}
