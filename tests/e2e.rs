//! End-to-end tests on `tests/fixtures/debug_reads.pod5`: the real signal → events → seeds → overlap
//! path on two high-quality R10.4.1 reads. The file is not distributed (`tests/fixtures/README.md`);
//! without it these tests are reported as ignored (`build.rs`).

use std::path::Path;

use haricot::params::Params;
use haricot::{pipeline, signal};

const DEBUG_POD5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/debug_reads.pod5"
);

fn fixture() -> Vec<signal::RawRead> {
    let mut reads = Vec::new();
    signal::for_each_raw_read(Path::new(DEBUG_POD5), None, |r| {
        reads.push(r);
        Ok(())
    })
    .expect("reading tests/fixtures/debug_reads.pod5");
    reads
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
fn production_detection_finds_the_debug_overlap() {
    // The two debug reads are a known same-strand containment: read 0 (14,633 retained events, 77
    // bins of 192) lies inside read 1 (74,346). The production detection (`haricot-detect
    // --production`: the Hough overlapper, the pre-filter, the null gate, the canonical pair frame,
    // bin 192) finds it as one overlap that spans read 0, at the pair's event-rate slope ≈ 0.895.
    let reads = fixture();
    let params = pipeline::fbg_build::parse(["--production".to_string()].into_iter())
        .unwrap()
        .params;
    let sketches = pipeline::sketch_reads(&reads, &params);
    let lens: Vec<u32> = sketches.iter().map(|s| s.n_coord_events).collect();
    assert_eq!(lens, [14_633, 74_346]);
    let ovl = pipeline::overlaps(&sketches, &params);
    assert_eq!(ovl.len(), 1, "the single known overlap");
    let o = &ovl[0];
    assert_eq!((o.a, o.b), (0, 1));
    assert_eq!(o.a_bin_span, (0, 75), "read 0's bins, end to end");
    assert_eq!(o.b_bin_span, (277, 362));
    assert_eq!(o.aligned_bins, 76);
    assert!((o.slope - 0.895).abs() < 0.01, "slope {}", o.slope);
    assert!(o.z > 40.0, "z {} (the null gate keeps z >= 4)", o.z);
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
