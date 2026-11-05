//! Cross-implementation check: the Rust e-mer hashing must reproduce
//! RawHash2's reference output bit-for-bit.
//!
//! The fixtures under `tests/rawhash_compat/` are produced once by
//! `regenerate.sh`, which compiles a small C harness against RawHash2's own
//! `rsketch.c` and dumps `pos<TAB>hash` lines for a deterministic random
//! signal (see that directory's README). Here we feed the *same* signal
//! through the Rust [`sketch`] and assert the emitted `(pos, hash)` sequence is
//! identical — i.e. quantization, the sig-diff collapse, e-mer packing, and
//! `hash64` all agree with the C reference.
//!
//! The golden's `pos` is RawHash2's `f_pos`: the index in the **raw** segmented
//! stream (`rsketch.c:214-218`). We compare it against the seed's raw-stream index (the
//! retained stream's `event_pos` at the seed's `pos`; a seed does not store it), which
//! is that same quantity. [`Seed::pos`] is deliberately *not* compared: it is our
//! coordinate, the index in the merged (retained) stream, and it intentionally
//! differs from RawHash2 — see [`emer::sketch::retained_quantized`] for why. So
//! this test pins hashing fidelity (quantize → sig-diff → pack → `hash64`) and the
//! set of events surviving the filter, which is all it ever meant to pin.
//!
//! All configurations use the dense sketch (`w = 0`) so seeds come out in event
//! order, giving an exact position-by-position comparison. (Minimizer selection
//! is deliberately not compared here: `src/emer/sketch.rs` documents an
//! intentional deviation from RawHash2's duplicate-handling rotation.)
//!
//! Tests skip gracefully if the fixtures are absent, matching `tests/e2e.rs`.

use std::path::Path;

use haricot::emer::sketch;
use haricot::params::{Params, Strand};

const DIR: &str = "tests/rawhash_compat";

/// One golden configuration. Mirrors a `emit` line in `regenerate.sh`; the
/// non-listed quant fields (`q`, `n_buckets`, `fine_*`) match the R10.4.1
/// preset that the harness was invoked with.
struct Config {
    /// Golden file stem under [`DIR`] (without `.tsv`).
    name: &'static str,
    e: u32,
    diff: i32,
}

const CONFIGS: &[Config] = &[
    Config {
        name: "e8_q3_diff0",
        e: 8,
        diff: 0,
    }, // R10.4.1 default
    Config {
        name: "e8_q3_nofilter",
        e: 8,
        diff: -1,
    }, // sig-diff disabled
    Config {
        name: "e4_q3_nofilter",
        e: 4,
        diff: -1,
    }, // shorter e-mer
    Config {
        name: "e6_q3_diff0",
        e: 6,
        diff: 0,
    }, // mid e + collapse
];

fn read_signal() -> Option<Vec<f32>> {
    let path = format!("{DIR}/signal.txt");
    if !Path::new(&path).exists() {
        eprintln!("skipping hash-compat test: {path} not present");
        return None;
    }
    let text = std::fs::read_to_string(&path).expect("reading signal fixture");
    Some(
        text.lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.trim().parse::<f32>().expect("signal value is an f32"))
            .collect(),
    )
}

/// Parse a golden `pos<TAB>hash` dump.
fn read_golden(name: &str) -> Vec<(u32, u64)> {
    let path = format!("{DIR}/{name}.tsv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("reading {path}: {e}"));
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let (pos, hash) = l
                .split_once('\t')
                .unwrap_or_else(|| panic!("malformed golden line: {l:?}"));
            (
                pos.parse::<u32>().expect("pos is u32"),
                hash.parse::<u64>().expect("hash is u64"),
            )
        })
        .collect()
}

fn params_for(cfg: &Config) -> Params {
    // Start from the preset (q=3, n_buckets=8, fine_min=-2, fine_max=2,
    // fine_range=0.75, k=9) and override only what the harness varied.
    let mut p = Params::r10_4_1();
    p.sketch.e = cfg.e;
    p.sketch.w = 0; // dense sketch -> event order, exact comparison
    p.sketch.diff = cfg.diff;
    p
}

#[test]
fn rust_hashes_match_rawhash2() {
    let Some(events) = read_signal() else { return };

    for cfg in CONFIGS {
        let params = params_for(cfg);
        let seeds = sketch::sketch(&events, /*read_id=*/ 0, Strand::Forward, &params);
        // a seed's raw-stream index: the retained stream's at the seed's coordinate
        let retained = sketch::retained_quantized(&events, &params);
        let golden = read_golden(cfg.name);

        assert_eq!(
            seeds.len(),
            golden.len(),
            "[{}] seed count differs: rust={} rawhash2={}",
            cfg.name,
            seeds.len(),
            golden.len()
        );

        for (i, (seed, &(gpos, ghash))) in seeds.iter().zip(golden.iter()).enumerate() {
            // `event_pos` (raw stream index) is RawHash2's `f_pos`; `seed.pos` is our
            // merged-stream coordinate and is expected to differ. See the module doc.
            let event_pos = retained[seed.pos as usize].event_pos;
            assert_eq!(
                (event_pos, seed.hash()),
                (gpos, ghash),
                "[{}] mismatch at seed {i}: rust=(event_pos={}, hash={}) rawhash2=(pos={}, hash={})",
                cfg.name,
                event_pos,
                seed.hash(),
                gpos,
                ghash
            );
        }
    }
}
