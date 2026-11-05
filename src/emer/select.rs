//! Indel-robust seed *selection* — syncmers and strobemers.
//!
//! Motivation: raw signal can segment the same DNA
//! into a **different number of events** on different reads, which acts like
//! pervasive insertion/deletion of events. How much this bites is unsettled — the
//! debug pair's diagonal sits at slope ≈1.3 (a ~30% mismatch), but over the long
//! fixture's 741 true pairs the median is 0.995, so treat ≈1.3 as a property of
//! that pair rather than of the chemistry. Classic `(w, e)` window minimizers
//! pick a per-window arg-min, which can shift to a different event whenever the
//! window contents change, so two reads over the same DNA need not select the same
//! representatives. This is a robustness margin, not a rescue from zero: measured
//! on the debug pair, `w=2` minimizers still yield ~20k anchors, and a bin-DP at the
//! pair's own slope recovers ~0.99 of the contained read (`bin=128` events,
//! `kmax=16`). The pipeline uses `w=2` minimizers.
//!
//! **Syncmers** and **strobemers** select by the seed's own content instead of a
//! neighbourhood window, so the choice is stable as long as the local sequence is:
//!   - an **open syncmer** keeps an e-mer iff the minimum *s*-mer inside it sits at
//!     a fixed offset `t` — a property of the e-mer alone (Edgar 2021);
//!   - a **randstrobe** links a first *s*-mer to a second one chosen by content
//!     from a downstream window, so a single indel between the strobes is absorbed
//!     rather than breaking the seed (Sahlin 2021).
//!
//! All selectors run on the retained quantized-event stream
//! ([`crate::emer::sketch::retained_quantized`]) and emit [`Emer`]s. Syncmer
//! hashes equal the dense e-mer hash at that position (they *are* dense e-mers,
//! sub-sampled), so anchoring is unchanged; strobemers emit a combined hash.

use crate::emer::sketch::{
    Emer, HASH_MASK, QEvent, dense_emers_from_retained, emer_hash, hash64, minimizers,
};
use crate::params::Params;

/// Which seeding scheme to use when sketching a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SketchMethod {
    /// Every e-mer (no sub-sampling) — the sensitivity ceiling.
    Dense,
    /// Classic `(w, e)` window minimizers (the current default; the baseline).
    Minimizer,
    /// Open syncmer: keep an e-mer iff its minimum s-mer is at offset `t`.
    OpenSyncmer,
    /// Closed syncmer: keep iff the minimum s-mer is at the first or last offset.
    ClosedSyncmer,
    /// Order-2 randstrobe (two linked s-mers).
    Strobemer,
}

impl SketchMethod {
    pub const ALL: [SketchMethod; 5] = [
        SketchMethod::Dense,
        SketchMethod::Minimizer,
        SketchMethod::OpenSyncmer,
        SketchMethod::ClosedSyncmer,
        SketchMethod::Strobemer,
    ];

    pub fn label(self) -> &'static str {
        match self {
            SketchMethod::Dense => "dense",
            SketchMethod::Minimizer => "minimizer",
            SketchMethod::OpenSyncmer => "open-syncmer",
            SketchMethod::ClosedSyncmer => "closed-syncmer",
            SketchMethod::Strobemer => "strobemer",
        }
    }

    /// Parse a method name (accepts `-`/`_` and a few aliases).
    pub fn parse(s: &str) -> Result<SketchMethod, String> {
        match s.to_ascii_lowercase().replace('_', "-").as_str() {
            "dense" => Ok(SketchMethod::Dense),
            "minimizer" | "min" => Ok(SketchMethod::Minimizer),
            "open-syncmer" | "syncmer" | "open" => Ok(SketchMethod::OpenSyncmer),
            "closed-syncmer" | "closed" => Ok(SketchMethod::ClosedSyncmer),
            "strobemer" | "strobe" | "randstrobe" => Ok(SketchMethod::Strobemer),
            other => Err(format!(
                "unknown sketch method {other:?} (expected one of: {})",
                Self::ALL.map(|m| m.label()).join(", ")
            )),
        }
    }
}

/// Tuning for the selectors, kept separate from [`crate::params::Params`] while
/// they are not in any verified preset (the pipeline uses minimizers).
#[derive(Debug, Clone, Copy)]
pub struct SelectCfg {
    /// e-mer length in events (matches `Params::sketch.e`).
    pub e: usize,
    /// bits per quantized event (matches `Params::quant.q`).
    pub q: u32,
    /// Minimizer window in e-mers (for [`SketchMethod::Minimizer`]).
    pub w: usize,
    /// s-mer length for syncmers/strobemers (`s < e`).
    pub s: usize,
    /// Target offset of the min s-mer for an open syncmer (`0..=e-s`).
    pub t: usize,
    /// Strobe link window `[w_min, w_max]` in s-mer steps (for strobemers).
    pub w_min: usize,
    pub w_max: usize,
}

impl SelectCfg {
    /// Sensible defaults derived from `(e, q)`: `s = e/2`, open-syncmer target at
    /// the centre, strobe window `[1, max(e, 2)]`.
    pub fn from_eq(e: usize, q: u32) -> Self {
        let s = (e / 2).max(1);
        SelectCfg {
            e,
            q,
            w: 2,
            s,
            t: (e - s) / 2,
            w_min: 1,
            w_max: e.max(2),
        }
    }

    /// Build from a [`Params`], honouring `sketch.{e,w,s}` and `quant.q`. `s` is
    /// clamped to `1..=e-1` so syncmers/strobemers have room inside an e-mer.
    pub fn from_params(p: &Params) -> Self {
        let e = p.sketch.e as usize;
        let mut cfg = SelectCfg::from_eq(e, p.quant.q);
        cfg.w = p.sketch.w as usize;
        cfg.s = (p.sketch.s as usize).clamp(1, e.saturating_sub(1).max(1));
        cfg.t = (e - cfg.s) / 2;
        cfg.w_max = e.max(2);
        cfg
    }
}

/// Hash of every `s`-mer in the quantized stream: `out[i]` covers `buckets[i..i+s]`.
fn smer_hashes(buckets: &[u32], s: usize, q: u32) -> Vec<u64> {
    if s == 0 || buckets.len() < s {
        return Vec::new();
    }
    (0..=buckets.len() - s)
        .map(|i| emer_hash(&buckets[i..i + s], q))
        .collect()
}

/// Open/closed syncmers over the retained stream.
///
/// An e-mer at retained index `i` contains the s-mers at offsets `0..=e-s`. We
/// take the offset of the minimum-hash s-mer (ties → smallest offset, so the rule
/// is identical on both reads). Open: keep iff that offset == `t`. Closed: keep
/// iff it is `0` or `e-s` (the endpoints), which is denser but symmetric.
pub fn syncmers(retained: &[QEvent], cfg: &SelectCfg, closed: bool) -> Vec<Emer> {
    let (e, s, q) = (cfg.e, cfg.s, cfg.q);
    if s >= e || s == 0 || retained.len() < e {
        return Vec::new();
    }
    let buckets: Vec<u32> = retained.iter().map(|r| r.bucket).collect();
    let sh = smer_hashes(&buckets, s, q);
    let span = e - s; // last s-mer offset inside an e-mer
    let t = cfg.t.min(span);
    let mut out = Vec::new();
    for i in 0..=retained.len() - e {
        // argmin over s-mer offsets 0..=span (first-wins tie-break).
        let mut best_off = 0usize;
        let mut best_h = sh[i];
        for off in 1..=span {
            if sh[i + off] < best_h {
                best_h = sh[i + off];
                best_off = off;
            }
        }
        let keep = if closed {
            best_off == 0 || best_off == span
        } else {
            best_off == t
        };
        if keep {
            out.push(Emer {
                hash: emer_hash(&buckets[i..i + e], q),
                pos: retained[i].pos,
                event_pos: retained[i].event_pos,
            });
        }
    }
    out
}

/// Order-2 randstrobes over the retained stream (dense in start position).
///
/// Strobe 1 is the s-mer at `i`; strobe 2 is the s-mer `j` in the window
/// `[i+w_min, i+w_max]` minimizing the link `(h1 + h_j) mod 2^32` (Sahlin's
/// randstrobe rule). The seed hash mixes both strobes; `pos` is strobe 1's start,
/// so it is comparable to dense/syncmer positions. Start positions near the end
/// of the read with no room for the window are dropped.
pub fn randstrobes(retained: &[QEvent], cfg: &SelectCfg) -> Vec<Emer> {
    let (s, q) = (cfg.s, cfg.q);
    let buckets: Vec<u32> = retained.iter().map(|r| r.bucket).collect();
    let sh = smer_hashes(&buckets, s, q);
    let m = sh.len();
    if m == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for i in 0..m {
        let lo = i + cfg.w_min;
        let hi = (i + cfg.w_max).min(m - 1);
        if lo > hi || lo >= m {
            continue; // no room for a second strobe
        }
        let h1 = sh[i];
        let mut best_j = lo;
        let mut best_link = u64::MAX;
        for (j, &hj) in sh.iter().enumerate().take(hi + 1).skip(lo) {
            let link = h1.wrapping_add(hj) & HASH_MASK;
            if link < best_link {
                best_link = link;
                best_j = j;
            }
        }
        let h2 = sh[best_j];
        // Order-sensitive mix so (h1,h2) != (h2,h1); stays in 32 bits.
        let combined = (h1 ^ h2.rotate_left(13)) & HASH_MASK;
        out.push(Emer {
            hash: hash64(combined, HASH_MASK),
            pos: retained[i].pos,
            event_pos: retained[i].event_pos,
        });
    }
    out
}

/// Dispatch: select seeds from a retained quantized stream by `method`.
pub fn select(retained: &[QEvent], method: SketchMethod, cfg: &SelectCfg) -> Vec<Emer> {
    match method {
        SketchMethod::Dense => dense_emers_from_retained(retained, cfg.e, cfg.q),
        SketchMethod::Minimizer => {
            let dense = dense_emers_from_retained(retained, cfg.e, cfg.q);
            minimizers(&dense, cfg.w)
        }
        SketchMethod::OpenSyncmer => syncmers(retained, cfg, false),
        SketchMethod::ClosedSyncmer => syncmers(retained, cfg, true),
        SketchMethod::Strobemer => randstrobes(retained, cfg),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a retained stream from explicit buckets. `pos` is the coordinate
    /// (index in this retained stream); `event_pos` is the raw index it came from —
    /// equal here only because the helper fabricates an already-merged stream.
    fn stream(buckets: &[u32]) -> Vec<QEvent> {
        buckets
            .iter()
            .enumerate()
            .map(|(i, &b)| QEvent {
                pos: i as u32,
                event_pos: i as u32,
                bucket: b,
            })
            .collect()
    }

    /// A deterministic pseudo-random bucket sequence over `0..n_buckets`.
    fn pseudo(len: usize, n_buckets: u32) -> Vec<u32> {
        let mut x: u64 = 0x9E3779B97F4A7C15;
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                (x % n_buckets as u64) as u32
            })
            .collect()
    }

    fn cfg() -> SelectCfg {
        SelectCfg::from_eq(8, 3)
    }

    #[test]
    fn syncmers_are_a_subset_of_dense_with_matching_hashes() {
        let buckets = pseudo(500, 8);
        let r = stream(&buckets);
        let dense = dense_emers_from_retained(&r, 8, 3);
        let sync = syncmers(&r, &cfg(), false);
        assert!(!sync.is_empty() && sync.len() < dense.len());
        use std::collections::HashMap;
        let dense_at: HashMap<u32, u64> = dense.iter().map(|e| (e.pos, e.hash)).collect();
        for s in &sync {
            // every syncmer coincides with the dense e-mer at the same position
            assert_eq!(dense_at.get(&s.pos), Some(&s.hash));
        }
    }

    #[test]
    fn open_syncmer_density_is_about_one_over_window() {
        // Expected fraction kept ≈ 1/(e-s+1) for a uniform-ish stream.
        let buckets = pseudo(20_000, 8);
        let r = stream(&buckets);
        let c = cfg();
        let sync = syncmers(&r, &c, false);
        let dense = dense_emers_from_retained(&r, c.e, c.q);
        let frac = sync.len() as f64 / dense.len() as f64;
        let expected = 1.0 / (c.e - c.s + 1) as f64;
        assert!(
            (frac - expected).abs() < 0.5 * expected,
            "syncmer fraction {frac:.3} not near 1/(e-s+1)={expected:.3}"
        );
    }

    /// Under synthetically inserted events, syncmers retain more shared seed
    /// *hashes* between the original and perturbed stream than minimizers do —
    /// the robustness property these selectors exist for.
    #[test]
    fn syncmers_survive_an_indel_better_than_minimizers() {
        use std::collections::HashSet;
        let base = pseudo(4000, 8);
        // Insert a spurious event every ~13 positions: a ~7% event-rate inflation,
        // milder than the real 30% but enough to show the effect.
        let mut perturbed = Vec::with_capacity(base.len() + base.len() / 13);
        for (i, &b) in base.iter().enumerate() {
            perturbed.push(b);
            if i % 13 == 7 {
                perturbed.push((b + 3) % 8); // an inserted, different event
            }
        }
        let (r0, r1) = (stream(&base), stream(&perturbed));
        let c = cfg();

        let shared = |a: &[Emer], b: &[Emer]| -> usize {
            let sa: HashSet<u64> = a.iter().map(|e| e.hash).collect();
            let sb: HashSet<u64> = b.iter().map(|e| e.hash).collect();
            sa.intersection(&sb).count()
        };

        let min0 = select(&r0, SketchMethod::Minimizer, &c);
        let min1 = select(&r1, SketchMethod::Minimizer, &c);
        let syn0 = select(&r0, SketchMethod::OpenSyncmer, &c);
        let syn1 = select(&r1, SketchMethod::OpenSyncmer, &c);

        // Normalize by seed count so the comparison is at matched-ish density.
        let min_rate = shared(&min0, &min1) as f64 / min0.len().min(min1.len()) as f64;
        let syn_rate = shared(&syn0, &syn1) as f64 / syn0.len().min(syn1.len()) as f64;
        assert!(
            syn_rate > min_rate,
            "expected syncmers to share more under indels: syncmer {syn_rate:.3} vs minimizer {min_rate:.3}"
        );
    }

    #[test]
    fn randstrobes_are_deterministic_and_nonempty() {
        let buckets = pseudo(1000, 8);
        let r = stream(&buckets);
        let a = randstrobes(&r, &cfg());
        let b = randstrobes(&r, &cfg());
        assert!(!a.is_empty());
        assert_eq!(a.len(), b.len());
        assert!(
            a.iter()
                .zip(&b)
                .all(|(x, y)| x.hash == y.hash && x.pos == y.pos)
        );
    }

    #[test]
    fn short_streams_yield_no_seeds() {
        let r = stream(&[1, 2, 3]); // shorter than e=8
        let c = cfg();
        assert!(syncmers(&r, &c, false).is_empty());
        assert!(select(&r, SketchMethod::OpenSyncmer, &c).is_empty());
    }
}
