//! `haricot-ovl-filter`'s filters, dump, depth estimate and read model, shared with `haricot-layout`
//! (`#[path]`-included like `cache.rs`): the same code on a cache read from disk or on overlaps
//! held in memory.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result};
use haricot::overlap::Overlap;
use haricot::overlap::read_model::{ReadModelParams, read_model};

/// Read name -> read number.
pub fn name_index(names: &[String]) -> HashMap<&str, u32> {
    names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect()
}

/// `--drop-reads`: one read NAME per line; the matched reads' numbers, and how many lines matched
/// and did not.
pub fn parse_drop_reads(txt: &str, idx: &HashMap<&str, u32>) -> (HashSet<u32>, usize, usize) {
    let mut drop_reads: HashSet<u32> = HashSet::new();
    let (mut hit, mut miss) = (0usize, 0usize);
    for line in txt.lines() {
        let n = line.trim();
        if n.is_empty() {
            continue;
        }
        match idx.get(n) {
            Some(&i) => {
                drop_reads.insert(i);
                hit += 1;
            }
            None => miss += 1,
        }
    }
    (drop_reads, hit, miss)
}

/// `--drop-pairs`: `name_a<TAB>name_b` per line, in either order; the matched pairs as (smaller,
/// larger) read numbers, and how many lines matched and did not.
pub fn parse_drop_pairs(
    txt: &str,
    idx: &HashMap<&str, u32>,
) -> (HashSet<(u32, u32)>, usize, usize) {
    let mut drop_pairs: HashSet<(u32, u32)> = HashSet::new();
    let (mut hit, mut miss) = (0usize, 0usize);
    for line in txt.lines() {
        let mut f = line.split('\t');
        let (Some(x), Some(y)) = (f.next(), f.next()) else {
            continue;
        };
        match (idx.get(x.trim()), idx.get(y.trim())) {
            (Some(&i), Some(&j)) => {
                drop_pairs.insert(if i < j { (i, j) } else { (j, i) });
                hit += 1;
            }
            _ => miss += 1,
        }
    }
    (drop_pairs, hit, miss)
}

/// The subtractive filters, in the order they are tested.
pub struct Filters {
    pub drop_reads: HashSet<u32>,
    pub drop_pairs: HashSet<(u32, u32)>,
    /// `Overlap::z` is NaN for a synthesized overlap, which means "not measured" and must pass
    /// every gate. Dropping on `z < cut` keeps NaN (every comparison with it is false); keeping on
    /// `z >= cut` would silently delete those.
    pub min_z: f32,
    /// wtdbg2's D: reject iff BOTH reads overhang on the left AND both on the right (-1: off).
    pub max_oh: i64,
    pub scores: HashMap<(u32, u32), f64>,
    pub min_score: f64,
}

impl Default for Filters {
    fn default() -> Self {
        Filters {
            drop_reads: HashSet::new(),
            drop_pairs: HashSet::new(),
            min_z: f32::NEG_INFINITY,
            max_oh: -1,
            scores: HashMap::new(),
            min_score: f64::NEG_INFINITY,
        }
    }
}

/// The overlaps each filter dropped (each counted by the first filter that drops it).
#[derive(Default)]
pub struct Tally {
    pub by_read: usize,
    pub by_pair: usize,
    pub by_z: usize,
    pub by_oh: usize,
    pub by_score: usize,
}

impl Filters {
    /// Whether `o` survives every filter; a dropped overlap is counted in `t` by its reason.
    pub fn keep(&self, o: &Overlap, counts: &[u32], t: &mut Tally) -> bool {
        if self.drop_reads.contains(&o.a) || self.drop_reads.contains(&o.b) {
            t.by_read += 1;
            return false;
        }
        if self.drop_pairs.contains(&(o.a, o.b)) {
            t.by_pair += 1;
            return false;
        }
        if o.z < self.min_z {
            t.by_z += 1;
            return false;
        }
        if self.max_oh >= 0 {
            let bin = o.bin_events.max(1) as i64;
            let a_bins = counts[o.a as usize] as i64 / bin;
            let b_bins = counts[o.b as usize] as i64 / bin;
            let left = (o.a_bin_span.0 as i64).min(o.b_bin_span.0 as i64);
            let right = (a_bins - o.a_bin_span.1 as i64).min(b_bins - o.b_bin_span.1 as i64);
            if left.max(right) > self.max_oh {
                t.by_oh += 1;
                return false;
            }
        }
        if !self.scores.is_empty()
            && self.scores.get(&(o.a, o.b)).copied().unwrap_or(f64::NAN) < self.min_score
        {
            t.by_score += 1;
            return false;
        }
        true
    }
}

impl Tally {
    /// The report line: `before -> after overlaps (dropped ...)`.
    pub fn line(&self, before: usize, after: usize) -> String {
        format!(
            "{before} -> {after} overlaps (dropped {} by read, {} by pair, {} by z, \
             {} by overhang, {} by score)",
            self.by_read, self.by_pair, self.by_z, self.by_oh, self.by_score
        )
    }
}

/// `--dump`: the overlaps as a TSV keyed by read NAME, with the aligned spans and each read's event
/// count.
pub fn dump_tsv(ovl: &[Overlap], counts: &[u32], names: &[String]) -> String {
    let mut s = String::from(
        "a\tb\tz\tn_anchors\tmatched\tcnt\taligned_bins\tgap\tbin_events\tslope\
         \ta_lo\ta_hi\tb_lo\tb_hi\ta_ev\tb_ev\n",
    );
    for o in ovl {
        s.push_str(&format!(
            "{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            names[o.a as usize],
            names[o.b as usize],
            o.z,
            o.n_anchors,
            o.matched,
            o.cnt,
            o.aligned_bins,
            o.gap,
            o.bin_events,
            o.slope,
            o.a_bin_span.0,
            o.a_bin_span.1,
            o.b_bin_span.0,
            o.b_bin_span.1,
            counts[o.a as usize],
            counts[o.b as usize]
        ));
    }
    s
}

/// `--read-model`: cut reads at their junctions and keep each read's best-supported segment
/// (`--read-model-tsv` writes one row per cut read to `tsv`; `--agree REF.bin`, `agree`, counts the
/// dropped pairs a reference cache also lacks).
pub fn read_model_apply(
    ovl: Vec<Overlap>,
    counts: &[u32],
    names: &[String],
    tsv: Option<&str>,
    agree: Option<&str>,
) -> Result<Vec<Overlap>> {
    let rm = read_model(&ovl, counts, &ReadModelParams::default());
    let crossing: usize = rm.clips.iter().map(|c| c.crossing).sum();
    eprintln!(
        "read_model: {} reads with overlaps; {} clipped ({} at >= 2 junctions), {} junctions, \
         {} protected by containment; {} of {} overlaps dropped ({} across a junction)",
        rm.reads,
        rm.clips.len(),
        rm.multi(),
        rm.junctions(),
        rm.contained_protected,
        rm.dropped(),
        ovl.len(),
        crossing
    );
    if let Some(p) = tsv {
        let mut s = String::from(
            "read\tn_bins\tn_junctions\tjunctions\tkept\tkept_lo\tkept_hi\tdropped\tcrossing\n",
        );
        for c in &rm.clips {
            let j: Vec<String> = c
                .junctions
                .iter()
                .map(|(a, b)| format!("{a}-{b}"))
                .collect();
            let hi = (c.n_bins as i64 - 1).max(0);
            s.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                names[c.read as usize],
                c.n_bins,
                c.junctions.len(),
                j.join(","),
                c.kept,
                c.kept_bins.0.clamp(0, hi),
                c.kept_bins.1.clamp(0, hi),
                c.dropped,
                c.crossing
            ));
        }
        std::fs::write(p, s).with_context(|| format!("writing {p}"))?;
    }
    if let Some(p) = agree {
        let bytes = std::fs::read(p).with_context(|| format!("reading {p}"))?;
        let (_, rnames, rovl): (Vec<u32>, Vec<String>, Vec<Overlap>) =
            bincode::deserialize(&bytes).with_context(|| format!("decoding {p}"))?;
        let key = |x: &str, y: &str| -> (String, String) {
            if x < y {
                (x.to_string(), y.to_string())
            } else {
                (y.to_string(), x.to_string())
            }
        };
        let in_ref: HashSet<(String, String)> = rovl
            .iter()
            .map(|o| key(&rnames[o.a as usize], &rnames[o.b as usize]))
            .collect();
        let pair = |o: &Overlap| key(&names[o.a as usize], &names[o.b as usize]);
        let ref_dropped = ovl.iter().filter(|o| !in_ref.contains(&pair(o))).count();
        let ours = ovl.iter().zip(&rm.keep).filter(|(_, k)| !**k);
        let (mut n_ours, mut both) = (0usize, 0usize);
        for (o, _) in ours {
            n_ours += 1;
            both += usize::from(!in_ref.contains(&pair(o)));
        }
        eprintln!(
            "agree: read model dropped {n_ours} pairs; {p} lacks {ref_dropped} of the input's \
             pairs; both {both} ({:.1}% of ours, {:.1}% of the reference's)",
            100.0 * both as f64 / n_ours.max(1) as f64,
            100.0 * both as f64 / ref_dropped.max(1) as f64
        );
    }
    Ok(ovl
        .into_iter()
        .zip(rm.keep)
        .filter_map(|(o, k)| k.then_some(o))
        .collect())
}

/// numpy's median: the middle value, or the mean of the two middle values.
pub fn np_median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// The label-free depth estimate: per read (each read in an overlap, with
/// `counts[r] / bin_events + 1` bins), 1 + the overlaps whose span covers each bin; each read's
/// median; the median over reads. The bin size is the first overlap's (192 without overlaps).
pub fn depth_estimate(ovl: &[&Overlap], counts: &[u32]) -> f64 {
    let be = ovl.first().map_or(192, |o| o.bin_events);
    let mut diff: HashMap<u32, Vec<i64>> = HashMap::new();
    for o in ovl {
        for (r, (lo, hi)) in [(o.a, o.a_bin_span), (o.b, o.b_bin_span)] {
            let d = diff
                .entry(r)
                .or_insert_with(|| vec![0; (counts[r as usize] / be + 1) as usize + 1]);
            let last = d.len() - 1;
            d[lo as usize] += 1;
            d[(hi as usize + 1).min(last)] -= 1;
        }
    }
    let mut med: Vec<f64> = diff
        .values()
        .map(|d| {
            let mut run = 0i64;
            let mut cov: Vec<f64> = d[..d.len() - 1]
                .iter()
                .map(|&x| {
                    run += x;
                    (1 + run) as f64
                })
                .collect();
            np_median(&mut cov)
        })
        .collect();
    if med.is_empty() {
        return f64::NAN;
    }
    np_median(&mut med)
}

/// `x` as Python's `f"{x:.1f}"`: the exact value rounded to one decimal, ties to even. A median
/// of medians of bin counts is a multiple of 1/4, where this is exact integer arithmetic.
pub fn py_fixed1(x: f64) -> String {
    let q = x * 4.0;
    if !(q.is_finite() && q.fract() == 0.0 && q.abs() < 1e15) {
        return format!("{x:.1}");
    }
    // tenths = 10 x = 2.5 q: an integer, or a tie half-way between two (q odd).
    let q = q as i64;
    let (t, tie) = ((10 * q).div_euclid(4), (10 * q).rem_euclid(4) == 2);
    let t = if tie && t % 2 != 0 { t + 1 } else { t };
    let sign = if t < 0 { "-" } else { "" };
    format!("{sign}{}.{}", t.abs() / 10, t.abs() % 10)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// As Python prints `f"{x:.1f}"` (`python3 -c 'print(f"{12.25:.1f}")'` -> 12.2).
    #[test]
    fn one_decimal_as_python() {
        for (x, want) in [
            (12.25, "12.2"),
            (12.75, "12.8"),
            (12.0, "12.0"),
            (32.0, "32.0"),
            (6.5, "6.5"),
            (0.25, "0.2"),
            (0.75, "0.8"),
            (19.25, "19.2"),
            (8.5, "8.5"),
            (2.0 / 3.0, "0.7"),
        ] {
            assert_eq!(py_fixed1(x), want, "{x}");
        }
    }
}
