//! The keep cut of `haricot-layout`, on overlaps held in memory: the forests' drop lists, the depth
//! estimate and the depth-adaptive keep rule that picks the stage-2 keep level from it (Algorithm S2,
//! lines 21-24), and the candidate dump.

use std::collections::{HashMap, HashSet};

use crate::overlap::Overlap;

/// The keep rule's name: `haricot-layout --keep-level`, and `keep_rule=` in the first line
/// `haricot assemble` prints.
pub const KEEP_RULE: &str = "depth-adaptive";

/// The keep rule's earlier name, still accepted by `haricot-layout --keep-level`.
pub const KEEP_RULE_ALIAS: &str = "policy7";

/// The depth-adaptive keep rule: the depth estimate as printed (depth.txt, [`py_fixed1`] of
/// [`depth_estimate`] over the overlaps stage 2 keeps at 0.995) picks the stage-2 keep level,
/// stricter at higher depth: below 8.5x 0.995, below 19x 0.990, else 0.980. `NaN` (no overlaps)
/// reads as 0: 0.995.
pub fn depth_adaptive_level(depth: &str) -> &'static str {
    let d: f64 = depth.trim().parse().unwrap_or(0.0);
    if d >= 19.0 {
        "0.980"
    } else if d >= 8.5 {
        "0.990"
    } else {
        "0.995"
    }
}

/// Read name -> read number.
pub fn name_index(names: &[String]) -> HashMap<&str, u32> {
    names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect()
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

/// The subtractive filter: a drop list of pairs.
#[derive(Default)]
pub struct Filters {
    pub drop_pairs: HashSet<(u32, u32)>,
}

/// The overlaps the filter dropped.
#[derive(Default)]
pub struct Tally {
    pub by_pair: usize,
}

impl Filters {
    /// Whether `o` survives the filter; a dropped overlap is counted in `t`.
    pub fn keep(&self, o: &Overlap, t: &mut Tally) -> bool {
        if self.drop_pairs.contains(&(o.a, o.b)) {
            t.by_pair += 1;
            return false;
        }
        true
    }
}

impl Tally {
    /// The report line: `before -> after overlaps (dropped ...)`.
    pub fn line(&self, before: usize, after: usize) -> String {
        format!(
            "{before} -> {after} overlaps (dropped {} by pair)",
            self.by_pair
        )
    }
}

/// The candidate dump (W/cand.dump.tsv): the overlaps as a TSV keyed by read NAME, with the aligned
/// spans and each read's event count.
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

    /// The rule on the printed depth.
    #[test]
    fn keep_level_follows_the_printed_depth() {
        assert_eq!(depth_adaptive_level("8.2"), "0.995");
        assert_eq!(depth_adaptive_level("8.5\n"), "0.990");
        assert_eq!(depth_adaptive_level("18.8"), "0.990");
        assert_eq!(depth_adaptive_level("19.0"), "0.980");
        assert_eq!(depth_adaptive_level("NaN"), "0.995");
    }
}
