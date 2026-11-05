//! Score one dataset's candidates with saved forests (`M/TAG.rf<seed>.rfb`, `M/TAG.rfmeta`, written by
//! `training/export_forest.py`), and write the score table and drop lists.
//!
//! Loading a joblib forest takes 3-8 s (zlib + unpickling); the compact `.rfb` export loads in a
//! fraction of that. The probabilities are sklearn's:
//!
//!   * the features are the TSV columns named in `TAG.rfmeta`, parsed as `f64` and cast to `f32`
//!     (as training reads them: pandas' `read_csv` then `to_numpy(float32)`), then `NaN -> -1`,
//!     `+inf -> 1e6`, `-inf -> -1e6` (`np.nan_to_num`); a forest column in no input table is an
//!     error (it would read -1 in every row), a column `NaN` in every row a warning on stderr; rows
//!     are parsed in parallel;
//!   * a tree's probability is its leaf's class-1 fraction (sklearn >= 1.4 stores fractions),
//!     reached by `x <= t` with `t` the largest `f32` at or under sklearn's `f64` threshold: for an
//!     `f32` input the same test as sklearn's;
//!   * a forest's probability is its trees' mean, summed in tree order; the ensemble's, the
//!     forests' mean, summed in seed order (as `training/train_stage1.py`'s `fit_forests`). sklearn
//!     sums a forest's trees in thread-completion order, so its own scores differ run to run by
//!     ~1 ulp; these are one such order, deterministic.
//!
//! Outputs: `OUT/TAG__K.scores.tsv` (`name_a name_b cls score`, the score in its shortest round-trip
//! digits, as Python prints a float) and `OUT/TAG__K__keep<k>.drop.tsv` per frozen threshold
//! (`score >= t` kept; the drop list holds every candidate not kept, stage-1 rejects included at
//! stage 2).
//!
//! `--accept A.tsv --read-names R.tsv [--accept-th 0.5]` also writes the scored pairs at
//! `score >= th` as read indices (`R.tsv`: `index<TAB>name`), lower first, sorted: the contig-end
//! rescue's accepted pool pairs.
//!
//! Usage: haricot-forest --stage 1 --model-dir M --tag all5 --test K --features-dir FD --out-dir OUT
//!        haricot-forest --stage 2 --model-dir M --tag s2all5 --test K --features-dir FD --s2-dir SD
//!                 --out-dir OUT
//! An unknown flag, a flag without its value and a value that does not parse are errors
//! (`cliargs.rs`).
use std::io::Write;

use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

#[path = "../cliargs.rs"]
mod cliargs;

/// The flags (each takes a value; `cliargs.rs`).
const VALUED: &[&str] = &[
    "--stage",
    "--model-dir",
    "--tag",
    "--test",
    "--features-dir",
    "--out-dir",
    "--s2-dir",
    "--accept",
    "--read-names",
    "--accept-th",
];
const USAGE: &str = "usage: haricot-forest --stage 1|2 --model-dir M --tag T --test K --features-dir FD \
                     [--s2-dir SD] --out-dir OUT [--accept A.tsv --read-names R.tsv [--accept-th 0.5]]";

/// The forest features that none of the input tables (`heads`: each table's columns) has. Such a
/// feature would reach the forest as -1 in every row (`np.nan_to_num` of an absent column): an
/// error, not a default. A feature in two tables counts as present (pandas' merge suffixes it, as
/// in training).
fn absent_features<'a>(feats: &'a [String], heads: &[&FxHashMap<String, usize>]) -> Vec<&'a str> {
    feats
        .iter()
        .filter(|f| heads.iter().all(|h| !h.contains_key(f.as_str())))
        .map(|f| f.as_str())
        .collect()
}

/// One tree node, as `training/export_forest.py` writes it.
#[derive(Clone, Copy)]
struct Node {
    /// Feature index, or `LEAF`.
    feature: u16,
    threshold: f32,
    /// Right child (internal) or value index (leaf); the left child is the next node.
    x: u32,
}

const LEAF: u16 = 0xFFFF;

/// A scored row's `(name_a, name_b, cls)`.
type RowName<'a> = (&'a str, &'a str, &'a str);

struct Forest {
    n_features: usize,
    /// `(first node, node count)` per tree, into `nodes`.
    trees: Vec<(usize, usize)>,
    nodes: Vec<Node>,
    values: Vec<f64>,
}

fn u32_at(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

impl Forest {
    fn read(path: &str) -> Result<Forest> {
        let b = std::fs::read(path).with_context(|| format!("reading {path}"))?;
        ensure!(
            b.len() >= 16 && &b[..4] == b"RFB1",
            "{path}: not an RFB1 forest"
        );
        let (n_features, n_trees, n_values) = (
            u32_at(&b, 4) as usize,
            u32_at(&b, 8) as usize,
            u32_at(&b, 12) as usize,
        );
        let mut at = 16;
        let values: Vec<f64> = b[at..at + 8 * n_values]
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        at += 8 * n_values;
        let (mut trees, mut nodes) = (Vec::with_capacity(n_trees), Vec::new());
        for _ in 0..n_trees {
            let n = u32_at(&b, at) as usize;
            at += 4;
            trees.push((nodes.len(), n));
            nodes.extend(b[at..at + 12 * n].chunks_exact(12).map(|c| Node {
                feature: u16::from_le_bytes([c[0], c[1]]),
                threshold: f32::from_le_bytes(c[4..8].try_into().unwrap()),
                x: u32::from_le_bytes(c[8..12].try_into().unwrap()),
            }));
            at += 12 * n;
        }
        ensure!(at == b.len(), "{path}: {} trailing bytes", b.len() - at);
        Ok(Forest {
            n_features,
            trees,
            nodes,
            values,
        })
    }

    /// Add each tree's probability for every row of `x` (`n_features` values per row) to
    /// `sum`, tree by tree.
    fn accumulate(&self, x: &[f32], sum: &mut [f64]) {
        let f = self.n_features;
        for &(start, len) in &self.trees {
            let nodes = &self.nodes[start..start + len];
            for (row, s) in x.chunks_exact(f).zip(sum.iter_mut()) {
                let mut k = 0usize;
                loop {
                    let nd = nodes[k];
                    if nd.feature == LEAF {
                        *s += self.values[nd.x as usize];
                        break;
                    }
                    k = if row[nd.feature as usize] <= nd.threshold {
                        k + 1
                    } else {
                        nd.x as usize
                    };
                }
            }
        }
    }
}

/// A TSV as rows of fields, with its header.
struct Table {
    cols: FxHashMap<String, usize>,
    rows: Vec<Vec<String>>,
}

impl Table {
    /// The rows are split in parallel (rayon's `par_lines`: `str::lines`' lines, in order).
    fn read(path: &str) -> Result<Table> {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
        let (head, body) = text.split_once('\n').unwrap_or((text.as_str(), ""));
        let head = head.strip_suffix('\r').unwrap_or(head);
        ensure!(!text.is_empty(), "{path} is empty");
        let cols = head
            .split('\t')
            .enumerate()
            .map(|(i, c)| (c.to_string(), i))
            .collect();
        let rows = body
            .par_lines()
            .map(|l| l.split('\t').map(str::to_string).collect())
            .collect();
        Ok(Table { cols, rows })
    }

    fn col(&self, name: &str) -> Result<usize> {
        self.cols
            .get(name)
            .copied()
            .with_context(|| format!("no column {name}"))
    }
}

/// A feature field as pandas reads it: a float (`NaN`, `inf` spelled as Rust and pandas both read
/// them), cast to `f32`; an absent column `NaN`. [`nan_to_num`] then makes it the forest's input.
fn raw_feature(field: Option<&str>) -> Result<f32> {
    Ok(match field {
        None => f32::NAN,
        Some(s) => s
            .parse::<f64>()
            .with_context(|| format!("unreadable feature value {s:?}"))? as f32,
    })
}

/// The forest features that read `NaN` in every one of `rows` rows (`nan[j]` = the rows where
/// feature `j` is `NaN`): a column computed as `NaN` (the signal-DTW columns under `haricot-pair-features
/// --no-dtw`; a column absent from every table is an error first, [`absent_features`]), or a stage-2
/// feature that is in both tables (pandas' merge suffixes it, so it reads `NaN` as in training).
/// [`nan_to_num`] turns them into -1 for every row, which no training row of a forest that uses the
/// column looked like.
fn all_nan_features<'a>(feats: &'a [String], nan: &[usize], rows: usize) -> Vec<&'a str> {
    if rows == 0 {
        return Vec::new();
    }
    feats
        .iter()
        .zip(nan)
        .filter(|&(_, &n)| n == rows)
        .map(|(f, _)| f.as_str())
        .collect()
}

/// `np.nan_to_num(nan=-1, posinf=1e6, neginf=-1e6)`.
fn nan_to_num(v: f32) -> f32 {
    if v.is_nan() {
        -1.0
    } else if v == f32::INFINITY {
        1e6
    } else if v == f32::NEG_INFINITY {
        -1e6
    } else {
        v
    }
}

/// `x` as Python's `repr(float)` (and numpy's `str`) writes it: the shortest digits that read
/// back to `x`, fixed notation for decimal exponents in `[-5, 15]` and scientific otherwise
/// (`1e-05`, `1.5e+16`), a trailing `.0` on integral values.
fn py_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.into();
    }
    let sign = if x.is_sign_negative() { "-" } else { "" };
    if x == 0.0 {
        return format!("{sign}0.0");
    }
    // Rust's `{:e}` is the shortest round-trip form, `d.ddde<exp>`.
    let e = format!("{:e}", x.abs());
    let (mant, exp) = e.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let decpt = exp + 1; // digits before the decimal point
    let nd = digits.len() as i32;
    if -4 < decpt && decpt <= 16 {
        if decpt <= 0 {
            format!("{sign}0.{}{digits}", "0".repeat((-decpt) as usize))
        } else if decpt >= nd {
            format!("{sign}{digits}{}.0", "0".repeat((decpt - nd) as usize))
        } else {
            let (a, b) = digits.split_at(decpt as usize);
            format!("{sign}{a}.{b}")
        }
    } else {
        let (a, b) = digits.split_at(1);
        let frac = if b.is_empty() {
            String::new()
        } else {
            format!(".{b}")
        };
        let es = if exp < 0 { '-' } else { '+' };
        format!("{sign}{a}{frac}e{es}{:02}", exp.abs())
    }
}

fn main() -> Result<()> {
    let t0 = std::time::Instant::now();
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = cliargs::Cli::parse(&argv, VALUED, &[], USAGE)?;
    args.positional(0, USAGE)?;
    let stage: u32 = args.opt("--stage").context("--stage 1|2")?.parse()?;
    let md = args.opt("--model-dir").context("--model-dir")?;
    let tag = args.opt("--tag").context("--tag")?;
    let key = args.opt("--test").context("--test")?;
    let fd = args.opt("--features-dir").context("--features-dir")?;
    let out = args.opt("--out-dir").context("--out-dir")?;
    ensure!(!key.contains(','), "one dataset per call");
    std::fs::create_dir_all(&out)?;

    // features, seeds, thresholds (descending by keep level)
    let meta = std::fs::read_to_string(format!("{md}/{tag}.rfmeta"))
        .with_context(|| format!("reading {md}/{tag}.rfmeta"))?;
    let (mut feats, mut seeds, mut thr) = (Vec::new(), Vec::new(), Vec::new());
    for line in meta.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "features" => feats = f[1..].iter().map(|s| s.to_string()).collect(),
            "seeds" => seeds = f[1..].iter().map(|s| s.to_string()).collect(),
            "threshold" => thr.push((f[1].parse::<f64>()?, f[2].parse::<f64>()?)),
            other => bail!("unknown rfmeta line {other}"),
        }
    }
    thr.sort_by(|a, b| b.0.total_cmp(&a.0));
    let forests: Vec<Forest> = seeds
        .par_iter()
        .map(|sd| Forest::read(&format!("{md}/{tag}.rf{sd}.rfb")))
        .collect::<Result<_>>()?;
    for f in &forests {
        ensure!(
            f.n_features == feats.len(),
            "forest and rfmeta disagree on features"
        );
    }

    // The scored rows `d` (name_a, name_b, cls, features) and the full candidate set.
    let ft = Table::read(&format!("{fd}/{key}.tsv"))?;
    let (fa, fb, fc) = (ft.col("name_a")?, ft.col("name_b")?, ft.col("cls")?);
    let full: Vec<(&str, &str, &str)> = ft
        .rows
        .iter()
        .map(|r| (r[fa].as_str(), r[fb].as_str(), r[fc].as_str()))
        .collect();
    let nf = feats.len();
    ensure!(stage == 1 || stage == 2, "--stage 1|2");
    let st = if stage == 2 {
        let sd = args.opt("--s2-dir").context("--s2-dir")?;
        Some(Table::read(&format!("{sd}/{key}.s2.tsv"))?)
    } else {
        None
    };
    // Every column the forest reads must exist in an input table (stage 2: the stage-2 table or the
    // stage-1 features); a renamed or never-written column would otherwise read -1 in every row.
    let mut heads = vec![&ft.cols];
    heads.extend(st.as_ref().map(|t| &t.cols));
    let absent = absent_features(&feats, &heads);
    ensure!(
        absent.is_empty(),
        "{md}/{tag}.rfmeta: {} forest feature(s) in no input table: {}",
        absent.len(),
        absent.join(",")
    );
    // The rows are parsed in parallel into their raw feature values (`NaN` kept), row-major in row
    // order; `no_s1`: the stage-2 rows without a stage-1 row. The `NaN` count per feature and
    // `np.nan_to_num` follow, over the whole matrix.
    let (names, mut x, no_s1): (Vec<RowName<'_>>, Vec<f32>, usize) = if stage == 1 {
        let idx: Vec<Option<usize>> = feats.iter().map(|c| ft.cols.get(c).copied()).collect();
        let x: Vec<Vec<f32>> = ft
            .rows
            .par_iter()
            .map(|r| {
                idx.iter()
                    .map(|i| raw_feature(i.map(|i| r[i].as_str())))
                    .collect::<Result<Vec<f32>>>()
            })
            .collect::<Result<_>>()?;
        (full.clone(), x.concat(), 0)
    } else {
        let st = st.as_ref().unwrap();
        let (sa, sb) = (st.col("name_a")?, st.col("name_b")?);
        // pandas' left merge on the names: the stage-2 rows in order, each with its stage-1
        // row's columns. A column in both tables would come out suffixed (`_x`, `_y`), so a
        // feature by that name reads as absent.
        let at: FxHashMap<(&str, &str), usize> = ft
            .rows
            .iter()
            .enumerate()
            .map(|(i, r)| ((r[fa].as_str(), r[fb].as_str()), i))
            .collect();
        let rows: Vec<(RowName<'_>, bool, Vec<f32>)> = st
            .rows
            .par_iter()
            .map(|r| {
                let (a, b) = (r[sa].as_str(), r[sb].as_str());
                let f = at.get(&(a, b)).copied();
                let name = (a, b, f.map_or("", |i| ft.rows[i][fc].as_str()));
                let x = feats
                    .iter()
                    .map(|c| {
                        let v = match (st.cols.get(c), ft.cols.get(c)) {
                            (Some(&i), None) => Some(r[i].as_str()),
                            (None, Some(&j)) => f.map(|i| ft.rows[i][j].as_str()),
                            _ => None,
                        };
                        raw_feature(v)
                    })
                    .collect::<Result<Vec<f32>>>()?;
                Ok::<_, anyhow::Error>((name, f.is_none(), x))
            })
            .collect::<Result<_>>()?;
        let no_s1 = rows.iter().filter(|r| r.1).count();
        let mut names = Vec::with_capacity(rows.len());
        let mut x = Vec::with_capacity(rows.len() * nf);
        for (name, _, v) in rows {
            names.push(name);
            x.extend(v);
        }
        (names, x, no_s1)
    };
    // Rows where each feature reads NaN, then np.nan_to_num's replacements (the forest's input).
    let nan: Vec<usize> = x
        .par_chunks(nf.max(1))
        .fold(
            || vec![0usize; nf],
            |mut acc, row| {
                for (a, v) in acc.iter_mut().zip(row) {
                    *a += v.is_nan() as usize;
                }
                acc
            },
        )
        .reduce(
            || vec![0usize; nf],
            |mut a, b| {
                for (p, q) in a.iter_mut().zip(b) {
                    *p += q;
                }
                a
            },
        );
    x.par_iter_mut().for_each(|v| *v = nan_to_num(*v));
    // A forest column that is NaN in every row (computed as NaN: e.g. DTW-trained forests over
    // `haricot-pair-features --no-dtw` tables) reaches the forest as -1 everywhere. Reported, never acted on
    // (the scores are those of the forests).
    let flat = all_nan_features(&feats, &nan, names.len());
    if !flat.is_empty() {
        eprintln!(
            "haricot-forest WARNING [{tag} {key}]: {} of {nf} features NaN (-1) in all {} rows: {}",
            flat.len(),
            names.len(),
            flat.join(",")
        );
    }
    if no_s1 > 0 {
        eprintln!(
            "haricot-forest WARNING [{tag} {key}]: {no_s1} stage-2 rows lack a stage-1 row (features -1)"
        );
    }

    let t_read = t0.elapsed().as_secs_f64();
    // Predict: rows in chunks across threads; within a chunk, forest by forest, tree by tree.
    const CHUNK: usize = 1024;
    let n = names.len();
    let score: Vec<f64> = x
        .par_chunks(CHUNK * nf)
        .flat_map_iter(|xc| {
            let rows = xc.len() / nf;
            let mut ens = vec![0.0f64; rows];
            for f in &forests {
                let mut sum = vec![0.0f64; rows];
                f.accumulate(xc, &mut sum);
                for (e, s) in ens.iter_mut().zip(&sum) {
                    *e += s / f.trees.len() as f64;
                }
            }
            ens.into_iter().map(|e| e / forests.len() as f64)
        })
        .collect();
    assert_eq!(score.len(), n);
    let t_pred = t0.elapsed().as_secs_f64();

    // pandas reads a column holding nothing but true/false spellings as booleans: it writes
    // them back as `True` / `False`, and no value then equals the string "true".
    const TRUE: [&str; 3] = ["true", "True", "TRUE"];
    const FALSE: [&str; 3] = ["false", "False", "FALSE"];
    let boolean = !full.is_empty()
        && full
            .iter()
            .all(|r| TRUE.contains(&r.2) || FALSE.contains(&r.2));
    let cls_out = |c: &str| -> String {
        match (boolean, TRUE.contains(&c), FALSE.contains(&c)) {
            (true, true, _) => "True".into(),
            (true, _, true) => "False".into(),
            _ => c.into(),
        }
    };
    let is = |c: &str, v: &str| !boolean && c == v;
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!(
        "{out}/{tag}__{key}.scores.tsv"
    ))?);
    writeln!(w, "name_a\tname_b\tcls\tscore")?;
    let lines: Vec<String> = names
        .par_iter()
        .zip(&score)
        .map(|((a, b, c), s)| format!("{a}\t{b}\t{}\t{}\n", cls_out(c), py_float(*s)))
        .collect();
    for l in &lines {
        w.write_all(l.as_bytes())?;
    }
    w.flush()?;
    drop(lines);
    let n_true = full.iter().filter(|r| is(r.2, "true")).count();
    // labelled rows (training and truth-scored sets); new data has none, and its lines say so
    let labelled = full
        .iter()
        .any(|r| is(r.2, "true") || is(r.2, "false") || is(r.2, "misplaced"));
    println!(
        "[{tag}] APPLY {key}: {} scored of {} candidates{}",
        n,
        full.len(),
        if labelled {
            format!(" ({n_true} true)")
        } else {
            " (no labels)".to_string()
        }
    );
    for &(k, t) in &thr {
        let kept: FxHashSet<(&str, &str)> = names
            .iter()
            .zip(&score)
            .filter(|(_, s)| **s >= t)
            .map(|(r, _)| (r.0, r.1))
            .collect();
        let (mut kt, mut kn) = (0usize, 0usize);
        for (r, s) in names.iter().zip(&score) {
            if *s >= t {
                kt += is(r.2, "true") as usize;
                kn += (is(r.2, "false") || is(r.2, "misplaced")) as usize;
            }
        }
        if labelled {
            println!(
                "   keep {k:.3}: true kept {kt}/{n_true} ({:.3})  negatives kept {kn}  precision {:.4}",
                kt as f64 / n_true.max(1) as f64,
                kt as f64 / (kt + kn).max(1) as f64
            );
        } else {
            println!("   keep {k:.3}: {} kept of {n}", kept.len());
        }
        let mut w = std::io::BufWriter::new(std::fs::File::create(format!(
            "{out}/{tag}__{key}__keep{k:.3}.drop.tsv"
        ))?);
        for r in full.iter().filter(|r| !kept.contains(&(r.0, r.1))) {
            writeln!(w, "{}\t{}", r.0, r.1)?;
        }
        w.flush()?;
    }
    if let Some(acc) = args.opt("--accept") {
        let rn = args
            .opt("--read-names")
            .context("--accept needs --read-names")?;
        let th: f64 = args.num("--accept-th", 0.5)?;
        let text = std::fs::read_to_string(&rn).with_context(|| format!("reading {rn}"))?;
        let mut rid: FxHashMap<&str, u64> = FxHashMap::default();
        for line in text.lines() {
            let (i, n) = line
                .split_once('\t')
                .with_context(|| format!("{rn}: {line:?}"))?;
            rid.insert(n, i.parse()?);
        }
        let id = |n: &str| {
            rid.get(n)
                .copied()
                .with_context(|| format!("{n} not in {rn}"))
        };
        let mut kept: Vec<(u64, u64)> = Vec::new();
        for ((a, b, _), s) in names.iter().zip(&score) {
            if *s >= th {
                let (x, y) = (id(a)?, id(b)?);
                kept.push((x.min(y), x.max(y)));
            }
        }
        kept.sort_unstable();
        let mut w = std::io::BufWriter::new(std::fs::File::create(&acc)?);
        for (i, j) in &kept {
            writeln!(w, "{i}\t{j}")?;
        }
        w.flush()?;
        println!("pool {n} scored; {} at score >= {th}", kept.len());
    }
    eprintln!(
        "[{tag}] PHASE read + parse {t_read:.2} s, predict {:.2} s ({n} rows x {} forests x {} trees), write {:.2} s",
        t_pred - t_read,
        forests.len(),
        forests.first().map_or(0, |f| f.trees.len()),
        t0.elapsed().as_secs_f64() - t_pred
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The calls in the tree parse: `haricot assemble`'s stage 1, stage 2 and rescue forest.
    #[test]
    fn the_calls_in_the_tree_parse() {
        let cli = |l: &str| {
            let v: Vec<String> = l.split_whitespace().map(String::from).collect();
            cliargs::Cli::parse(&v, VALUED, &[], USAGE)
        };
        for l in [
            "--stage 1 --model-dir M --tag t --test k --features-dir F --out-dir O",
            "--stage 2 --model-dir M --tag s2t --test k --features-dir F --s2-dir S --out-dir O",
            "--stage 1 --model-dir M --tag t --test k --features-dir F --out-dir O --accept A \
             --read-names R",
            "--stage 1 --model-dir M --tag t --test k --features-dir F --out-dir O --accept A \
             --read-names R --accept-th 0.5",
        ] {
            let a = cli(l).unwrap();
            a.positional(0, "").unwrap();
            assert_eq!(a.num("--accept-th", 0.5).unwrap(), 0.5);
        }
        assert!(cli("--stage 1 --model-dir M --tags t").is_err());
        assert!(cli("--stage 1 x").unwrap().positional(0, "").is_err());
    }

    /// The forest's input is `np.nan_to_num`'s, and a column NaN in every row is reported (absent
    /// columns and `--no-dtw` columns read NaN): only those, and nothing over an empty table.
    #[test]
    fn nan_features_are_minus_one_and_all_nan_columns_are_named() {
        let v: Vec<f32> = ["nan", "inf", "-inf", "0.25"]
            .iter()
            .map(|s| nan_to_num(raw_feature(Some(s)).unwrap()))
            .collect();
        assert_eq!(v, vec![-1.0, 1e6, -1e6, 0.25]);
        assert!(raw_feature(None).unwrap().is_nan());
        assert!(raw_feature(Some("x")).is_err());
        let feats: Vec<String> = ["z", "dtw_mean", "tri_frac"].map(String::from).to_vec();
        assert_eq!(all_nan_features(&feats, &[0, 7, 3], 7), vec!["dtw_mean"]);
        assert!(all_nan_features(&feats, &[0, 0, 0], 0).is_empty());
    }

    /// Python's `repr` of these floats, as printed by `python3 -c 'print(repr(x))'`.
    #[test]
    fn floats_print_as_python_prints_them() {
        for (x, want) in [
            (0.0, "0.0"),
            (1.0, "1.0"),
            (0.5, "0.5"),
            (0.9800781849860145, "0.9800781849860145"),
            (0.1, "0.1"),
            (1e-05, "1e-05"),
            (1.25e-05, "1.25e-05"),
            (0.0001, "0.0001"),
            (0.00012345, "0.00012345"),
            (123456789012345.0, "123456789012345.0"),
            (1e16, "1e+16"),
            (1.5e16, "1.5e+16"),
            (1234567890123456.0, "1234567890123456.0"),
            (2.0 / 3.0, "0.6666666666666666"),
            (-0.25, "-0.25"),
        ] {
            assert_eq!(py_float(x), want, "{x}");
        }
    }

    /// A forest feature no input table has is reported (stage 2: neither table); one in either
    /// table is not.
    #[test]
    fn absent_forest_features_are_named() {
        let h = |cols: &[&str]| -> FxHashMap<String, usize> {
            cols.iter()
                .enumerate()
                .map(|(i, c)| (c.to_string(), i))
                .collect()
        };
        let feats: Vec<String> = ["z", "s1", "dtw_ratio"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let (s1t, s2t) = (
            h(&["name_a", "name_b", "cls", "z", "dtw_ratio"]),
            h(&["name_a", "name_b", "s1"]),
        );
        assert_eq!(absent_features(&feats, &[&s1t]), vec!["s1"]);
        assert!(absent_features(&feats, &[&s1t, &s2t]).is_empty());
        assert_eq!(
            absent_features(&feats, &[&h(&["name_a", "z"])]),
            vec!["s1", "dtw_ratio"]
        );
    }
}
