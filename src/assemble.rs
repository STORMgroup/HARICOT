//! `haricot assemble`: the full assembler, reads to layout, as run for the manuscript.
//!
//! Each step is one pipeline program (`pipeline/`, built next to `haricot` in `target/release`), run as a child
//! process with its arguments fixed here; the programs' defaults are the manuscript configuration. Per step, one row
//! of `OUT/time.tsv` (wall s, max RSS KB, user s, sys s, exit) and a log in `OUT/work/`.
//!
//! Steps: candidate detection (`haricot-detect --production`) -> stage-1 features (`haricot-pair-features`) and
//! forest (`haricot-forest`) -> stage-2 graph features (`haricot-graph-features`) and forest -> keep cut, graph,
//! layout, end bridges (`haricot-layout layout`) -> contig-end rescue (`haricot-rescue`, then features and the
//! rescue forest on its pool) -> the rescue's admission and the final build (`haricot-layout final`).

use std::ffi::{OsStr, OsString};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

use anyhow::{Context, Result, bail, ensure};
use sha2::{Digest, Sha256};

/// The forests of the manuscript runs, as published: file name -> SHA-256.
const PUBLISHED_FORESTS: &[(&str, &str)] = &[
    (
        "stage1.rf0.rfb",
        "2640cfe36e052d2017f10cd6df2bdfb5ef77686a3e889f685f86d7727ac26835",
    ),
    (
        "stage1.rf1.rfb",
        "5c237508c8550a2293d1626cb16b7e9894f9bf1e1867c46ace3be34c981f3f70",
    ),
    (
        "stage1.rf2.rfb",
        "26bb9f26d8603903da5d17440ea0621d63459ddafdbd6a21ed40020c8620b994",
    ),
    (
        "stage1.rf3.rfb",
        "5c8a3e8a98c22a8933e9a26f04d0b655d54189e9a8e30823cdadef894b602fd2",
    ),
    (
        "stage1.rf4.rfb",
        "c9f00c06a6af671b06ca628333efc69c2d21dcf5d9f94726e730ae3d86b83761",
    ),
    (
        "stage1.rfmeta",
        "f9570e80b7cd7cfde4d6cb772db9a9958b14e40bd5d0df610cb7371a7e8fdd62",
    ),
    (
        "stage1.thresholds.json",
        "88bcadc094a4d57d550563c27f2f4ae95590f69dfa8f7470f2d6573f3019f72e",
    ),
    (
        "stage2.rf0.rfb",
        "9ba8002beedd65a250aece4a91390dc8e1b8e1e81a175fda59b3de1ed5f9d93d",
    ),
    (
        "stage2.rf1.rfb",
        "1e7e3afb0e5e8463a7a02822dd71cb356198cb73e0cf2dae9ea27f5b4798ca52",
    ),
    (
        "stage2.rf2.rfb",
        "a8fbccd84af2f8b7f1d778e8c63a58f3cc5b69e4df166e3f2d0b068776e78bbd",
    ),
    (
        "stage2.rf3.rfb",
        "ccb8374b59f1b31e55a427966d4c3f612a5bb999616e6bc8f976c6f1551bae9f",
    ),
    (
        "stage2.rf4.rfb",
        "4aab96c9558efff9d2d9a876c61ac21cfd3366474710add7bc00ea2fdd622959",
    ),
    (
        "stage2.rfmeta",
        "405a4d8bb837e92fd08e9f1dc627d8ab48761e6bdcd8648e76783e87740baa09",
    ),
    (
        "stage2.thresholds.json",
        "8a3848e44432b5bc7a42f8ad3744c6ca2bd88d909eaa0306b3a94d70cee7990b",
    ),
    (
        "rescue.rf0.rfb",
        "f6237b400be4a7c91907bc9cca61cef4f19c3b045696d060efe2e44f006b3f73",
    ),
    (
        "rescue.rf1.rfb",
        "c72cb521808e85f6444cfbcc090eef2a6d9e14648eddc64ddf967a10ec0dcce0",
    ),
    (
        "rescue.rf2.rfb",
        "b06b5623321af720a61df7ccb89e0f39015b50c7c30d7ae7b7e2f3f4813b6a9d",
    ),
    (
        "rescue.rf3.rfb",
        "e5816d6927338c9571f7aadadd34deaf762466598f2881cce3d9118d7228240a",
    ),
    (
        "rescue.rf4.rfb",
        "5057b2e8ed7063dea5106278dbedff2570aee41ac73a54a8a44b0919ac14f2fb",
    ),
    (
        "rescue.rfmeta",
        "daed11d4c65877bcda610d18336c4ecb621bc1636d904cc684b31aede27aff2f",
    ),
    (
        "rescue.thresholds.json",
        "9ec335d41749f33e311b8ceeafdc3dbd5492f849421becb1ef658652c49b2a70",
    ),
];

/// The pipeline programs `assemble` runs.
const PROGRAMS: [&str; 6] = [
    "haricot-detect",
    "haricot-pair-features",
    "haricot-graph-features",
    "haricot-forest",
    "haricot-layout",
    "haricot-rescue",
];

/// The stage-2 keep cut: `policy7` (the depth of the stage-2 keep-0.995 set picks the level; `haricot-layout
/// --keep-level`).
const KEEP_RULE: &str = "policy7";

/// The models directory a build uses when `--models` is not given: `models/` in the source tree it was built from.
pub const DEFAULT_MODELS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/models");

/// What `haricot assemble` takes from the command line.
#[derive(Debug)]
pub struct Options {
    pub reads: PathBuf,
    pub out: PathBuf,
    pub models: PathBuf,
    pub threads: Option<usize>,
    pub keep_work: bool,
    pub force: bool,
    pub bin_dir: Option<PathBuf>,
}

/// Where a step's stdout or stderr goes.
#[derive(Clone, Copy, PartialEq)]
enum Sink {
    /// a log file in the work directory
    Log(&'static str),
    /// the terminal (this process's stream)
    Inherit,
}

struct Runner {
    bin: PathBuf,
    work: PathBuf,
    time_tsv: PathBuf,
    threads: Option<usize>,
}

impl Runner {
    fn open_log(&self, name: &str) -> Result<File> {
        let p = self.work.join(name);
        File::create(&p).with_context(|| format!("creating {}", p.display()))
    }

    /// Run one step; append its row to time.tsv; stop on failure.
    fn step(&self, name: &str, prog: &str, args: &[OsString], out: Sink, err: Sink) -> Result<()> {
        let mut cmd = Command::new(self.bin.join(prog));
        cmd.args(args).stdin(Stdio::null());
        if let Some(t) = self.threads {
            cmd.env("RAYON_NUM_THREADS", t.to_string());
        }
        let shared = match (out, err) {
            (Sink::Log(a), Sink::Log(b)) if a == b => Some(self.open_log(a)?),
            _ => None,
        };
        match (&shared, out) {
            (Some(f), _) => {
                cmd.stdout(f.try_clone()?);
            }
            (None, Sink::Log(n)) => {
                cmd.stdout(self.open_log(n)?);
            }
            (None, Sink::Inherit) => {
                cmd.stdout(Stdio::inherit());
            }
        }
        match (&shared, err) {
            (Some(f), _) => {
                cmd.stderr(f.try_clone()?);
            }
            (None, Sink::Log(n)) => {
                cmd.stderr(self.open_log(n)?);
            }
            (None, Sink::Inherit) => {
                cmd.stderr(Stdio::inherit());
            }
        }
        drop(shared);
        let t0 = Instant::now();
        let child = cmd
            .spawn()
            .with_context(|| format!("starting {}", self.bin.join(prog).display()))?;
        let usage = wait_usage(child.id())?;
        let wall = t0.elapsed().as_secs_f64();
        let mut tl = OpenOptions::new().append(true).open(&self.time_tsv)?;
        writeln!(
            tl,
            "{name}\t{wall:.2}\t{}\t{:.2}\t{:.2}\t{}",
            usage.max_rss_kb, usage.user_s, usage.sys_s, usage.exit
        )?;
        if usage.exit != 0 {
            let log = match (out, err) {
                (Sink::Log(n), _) | (_, Sink::Log(n)) => {
                    format!(": {}", self.work.join(n).display())
                }
                _ => String::new(),
            };
            bail!("step {name} ({prog}) failed with exit {}{log}", usage.exit);
        }
        Ok(())
    }
}

struct Usage {
    max_rss_kb: i64,
    user_s: f64,
    sys_s: f64,
    /// the exit code; 128 + the signal for a killed process
    exit: i32,
}

/// Wait for the child `pid` and return its own resource use (wait4).
fn wait_usage(pid: u32) -> Result<Usage> {
    let mut status: libc::c_int = 0;
    // SAFETY: rusage is plain data; wait4 fills it for the child we spawned.
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    loop {
        // SAFETY: valid pointers to local storage; pid is our own child.
        let r = unsafe { libc::wait4(pid as libc::pid_t, &mut status, 0, &mut ru) };
        if r == pid as libc::pid_t {
            break;
        }
        let e = std::io::Error::last_os_error();
        if e.kind() != std::io::ErrorKind::Interrupted {
            return Err(e).context("waiting for a pipeline step");
        }
    }
    let exit = if libc::WIFEXITED(status) {
        libc::WEXITSTATUS(status)
    } else if libc::WIFSIGNALED(status) {
        128 + libc::WTERMSIG(status)
    } else {
        255
    };
    let secs = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 / 1e6;
    // Linux reports ru_maxrss in KB, macOS in bytes.
    let max_rss_kb = if cfg!(target_os = "macos") {
        ru.ru_maxrss as i64 / 1024
    } else {
        ru.ru_maxrss as i64
    };
    Ok(Usage {
        max_rss_kb,
        user_s: secs(ru.ru_utime),
        sys_s: secs(ru.ru_stime),
        exit,
    })
}

fn sha256_file(p: &Path) -> Result<String> {
    let mut f = File::open(p).with_context(|| format!("reading {}", p.display()))?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(h.finalize().iter().map(|b| format!("{b:02x}")).collect())
}

/// The forests in `dir`: checks, the stage-1 keep-0.995 threshold, the SHA-256 of every file read.
struct Forests {
    t1: String,
    files: Vec<(String, String)>,
    published: bool,
}

/// The seeds an `.rfmeta` lists (`seeds\t0\t1...`): the forests `haricot-forest` reads.
fn rfmeta_seeds(meta: &str) -> Vec<String> {
    meta.lines()
        .find_map(|l| l.strip_prefix("seeds\t"))
        .map(|s| s.split('\t').map(str::to_string).collect())
        .unwrap_or_default()
}

/// The stage-1 keep-0.995 cut as the `.rfmeta` writes it (`threshold\t0.995\t<value>`), compared as text.
fn rfmeta_t1(meta: &str) -> Option<String> {
    meta.lines()
        .find_map(|l| l.strip_prefix("threshold\t0.995\t"))
        .map(str::to_string)
}

/// The keep-0.995 value in a `.thresholds.json` (`"0.995": <number>`), as text.
fn json_t1(json: &str) -> Option<String> {
    let i = json.find("\"0.995\":")? + "\"0.995\":".len();
    let v: String = json[i..]
        .trim_start()
        .chars()
        .take_while(|c| c.is_ascii_digit() || "+-.eE".contains(*c))
        .collect();
    (!v.is_empty()).then_some(v)
}

impl Forests {
    fn check(dir: &Path) -> Result<Forests> {
        let hint = "unpack haricot-forests-v1.tar.gz from https://doi.org/10.5281/zenodo.23963126 \
                    in the checkout (README, Forests) or pass --models DIR";
        ensure!(
            dir.is_dir(),
            "no forests directory {}: {hint}",
            dir.display()
        );
        let mut files = Vec::new();
        let mut t1 = None;
        for prefix in ["stage1", "stage2", "rescue"] {
            let mp = dir.join(format!("{prefix}.rfmeta"));
            ensure!(mp.is_file(), "no {}: {hint}", mp.display());
            let meta = fs::read_to_string(&mp)?;
            let seeds = rfmeta_seeds(&meta);
            ensure!(!seeds.is_empty(), "{}: no seeds line", mp.display());
            let mut names = vec![format!("{prefix}.rfmeta")];
            names.extend(seeds.iter().map(|s| format!("{prefix}.rf{s}.rfb")));
            let jp = dir.join(format!("{prefix}.thresholds.json"));
            if prefix == "stage1" {
                let v = rfmeta_t1(&meta)
                    .with_context(|| format!("{}: no keep-0.995 threshold", mp.display()))?;
                if jp.is_file() {
                    let j = json_t1(&fs::read_to_string(&jp)?);
                    ensure!(
                        j.as_deref() == Some(v.as_str()),
                        "the keep-0.995 threshold differs: {} {:?}, {} {v:?}",
                        jp.display(),
                        j,
                        mp.display()
                    );
                }
                t1 = Some(v);
            }
            if jp.is_file() {
                names.push(format!("{prefix}.thresholds.json"));
            }
            for n in names {
                let p = dir.join(&n);
                ensure!(p.is_file(), "no {}: {hint}", p.display());
                files.push((n, sha256_file(&p)?));
            }
        }
        let published = files.len() == PUBLISHED_FORESTS.len()
            && files
                .iter()
                .all(|(n, h)| PUBLISHED_FORESTS.iter().any(|(pn, ph)| pn == n && ph == h));
        Ok(Forests {
            t1: t1.expect("stage1 checked"),
            files,
            published,
        })
    }
}

/// Output of `git` in the source tree this binary was built from (provenance line only).
fn git(args: &[&str]) -> Option<String> {
    let o = Command::new("git")
        .arg("-C")
        .arg(env!("CARGO_MANIFEST_DIR"))
        .args(args)
        .stderr(Stdio::null())
        .output()
        .ok()?;
    o.status
        .success()
        .then(|| String::from_utf8_lossy(&o.stdout).trim().to_string())
}

fn hostname() -> String {
    let mut buf = [0u8; 256];
    // SAFETY: buf is writable for its length.
    let r = unsafe { libc::gethostname(buf.as_mut_ptr() as *mut libc::c_char, buf.len()) };
    if r != 0 {
        return "unknown".into();
    }
    let n = buf.iter().position(|&b| b == 0).unwrap_or(buf.len());
    String::from_utf8_lossy(&buf[..n]).into_owned()
}

/// Lines of a log containing any of `pats` (the last `max` of them).
fn log_lines(p: &Path, pats: &[&str], max: usize) -> Vec<String> {
    let Ok(f) = File::open(p) else {
        return Vec::new();
    };
    let v: Vec<String> = BufReader::new(f)
        .lines()
        .map_while(Result::ok)
        .filter(|l| pats.iter().any(|p| l.contains(p)))
        .collect();
    v[v.len().saturating_sub(max)..].to_vec()
}

fn os(parts: &[&dyn AsRef<OsStr>]) -> Vec<OsString> {
    parts.iter().map(|p| p.as_ref().to_os_string()).collect()
}

/// An existing OUT is replaced only with --force and only if it holds an earlier run's `time.tsv`.
fn prepare_out(out: &Path, force: bool) -> Result<()> {
    if out.exists() {
        let empty = out.is_dir() && fs::read_dir(out)?.next().is_none();
        if !empty {
            ensure!(
                out.join("time.tsv").is_file(),
                "{} exists and is not an earlier haricot assemble output: not replacing it",
                out.display()
            );
            ensure!(
                force,
                "{} holds an earlier run: pass --force to replace it",
                out.display()
            );
        }
        fs::remove_dir_all(out).with_context(|| format!("removing {}", out.display()))?;
    }
    Ok(())
}

pub fn run(o: &Options) -> Result<()> {
    ensure!(o.reads.is_file(), "no reads file {}", o.reads.display());
    let name = o
        .reads
        .file_name()
        .context("reads path has no file name")?
        .to_string_lossy()
        .into_owned();
    let key = name.strip_suffix(".pod5").unwrap_or(&name).to_string();
    let bin = match &o.bin_dir {
        Some(b) => b.clone(),
        None => std::env::current_exe()?
            .parent()
            .context("executable directory")?
            .to_path_buf(),
    };
    for p in PROGRAMS {
        ensure!(
            bin.join(p).is_file(),
            "no {} next to haricot: build the workspace (cargo build --release) and run target/release/haricot",
            bin.join(p).display()
        );
    }
    if let Some(t) = o.threads {
        ensure!(t >= 1, "--threads must be at least 1");
    }
    let forests = Forests::check(&o.models)?;

    prepare_out(&o.out, o.force)?;
    let w = o.out.join("work");
    for d in ["feat", "m", "s2", "feat9", "s9"] {
        fs::create_dir_all(w.join(d))?;
    }
    let time_tsv = o.out.join("time.tsv");
    File::create(&time_tsv)?;
    let mut sums = File::create(w.join("forests.sha256"))?;
    for (n, h) in &forests.files {
        writeln!(sums, "{h}  {n}")?;
    }

    let mut bsum = Sha256::new();
    for p in PROGRAMS {
        bsum.update(sha256_file(&bin.join(p))?.as_bytes());
    }
    let bsum: String = bsum
        .finalize()
        .iter()
        .take(6)
        .map(|b| format!("{b:02x}"))
        .collect();
    let threads = o
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    println!(
        "haricot assemble {} commit={} dirty={} reads={} key={key} out={} models={} forests={} keep_rule={KEEP_RULE} \
         programs={} programs_sha256={bsum} threads={threads} host={}",
        env!("CARGO_PKG_VERSION"),
        git(&["rev-parse", "--short=7", "HEAD"]).unwrap_or_else(|| "unknown".into()),
        git(&["status", "--porcelain", "--untracked-files=no"])
            .map_or("unknown", |s| if s.is_empty() { "0" } else { "1" }),
        o.reads.display(),
        o.out.display(),
        o.models.display(),
        if forests.published {
            "published"
        } else {
            "CUSTOM (checksums differ from the published forests)"
        },
        bin.display(),
        hostname(),
    );

    let r = Runner {
        bin,
        work: w.clone(),
        time_tsv: time_tsv.clone(),
        threads: o.threads,
    };
    let (reads, m) = (&o.reads, &o.models);
    let ovl = w.join("ovl.bin");
    let occ = w.join("occ.bin");
    let keep: &[&dyn AsRef<OsStr>] = if o.keep_work {
        &[&"--keep-intermediates"]
    } else {
        &[]
    };
    let s1_scores = w.join(format!("m/stage1__{key}.scores.tsv"));
    let both = |n| (Sink::Log(n), Sink::Log(n));

    // 1. candidates: Hough detection, the partner pre-filter, the relative repeat-mask cap (--production)
    let (a, b) = both("detect.log");
    r.step(
        "detect",
        "haricot-detect",
        &os(&[
            reads,
            &"--production",
            &"--outdir",
            &w.join("p2f"),
            &"--save-overlaps",
            &ovl,
            &"--no-graph",
        ]),
        a,
        b,
    )?;
    for l in log_lines(
        &w.join("detect.log"),
        &["repeat-mask cap", "overlaps:", "prefilter:"],
        3,
    ) {
        println!("{l}");
    }
    fs::remove_dir_all(w.join("p2f")).ok();

    // 2. stage-1 features (and the seed-occurrence table the rescue reads), stage-1 forest, stage-2 graph features,
    //    stage-2 forest
    let (a, b) = both("pair_features.log");
    r.step(
        "pair_features",
        "haricot-pair-features",
        &os(&[
            reads,
            &ovl,
            &"--key",
            &key,
            &"--out",
            &w.join(format!("feat/{key}.tsv")),
            &"--seed-occ-rel",
            &"17",
            &"--save-occ",
            &occ,
        ]),
        a,
        b,
    )?;
    r.step(
        "stage1",
        "haricot-forest",
        &os(&[
            &"--stage",
            &"1",
            &"--model-dir",
            m,
            &"--tag",
            &"stage1",
            &"--test",
            &key,
            &"--features-dir",
            &w.join("feat"),
            &"--out-dir",
            &w.join("m"),
        ]),
        Sink::Log("stage1.log"),
        Sink::Inherit,
    )?;
    let (a, b) = both("graph_features.log");
    r.step(
        "graph_features",
        "haricot-graph-features",
        &os(&[
            &ovl,
            &s1_scores,
            &"--t1",
            &forests.t1,
            &"--out",
            &w.join(format!("s2/{key}.s2.tsv")),
        ]),
        a,
        b,
    )?;
    r.step(
        "stage2",
        "haricot-forest",
        &os(&[
            &"--stage",
            &"2",
            &"--model-dir",
            m,
            &"--tag",
            &"stage2",
            &"--test",
            &key,
            &"--features-dir",
            &w.join("feat"),
            &"--s2-dir",
            &w.join("s2"),
            &"--out-dir",
            &w.join("m"),
        ]),
        Sink::Log("stage2.log"),
        Sink::Inherit,
    )?;

    // 3. the keep cut, its build, the end bridges and their build
    let mut args = os(&[
        &"layout",
        &ovl,
        &"--work",
        &w,
        &"--key",
        &key,
        &"--s2-drop",
        &w.join(format!("m/stage2__{key}__keep")),
        &"--scores",
        &s1_scores,
        &"--keep-level",
        &KEEP_RULE,
    ]);
    args.extend(os(keep));
    r.step(
        "layout",
        "haricot-layout",
        &args,
        Sink::Inherit,
        Sink::Log("layout.log"),
    )?;

    // 4. contig-end rescue: end-read pairs the detection never proposed; the rescue forest over its pool
    let (a, b) = both("rescue.log");
    r.step(
        "rescue",
        "haricot-rescue",
        &os(&[
            &"rescue",
            reads,
            &w.join("eb"),
            &w.join("thr_eb.tsv"),
            &ovl,
            &"--out-dir",
            &w.join("q9c"),
            &"--dovetails-only",
            &"--occ",
            &occ,
        ]),
        a,
        b,
    )?;
    for l in log_lines(&w.join("rescue.log"), &["contig ends", "pool:"], usize::MAX) {
        println!("{l}");
    }
    File::create(w.join("acc.tsv"))?;
    if fs::metadata(w.join("q9c/pool.tsv"))
        .map(|m| m.len() > 0)
        .unwrap_or(false)
    {
        let (a, b) = both("rescue_features.log");
        r.step(
            "rescue_features",
            "haricot-pair-features",
            &os(&[
                reads,
                &w.join("q9c/cand_pool.bin"),
                &"--key",
                &key,
                &"--only",
                &w.join("q9c/pool.tsv"),
                &"--occ",
                &occ,
                &"--out",
                &w.join(format!("feat9/{key}.tsv")),
            ]),
            a,
            b,
        )?;
        r.step(
            "rescue_forest",
            "haricot-forest",
            &os(&[
                &"--stage",
                &"1",
                &"--model-dir",
                m,
                &"--tag",
                &"rescue",
                &"--test",
                &key,
                &"--features-dir",
                &w.join("feat9"),
                &"--out-dir",
                &w.join("s9"),
                &"--accept",
                &w.join("acc.tsv"),
                &"--read-names",
                &w.join("eb/read_names.tsv"),
            ]),
            Sink::Log("rescue_forest.log"),
            Sink::Inherit,
        )?;
    }

    // 5. the rescue's admission and the final build into OUT
    let mut args = os(&[&"final", &ovl, &"--work", &w, &"--out", &o.out]);
    args.extend(os(keep));
    r.step(
        "final",
        "haricot-layout",
        &args,
        Sink::Inherit,
        Sink::Log("final.log"),
    )?;

    if !o.keep_work {
        for e in fs::read_dir(&w)? {
            let p = e?.path();
            let n = p
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();
            if p.is_file() && (n.ends_with(".bin") || n.ends_with(".dump.tsv")) {
                fs::remove_file(&p)?;
            }
        }
        if let Ok(rd) = fs::read_dir(w.join("q9c")) {
            for e in rd {
                let p = e?.path();
                if p.extension().is_some_and(|x| x == "bin") {
                    fs::remove_file(&p)?;
                }
            }
        }
        for d in ["feat", "feat9", "s2", "pol7", "eb"] {
            let p = w.join(d);
            if p.exists() {
                fs::remove_dir_all(&p)?;
            }
        }
    }

    let (mut wall, mut cpu, mut rss, mut n) = (0.0, 0.0, 0i64, 0);
    for l in fs::read_to_string(&time_tsv)?.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() == 6 {
            wall += f[1].parse::<f64>().unwrap_or(0.0);
            cpu += f[3].parse::<f64>().unwrap_or(0.0) + f[4].parse::<f64>().unwrap_or(0.0);
            rss = rss.max(f[2].parse().unwrap_or(0));
            n += 1;
        }
    }
    println!(
        "TOTAL wall {wall:.1} s  cpu {cpu:.1} s  peak RSS {:.2} GB  steps {n}",
        rss as f64 / 1048576.0
    );
    println!("output: {}", o.out.join("unitigs.gfa").display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thresholds_read_as_text() {
        let meta = "features\tz\tslope\nseeds\t0\t1\t2\nthreshold\t0.995\t0.43939824982790304\nthreshold\t0.99\t0.5\n";
        assert_eq!(rfmeta_seeds(meta), vec!["0", "1", "2"]);
        assert_eq!(rfmeta_t1(meta).as_deref(), Some("0.43939824982790304"));
        let json =
            "{\n \"thresholds\": {\n  \"0.995\": 0.43939824982790304,\n  \"0.99\": 0.5\n }\n}";
        assert_eq!(json_t1(json).as_deref(), Some("0.43939824982790304"));
        assert_eq!(json_t1("{}"), None);
    }

    #[test]
    fn published_table_is_complete() {
        assert_eq!(PUBLISHED_FORESTS.len(), 21);
        for p in ["stage1", "stage2", "rescue"] {
            assert_eq!(
                PUBLISHED_FORESTS
                    .iter()
                    .filter(|(n, _)| n.starts_with(p))
                    .count(),
                7
            );
        }
        assert!(PUBLISHED_FORESTS.iter().all(|(_, h)| h.len() == 64));
    }

    #[test]
    fn out_dir_is_never_replaced_unless_it_is_a_run() {
        let d = std::env::temp_dir().join(format!("haricot_assemble_test_{}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        fs::write(d.join("x.txt"), "x").unwrap();
        assert!(prepare_out(&d, true).is_err());
        assert!(d.join("x.txt").exists());
        fs::write(d.join("time.tsv"), "").unwrap();
        assert!(prepare_out(&d, false).is_err());
        assert!(prepare_out(&d, true).is_ok());
        assert!(!d.exists());
    }
}
