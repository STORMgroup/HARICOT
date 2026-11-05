#!/usr/bin/env python3
"""Training is deterministic and invariant to read names and row order: train one small fold three
times, twice from the labelled tables as written (determinism: every file the two runs write byte-identical) and once
from a copy with every read renamed (a seeded hash) and every table's rows shuffled (seeded), and compare.

Each run: train_stage1.py (stage 1) -> haricot-graph-features on the caches with its out-of-bag scores ->
train_stage2.py (stage 2) -> export_forest.py, with --cap-per-set and --exclude-reads (the list renamed in the second run). The second run's
stage-2 tables are the first run's, renamed and shuffled: that is what haricot-graph-features writes for renamed caches when the
stage-1 out-of-bag scores agree, which is checked (the scores per pair, names mapped back).
PASS = (determinism) the two as-written runs' score, drop and stage-2 tables, forests, rfmeta and thresholds files are
byte-identical; (invariance) every score file of the renamed run (stage-1 out-of-bag and test scores, stage-2 test
scores) agrees pair for pair, as text, with the as-written run's (names mapped back), and for the stage-1 and stage-2
tags every TAG.rf<seed>.rfb and TAG.rfmeta is byte-identical and TAG.thresholds.json agrees but for the list's hash;
and train_stage1.forest_proba equals sklearn's single-threaded predict_proba exactly on the test table. Exit 0 = PASS.

--det is train.sh's WORK/det (KEY.ovl.bin per set).
Usage: train_invariance.py --feat DIR --det DIR --bin DIR --exclude-reads FILE --train K1,K2 --test K --out DIR
                           [--seeds 0,1] [--trees 60] [--cap 3000]
"""
import argparse
import filecmp
import hashlib
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SC = HERE


def run(cmd, log):
    with open(log, "a") as fh:
        fh.write("$ " + " ".join(cmd) + "\n")
        fh.flush()
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"FAILED ({r.returncode}): {' '.join(cmd)} (see {log})")


def new_name(x):
    return "rn" + hashlib.sha1(f"inv-7:{x}".encode()).hexdigest()[:20]


def rename_table(src, dst, seed):
    """Rename name_a / name_b of a TSV (text level: every other byte unchanged) and shuffle its rows (seeded)."""
    with open(src) as fh:
        head = fh.readline()
        rows = fh.read().splitlines()
    cols = head.rstrip("\n").split("\t")
    ia, ib = cols.index("name_a"), cols.index("name_b")
    out = []
    for ln in rows:
        f = ln.split("\t")
        f[ia], f[ib] = new_name(f[ia]), new_name(f[ib])
        out.append("\t".join(f))
    random.Random(seed).shuffle(out)
    with open(dst, "w") as fh:
        fh.write(head)
        fh.write("".join(x + "\n" for x in out))


def one(label, feat, s2src, excl, a, out):
    md = os.path.join(out, label)
    os.makedirs(md, exist_ok=True)
    log = os.path.join(out, f"{label}.log")
    common = ["--seeds", a.seeds, "--save-models", "--cap-per-set", str(a.cap), "--trees", str(a.trees),
              "--exclude-reads", excl]
    run([sys.executable, f"{SC}/train_stage1.py", "--train", a.train, "--test", a.test, "--tag", "inv",
         "--features-dir", feat, "--out-dir", md] + common, log)
    t1 = json.load(open(f"{md}/inv.thresholds.json"))["thresholds"]["0.995"]
    s2 = os.path.join(out, f"{label}_s2")
    os.makedirs(s2, exist_ok=True)
    if s2src is None:  # the caches: haricot-graph-features on this run's stage-1 scores
        for k in a.train.split(","):
            run([f"{a.bin}/haricot-graph-features", f"{a.det}/{k}.ovl.bin", f"{md}/inv__{k}.oob.tsv", "--t1", repr(t1),
                 "--out", f"{s2}/{k}.s2.tsv", "--symmetric"], log)
        run([f"{a.bin}/haricot-graph-features", f"{a.det}/{a.test}.ovl.bin", f"{md}/inv__{a.test}.scores.tsv", "--t1",
             repr(t1), "--out", f"{s2}/{a.test}.s2.tsv", "--symmetric"], log)
    else:
        for i, k in enumerate(a.train.split(",") + [a.test]):
            rename_table(f"{s2src}/{k}.s2.tsv", f"{s2}/{k}.s2.tsv", 100 + i)
    run([sys.executable, f"{SC}/train_stage2.py", "--train", a.train, "--test", a.test, "--tag", "s2inv",
         "--s2-dir", s2, "--features-dir", feat, "--out-dir", md] + common, log)
    for tag in ("inv", "s2inv"):
        run([sys.executable, f"{SC}/export_forest.py", "--model-dir", md, "--tag", tag], log)
    return md, s2


def score_map(path, back=None):
    """{(name_a, name_b): (cls, score text)} of a score file (name_a name_b cls score), names mapped through `back`."""
    out = {}
    with open(path) as fh:
        fh.readline()
        for ln in fh:
            a, b, c, s = ln.rstrip("\n").split("\t")
            if back:
                a, b = back[a], back[b]
            out[(a, b)] = (c, s)
    return out


def same_tree(d1, d2):
    """Every file of d1 (but joblib forests and logs) byte-identical in d2, and no file missing: (n, differing)."""
    n, bad = 0, []
    for f in sorted(os.listdir(d1)):
        if f.endswith((".joblib", ".log")) or not os.path.isfile(os.path.join(d1, f)):
            continue
        n += 1
        g = os.path.join(d2, f)
        if not (os.path.exists(g) and filecmp.cmp(os.path.join(d1, f), g, shallow=False)):
            bad.append(f)
    return n, bad


def sklearn_check(md, feat, test):
    """forest_proba against sklearn's own predict_proba with n_jobs=1 (the trees summed in order), on the test table."""
    import joblib
    import numpy as np
    sys.path.insert(0, SC)
    import train_stage1 as T
    j = json.load(open(f"{md}/inv.thresholds.json"))
    d = T.load(feat, test)
    X = np.nan_to_num(d[j["features"]].to_numpy(dtype=np.float32), nan=-1.0, posinf=1e6, neginf=-1e6)
    rf = joblib.load(f"{md}/inv.rf{j['seeds'][0]}.joblib")
    rf.n_jobs = 1
    return bool(np.array_equal(rf.predict_proba(X)[:, 1], T.forest_proba(rf, X))), len(X)


def main():
    ap = argparse.ArgumentParser()
    for x in ("--feat", "--det", "--bin", "--exclude-reads", "--train", "--test", "--out"):
        ap.add_argument(x, required=True)
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--trees", type=int, default=60)
    ap.add_argument("--cap", type=int, default=3000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    keys = a.train.split(",") + [a.test]
    # the renamed, shuffled copy of the labelled tables and of the excluded list
    fr = os.path.join(a.out, "feat_renamed")
    os.makedirs(fr, exist_ok=True)
    back = {}
    for i, k in enumerate(keys):
        rename_table(f"{a.feat}/{k}.tsv", f"{fr}/{k}.tsv", i)
        with open(f"{a.feat}/{k}.tsv") as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            ia, ib = cols.index("name_a"), cols.index("name_b")
            for ln in fh:
                f = ln.split("\t")
                back[new_name(f[ia])], back[new_name(f[ib])] = f[ia], f[ib]
    excl = [x.strip() for x in open(a.exclude_reads) if x.strip()]
    er = os.path.join(a.out, "excluded_renamed.txt")
    with open(er, "w") as fh:
        fh.write("".join(f"{new_name(x)}\n" for x in excl))
    named = sum(1 for x in excl if x in back.values())
    print(f"tables {keys}; {len(back):,} reads renamed; excluded list {len(excl)} reads, {named} of them in these tables")
    m1, s1 = one("as_written", a.feat, None, a.exclude_reads, a, a.out)
    m1b, s1b = one("as_written_2", a.feat, None, a.exclude_reads, a, a.out)
    m2, _ = one("renamed_shuffled", fr, s1, er, a, a.out)
    ok = True
    for d1, d2 in ((m1, m1b), (s1, s1b)):
        n, bad = same_tree(d1, d2)
        ok &= not bad
        print(f"  determinism, {os.path.basename(d1)} vs {os.path.basename(d2)}: {n} files, "
              + ("all byte-identical" if not bad else f"DIFFER: {bad[:8]}"))
    files = [f"inv__{k}.oob.tsv" for k in a.train.split(",")] + [f"inv__{a.test}.scores.tsv", f"s2inv__{a.test}.scores.tsv"]
    for f in files:
        o1, o2 = score_map(f"{m1}/{f}"), score_map(f"{m2}/{f}", back)
        same = o1 == o2
        ok &= same
        print(f"  renamed + shuffled, {f}: {len(o1):,} pairs, {'IDENTICAL' if same else 'DIFFER'}")
    eq, n = sklearn_check(m1, a.feat, a.test)
    ok &= eq
    print(f"  forest_proba vs sklearn predict_proba (n_jobs=1) on {a.test}: {n:,} rows, {'IDENTICAL' if eq else 'DIFFER'}")
    for tag in ("inv", "s2inv"):
        j1, j2 = (json.load(open(f"{m}/{tag}.thresholds.json")) for m in (m1, m2))
        for m in (j1, j2):
            m.get("excluded_reads", {}).pop("sha256", None)
        jsame = j1 == j2
        files = [f"{tag}.rf{s}.rfb" for s in j1["seeds"]] + [f"{tag}.rfmeta"]
        fsame = {f: filecmp.cmp(f"{m1}/{f}", f"{m2}/{f}", shallow=False) for f in files}
        ok &= jsame and all(fsame.values())
        print(f"  {tag}: " + ", ".join(f"{f} {'IDENTICAL' if v else 'DIFFER'}" for f, v in fsame.items())
              + f"; thresholds.json (but the list's hash) {'IDENTICAL' if jsame else 'DIFFER'}"
              + f"; excluded {j1.get('excluded_reads', {}).get('n')} / {j2.get('excluded_reads', {}).get('n')}")
    print(f"training determinism and invariance (two as-written runs; renamed + shuffled inputs): {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
