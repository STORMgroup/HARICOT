#!/usr/bin/env python3
"""Stage 2: re-decide the stage-1 ACCEPTED overlaps using what the accepted graph says.

Input per dataset: `haricot-graph-features` graph features (mutual support, offset consistency, region
partners, neighbour scores, rank) for every candidate stage 1 accepted at its frozen keep-0.995
threshold, merged with that candidate's stage-1 features and score. Labels are the exact ones
(`haricot-pair-features`). Trained on the training genomes (stage-1 scores there are OUT-OF-BAG, so stage 2
never sees in-sample scores), thresholds frozen on stage-2 out-of-bag scores at a fixed share of
the accepted true candidates kept, applied unchanged to the test datasets.

Writes per test dataset and operating point a drop list for `haricot-ovl-filter --drop-pairs` that removes
EVERY candidate not kept: the ones stage 1 rejected plus the ones stage 2 rejects, so it applies
to the permissive cache directly.

Training rows are drawn (the cap) and ordered (`canon_order`) by their values only, as in train_stage1.py, so renaming
the reads or reordering the tables gives the same forests; `--exclude-reads FILE` drops every row naming a listed read
before the cap and the bootstrap and records the list with the forest (train_stage1.record_excluded).

Usage: train_stage2.py --train K1,K2,... [--test K3,...] --tag TAG --features-dir DIR --s2-dir DIR --out-dir DIR
                       [--keep 0.995,0.99,0.98,0.95] [--seeds 0,1,2,3,4] [--save-models] [--cap-per-set N]
                       [--exclude-reads FILE]
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from train_stage1 import (LABELLED, canon_order, cap_mask, excluded_mask, fit_forests, read_excluded,
                       record_excluded, sub_u)

S2 = ["s1", "mut_a", "mut_b", "cons_a", "cons_b", "nbr_a", "nbr_b", "nb_score", "rank_a", "rank_b"]
S1 = ["z", "ov_frac_short", "oh_frac", "oh_left", "oh_right", "a_bins", "b_bins", "span_a",
      "hole_frac", "anc_per_bin", "ctx_abs", "ctx_pois", "deg_rel", "tri_frac", "tri_rel",
      "occ_rel_q90", "dtw_ratio"]


def load(fdir, sdir, key):
    f = pd.read_csv(os.path.join(fdir, f"{key}.tsv"), sep="\t", dtype={"cls": str})  # as train_stage1.load
    s = pd.read_csv(os.path.join(sdir, f"{key}.s2.tsv"), sep="\t")
    d = s.merge(f, on=["name_a", "name_b"], how="left")
    d["key"] = key
    for c in S1:
        if c not in d.columns:
            d[c] = np.nan
    return d, f


def stage2_selection(f, cap, true_only=False, excl=()):
    """The stage-1 rows of one set (its full feature table `f`, with `key`) that stage 2 may train on under
    `--cap-per-set CAP`: train_stage1.cap_mask's class-stratified selection. With `true_only` (`--cap-true-only`)
    the cap applies to the true rows only: the true rows are cap_mask's (the same rows),
    and every other labelled row (false, misplaced) is kept. Rows naming a read of `excl` are never selected."""
    ex = excluded_mask(f, list(excl))
    m = cap_mask(f, cap, exclude=ex)
    if true_only:
        m = m | (f["cls"].astype(str).isin(LABELLED) & (f["cls"].astype(str) != "true")).to_numpy()
    return m & ~ex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", default="", help="comma list of sets to score with the trained forests (printed)")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--keep", default="0.995,0.99,0.98,0.95")
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--s2-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--only-s2", action="store_true", help="graph features only (+ s1)")
    ap.add_argument("--extra-features", default="",
                    help="comma list of stage-1 feature columns to add (absent columns read as NaN)")
    ap.add_argument("--drop-features", default="", help="comma list of feature columns to leave out")
    ap.add_argument("--trees", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0, help="forest random_state")
    ap.add_argument("--seeds", default="",
                    help="comma list of random_states: an ENSEMBLE of forests whose scores are averaged "
                         "(train_stage1.fit_forests); overrides --seed")
    ap.add_argument("--save-models", action="store_true",
                    help="write each forest to OUT_DIR/TAG.rfSEED.joblib (export_forest.py exports them)")
    ap.add_argument("--subsample", type=float, default=1.0,
                    help="train on the rows train_stage1.py --subsample F trained stage 1 on (the same "
                         "draw, sub_u), whose s1 is out-of-bag. 1 (default): every labelled row")
    ap.add_argument("--cap-per-set", type=int, default=0,
                    help="train on the rows train_stage1.py --cap-per-set N trained stage 1 on (cap_mask over "
                         "each set's full table; their s1 is out-of-bag) among the stage-1 accepted ones, "
                         "so at most N per set. 0 (default): no cap")
    ap.add_argument("--cap-true-only", action="store_true",
                    help="with --cap-per-set: cap the true rows only (cap_mask's true rows); every false and "
                         "misplaced accepted row is kept. Off (default): the "
                         "class-stratified cap")
    ap.add_argument("--exclude-reads", default="",
                    help="a file of read ids (one per line): every row naming one is dropped from training, before the "
                         "cap and the bootstrap; the list and its SHA-256 are recorded with the forest (train_stage1.py)")
    a = ap.parse_args()
    if a.cap_true_only and not a.cap_per_set:
        raise SystemExit("--cap-true-only needs --cap-per-set")
    if a.cap_per_set and a.subsample < 1.0:
        raise SystemExit("--cap-per-set and --subsample are exclusive")
    extra = [f for f in a.extra_features.split(",") if f]
    S1.extend(f for f in extra if f not in S1)
    drop = set(a.drop_features.split(","))
    feats = [f for f in S2 + ([] if a.only_s2 else S1) if f not in drop]

    excl = read_excluded(a.exclude_reads)
    parts = []
    for k in a.train.split(","):
        d, f = load(a.features_dir, a.s2_dir, k)
        if a.cap_per_set:  # stage 1's capped rows of this set (out-of-bag s1), among the accepted ones
            f = f.assign(key=k)
            sel = f.loc[stage2_selection(f, a.cap_per_set, a.cap_true_only, excl), ["name_a", "name_b"]]
            d = d.merge(sel, on=["name_a", "name_b"], how="inner")
        parts.append(d)
    tr = pd.concat(parts, ignore_index=True)
    tr = tr[tr.cls.isin(LABELLED)]
    if a.exclude_reads:
        ex = excluded_mask(tr, excl)
        print(f"[{a.tag}] --exclude-reads {a.exclude_reads}: {len(excl)} reads; labelled accepted rows dropped per "
              f"set: " + (", ".join(f"{k} {n:,}" for k, n in tr[ex].groupby("key").size().items()) or "none"))
        tr = tr[~ex]
    assert not excluded_mask(tr, excl).any(), "a training row names an excluded read"
    if a.subsample < 1.0:
        tr = tr[sub_u(tr) < a.subsample]
    if a.cap_per_set:
        print(f"[{a.tag}] --cap-per-set {a.cap_per_set}{' --cap-true-only' if a.cap_true_only else ''}: "
              f"stage-2 training rows per set " + ", ".join(
            f"{k} {n:,}" for k, n in tr.groupby("key").size().items()))
    clean = lambda d: np.nan_to_num(d[feats].to_numpy(dtype=np.float32), nan=-1.0, posinf=1e6,
                                    neginf=-1e6)
    X = clean(tr)
    y = (tr.cls == "true").to_numpy().astype(int)
    o = canon_order(sub_u(tr), y, *X.T)  # the canonical row order (train_stage1.py)
    X, y = X[o], y[o]
    print(f"[{a.tag}] stage-2 train: {len(y):,} accepted candidates, {int((y == 0).sum()):,} negative")
    tests = {key: load(a.features_dir, a.s2_dir, key) for key in a.test.split(",") if key}
    seeds = [int(x) for x in a.seeds.split(",")] if a.seeds else [a.seed]
    oob, s_test, imp = fit_forests(X, y, seeds, [clean(d) for d, _ in tests.values()],
                                   save=os.path.join(a.out_dir, a.tag) if a.save_models else None,
                                   n_estimators=a.trees, min_samples_leaf=10)
    if len(seeds) > 1:
        print(f"[{a.tag}] ensemble of {len(seeds)} forests (seeds {a.seeds}), scores averaged")
    print(f"[{a.tag}] stage-2 out-of-bag AUC {roc_auc_score(y, oob):.4f}")
    keeps = [float(k) for k in a.keep.split(",")]
    thr = {k: float(np.quantile(oob[y == 1], 1 - k)) for k in keeps}
    for k in keeps:
        print(f"[{a.tag}]   keep {k}: threshold {thr[k]:.4f}; training negatives kept "
              f"{((oob >= thr[k]) & (y == 0)).sum() / max((y == 0).sum(), 1):.3f}")
    imp = sorted(zip(imp, feats), reverse=True)
    print(f"[{a.tag}] importance: " + ", ".join(f"{f} {v:.3f}" for v, f in imp[:14]))
    meta = {"thresholds": thr, "features": feats, "seeds": seeds}
    if a.exclude_reads:
        record_excluded(excl, a.out_dir, a.tag, meta)
    json.dump(meta, open(os.path.join(a.out_dir, f"{a.tag}.thresholds.json"), "w"), indent=1)

    for (key, (d, full)), s in zip(tests.items(), s_test):
        lab = d.cls.isin(["true", "false", "misplaced"]).to_numpy()
        yt = (d.cls[lab] == "true").to_numpy().astype(int)
        au = roc_auc_score(yt, s[lab]) if len(np.unique(yt)) == 2 else np.nan
        nt_all = int((full.cls == "true").sum())
        nn_all = int(full.cls.isin(["false", "misplaced"]).sum())
        print(f"\n[{a.tag}] TEST {key}: {len(d):,} stage-1 accepted of {len(full):,} candidates; "
              f"stage-2 AUC {au:.4f}")
        acc = set(zip(d.name_a, d.name_b))
        d.assign(score=s)[["name_a", "name_b", "cls", "score"]].to_csv(
            os.path.join(a.out_dir, f"{a.tag}__{key}.scores.tsv"), sep="\t", index=False)
        for k in keeps:
            keep = s >= thr[k]
            kt = int((keep & (d.cls == "true").to_numpy()).sum())
            kn = int((keep & d.cls.isin(["false", "misplaced"]).to_numpy()).sum())
            print(f"   keep {k}: true kept {kt:,}/{nt_all:,} ({kt / max(nt_all, 1):.3f} of ALL true)  "
                  f"negatives kept {kn:,}/{nn_all:,}  precision {kt / max(kt + kn, 1):.4f}")
            kept = set(zip(d.name_a[keep], d.name_b[keep]))
            drop = full[[(x, y) not in kept for x, y in zip(full.name_a, full.name_b)]]
            drop[["name_a", "name_b"]].to_csv(
                os.path.join(a.out_dir, f"{a.tag}__{key}__keep{k:.3f}.drop.tsv"),
                sep="\t", header=False, index=False)
        del acc


if __name__ == "__main__":
    main()
