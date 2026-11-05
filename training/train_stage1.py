#!/usr/bin/env python3
"""Stage 1: train the overlap classifier (random forests) on the training sets, freeze its thresholds, and score
the --test sets. training/train.sh runs it for the stage-1 and rescue forests.

Labels are the signal-exact oracle's (haricot-pair-features): `true` = in the strict oracle with a
correctly placed path; `false` = no placed intersection; `misplaced` = a true pair whose
detected path is off (harmful to the layout, so a negative); `short` = a real intersection the
strict oracle does not hold -- UNKNOWN, never trained on, scored like everything else.

Features are restricted to quantities comparable across genomes, depths and read lengths:
ratios, fractions, per-bin densities, bin-scale geometry (a bin is ~267 bp everywhere), and
statistics normalised by the dataset's own median. Raw counts that grow with coverage
(third-read counts, candidate degree) are excluded.

Operating points are FROZEN on the training sets: the score threshold that keeps a fixed
share (--keep) of the training sets' true candidates (out-of-bag scores, so the threshold is not
set on training fit). Each test dataset is then filtered at that same threshold.

Writes, per test dataset and operating point, a drop list (`name_a name_b`, for
`haricot-ovl-filter --drop-pairs`), and prints overlap-level recall/precision and AUC.

Training is a function of the rows' values, not of read names or row order: the per-set
cap draws rows by a hash of the set key and the row's canonical-frame values (`sub_u`, no read name), and the training
matrix is put in a canonical order (`canon_order`) before every fit, so renaming the reads or reordering a table gives
the same forests. `--exclude-reads FILE` drops every row naming a listed read from training (before the cap and the
bootstrap; the rows are still scored, out of sample) and records the list next to the forest (TAG.excluded_reads.txt,
its SHA-256 in TAG.thresholds.json).

Usage: train_stage1.py --train K1,K2,... [--test K3,...] --tag TAG --features-dir DIR --out-dir DIR
                    [--keep 0.995,0.99,0.98,0.95] [--seeds 0,1,2,3,4] [--save-models] [--cap-per-set N]
                    [--exclude-reads FILE] [--drop-features F1,F2,...]
"""
import argparse
import hashlib
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

FEATURES = [
    # geometry
    "z", "slope", "abs_log_slope", "a_bins", "b_bins", "span_a", "ov_frac_short",
    "oh_left", "oh_right", "oh_frac", "a_contained", "b_contained",
    # path shape
    "hole_frac", "max_hole", "resid_med", "resid_max", "resid_frac1",
    "anc_per_bin", "dens_cv", "bins_hit",
    # read context / repeat
    "high_occ_frac", "occ_rel_med", "occ_rel_q90", "ctx_rel_a", "ctx_rel_b", "ctx_abs",
    "deg_rel", "tri_frac", "tri_rel", "ctx_pois", "ctx_pois_abs",
    # signal-level DTW
    "dtw_mean", "dtw_decoy", "dtw_ratio", "dtw_good_frac", "dtw_r_q10", "dtw_r_q50",
    "dtw_r_q90", "dtw_r_max",
]


LABELLED = ["true", "false", "misplaced"]
# The row's identity for the draw: the set key and these canonical-frame values of the pair (haricot-pair-features),
# never its read names or its position in the table.
HASHC = ["z", "n_anchors", "aligned_bins", "a_lo", "a_hi", "b_lo", "b_hi", "line_s", "line_c", "dtw_mean", "dtw_decoy"]
# The rule as recorded in every forest's thresholds.json (excluded_reads.rule): kept verbatim, so a retrained forest's
# record is byte-identical to the published one (the text names the paths and script of the original training run,
# not files of this tree). The list itself: training/excluded_reads.py.
EXCL_RULE = ("reads whose id is in the base sets (data/signal_oracle/reads/<g>/ids.txt) of two or more genomes of the "
             "training pool, and the exact-cut fragments of those reads (data/signal_oracle/frag/*.frag.tsv); "
             "scripts/train/excluded_reads.sh")


def load(fdir, key):
    # `cls` as text: pandas parses a read chunk holding only true / false as booleans, and the rows of
    # such a chunk would then fall out of the labelled set
    d = pd.read_csv(os.path.join(fdir, f"{key}.tsv"), sep="\t", dtype={"cls": str})
    d["key"] = key
    for f in FEATURES:  # a table without a feature column reads it as NaN
        if f not in d.columns:
            d[f] = np.nan
    return d


def sub_u(d):
    """A deterministic uniform draw in [0, 1) per row, from the set key and the row's canonical-frame values (HASHC;
    no read name, no row position): `--subsample` keeps the labelled rows with a draw below F, `cap_mask` the rows
    with the smallest draws, the same rows in both stages (train_stage2.py)."""
    miss = [c for c in HASHC if c not in d.columns]
    if miss:
        raise SystemExit(f"the row draw needs the columns {miss} (haricot-pair-features' table)")
    hd = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").astype(np.float64).to_numpy() for c in HASHC})
    hd.insert(0, "key", d["key"].astype(str).to_numpy())
    h = pd.util.hash_pandas_object(hd, index=False).to_numpy()
    return (h >> np.uint64(11)).astype(np.float64) / float(1 << 53)


def tie_cols(d):
    """The tie-break of the draw: the row's feature and HASHC values (NaN as -1), as float64 columns."""
    cols = [c for c in dict.fromkeys(HASHC + FEATURES) if c in d.columns]
    return np.nan_to_num(np.column_stack([pd.to_numeric(d[c], errors="coerce").to_numpy(np.float64) for c in cols]),
                         nan=-1.0)


def canon_order(u, *cols):
    """The canonical row order: by the draw `u`, ties by `cols` in turn (values only), so no read name and no row
    position decides it. Rows equal in all of them are interchangeable."""
    return np.lexsort([np.asarray(c) for c in reversed(cols)] + [np.asarray(u)])


def read_excluded(path):
    """The read ids of an --exclude-reads file (one per line), as a sorted list."""
    return sorted({x.strip() for x in open(path) if x.strip()}) if path else []


def excluded_mask(d, excl):
    """Rows naming an excluded read (either side)."""
    if not excl:
        return np.zeros(len(d), dtype=bool)
    s = set(excl)
    return (d["name_a"].astype(str).isin(s) | d["name_b"].astype(str).isin(s)).to_numpy()


def record_excluded(excl, out_dir, tag, meta):
    """Write the list next to the forest (OUT_DIR/TAG.excluded_reads.txt, sorted, one id per line) and its SHA-256,
    length and rule into `meta` (TAG.thresholds.json)."""
    p = os.path.join(out_dir, f"{tag}.excluded_reads.txt")
    b = "".join(f"{x}\n" for x in excl).encode()
    with open(p, "wb") as fh:
        fh.write(b)
    meta["excluded_reads"] = {"sha256": hashlib.sha256(b).hexdigest(), "n": len(excl), "rule": EXCL_RULE,
                              "file": os.path.basename(p)}


def cap_mask(d, cap, exclude=None):
    """`--cap-per-set`: per training set (`key`), at most `cap` labelled rows, class-stratified: each class
    keeps round(cap x its share of the set's labelled rows) rows, those with the smallest draw `sub_u`. A set
    with at most `cap` labelled rows keeps them all. Returns a boolean mask over `d`'s rows (unlabelled rows
    are never selected). Per set, it depends only on that set's rows, so train_stage2.py recomputes the same
    selection from one set's table."""
    keys, cls = d["key"].astype(str).to_numpy(), d["cls"].astype(str).to_numpy()
    lab = np.isin(cls, LABELLED)
    if exclude is not None:  # --exclude-reads: never trained on, so not counted in the cap either
        lab = lab & ~exclude
    u = sub_u(d)
    keep = np.zeros(len(d), dtype=bool)
    for k in np.unique(keys[lab]):
        idx = np.flatnonzero(lab & (keys == k))
        if len(idx) <= cap:
            keep[idx] = True
            continue
        for c in np.unique(cls[idx]):
            ic = idx[cls[idx] == c]
            m = int(round(cap * len(ic) / len(idx)))
            keep[ic[canon_order(u[ic], *tie_cols(d.iloc[ic]).T)[:m]]] = True
    return keep


def forest_proba(rf, X, chunk=32768):
    """The forest's class-1 probability: sklearn's predict_proba (the mean of the trees' probabilities), with the sum
    taken in `rf.estimators_` order for every row. sklearn with n_jobs > 1 adds the trees' predictions under a lock in
    thread completion order (sklearn/ensemble/_forest.py `_accumulate_prediction`), so the last bits of a score varied
    from run to run. Rows are scored in chunks in parallel threads; within a row the order
    is fixed, so the result is sklearn's single-threaded one."""
    from concurrent.futures import ThreadPoolExecutor
    X = np.ascontiguousarray(X, dtype=np.float32)
    out = np.zeros(len(X))
    if not len(X):
        return out

    def part(lo):
        hi = min(lo + chunk, len(X))
        acc = np.zeros(hi - lo)
        for e in rf.estimators_:
            acc += e.predict_proba(X[lo:hi], check_input=False)[:, 1]
        out[lo:hi] = acc / len(rf.estimators_)

    with ThreadPoolExecutor(max_workers=len(os.sched_getaffinity(0))) as ex:
        list(ex.map(part, range(0, len(X), chunk)))
    return out


def fit_forests(X, y, seeds, score_sets, save=None, **kw):
    """One forest per seed, trained one at a time so only one is ever in memory; returns the mean
    out-of-bag score of the training rows, the mean score of each matrix in `score_sets`, and the
    mean feature importances. One seed is exactly the single forest. With `save` (a path prefix)
    each forest is also written to `<save>.rf<seed>.joblib`, which export_forest.py exports for haricot-forest.
    Scores are `forest_proba`'s (a fixed summation order: run to run identical); the out-of-bag ones sklearn's
    oob_decision_function_ (computed sequentially)."""
    oob, preds, imp = np.zeros(len(y)), [np.zeros(len(m)) for m in score_sets], 0.0
    for sd in seeds:
        rf = RandomForestClassifier(max_features="sqrt", n_jobs=-1, oob_score=True, random_state=sd, **kw)
        rf.fit(X, y)
        if save:
            import joblib
            joblib.dump(rf, f"{save}.rf{sd}.joblib", compress=3)
        oob += rf.oob_decision_function_[:, 1]
        for p, m in zip(preds, score_sets):
            if len(m):
                p += forest_proba(rf, m)
        imp = imp + rf.feature_importances_
        del rf
    n = len(seeds)
    return oob / n, [p / n for p in preds], imp / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", default="", help="comma list of sets to score with the trained forests (printed)")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--keep", default="0.995,0.99,0.98,0.95")
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--drop-features", default="")
    ap.add_argument("--extra-features", default="",
                    help="comma list of feature columns beyond FEATURES (absent columns read as NaN)")
    ap.add_argument("--trees", type=int, default=400)
    ap.add_argument("--leaf", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0, help="forest random_state")
    ap.add_argument("--seeds", default="",
                    help="comma list of random_states: an ENSEMBLE of that many forests, trained one at a "
                         "time, whose scores (out-of-bag and test) are averaged; overrides --seed")
    ap.add_argument("--save-models", action="store_true",
                    help="write each forest to OUT_DIR/TAG.rfSEED.joblib (export_forest.py exports them)")
    ap.add_argument("--subsample", type=float, default=1.0,
                    help="train on a deterministic share F of the labelled rows (sub_u); the other "
                         "labelled rows are scored like the unknown ones, out of sample, so every "
                         "training row still gets an out-of-sample score; thresholds from the sampled "
                         "rows' out-of-bag scores. 1 (default): every labelled row, as before")
    ap.add_argument("--cap-per-set", type=int, default=0,
                    help="train on at most N labelled rows per training set, class-stratified (cap_mask); "
                         "the other labelled rows are scored like the unknown ones, out of sample, as with "
                         "--subsample. 0 (default): no cap")
    ap.add_argument("--exclude-reads", default="",
                    help="a file of read ids (one per line): every row naming one is dropped from training, before the "
                         "cap and the bootstrap (it is still scored, out of sample); the list and its SHA-256 are "
                         "recorded with the forest (TAG.excluded_reads.txt, TAG.thresholds.json)")
    a = ap.parse_args()
    if a.cap_per_set and a.subsample < 1.0:
        raise SystemExit("--cap-per-set and --subsample are exclusive")
    extra = [f for f in a.extra_features.split(",") if f]
    feats = [f for f in FEATURES + extra if f not in set(a.drop_features.split(","))]
    FEATURES.extend(f for f in extra if f not in FEATURES)  # so load() fills absent ones
    os.makedirs(a.out_dir, exist_ok=True)

    clean = lambda d: np.nan_to_num(d[feats].to_numpy(dtype=np.float32), nan=-1.0, posinf=1e6,
                                    neginf=-1e6)
    full = pd.concat([load(a.features_dir, k) for k in a.train.split(",")], ignore_index=True)
    excl = read_excluded(a.exclude_reads)
    ex = excluded_mask(full, excl)
    if a.exclude_reads:
        lab_ex = ex & full.cls.isin(LABELLED).to_numpy()
        print(f"[{a.tag}] --exclude-reads {a.exclude_reads}: {len(excl)} reads; labelled rows dropped from training "
              f"per set: " + (", ".join(f"{k} {n:,}" for k, n in full[lab_ex].groupby("key").size().items()) or "none"))
    lab_all = full.cls.isin(LABELLED).to_numpy() & ~ex
    if a.subsample < 1.0:  # the trained rows; the rest of the labelled rows are scored like unknowns
        lab_all = lab_all & (sub_u(full) < a.subsample)
    if a.cap_per_set:
        lab_all = cap_mask(full, a.cap_per_set, exclude=ex)
        print(f"[{a.tag}] --cap-per-set {a.cap_per_set}: training rows per set " + ", ".join(
            f"{k} {n:,}" for k, n in full[lab_all].groupby("key").size().items()))
    assert not (lab_all & ex).any(), "a training row names an excluded read"
    tr = full[lab_all]
    X = clean(tr)
    y = (tr.cls == "true").to_numpy().astype(int)
    # the canonical row order (no read name, no table position): the forests see the same matrix whatever the names
    # and the order of the input rows
    o = canon_order(sub_u(tr), y, *X.T)
    X, y = X[o], y[o]
    print(f"[{a.tag}] train {a.train}: {len(y):,} candidates, {int((y == 0).sum()):,} negative, "
          f"{len(feats)} features")
    tests = {key: load(a.features_dir, key) for key in a.test.split(",") if key}
    seeds = [int(x) for x in a.seeds.split(",")] if a.seeds else [a.seed]
    oob_c, (s_unl, *s_test), imp = fit_forests(
        X, y, seeds, [clean(full[~lab_all])] + [clean(d) for d in tests.values()],
        save=os.path.join(a.out_dir, a.tag) if a.save_models else None,
        n_estimators=a.trees, min_samples_leaf=a.leaf)
    oob = np.empty(len(oob_c))
    oob[o] = oob_c  # back to the table's row order
    y = y[np.argsort(o)]
    if len(seeds) > 1:
        print(f"[{a.tag}] ensemble of {len(seeds)} forests (seeds {a.seeds}), scores averaged")
    print(f"[{a.tag}] out-of-bag AUC {roc_auc_score(y, oob):.4f}")
    keeps = [float(x) for x in a.keep.split(",")]
    thr = {k: float(np.quantile(oob[y == 1], 1 - k)) for k in keeps}
    for k in keeps:
        fr = ((oob >= thr[k]) & (y == 0)).sum() / max((y == 0).sum(), 1)
        print(f"[{a.tag}]   keep {k:.3f} of true -> threshold {thr[k]:.4f}; training negatives "
              f"kept {fr:.3f}")
    # Out-of-bag scores of the TRAINING candidates, for a second stage that must not see
    # in-sample scores; the unknown ("short") ones get plain predictions.
    s_full = np.empty(len(full))
    s_full[lab_all], s_full[~lab_all] = oob, s_unl
    for key in a.train.split(","):
        d = full[(full.key == key).to_numpy()].assign(score=s_full[(full.key == key).to_numpy()])
        d[["name_a", "name_b", "cls", "score"]].to_csv(
            os.path.join(a.out_dir, f"{a.tag}__{key}.oob.tsv"), sep="\t", index=False)
    import json
    meta = {"thresholds": thr, "train": a.train.split(","), "features": feats, "seeds": seeds}
    if a.exclude_reads:
        record_excluded(excl, a.out_dir, a.tag, meta)
    json.dump(meta, open(os.path.join(a.out_dir, f"{a.tag}.thresholds.json"), "w"), indent=1)
    imp = sorted(zip(imp, feats), reverse=True)
    print(f"[{a.tag}] importance: " + ", ".join(f"{f} {v:.3f}" for v, f in imp[:14]))

    for (key, d), s in zip(tests.items(), s_test):
        d["score"] = s
        lab = d.cls.isin(["true", "false", "misplaced"])
        yt = (d.cls[lab] == "true").to_numpy().astype(int)
        au = roc_auc_score(yt, s[lab.to_numpy()]) if len(np.unique(yt)) == 2 else np.nan
        n_true = int((d.cls == "true").sum())
        n_neg = int(d.cls.isin(["false", "misplaced"]).sum())
        print(f"\n[{a.tag}] TEST {key}: {len(d):,} candidates ({n_true:,} true, {n_neg:,} neg, "
              f"{int((d.cls == 'short').sum()):,} unknown); AUC {au:.4f}")
        for k in keeps:
            keep = s >= thr[k]
            kt = int((keep & (d.cls == "true")).sum())
            kn = int((keep & d.cls.isin(["false", "misplaced"])).sum())
            print(f"   keep {k:.3f}: true kept {kt:,}/{n_true:,} ({kt / max(n_true, 1):.3f})  "
                  f"negatives kept {kn:,}/{n_neg:,} ({kn / max(n_neg, 1):.3f})  "
                  f"precision {kt / max(kt + kn, 1):.4f}  dropped {int((~keep).sum()):,}")
            out = os.path.join(a.out_dir, f"{a.tag}__{key}__keep{k:.3f}.drop.tsv")
            d.loc[~keep, ["name_a", "name_b"]].to_csv(out, sep="\t", header=False, index=False)
        d[["name_a", "name_b", "cls", "score"]].to_csv(
            os.path.join(a.out_dir, f"{a.tag}__{key}.scores.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
