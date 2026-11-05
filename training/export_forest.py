#!/usr/bin/env python3
"""Export saved forests (train_stage1.py / train_stage2.py --save-models) for `haricot-forest`.

`MODEL_DIR/TAG.rf<seed>.joblib`, the seeds and features of `MODEL_DIR/TAG.thresholds.json`
-> `OUT_DIR/TAG.rf<seed>.rfb` per seed and `OUT_DIR/TAG.rfmeta` (features, seeds, thresholds).
Loading a 400-tree joblib forest takes 3-8 s (zlib + unpickling) to predict in under a second;
`haricot-forest` reads the compact form in a fraction of that and predicts the same probabilities.

`.rfb` (little-endian): b"RFB1", n_features u32, n_trees u32, n_values u32, values f64[n_values]
(the distinct leaf probabilities of class 1), then per tree: n_nodes u32 and its nodes in
preorder (a node's left child is the next node), each 12 bytes:
  internal: feature u16, 0u16, threshold f32, right-child index u32
  leaf:     0xFFFF u16, 0u16, 0f32, value index u32
The threshold is the largest f32 <= sklearn's float64 threshold: sklearn compares float32 inputs
`x <= t` in float64, and for any float32 `x` that is the same test as `x <= f32_down(t)`. Nodes
keep sklearn's numbering, which is already preorder (asserted).

Usage: export_forest.py --model-dir M --tag all5 [--out-dir M]
"""
import argparse
import json
import os
import struct

import joblib
import numpy as np


def export(rf, path):
    trees = [e.tree_ for e in rf.estimators_]
    vals = np.unique(np.concatenate([t.value[t.children_left == -1, 0, 1] for t in trees]))
    with open(path, "wb") as fh:
        fh.write(b"RFB1")
        fh.write(struct.pack("<III", rf.n_features_in_, len(trees), len(vals)))
        fh.write(vals.astype("<f8").tobytes())
        for t in trees:
            left, right, n = t.children_left, t.children_right, t.node_count
            leaf = left == -1
            inner = np.flatnonzero(~leaf)
            # sklearn's depth-first builder numbers nodes in preorder (a left child follows its
            # parent); the format relies on it.
            assert np.all(left[inner] == inner + 1), "nodes not in preorder"
            thr = t.threshold.astype(np.float32)
            down = thr.astype(np.float64) > t.threshold
            thr[down] = np.nextafter(thr[down], np.float32(-np.inf))
            assert np.all(thr[inner].astype(np.float64) <= t.threshold[inner])
            up = np.nextafter(thr[inner], np.float32(np.inf)).astype(np.float64)
            assert np.all(up > t.threshold[inner])
            rec = np.zeros(n, dtype=[("f", "<u2"), ("pad", "<u2"), ("t", "<f4"), ("x", "<u4")])
            rec["f"] = np.where(leaf, 0xFFFF, t.feature)
            rec["t"] = np.where(leaf, np.float32(0), thr)
            vi = np.searchsorted(vals, t.value[:, 0, 1])
            assert np.all(vals[vi[leaf]] == t.value[leaf, 0, 1])
            rec["x"] = np.where(leaf, vi, right)
            fh.write(struct.pack("<I", n))
            fh.write(rec.tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out-dir", default="")
    a = ap.parse_args()
    out = a.out_dir or a.model_dir
    os.makedirs(out, exist_ok=True)
    meta = json.load(open(os.path.join(a.model_dir, f"{a.tag}.thresholds.json")))
    for sd in meta["seeds"]:
        rf = joblib.load(os.path.join(a.model_dir, f"{a.tag}.rf{sd}.joblib"))
        assert list(rf.classes_) == [0, 1] and rf.n_outputs_ == 1
        assert rf.n_features_in_ == len(meta["features"])
        export(rf, os.path.join(out, f"{a.tag}.rf{sd}.rfb"))
        print(f"{a.tag}.rf{sd}: {len(rf.estimators_)} trees", flush=True)
    with open(os.path.join(out, f"{a.tag}.rfmeta"), "w") as fh:
        fh.write("features\t" + "\t".join(meta["features"]) + "\n")
        fh.write("seeds\t" + "\t".join(str(s) for s in meta["seeds"]) + "\n")
        for k, v in sorted(meta["thresholds"].items(), key=lambda kv: float(kv[0]), reverse=True):
            fh.write(f"threshold\t{k}\t{float(v)!r}\n")


if __name__ == "__main__":
    main()
