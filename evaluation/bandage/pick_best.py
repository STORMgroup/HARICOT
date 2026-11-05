#!/usr/bin/env python3
"""Pick the least-folded Bandage render from several random-seed candidates.

Bandage's force-directed (FMMM) layout is random per run and, at low iteration
counts, folds long cycles onto themselves (a circular genome drawn as a triangle /
figure-8). With --iter 4 the layouts converge cleanly; this just picks the roundest
of a few seeds. Score = ink pixels / convex-hull area (LOWER = ink spread over a
larger hull = a thin, open ring rather than a folded blob). Copies the winner to
<out>. Falls back to the first candidate if scoring fails.

Usage: pick_best.py <out.png> <cand1.png> <cand2.png> ...
"""
import shutil
import sys

import numpy as np
from PIL import Image

out, cands = sys.argv[1], sys.argv[2:]


def score(p):
    try:
        a = np.asarray(Image.open(p).convert("RGB"))
        ink = a.sum(2) < 720          # non-white pixels
        ys, xs = np.nonzero(ink)
        if len(xs) < 50:
            return 1e9
        from scipy.spatial import ConvexHull
        pts = np.column_stack([xs, ys])[::7]
        return ink.sum() / ConvexHull(pts).volume
    except Exception:
        return 1e9


best = min(cands, key=score)
shutil.copy(best, out)
print(f"  picked {best.split('/')[-1]} of {len(cands)}  (density {score(best):.4f})")
