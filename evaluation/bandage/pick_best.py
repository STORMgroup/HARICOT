#!/usr/bin/env python3
"""Pick the least-folded Bandage render from several random layouts of one graph.

Bandage's force-directed layout is random per run and can fold long cycles onto themselves (a circular genome drawn
as a triangle or a figure-8). This picks the roundest of a few layouts. Score = ink pixels / convex-hull area
(lower = ink spread over a larger hull = a thin, open ring rather than a folded blob). Copies the winner to OUT.png;
takes the first candidate if no layout can be scored. panels.sh runs it.

Usage: pick_best.py OUT.png CAND1.png [CAND2.png ...]
"""
import argparse
import shutil

import numpy as np
from PIL import Image


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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="the chosen render (written)")
    ap.add_argument("cands", nargs="+", help="the renders of the same graph")
    a = ap.parse_args()
    best = min(a.cands, key=score)
    shutil.copy(best, a.out)
    print(f"  picked {best.split('/')[-1]} of {len(a.cands)}  (density {score(best):.4f})")


if __name__ == "__main__":
    main()
