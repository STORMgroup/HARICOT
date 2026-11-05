#!/usr/bin/env python3
"""Bandage colour table for a sequence graph (Flye's assembly_graph.gfa, Rawsamble's asm.gfa, a FASTA as GFA).

A segment's position is the reference midpoint of its longest primary alignment block (minimap2
PAF of the segment sequences to the reference). Hue ramps 0 -> 0.85 along --contig, the same scale
for every tool of a dataset.
Segments placed on another contig are dark grey; unplaced segments light grey.

Usage: color_gfa.py SEGMENTS.paf --gfa G.gfa --fai contigs.tsv --contig NAME --out colors.csv
"""
import argparse
import colorsys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf")
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    length = next(int(l.split("\t")[1]) for l in open(a.fai) if l.split("\t")[0] == a.contig)
    best = {}
    for line in open(a.paf):
        f = line.rstrip("\n").split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        ql = int(f[3]) - int(f[2])
        if ql > best.get(f[0], (0,))[0]:
            best[f[0]] = (ql, f[5], (int(f[7]) + int(f[8])) // 2)
    segs = [l.split("\t")[1] for l in open(a.gfa) if l.startswith("S\t")]
    n_on = n_off = n_none = 0
    with open(a.out, "w") as fh:
        fh.write("Node,Colour\n")
        for s in segs:
            b = best.get(s)
            if b is None:
                col, n_none = "#cccccc", n_none + 1
            elif b[1] != a.contig:
                col, n_off = "#555555", n_off + 1
            else:
                h = 0.85 * min(max(b[2] / length, 0.0), 1.0)
                rgb = colorsys.hsv_to_rgb(h, 0.85, 0.95)
                col = "#%02x%02x%02x" % tuple(int(round(255 * c)) for c in rgb)
                n_on += 1
            fh.write(f"{s},{col}\n")
    print(f"{a.out}: {n_on} on {a.contig}, {n_off} elsewhere, {n_none} unplaced")


if __name__ == "__main__":
    main()
