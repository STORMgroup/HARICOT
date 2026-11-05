#!/usr/bin/env python3
"""Adjacencies of the SAMPLE's genome that the reference does not have, from split reads.

A read whose consecutive alignment segments land far apart (or on opposite strands) spans
a junction. Supported by one read it is a chimera; supported by several independent reads
it is the sample's own structure relative to the reference -- an inversion, a
translocation, the origin of a circular chromosome -- and an assembly that follows it is
RIGHT. Example: E. coli ATCC 25922 reads against NZ_CP009072.1 span a ~1.6 Mb inversion
(535 kb <-> 2,145 kb) with ~41 reads in all four orientations; chimeric junctions have 1 read.

Junctions are clustered at --res bp on both sides (orientation-free), and those with at
least --min-support reads are written as `pos_a pos_b support` in GLOBAL coordinates
(event_ref.py's axis).

Usage: junctions.py SEGS.tsv --out junctions.tsv [--min-support 3] [--res 20000] [--near 50000]
"""
import argparse
import collections


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("segs")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-support", type=int, default=3)
    ap.add_argument("--res", type=int, default=20000)
    ap.add_argument("--near", type=int, default=50000,
                    help="consecutive segments closer than this, same strand, are not a junction")
    a = ap.parse_args()

    segs = collections.defaultdict(list)
    for line in list(open(a.segs))[1:]:
        f = line.rstrip("\n").split("\t")
        segs[f[1]].append((int(f[7]), f[2], f[5], int(f[9]), int(f[10])))
    J = collections.defaultdict(set)
    for r, v in segs.items():
        v = sorted(x for x in v if x[1] == "mm2") or sorted(v)
        for x, y in zip(v, v[1:]):
            end_x = x[4] if x[2] == "+" else x[3]
            beg_y = y[3] if y[2] == "+" else y[4]
            if x[2] == y[2] and abs(end_x - beg_y) < a.near:
                continue
            lo, hi = sorted((end_x, beg_y))
            J[(lo // a.res, hi // a.res)].add(r)
    # Merge neighbouring grid cells so a junction on a cell boundary is not split.
    keys = sorted(J, key=lambda k: -len(J[k]))
    used, out = set(), []
    for k in keys:
        if k in used:
            continue
        reads = set()
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                kk = (k[0] + da, k[1] + db)
                if kk in J and kk not in used:
                    reads |= J[kk]
                    used.add(kk)
        out.append(((k[0] + 0.5) * a.res, (k[1] + 0.5) * a.res, len(reads)))
    kept = [x for x in out if x[2] >= a.min_support]
    with open(a.out, "w") as fh:
        fh.write("pos_a\tpos_b\tsupport\n")
        for pa, pb, n in kept:
            fh.write(f"{int(pa)}\t{int(pb)}\t{n}\n")
    single = sum(1 for x in out if x[2] == 1)
    print(f"{a.segs}: {len(out)} junction clusters, {single} single-read (chimeras), "
          f"{len(kept)} with >= {a.min_support} reads -> {a.out}")
    for pa, pb, n in sorted(kept, key=lambda x: -x[2])[:10]:
        print(f"  {n:4d} reads  {int(pa):>14,} <-> {int(pb):>14,}")


if __name__ == "__main__":
    main()
