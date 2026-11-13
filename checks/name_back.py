#!/usr/bin/env python3
"""A haricot assemble run on renamed reads (pod5_variant.py --mode reno / renr), with its read ids mapped back to the
source pod5's, so sym_compare.py reads it as a run on the source reads. invariance.sh runs it.

Copies RUN to OUT: every .tsv / .gfa / .txt file with each tab-separated field that is a renamed id replaced by its
source id; every other file except the .bin caches (the comparator reads their canonical dumps) as it is. RUN is only
read.

Usage: name_back.py RUN MAP OUT     (MAP: `dst_id<TAB>src_id` per read, pod5_variant.py's)
"""
import argparse
import os
import shutil
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="the run on the renamed reads")
    ap.add_argument("map", help="pod5_variant.py's --map")
    ap.add_argument("out", help="the run with the source ids (written; replaced if it exists)")
    a = ap.parse_args()
    run, mp, out = a.run, a.map, a.out
    back = {}
    with open(mp) as f:
        for line in f:
            x, y = line.rstrip("\n").split("\t")
            back[x] = y
    if len(set(back.values())) != len(back):
        sys.exit(f"{mp}: not one-to-one (a dup-mode map); name_back maps renamed runs only")
    if os.path.exists(out):
        shutil.rmtree(out)
    n_files = n_fields = 0
    for d, _, fs in os.walk(run):
        rel = os.path.relpath(d, run)
        os.makedirs(os.path.join(out, rel), exist_ok=True)
        for fn in fs:
            src, dst = os.path.join(d, fn), os.path.join(out, rel, fn)
            if fn.endswith(".bin"):
                continue
            if not fn.endswith((".tsv", ".gfa", ".txt")):
                shutil.copyfile(src, dst)
                continue
            n_files += 1
            with open(src) as fi, open(dst, "w") as fo:
                for line in fi:
                    end = "\n" if line.endswith("\n") else ""
                    fields = line[: len(line) - len(end)].split("\t")
                    for i, x in enumerate(fields):
                        y = back.get(x)
                        if y is not None:
                            fields[i] = y
                            n_fields += 1
                    fo.write("\t".join(fields) + end)
    print(f"{run} -> {out}: {n_files} text files, {n_fields} read ids mapped back")


if __name__ == "__main__":
    main()
