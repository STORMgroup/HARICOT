#!/usr/bin/env python3
"""The reads training never uses: every read whose id is in the base sets of two or more genomes of the training
table (a base set: the row whose key is its genome), and the exact-cut fragments of those reads (their own ids, from
the table's fragment tables), which carry the same signal under new ids. One id per line, sorted.

Usage: excluded_reads.py SETS.tsv OUT.txt   (SETS.tsv: train.sh's table)
"""
import csv
import sys

import pod5


def main():
    sets, out = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(sets), delimiter="\t"))
    genomes = {}
    for r in rows:
        if r["key"] == r["genome"]:
            with pod5.Reader(r["pod5"]) as rd:
                genomes[r["genome"]] = {str(x) for x in rd.read_ids}
    count = {}
    for ids in genomes.values():
        for x in ids:
            count[x] = count.get(x, 0) + 1
    base = {x for x, n in count.items() if n >= 2}
    frag = set()
    for f in sorted({r["fragments"] for r in rows if r.get("fragments", "-") not in ("", "-")}):
        with open(f) as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if c[1] in base:
                    frag.add(c[0])
    with open(out, "w") as w:
        w.write("".join(x + "\n" for x in sorted(base | frag)))
    print(f"excluded reads: {len(base)} in two or more base sets of {sorted(genomes)} + {len(frag)} fragments of them "
          f"= {len(base | frag)} ids -> {out}")


if __name__ == "__main__":
    main()
