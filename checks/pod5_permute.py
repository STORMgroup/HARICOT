#!/usr/bin/env python3
"""Write a NEW pod5 holding exactly the reads of SRC in another order. SRC is only read, never modified.

The pipeline numbers reads in pod5 iteration order (`pipeline::sketch_pod5_streaming`), so the order of
the written file is the order every step sees.

  --order identity   SRC's own order (writes only the order list unless --write is given)
  --order reverse    the reverse of SRC's order: every pair's lower-id read becomes the higher one
  --order shuffle    a seeded uniform shuffle (Python's random.Random(SEED).shuffle)

Always writes ORDER (one read id per line, in the new order). With a DST, writes the new pod5 and then
re-reads it to check that its read ids are exactly ORDER, in order (exit 1 otherwise).

Usage: pod5_permute.py SRC --order identity|reverse|shuffle [--seed N] --order-out ORDER [--dst DST]
"""
import argparse
import os
import random
import sys

import pod5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--order", required=True, choices=["identity", "reverse", "shuffle"])
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--order-out", required=True)
    ap.add_argument("--dst", default=None)
    a = ap.parse_args()
    if a.order == "shuffle" and a.seed is None:
        sys.exit("--order shuffle needs --seed")

    with pod5.Reader(a.src) as r:
        ids = [str(x) for x in r.read_ids]
    if len(set(ids)) != len(ids):
        sys.exit(f"{a.src}: duplicate read ids")
    order = list(ids)
    if a.order == "reverse":
        order.reverse()
    elif a.order == "shuffle":
        random.Random(a.seed).shuffle(order)
    with open(a.order_out, "w") as f:
        f.write("".join(f"{x}\n" for x in order))
    print(f"{a.src}: {len(ids)} reads, order {a.order}" + (f" seed {a.seed}" if a.seed is not None else "")
          + f" -> {a.order_out}")
    if a.dst is None:
        return

    # Every read is held in memory once (signal decompressed by to_read()), then written in the new order.
    # Bacterial 30x sets are <= ~5 GB of samples.
    reads = {}
    with pod5.Reader(a.src) as r:
        for rec in r.reads():
            reads[str(rec.read_id)] = rec.to_read()
    if len(reads) != len(ids):
        sys.exit(f"{a.src}: read {len(reads)} records for {len(ids)} ids")
    if os.path.abspath(a.dst) == os.path.abspath(a.src):
        sys.exit("--dst must differ from SRC")
    if os.path.exists(a.dst):
        os.remove(a.dst)  # a previous attempt's output (never SRC, checked above)
    with pod5.Writer(a.dst) as w:
        for x in order:
            w.add_read(reads.pop(x))
    with pod5.Reader(a.dst) as r:
        got = [str(x) for x in r.read_ids]
    if got != order:
        sys.exit(f"{a.dst}: written read order differs from {a.order_out}")
    print(f"{a.dst}: {len(got)} reads written and re-read in the intended order")


if __name__ == "__main__":
    main()
