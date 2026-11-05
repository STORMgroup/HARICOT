#!/usr/bin/env python3
"""Metamorphic inputs: a NEW pod5 with SRC's
reads changed in one way, and the map from its read ids back to SRC's. SRC is only read.

  --mode unc   the same reads in the same order, the signal stored uncompressed (another file layout)
  --mode reno  every read id replaced by a new random UUID, the ORDER of the ids preserved (the sorted new ids go to
               the sorted old ids): every read-id rank, so every canonical tie-break, is unchanged
  --mode renr  every read id replaced by a new random UUID (seeded): the read-id ranks reshuffled; file order kept
  --mode dupn  SRC, then one more copy of one read (the median-length read; ties by id) under a new id
  --mode dups  SRC, then one more copy of that read under the SAME id (an invalid input)

The order of the reads in DST is SRC's (dup modes: the copy last), so every read keeps its number (reads are numbered
in file order).
Writes DST and MAP (`dst_id<TAB>src_id`, one line per read of DST, in DST's order) and prints the copied read.

Usage: pod5_variant.py SRC --mode unc|reno|renr|dupn|dups [--seed N] --dst DST --map MAP
"""
import argparse
import copy
import dataclasses
import os
import random
import sys
import uuid

import pod5


def new_ids(rng, n):
    """n distinct random version-4 UUIDs."""
    out = set()
    while len(out) < n:
        out.add(uuid.UUID(int=rng.getrandbits(128), version=4))
    return list(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--mode", required=True, choices=["unc", "reno", "renr", "dupn", "dups"])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dst", required=True)
    ap.add_argument("--map", required=True)
    a = ap.parse_args()
    if os.path.abspath(a.dst) == os.path.abspath(a.src):
        sys.exit("--dst must differ from SRC")
    if os.path.exists(a.dst):
        os.remove(a.dst)  # a previous attempt's output (never SRC, checked above)
    rng = random.Random(a.seed)

    with pod5.Reader(a.src) as r:
        meta = [(str(rec.read_id), rec.num_samples) for rec in r.reads()]
    ids = [m[0] for m in meta]
    if len(set(ids)) != len(ids):
        sys.exit(f"{a.src}: duplicate read ids")
    rename = {x: x for x in ids}
    if a.mode == "reno":
        fresh = sorted(str(u) for u in new_ids(rng, len(ids)))
        rename = dict(zip(sorted(ids), fresh))
    elif a.mode == "renr":
        fresh = [str(u) for u in new_ids(rng, len(ids))]
        rename = dict(zip(ids, fresh))
    dup = None
    if a.mode in ("dupn", "dups"):
        by_len = sorted(meta, key=lambda m: (m[1], m[0]))
        dup = by_len[len(by_len) // 2][0]
    dup_new = str(new_ids(rng, 1)[0]) if a.mode == "dupn" else dup

    kind = pod5.SignalType.UncompressedSignal if a.mode == "unc" else pod5.SignalType.VbzSignal
    rows, held = [], None
    with pod5.Reader(a.src) as r, pod5.Writer(a.dst, signal_compression_type=kind) as w:
        for rec in r.reads():
            rd = rec.to_read()
            src_id = str(rd.read_id)
            if src_id == dup:
                held = copy.deepcopy(rd)
            rd = dataclasses.replace(rd, read_id=uuid.UUID(rename[src_id]))
            w.add_read(rd)
            rows.append((rename[src_id], src_id))
        if dup is not None:
            w.add_read(dataclasses.replace(held, read_id=uuid.UUID(dup_new)))
            rows.append((dup_new, dup))
    with open(a.map, "w") as f:
        f.write("".join(f"{x}\t{y}\n" for x, y in rows))
    with pod5.Reader(a.dst) as r:
        got = [str(x) for x in r.read_ids]
    if got != [x for x, _ in rows]:
        sys.exit(f"{a.dst}: written read ids differ from {a.map}")
    print(f"{a.src}: {len(ids)} reads -> {a.dst} ({len(got)} reads, mode {a.mode}, seed {a.seed})"
          + (f"; copied read {dup} ({dict(meta)[dup]} samples) as {dup_new}" if dup else ""))


if __name__ == "__main__":
    main()
