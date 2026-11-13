#!/usr/bin/env python3
"""A depth subset of a POD5, chosen by a salted read-id hash.

Keeps the reads with md5("depth:" + id)[:8] / 0xFFFFFFFF < COV / BASE. Order: the first kept read of every run first
(the pipeline's POD5 reader needs every run in the first record batch), then the other kept reads in the source's
order. truth/fragments.sh runs it on the fragment sets; datasets/README.md explains its use for the depth subsets.

Usage: python3 -I subset.py SRC.pod5 DST.pod5 COV BASE
"""
import argparse
import hashlib
import os
import sys

import pod5


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", help="the read set")
    ap.add_argument("dst", help="the subset's POD5 (written)")
    ap.add_argument("cov", type=float, help="the subset's depth")
    ap.add_argument("base", type=float, help="the read set's depth")
    a = ap.parse_args()
    src, dst, cov, base = a.src, a.dst, a.cov, a.base
    frac = cov / base
    kept, lead, n = [], {}, 0
    with pod5.Reader(src) as r:
        for rec in r.reads():
            n += 1
            rid = str(rec.read_id)
            if int(hashlib.md5(("depth:" + rid).encode()).hexdigest()[:8], 16) / 0xFFFFFFFF < frac:
                kept.append(rid)
                lead.setdefault(rec.run_info.acquisition_id, rid)
    leads = set(lead.values())
    tmp = dst + ".tmp.pod5"
    if os.path.exists(tmp):
        os.remove(tmp)
    k = 0
    with pod5.Reader(src) as r, pod5.Writer(tmp) as w:
        for rec in r.reads(selection=list(leads)):
            w.add_read(rec.to_read())
            k += 1
        rest = [x for x in kept if x not in leads]
        if rest:
            for rec in r.reads(selection=rest):
                w.add_read(rec.to_read())
                k += 1
    if k != len(kept):
        sys.exit(f"wrote {k} reads, expected {len(kept)}")
    os.replace(tmp, dst)
    print(f"{dst}: kept {k} of {n} reads ({frac:.3f}), {len(leads)} runs first", file=sys.stderr)


if __name__ == "__main__":
    main()
