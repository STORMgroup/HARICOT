#!/usr/bin/env python3
"""A depth subset of a POD5 by a salted read-id hash:
reads with md5("depth:" + id)[:8] / 0xFFFFFFFF < COV / BASE. Order: the first kept read of every run first (the pipeline's
POD5 reader needs every run in the first record batch), then the other kept reads in the
source's order.
Usage: python -I subset.py SRC.pod5 DST.pod5 COV BASE
"""
import hashlib
import os
import sys

import pod5


def main():
    src, dst, cov, base = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
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
