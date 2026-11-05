#!/usr/bin/env python3
"""Write one read set's POD5 from the downloaded public files, and check it read by read.

  python3 build_set.py SET_DIR DOWNLOADS -o OUT.pod5 [--partial]

SET_DIR    a set folder of the read lists (read_lists/<role>/<set>/: read_ids.txt, checks.tsv, index.tsv)
DOWNLOADS  the folder download.py wrote: each read is read from DOWNLOADS/<source>/<file> as index.tsv names it
           (a .fast5 file as the .pod5 that download.py converted it to)

The reads are written in read_ids.txt order (it puts one read of every run first: the pipeline's POD5 reader needs
every run in the first record batch). Every read is checked against checks.tsv (samples, MD5 of the raw int16
signal, calibration offset and scale, sample rate), then OUT is read back and checked again: same ids in the same
order, same signal MD5s, every run in the first record batch. Exit status 1 on any missing read or mismatch.

--partial  skip the reads whose file is not downloaded (a test on part of a source); a read missing from a file
           that is present, or any mismatch, still fails.
Needs: pod5 (Python package).
"""
import argparse
import hashlib
import os
import sys
import uuid
import warnings
from collections import OrderedDict

import pod5

# pod5 0.3.x marks the scaling fields deprecated; they are still copied so each read keeps all its fields
warnings.filterwarnings("ignore", category=DeprecationWarning)


def load_set(d):
    order = [x.strip() for x in open(os.path.join(d, "read_ids.txt")) if x.strip()]
    checks, where = {}, {}
    with open(os.path.join(d, "checks.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if head != ["read_id", "samples", "signal_md5", "offset", "scale", "sample_rate"]:
            sys.exit(f"{d}/checks.tsv: unexpected header {head}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            checks[f[0]] = (int(f[1]), f[2], float(f[3]), float(f[4]), int(f[5]))
    with open(os.path.join(d, "index.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if head != ["read_id", "source", "file", "run", "samples"]:
            sys.exit(f"{d}/index.tsv: unexpected header {head}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            where[f[0]] = (f[1], f[2])
    if len(set(order)) != len(order):
        sys.exit(f"{d}/read_ids.txt: duplicate ids")
    if set(order) != set(checks) or set(order) != set(where):
        sys.exit(f"{d}: read_ids.txt, checks.tsv and index.tsv name different reads")
    return order, checks, where


def local_path(root, source, file):
    p = os.path.join(root, source, file)
    return p[:-len(".fast5")] + ".pod5" if p.endswith(".fast5") else p


def locate(path, wanted):
    """read id (16 bytes) -> (batch, row) for the wanted ids in one POD5"""
    loc = {}
    with pod5.Reader(path) as r:
        for bi in range(r.batch_count):
            col = r.read_table.get_batch(bi).column("read_id")
            buf = memoryview(col.buffers()[1])[col.offset * 16:(col.offset + len(col)) * 16]
            for row in range(len(col)):
                k = bytes(buf[row * 16:row * 16 + 16])
                if k in wanted:
                    loc[k] = (bi, row)
    return loc


class Readers:
    """open readers, at most MAX at a time (least recently used closed first)"""
    MAX = 64

    def __init__(self):
        self.open = OrderedDict()

    def get(self, path):
        if path in self.open:
            self.open.move_to_end(path)
        else:
            if len(self.open) >= self.MAX:
                self.open.popitem(last=False)[1].close()
            self.open[path] = pod5.Reader(path)
        return self.open[path]

    def close(self):
        for r in self.open.values():
            r.close()


def check(rec, sig, want):
    samples, md5, offset, scale, rate = want
    got = (len(sig), hashlib.md5(sig.tobytes()).hexdigest(), float(rec.calibration.offset),
           float(rec.calibration.scale), int(rec.run_info.sample_rate))
    return [f"{n} {g} != {w}" for n, g, w in zip(("samples", "signal_md5", "offset", "scale", "sample_rate"),
                                                  got, (samples, md5, offset, scale, rate)) if g != w]


def to_read(rec, sig):
    return pod5.Read(read_id=rec.read_id, pore=rec.pore, calibration=rec.calibration, median_before=rec.median_before,
                     end_reason=rec.end_reason, read_number=rec.read_number, run_info=rec.run_info,
                     start_sample=rec.start_sample, num_minknow_events=rec.num_minknow_events,
                     tracked_scaling=rec.tracked_scaling, predicted_scaling=rec.predicted_scaling,
                     num_reads_since_mux_change=rec.num_reads_since_mux_change,
                     time_since_mux_change=rec.time_since_mux_change, open_pore_level=rec.open_pore_level,
                     signal=sig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("set_dir")
    ap.add_argument("downloads")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--partial", action="store_true")
    a = ap.parse_args()

    order, checks, where = load_set(a.set_dir)
    name = os.path.basename(os.path.normpath(a.set_dir))
    by_path = {}
    for rid in order:
        by_path.setdefault(local_path(a.downloads, *where[rid]), []).append(rid)
    absent = [p for p in by_path if not os.path.isfile(p)]
    missing = [r for p in absent for r in by_path[p]]
    print(f"{name}: {len(order)} reads in {len(by_path)} files; {len(absent)} files not downloaded "
          f"({len(missing)} reads)", file=sys.stderr)
    if absent and not a.partial:
        print(f"{name}: not downloaded, e.g. {absent[:3]}", file=sys.stderr)
        sys.exit(1)

    loc, lost = {}, []
    for p, rids in by_path.items():
        if p in absent:
            continue
        found = locate(p, {uuid.UUID(r).bytes for r in rids})
        for r in rids:
            k = uuid.UUID(r).bytes
            if k in found:
                loc[r] = (p,) + found[k]
            else:
                lost.append(f"{r} ({p})")
    if lost:
        print(f"{name}: {len(lost)} reads are not in the file index.tsv names, e.g.", *lost[:5], sep="\n  ",
              file=sys.stderr)
        sys.exit(1)

    readers, bad, written = Readers(), [], []
    tmp = a.out + ".tmp.pod5"
    if os.path.exists(tmp):
        os.remove(tmp)
    with pod5.Writer(tmp) as w:
        for rid in order:
            if rid not in loc:
                continue
            p, bi, row = loc[rid]
            rec = readers.get(p).get_batch(bi).get_read(row)
            sig = rec.signal
            errs = check(rec, sig, checks[rid])
            if errs:
                bad.append(f"{rid}: {'; '.join(errs)} ({p})")
            else:
                w.add_read(to_read(rec, sig))
                written.append(rid)
    readers.close()
    if bad:
        print(f"{name}: {len(bad)} reads differ from checks.tsv:", *bad[:10], sep="\n  ", file=sys.stderr)
        os.remove(tmp)
        sys.exit(1)

    # read the output back
    got, runs, first = [], set(), set()
    with pod5.Reader(tmp) as r:
        for bi in range(r.batch_count):
            for rec in r.get_batch(bi).reads():
                rid = str(rec.read_id)
                got.append(rid)
                runs.add(rec.run_info.acquisition_id)
                if bi == 0:
                    first.add(rec.run_info.acquisition_id)
                if check(rec, rec.signal, checks[rid]):
                    bad.append(rid)
    errs = []
    if got != written:
        errs.append("the output's read order differs from read_ids.txt")
    if bad:
        errs.append(f"{len(bad)} reads of the output differ from checks.tsv")
    if runs != first and not a.partial:
        errs.append(f"{len(runs - first)} runs have no read in the first record batch")
    if errs:
        print(f"{name}: " + "; ".join(errs), file=sys.stderr)
        os.remove(tmp)
        sys.exit(1)
    os.replace(tmp, a.out)
    note = f"; {len(missing)} not downloaded (--partial)" if missing else ""
    print(f"{name}: wrote {len(written)} reads to {a.out}, every one matching checks.tsv, in read_ids.txt order, "
          f"{len(runs)} runs ({len(first)} in the first batch){note}", file=sys.stderr)


if __name__ == "__main__":
    main()
