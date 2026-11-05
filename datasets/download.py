#!/usr/bin/env python3
"""Download the public raw-signal files of HARICOT's read sets and check them against the archives' checksums.

  python3 download.py OUT --set SET_DIR [--set SET_DIR ...]     the files these read sets need
  python3 download.py OUT --source SOURCE [--files FILE ...]    a whole source, or some of its files
  --list   print the files and their size, download nothing

SET_DIR  a set folder of the read lists (read_lists/<role>/<set>/); its index.tsv names each read's source and file
SOURCE   a source of sources.tsv; FILE a file of that source in files.tsv (path or base name)
OUT      files land in OUT/<source>/<file>; give OUT to build_set.py

By source type (sources.tsv type):
  tar  (figshare)  whole source: the archive is streamed, its published MD5 checked, every member extracted.
                   Some files (--set, --files): each member is fetched by HTTP range at its offset in files.tsv;
                   no per-member checksum is published (build_set.py checks every read).
  tgz  (ENA)       the archive is streamed to its end, its published MD5 checked, the wanted members extracted.
                   With --files the stream stops after the last wanted member and the MD5 is not checked.
  s3   (AWS)       one object per file; size and S3 ETag (MD5, or MD5 of the part MD5s) checked against files.tsv.
fast5 sources (sources.tsv format) are converted to POD5 with `pod5 convert fast5`; the fast5 is then removed.
Files present already are kept (a partial download carries a .part suffix and is resumed or redone).
Needs: Python 3.8+, curl; pod5 (command, from the pod5 Python package) for fast5 sources.
"""
import argparse
import csv
import hashlib
import math
import os
import shutil
import subprocess
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
MIB = 1 << 20


def table(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def safe_join(root, name):
    p = os.path.normpath(os.path.join(root, name))
    if not p.startswith(os.path.normpath(root) + os.sep):
        sys.exit(f"unsafe path {name}")
    return p


def curl(url, out, rng=None):
    cmd = ["curl", "-sS", "-L", "--fail", "--retry", "8", "--retry-delay", "10", "--retry-all-errors",
           "--connect-timeout", "60", "-o", out]
    cmd += ["-r", rng] if rng else ["-C", "-"]
    subprocess.run(cmd + [url], check=True)


def md5_file(path, block=8 * MIB):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def s3_etag_ok(path, etag):
    """an S3 ETag is the object's MD5, or for a multipart upload the MD5 of the parts' MD5s and '-N';
    the part size is not published: try whole MiB sizes that give N parts (8 and 16 MiB first)"""
    if "-" not in etag:
        return md5_file(path) == etag
    size, n = os.path.getsize(path), int(etag.split("-")[1])
    lo = max(1, math.ceil(size / n / MIB))
    sizes = sorted((p for p in range(lo, lo + 64) if math.ceil(size / (p * MIB)) == n),
                   key=lambda p: (p not in (8, 16), p))
    for p in sizes:
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for part in iter(lambda: fh.read(p * MIB), b""):
                h.update(hashlib.md5(part).digest())
        if f"{h.hexdigest()}-{n}" == etag:
            return True
    return False


def fetch_s3(src, row, root):
    dst = safe_join(root, row["file"])
    fast5 = src["format"] == "fast5"
    final = dst[:-len(".fast5")] + ".pod5" if fast5 else dst
    if os.path.exists(final):
        return "present"
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if not os.path.exists(dst):
        curl(src["url"] + row["file"], dst + ".part")
        got = os.path.getsize(dst + ".part")
        if got != int(row["bytes"]):
            sys.exit(f"{row['file']}: {got} bytes, expected {row['bytes']}")
        if not s3_etag_ok(dst + ".part", row["etag"]):
            sys.exit(f"{row['file']}: S3 ETag mismatch (expected {row['etag']})")
        os.replace(dst + ".part", dst)
    if not fast5:
        return "size and ETag ok"
    if os.path.exists(final + ".part.pod5"):
        os.remove(final + ".part.pod5")
    subprocess.run(["pod5", "convert", "fast5", "--strict", "-t", "2", "-o", final + ".part.pod5", dst],
                   check=True, stdout=subprocess.DEVNULL)
    os.replace(final + ".part.pod5", final)
    os.remove(dst)
    return "size and ETag ok, converted to POD5"


def fetch_member(src, row, root):
    """one member of a plain tar by HTTP range; the tar header before it must name it"""
    dst = safe_join(root, row["file"])
    if os.path.exists(dst):
        return "present"
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    off, n = int(row["offset"]), int(row["bytes"])
    # the download link redirects to a storage URL that expires within a minute: curl follows it on every request
    curl(src["url"], dst + ".hdr", f"{off - 512}-{off - 1}")
    with open(dst + ".hdr", "rb") as fh:
        h = fh.read()
    os.remove(dst + ".hdr")
    name, prefix = h[:100].rstrip(b"\0").decode(), h[345:500].rstrip(b"\0").decode()
    name = prefix + "/" + name if prefix else name
    if len(h) != 512 or not row["file"].startswith(name) or h[156:157] not in (b"0", b"\0"):
        sys.exit(f"{row['file']}: the tar header at byte {off - 512} names {name!r}")
    curl(src["url"], dst + ".part", f"{off}-{off + n - 1}")
    if os.path.getsize(dst + ".part") != n:
        sys.exit(f"{row['file']}: {os.path.getsize(dst + '.part')} bytes, expected {n}")
    os.replace(dst + ".part", dst)
    return "by range"


class HashingReader:
    def __init__(self, fh):
        self.fh, self.md5, self.n = fh, hashlib.md5(), 0

    def read(self, size=-1):
        b = self.fh.read(size)
        self.md5.update(b)
        self.n += len(b)
        return b


def stream_archive(src, root, wanted, stop_early):
    """stream a tar or tar.gz; extract the members whose base name is wanted (all if None); check the archive MD5"""
    proc = subprocess.Popen(["curl", "-sS", "-L", "--fail", "--retry", "8", "--retry-all-errors", src["url"]],
                            stdout=subprocess.PIPE)
    hr = HashingReader(proc.stdout)
    left = None if wanted is None else set(wanted)
    made, stopped = [], False
    with tarfile.open(fileobj=hr, mode="r|gz" if src["type"] == "tgz" else "r|") as tf:
        for m in tf:
            if not m.isfile() or (left is not None and m.name not in left):
                continue
            dst = safe_join(root, m.name)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with tf.extractfile(m) as fi, open(dst + ".part", "wb") as fo:
                shutil.copyfileobj(fi, fo, 8 * MIB)
            os.replace(dst + ".part", dst)
            made.append(dst)
            print(f"  {m.name} ({m.size} bytes)", file=sys.stderr, flush=True)
            if left is not None:
                left.discard(m.name)
                if not left and stop_early:
                    stopped = True
                    break
    if stopped:
        proc.kill()
        proc.wait()
        print(f"{src['source']}: stopped after the wanted members; archive MD5 not checked", file=sys.stderr)
    else:
        while hr.read(8 * MIB):
            pass
        if proc.wait() != 0:
            sys.exit(f"{src['source']}: download failed (curl exit {proc.returncode})")
        if hr.md5.hexdigest() != src["md5"]:
            for p in made:
                os.remove(p)
            sys.exit(f"{src['source']}: archive MD5 {hr.md5.hexdigest()} != published {src['md5']}; "
                     "extracted files removed")
        print(f"{src['source']}: archive MD5 {hr.md5.hexdigest()} = published ({hr.n} bytes)", file=sys.stderr)
    if left:
        sys.exit(f"{src['source']}: {len(left)} wanted members not in the archive, e.g. {sorted(left)[:3]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out")
    ap.add_argument("--set", action="append", default=[], dest="sets")
    ap.add_argument("--source")
    ap.add_argument("--files", nargs="+")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if bool(a.sets) == bool(a.source) or (a.files and not a.source):
        ap.error("give --set SET_DIR ... or --source SOURCE [--files FILE ...]")
    sources = {r["source"]: r for r in table(os.path.join(HERE, "sources.tsv"))}
    files = table(os.path.join(HERE, "files.tsv"))

    wanted = {}  # source -> rows of files.tsv
    if a.source:
        if a.source not in sources:
            sys.exit(f"unknown source {a.source}; sources: {' '.join(sources)}")
        rows = [r for r in files if r["source"] == a.source]
        if a.files:
            names = set(a.files)
            rows = [r for r in rows if r["file"] in names or os.path.basename(r["file"]) in names]
            if len(rows) != len(names):
                sys.exit(f"{a.source}: {len(names) - len(rows)} of the --files are not in files.tsv")
        wanted[a.source] = rows
    else:
        row_of = {(r["source"], r["file"]): r for r in files}
        for d in a.sets:
            with open(os.path.join(d, "index.tsv")) as fh:
                if fh.readline().rstrip("\n").split("\t")[:3] != ["read_id", "source", "file"]:
                    sys.exit(f"{d}/index.tsv: not a public index (read_id source file run samples)")
                need = {tuple(line.split("\t")[1:3]) for line in fh}
            miss = sorted(x for x in need if x not in row_of)
            if miss:
                sys.exit(f"{d}: {len(miss)} files of index.tsv are not in files.tsv, e.g. {miss[:3]}")
            for x in need:
                wanted.setdefault(x[0], {})[x[1]] = row_of[x]
        wanted = {s: sorted(v.values(), key=lambda r: r["file"]) for s, v in wanted.items()}

    total = 0
    for s, rows in wanted.items():
        b = sum(int(r["bytes"]) for r in rows)
        total += b
        print(f"{s}: {len(rows)} files, {b / 1e9:.1f} GB", file=sys.stderr)
        if a.list:
            for r in rows:
                print(s, r["file"], r["bytes"], sep="\t")
    print(f"total {total / 1e9:.1f} GB", file=sys.stderr)
    if a.list:
        return

    for s, rows in wanted.items():
        src, root = sources[s], os.path.join(a.out, s)
        os.makedirs(root, exist_ok=True)
        whole = a.source is not None and not a.files
        if src["type"] == "tgz" or (src["type"] == "tar" and whole):
            stream_archive(src, root, None if whole else [r["file"] for r in rows], stop_early=bool(a.files))
        elif src["type"] == "tar":
            for r in rows:
                print(f"  {r['file']}: {fetch_member(src, r, root)}", file=sys.stderr, flush=True)
        elif src["type"] == "s3":
            for r in rows:
                print(f"  {r['file']}: {fetch_s3(src, r, root)}", file=sys.stderr, flush=True)
        else:
            sys.exit(f"{s}: unknown source type {src['type']}")
        print(f"{s}: done -> {root}", file=sys.stderr)


if __name__ == "__main__":
    main()
