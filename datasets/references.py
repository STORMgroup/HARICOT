#!/usr/bin/env python3
"""Download the manuscript's reference genomes and check every sequence.

  python3 references.py OUT [--files FILE ...]

Writes OUT/<file> for each file of references.tsv (cj.fa, ..., athal.fa, dmel.fa, hg002v1.2_matY.fa, chr18_PAT.fa):
the sequences of reference_sequences.tsv, in its order, renamed as it says, 80 bases per line. Each sequence's length
and MD5 (of the bases as published, case kept) must equal the table's; otherwise that file is not written and the exit
status is 1. Sources are fetched once into OUT/.sources/ (HG002 v1.2: 1.8 GB; its sequences are held in memory, ~4 GB).
Needs: Python 3.8+, curl.
"""
import argparse
import csv
import gzip
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def table(p):
    with open(os.path.join(HERE, p), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fasta(path, want):
    """(name, bases) of the records whose header's first word is in want"""
    op = gzip.open if path.endswith(".gz") else open
    name, parts = None, []
    with op(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name in want:
                    yield name, b"".join(parts)
                name, parts = line[1:].split()[0].decode(), []
            elif name in want:
                parts.append(line.strip())
    if name in want:
        yield name, b"".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out")
    ap.add_argument("--files", nargs="+")
    a = ap.parse_args()
    url_of = {r["file"]: r["url"] for r in table("references.tsv")}
    seqs = table("reference_sequences.tsv")
    files = a.files or list(url_of)
    unknown = [f for f in files if f not in url_of]
    if unknown:
        sys.exit(f"unknown --files {unknown}; files: {' '.join(url_of)}")
    src_dir = os.path.join(a.out, ".sources")
    os.makedirs(src_dir, exist_ok=True)

    got = {}  # (url, source name) -> bases
    for url in dict.fromkeys(url_of[f] for f in files):
        local = os.path.join(src_dir, hashlib.md5(url.encode()).hexdigest()[:8] + "_" +
                             os.path.basename(url.split("?")[0]).replace("efetch.fcgi", "efetch.fa"))
        if not os.path.exists(local):
            subprocess.run(["curl", "-sS", "-L", "--fail", "--retry", "5", "-o", local + ".part", url], check=True)
            os.replace(local + ".part", local)
        want = {r["source_name"] for r in seqs if url_of[r["file"]] == url and r["file"] in files}
        for name, s in fasta(local, want):
            got[(url, name)] = s
        print(f"{os.path.basename(local)}: {sum(1 for u, n in got if u == url)} of {len(want)} sequences", file=sys.stderr)

    bad = []
    for f in files:
        out, errs = [], []
        for r in (r for r in seqs if r["file"] == f):
            s = got.get((url_of[f], r["source_name"]))
            if s is None:
                errs.append(f"{r['source_name']} not in {url_of[f]}")
                continue
            md5 = hashlib.md5(s).hexdigest()
            if len(s) != int(r["length"]) or md5 != r["md5"]:
                errs.append(f"{r['name']}: length {len(s)}, MD5 {md5}; expected {r['length']}, {r['md5']}")
            out.append((r["name"], s))
        if errs:
            bad += [f"{f}: {e}" for e in errs]
            continue
        p = os.path.join(a.out, f)
        with open(p + ".part", "wb") as w:
            for name, s in out:
                w.write(b">" + name.encode() + b"\n")
                for i in range(0, len(s), 80):
                    w.write(s[i:i + 80] + b"\n")
        os.replace(p + ".part", p)
        print(f"{f}: {len(out)} sequences, {sum(len(s) for _, s in out)} bases, each matching its length and MD5",
              file=sys.stderr)
    if bad:
        print("mismatches:", *bad, sep="\n  ", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
