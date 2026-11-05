"""Readers shared by the evaluation harness (evaluation only; nothing here feeds the assembler).

Formats (all little-endian):
  evref.bin (truth/event_ref.py), per read:
      u32 name_len, name, u32 n, n x i64 ref, n x i32 base, n x i32 seg, n x u8 status
      ref = global reference coordinate (contigs in header order, each offset by the previous ones'
      length + 10 Mb), -1 unplaced; base = read base of the retained event from the move table
      (reference-free); status 0 UNPLACED 1 CONF 2 SINGLE 3 LOWQ 4 CONFLICT
  evbase.bin (a base map without placements), per read:
      u32 name_len, name, u32 n, n x i32 base
  spans.bin (haricot-truth-events), per read:
      u32 name_len, name, u32 n_samples, u32 n, n x u32 start, n x u32 end
  unitig_paths.tsv (haricot-detect --dump-paths): contig, step, anchor_read, ev_beg, ev_end, dir,
      advance_events, coord_events
"""
import gzip
import os
import struct

import numpy as np

UNPLACED, CONF, SINGLE, LOWQ, CONFLICT = 0, 1, 2, 3, 4


def read_evref(path, want=None, fields=("ref", "base", "seg", "status")):
    """name -> dict of the requested arrays; reads not in `want` are seeked past."""
    out = {}
    with open(path, "rb") as fh:
        size = os.fstat(fh.fileno()).st_size
        off = 0
        while off < size:
            fh.seek(off)
            (nl,) = struct.unpack("<I", fh.read(4))
            name = fh.read(nl).decode()
            (n,) = struct.unpack("<I", fh.read(4))
            base = off + 4 + nl + 4
            off = base + 8 * n + 4 * n + 4 * n + n
            if want is not None and name not in want:
                continue
            rec = {"n": n}
            if "ref" in fields:
                fh.seek(base)
                rec["ref"] = np.frombuffer(fh.read(8 * n), dtype="<i8")
            if "base" in fields:
                fh.seek(base + 8 * n)
                rec["base"] = np.frombuffer(fh.read(4 * n), dtype="<i4")
            if "seg" in fields:
                fh.seek(base + 12 * n)
                rec["seg"] = np.frombuffer(fh.read(4 * n), dtype="<i4")
            if "status" in fields:
                fh.seek(base + 16 * n)
                rec["status"] = np.frombuffer(fh.read(n), dtype=np.uint8)
            out[name] = rec
    return out


def read_evbase(path, want=None):
    """name -> base array, from an evbase.bin (reference-free) OR an evref.bin (its base field).
    The format is told by the extension."""
    if path.endswith(".evref.bin"):
        return {k: v["base"] for k, v in read_evref(path, want, fields=("base",)).items()}
    out = {}
    with open(path, "rb") as fh:
        size = os.fstat(fh.fileno()).st_size
        off = 0
        while off < size:
            fh.seek(off)
            (nl,) = struct.unpack("<I", fh.read(4))
            name = fh.read(nl).decode()
            (n,) = struct.unpack("<I", fh.read(4))
            base = off + 4 + nl + 4
            off = base + 4 * n
            if want is None or name in want:
                out[name] = np.frombuffer(fh.read(4 * n), dtype="<i4")
    return out


def read_spans(path):
    """name -> (n_samples, start u32 array, end u32 array)."""
    buf = open(path, "rb").read()
    off, out = 0, {}
    while off < len(buf):
        nl = int.from_bytes(buf[off:off + 4], "little")
        off += 4
        name = buf[off:off + nl].decode()
        off += nl
        ns, n = struct.unpack_from("<II", buf, off)
        off += 8
        st = np.frombuffer(buf, dtype="<u4", count=n, offset=off)
        off += 4 * n
        en = np.frombuffer(buf, dtype="<u4", count=n, offset=off)
        off += 4 * n
        out[name] = (ns, st, en)
    return out


def read_paths(outdir):
    """unitig -> list of (anchor_read, ev_beg, ev_end, dir, advance_events, coord_events)."""
    steps = {}
    with open(f"{outdir}/unitig_paths.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            steps.setdefault(f[0], []).append(
                (f[2], int(f[3]), int(f[4]), f[5], int(f[6]), int(f[7]) if len(f) > 7 else -1))
    return steps


def read_fastq(path, want):
    """name -> sequence for the wanted reads (gzip or plain FASTQ)."""
    op = gzip.open if path.endswith(".gz") else open
    seqs = {}
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            nm = h[1:].split()[0]
            if want is None or nm in want:
                seqs[nm] = s
    return seqs


def contig_table(fai, coord_gap=10_000_000):
    """[(name, global start, length)] in .fai / contigs.tsv order (event_ref.py's axis)."""
    out, cur = [], 0
    for line in open(fai):
        nm, ln = line.split("\t")[:2]
        out.append((nm, cur, int(ln)))
        cur += int(ln) + coord_gap
    return out
