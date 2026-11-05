"""Shared by the evaluation tests (tests/eval/test_*.py): repository paths, running a repository script, and the
writer of synthetic truth files. Run the suite from anywhere (Python 3 with numpy):

    python -m unittest discover -s tests/eval

Every test writes its own small inputs (tables, a binary truth file, reads) into a temporary directory.
"""
import os
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CONF, SINGLE, LOWQ, CONFLICT = 1, 2, 3, 4  # truth/event_ref.py's event status codes

sys.path.insert(0, os.path.join(ROOT, "evaluation"))


def run(rel, *args, cwd=None, env=None, check=True):
    """Run the repository script `rel` with this interpreter, isolated (-I, as the scripts run them)."""
    r = subprocess.run([sys.executable, "-I", os.path.join(ROOT, rel), *map(str, args)], cwd=cwd, env=env,
                       capture_output=True, text=True)
    if check and r.returncode != 0:
        raise AssertionError(f"{rel} exited {r.returncode}\n--- stdout\n{r.stdout}\n--- stderr\n{r.stderr}")
    return r


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def read(path):
    with open(path) as fh:
        return fh.read()


def one_row(path):
    """A TSV of a header and one row (a scorer's --tsv on a new file) -> dict."""
    hdr, val = read(path).rstrip("\n").split("\n")
    return dict(zip(hdr.split("\t"), val.split("\t")))


def tsv_rows(path):
    lines = read(path).rstrip("\n").split("\n")
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, x.split("\t"))) for x in lines[1:]]


def write_evref(path, reads):
    """reads: {name: [(global reference position or -1, status), ...]} -> an evref.bin as truth/event_ref.py writes
    it (per read: u32 name_len, name, u32 n, n x i64 ref, n x i32 base, n x i32 seg, n x u8 status); base = the
    event index, seg 0."""
    with open(path, "wb") as fh:
        for name, ev in reads.items():
            nm, n = name.encode(), len(ev)
            fh.write(struct.pack("<I", len(nm)) + nm + struct.pack("<I", n))
            fh.write(struct.pack(f"<{n}q", *[r for r, _ in ev]))
            fh.write(struct.pack(f"<{n}i", *range(n)))
            fh.write(struct.pack(f"<{n}i", *([0] * n)))
            fh.write(bytes(s for _, s in ev))
