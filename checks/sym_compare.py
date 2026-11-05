#!/usr/bin/env python3
"""The invariance comparator (invariance.sh): two haricot assemble --keep-work runs on the
same reads in different pod5 orders (the reads numbered in each pod5's order) or under other read names (ids mapped
back), every output compared as a function of the read set alone.

  final / builds  the 8 graph and layout files (fbg.gfa fbg_nodes.tsv unitigs.gfa unitig_nodes.tsv unitig_paths.tsv
                  contigs.gfa contig_paths.tsv contig_nodes.tsv) of OUT/ and of the pol7 / eb builds: byte for byte
  canon           every candidate cache's canonical dump (`haricot-ovl-filter --dump-canonical`, written by invariance.sh into
                  work/canon/): rows keyed by the two read names in the pair's canonical frame, every column exact
                  except `slope`, which may differ by a relative 1e-12 (a stored slope is the inverse of the frame's
                  when the owner is not the query, so its canonical view is 1 / (1 / s))
  table           name-keyed TSVs with a header (stage-1 features, stage-1 / stage-2 scores, stage-2 features, the
                  rescue's features and scores): rows keyed by the unordered name pair, every other column exact as
                  printed (the features are computed in the canonical frame, so they do not depend on the owner)
  name_pairs      drop lists and the end-bridge list (name pairs): equal as sets of unordered pairs
  idx_pairs       pool.tsv, acc.tsv, admit.tsv (read indices): mapped through each run's read_names.tsv, equal as sets
  threads         the read-placement dumps (work/thr.tsv, work/thr_eb.tsv, OUT/threads.tsv): read numbers mapped
                  to names, equal as multisets of rows
  pool_detail     q9c/pool_detail.tsv: the header lines exact; rows with the two reads in name order (a row read the
                  other way round has its shape mirrored and its end pairs reversed), end pairs as a set; equal as sets
  dump            cand.dump.tsv (`haricot-ovl-filter --dump` of the candidates, owner first): rows in the canonical frame (query =
                  the longer read, equal lengths the smaller name), every column but `slope` exact (the canonical cache
                  dump of the same cache carries the slope)
  names           read_names.tsv: the same set of names (the numbering itself differs by construction)
  bytes           everything else (depth.txt, ...): byte for byte
Not compared: *.log, *.stdout, time.tsv (timings) and the *.bin caches themselves (compared through their dumps).

Prints one line per item that is not SAME and a summary line; writes every item to --out. Exit 0 when every item is
SAME, 1 otherwise. Reads only.

Usage: sym_compare.py A B --label L [--out OUT.tsv]
"""
import argparse
import filecmp
import math
import os
import sys
from collections import Counter

GRAPH = set("fbg.gfa fbg_nodes.tsv unitigs.gfa unitig_nodes.tsv unitig_paths.tsv contigs.gfa contig_paths.tsv "
            "contig_nodes.tsv".split())
SLOPE_REL = 1e-12
MIRROR = {"a_in_b": "b_in_a", "b_in_a": "a_in_b", "dove_ab": "dove_ba", "dove_ba": "dove_ab"}


def names_of(run):
    out = {}
    with open(os.path.join(run, "read_names.tsv")) as f:
        for line in f:
            i, n = line.rstrip("\n").split("\t")
            out[int(i)] = n
    return out


def lines(p):
    with open(p) as f:
        return [x.rstrip("\n") for x in f if x.strip()]


def read_tsv(p):
    ls = lines(p)
    if not ls:
        return [], []
    return ls[0].split("\t"), [x.split("\t") for x in ls[1:]]


def upair(x, y):
    return (x, y) if x < y else (y, x)


def setdiff(a, b, what="rows"):
    """SAME or DIFF with counts, for two collections compared as sets / multisets."""
    ca, cb = Counter(a), Counter(b)
    if ca == cb:
        return "SAME", f"{sum(ca.values())} {what}"
    oa, ob = ca - cb, cb - ca
    ex = next(iter(oa or ob))
    return "DIFF", (f"{sum(ca.values())} vs {sum(cb.values())} {what}; {sum(oa.values())} only in A, "
                    f"{sum(ob.values())} only in B; e.g. {str(ex)[:200]}")


def r_bytes(pa, pb, na, nb):
    return ("SAME", "bytes") if filecmp.cmp(pa, pb, shallow=False) else ("DIFF", "bytes differ")


def r_names(pa, pb, na, nb):
    return setdiff(na.values(), nb.values(), "names")


def r_name_pairs(pa, pb, na, nb):
    def get(p):
        return [upair(*x.split("\t")[:2]) for x in lines(p)]
    return setdiff(get(pa), get(pb), "pairs")


def r_idx_pairs(pa, pb, na, nb):
    def get(p, n):
        out = []
        for x in lines(p):
            f = x.split()
            out.append(upair(n[int(f[0])], n[int(f[1])]))
        return out
    return setdiff(get(pa, na), get(pb, nb), "pairs")


def r_threads(pa, pb, na, nb):
    def get(p, n):
        out = []
        for x in lines(p):
            f = x.split("\t")
            if f[0] in ("T", "C"):
                f[1] = n[int(f[1])]
            out.append(tuple(f))
        return out
    return setdiff(get(pa, na), get(pb, nb), "rows")


def r_pool_detail(pa, pb, na, nb):
    def get(p, n):
        ls = lines(p)
        head = [x for x in ls if x.startswith("#") or x.startswith("a\t")]
        rows = []
        for x in ls:
            if x.startswith("#") or x.startswith("a\t"):
                continue
            a, b, sc, shape, ends = x.split("\t")
            a, b = n[int(a)], n[int(b)]
            eps = [e for e in ends.split(";") if e]
            if a > b:
                a, b, shape = b, a, MIRROR.get(shape, shape)
                eps = ["-".join(reversed(e.split("-"))) for e in eps]
            rows.append((a, b, sc, shape, tuple(sorted(eps))))
        return head, rows
    ha, ra = get(pa, na)
    hb, rb = get(pb, nb)
    if ha != hb:
        return "DIFF", f"header lines differ: {ha} vs {hb}"
    return setdiff(ra, rb, "pool rows")


def r_dump(pa, pb, na, nb):
    def get(p):
        hdr, rows = read_tsv(p)
        out = []
        for r in rows:
            v = dict(zip(hdr, r))
            a_ev, b_ev = int(v["a_ev"]), int(v["b_ev"])
            q_first = a_ev > b_ev or (a_ev == b_ev and v["a"] < v["b"])
            if not q_first:
                for x, y in (("a", "b"), ("a_lo", "b_lo"), ("a_hi", "b_hi"), ("a_ev", "b_ev")):
                    v[x], v[y] = v[y], v[x]
            out.append(tuple(v[c] for c in hdr if c != "slope"))
        return out
    return setdiff(get(pa), get(pb), "overlaps")


def r_table(pa, pb, na, nb):
    ha, ra = read_tsv(pa)
    hb, rb = read_tsv(pb)
    if ha != hb:
        return "DIFF", "headers differ"
    ia, ib = ha.index("name_a"), ha.index("name_b")
    rest = [j for j in range(len(ha)) if j not in (ia, ib)]

    def keyed(rows):
        d = {}
        for r in rows:
            d.setdefault(upair(r[ia], r[ib]), []).append(tuple(r[j] for j in rest))
        return d
    A, B = keyed(ra), keyed(rb)
    oa, ob = set(A) - set(B), set(B) - set(A)
    diff = [k for k in set(A) & set(B) if sorted(A[k]) != sorted(B[k])]
    if not (oa or ob or diff):
        return "SAME", f"{len(ra)} rows"
    cols = Counter()
    for k in diff:
        x, y = sorted(A[k])[0], sorted(B[k])[0]
        cols.update(ha[rest[j]] for j in range(len(rest)) if x[j] != y[j])
    return "DIFF", (f"{len(ra)} vs {len(rb)} rows; {len(oa)} pairs only in A, {len(ob)} only in B, {len(diff)} differ"
                    f" (columns: {', '.join(f'{c} {n}' for c, n in cols.most_common(8))})")


def r_canon(pa, pb, na, nb):
    ha, ra = read_tsv(pa)
    hb, rb = read_tsv(pb)
    if ha != hb:
        return "DIFF", "headers differ"
    js = ha.index("slope")
    A = {(r[0], r[1]): r for r in ra}
    B = {(r[0], r[1]): r for r in rb}
    if len(A) != len(ra) or len(B) != len(rb):
        return "DIFF", "a pair appears twice"
    oa, ob = set(A) - set(B), set(B) - set(A)
    bad, worst, slope_off = Counter(), 0.0, 0
    for k in set(A) & set(B):
        x, y = A[k], B[k]
        for j in range(len(ha)):
            if j == js:
                sx, sy = float(x[j]), float(y[j])
                if sx != sy:
                    rel = abs(sx - sy) / max(abs(sx), abs(sy)) if not (math.isnan(sx) or math.isnan(sy)) else math.inf
                    worst = max(worst, rel)
                    if rel > SLOPE_REL:
                        bad["slope"] += 1
                    else:
                        slope_off += 1
            elif x[j] != y[j]:
                bad[ha[j]] += 1
    tail = f"; slope differs within {SLOPE_REL:g} on {slope_off} (max rel {worst:.2e})"
    if not (oa or ob or bad):
        return "SAME", f"{len(ra)} overlaps{tail}"
    return "DIFF", (f"{len(ra)} vs {len(rb)} overlaps; {len(oa)} only in A, {len(ob)} only in B; differing columns: "
                    f"{', '.join(f'{c} {n}' for c, n in bad.most_common(8))}{tail}")


def has_name_header(p):
    with open(p) as f:
        h = f.readline().rstrip("\n").split("\t")
    return "name_a" in h and "name_b" in h


def rule_of(rel, pa):
    base = os.path.basename(rel)
    if base.endswith(".log") or base.endswith(".stdout") or base == "time.tsv" or base.endswith(".bin"):
        return None, None
    if base == "read_names.tsv":
        return "names", r_names
    if base in GRAPH:
        return ("final" if "/" not in rel else "builds"), r_bytes
    if rel.startswith("work/canon/"):
        return "canon", r_canon
    if base in ("thr.tsv", "thr_eb.tsv", "threads.tsv"):
        return "threads", r_threads
    if rel == "work/q9c/pool_detail.tsv":
        return "pool_detail", r_pool_detail
    if rel in ("work/q9c/pool.tsv", "work/acc.tsv", "work/admit.tsv"):
        return "idx_pairs", r_idx_pairs
    if base.endswith(".drop.tsv"):
        return "name_pairs", r_name_pairs
    if rel == "work/cand.dump.tsv":
        return "dump", r_dump
    if base.endswith(".tsv") and os.path.getsize(pa) > 0 and has_name_header(pa):
        return "table", r_table
    return "bytes", r_bytes


def files(root):
    out = set()
    for d, _, fs in os.walk(root):
        for f in fs:
            out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    na, nb = names_of(a.a), names_of(a.b)
    fa, fb = files(a.a), files(a.b)
    rows = []
    for rel in sorted(fa | fb):
        pa, pb = os.path.join(a.a, rel), os.path.join(a.b, rel)
        p0 = pa if rel in fa else pb
        kind, fn = rule_of(rel, p0)
        if fn is None:
            continue
        if rel not in fa or rel not in fb:
            rows.append((rel, kind, "MISSING", f"only in {'A' if rel in fa else 'B'}"))
            continue
        try:
            st, det = fn(pa, pb, na, nb)
        except Exception as e:  # a parse failure is a difference to trace, not a pass
            st, det = "ERROR", f"{type(e).__name__}: {e}"
        rows.append((rel, kind, st, det))
    finals = [r for r in rows if r[1] == "final"]
    bad = [r for r in rows if r[2] != "SAME"]
    for r in bad:
        print(f"  {a.label} {r[2]:7} {r[1]:11} {r[0]}: {r[3]}")
    nf = sum(r[2] == "SAME" for r in finals)
    verdict = "PASS" if not bad and len(finals) == len(GRAPH) else "FAIL"
    print(f"E {a.label}: {verdict} -- {len(rows) - len(bad)} of {len(rows)} items SAME; final graph files "
          f"{nf} of {len(GRAPH)} byte-identical"
          + ("" if len(finals) == len(GRAPH) else f" ({len(GRAPH) - len(finals)} absent)"))
    if a.out:
        with open(a.out, "w") as f:
            f.write("label\tpath\trule\tstatus\tdetail\n")
            for r in rows:
                f.write(f"{a.label}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
