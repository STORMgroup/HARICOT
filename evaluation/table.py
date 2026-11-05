#!/usr/bin/env python3
"""The manuscript's metrics per read set and tool, from the runs (run_tools.sh, haricot assemble) and their scoring
(score.sh into RUN/score).

Per row: QUAST on the main contig, contigs >= 1 kb (contigs, largest contig, NG50, NGA50, misassemblies, genome
fraction; pat_NGA50: against the paternal copy, where scored), the graph's connected components (components.py),
elapsed / CPU time and peak memory summed over the tool's timed steps, failed steps, and Rawsamble's miniasm protocol
(PLAIN / DEVIATION / FAIL). Rows: ours_seq (HARICOT), flye, wtdbg2, wtd_lay (wtdbg2's layout without consensus; wtdbg2's
times), rawsamble, and wtdbg2_L0 where it ran.

Usage: table.py ROOT KEY [KEY ...] --out table.tsv     (ROOT/KEY: one set's runs, ROOT/KEY/score: its scoring)
"""
import argparse
import csv
import os
import re
import sys

NOTE = re.compile(r"^Command (exited with non-zero status|terminated by signal) \d+$")
SIGNAL = re.compile(r"^Command terminated by signal (\d+)$")
# label: (time files, graph tool, graph level)
TOOLS = [("ours_seq", ["ours/time.tsv"], "ours", "layout"), ("flye", ["flye.time.tsv"], "flye", "graph"),
         ("wtdbg2", ["wtdbg2.time.tsv"], "wtdbg2", "frg_graph"), ("wtd_lay", ["wtdbg2.time.tsv"], "", ""),
         ("rawsamble", ["rawsamble.time.tsv"], "rawsamble", "graph"),
         ("wtdbg2_L0", ["wtdbg2_L0.time.tsv"], "", "")]
QUAST = [("contigs", "# contigs"), ("largest", "Largest contig"), ("NG50", "NG50"), ("NGA50", "NGA50"),
         ("misassemblies", "# misassemblies"), ("genome_frac", "Genome fraction (%)")]
COLS = ["key", "tool"] + [c for c, _ in QUAST] + ["pat_NGA50", "components", "wall_s", "cpu_s", "rss_gb",
                                                   "failed_steps", "rawsamble_protocol"]


def quast_report(path):
    """QUAST's report.tsv of one assembly: metric -> value as written; {} without the file."""
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                out[p[0]] = p[1]
    return out


def time_totals(paths):
    """GNU time rows (step, wall, max RSS KB, user, sys, exit) of several files: wall and user + sys summed, the peak
    RSS, failed steps (exit not 0; a row after GNU time's "terminated by signal N" line, whose exit reads 0, too).
    None if a file is absent."""
    wall = cpu = 0.0
    rss = failed = 0
    for p in paths:
        if not os.path.exists(p):
            return None
        sig = None
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) == 6:
                wall += float(f[1])
                rss = max(rss, int(f[2]))
                cpu += float(f[3]) + float(f[4])
                failed += f[5].strip() != "0" or sig is not None
                sig = None
            elif NOTE.match(line.strip()):
                m = SIGNAL.match(line.strip())
                sig = m.group(1) if m else None
    return dict(wall=f"{wall:.0f}", cpu=f"{cpu:.0f}", rss=f"{rss / 1048576:.2f}", failed=failed)


def components(path):
    """components.tsv (score.sh): (tool, level) -> components."""
    if not os.path.exists(path):
        return {}
    return {(r["tool"], r["level"]): r["components"] for r in csv.DictReader(open(path), delimiter="\t")}


def rows(root, key):
    run, score = os.path.join(root, key), os.path.join(root, key, "score")
    comps = components(os.path.join(score, "components.tsv"))
    proto_file = os.path.join(run, "rawsamble", "miniasm_protocol.txt")
    proto = open(proto_file).read().split(":")[0].strip() if os.path.exists(proto_file) else ""
    for label, times, gtool, glevel in TOOLS:
        if label == "wtdbg2_L0" and not os.path.isdir(os.path.join(run, "wtdbg2_L0")):
            continue
        q = quast_report(os.path.join(score, f"quast_{label}", "report.tsv"))
        qp = quast_report(os.path.join(score, f"quast_{label}_pat", "report.tsv"))
        t = time_totals([os.path.join(run, x) for x in times])
        r = dict(key=key, tool=label, pat_NGA50=qp.get("NGA50", ""), components=comps.get((gtool, glevel), ""),
                 wall_s=t["wall"] if t else "", cpu_s=t["cpu"] if t else "", rss_gb=t["rss"] if t else "",
                 failed_steps=t["failed"] if t else "", rawsamble_protocol=proto if label == "rawsamble" else "")
        r.update({c: q.get(m, "") for c, m in QUAST})
        yield r


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root")
    ap.add_argument("keys", nargs="+")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = [r for k in a.keys for r in rows(a.root, k)]
    with open(a.out, "w") as w:
        w.write("\t".join(COLS) + "\n")
        for r in out:
            w.write("\t".join(str(r[c]) for c in COLS) + "\n")
    print(f"{len(out)} rows -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
