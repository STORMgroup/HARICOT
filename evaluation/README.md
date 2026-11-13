# Evaluation

These scripts run the other assemblers on a read set, score every assembly with QUAST and build the manuscript's
results table and appendix figure. `reproduce/` runs them on all 56 test sets.

HARICOT's output is a layout without sequence. For QUAST, `anchor_seq.py` writes each unitig as the basecalled bases
of its anchor reads, with no consensus. `wtdbg2_lay.py` writes wtdbg2's layout the same way, as a control for what the
missing consensus costs (Supplementary Table S11).

## Files

| file | contents |
|---|---|
| `run_tools.sh` | Rawsamble on the raw signal, Flye and wtdbg2 on the basecalls of one read set; every step timed with GNU time |
| `score.sh` | HARICOT's evaluation sequence, wtdbg2's layout without consensus, QUAST on every assembly, the graphs' connected components |
| `anchor_seq.py` | HARICOT's evaluation sequence: each layout step's anchor read over the events the step advances |
| `wtdbg2_lay.py` | wtdbg2's contig layout without its consensus: one read per step, cut at the layout's offsets |
| `components.py` | connected components of an assembly graph (GFA, or wtdbg2's `.dot`) |
| `filter_ovl.awk` | Rawsamble's overlaps reduced to the best 25, 15 or 10 per read, when miniasm assembles nothing from all of them |
| `table.py` | one row per set and tool: QUAST's values, components, time and peak memory |
| `bandage/panels.sh` | Bandage panels of one set's graphs, each segment coloured by its position on the reference |
| `bandage/figure.py` | the appendix figure from the panels of every dataset |
| `lib/evio.py` | readers of the truth's binary files (`truth/`) and of `unitig_paths.tsv` |

## Usage

One read set, after `haricot assemble` (`RUN/ours`) and the truth (`truth/truth.sh`, which also writes the set's
basecalls):

```bash
evaluation/run_tools.sh SET.pod5 TRUTH/reads.fq.gz RUN --genome-size BP --rawhash2-dir RAWHASH2
evaluation/score.sh RUN TRUTH REF.fa RUN/score --main CONTIG [--large] [--paternal PAT.fa]
python3 evaluation/table.py ROOT SET [SET ...] --out table.tsv
```

Each script documents its inputs, outputs and the tools it needs in its header.

- Each tool gets exactly the set's reads. Basecalling is not timed.
- QUAST scores the main contig only, on contigs of at least 1 kb, with `--large` for the eukaryotes. It aligns with
  the minimap2 of its own environment (2.30); another version changes the reports.
- Rawsamble: when miniasm produces no assembly from RawHash2's overlaps, it is rerun on the overlaps filtered by
  `filter_ovl.awk`. `table.py`'s `rawsamble_protocol` records `PLAIN`, `DEVIATION` (the filter was needed) or `FAIL`.
- Flye and wtdbg2 vary slightly between runs. `--tools wtdbg2_L0` runs wtdbg2 again with `-L 0`, which loads the same
  reads on sets with a 5 kb read floor; the manuscript's runs include it.

## `table.tsv`

One row per set and tool: `ours_seq` (HARICOT), `flye`, `wtdbg2`, `wtdbg2_L0`, `wtd_lay` (wtdbg2's layout without its
consensus) and `rawsamble`. Columns:

- QUAST's contigs, largest contig, NG50, NGA50, misassemblies and genome fraction; `pat_NGA50` against the paternal
  chr. 18 for HG002;
- the connected components of the assembly graph;
- wall and CPU time summed over the tool's steps, the peak memory of its largest step, and the failed steps;
- `rawsamble_protocol`.
