# Reproducing the manuscript

These scripts rebuild the manuscript's results from public data: the 56 test sets of 14 datasets, their assemblies by
HARICOT, Rawsamble, Flye and wtdbg2, the QUAST scores, the results table, the appendix figure and the forests. Run them
from the root of the checkout. Each script prints its usage with `-h`.

## Requirements

- Linux x86-64, bash, git, curl, GNU time at `/usr/bin/time`, conda, [rustup](https://rustup.rs), a C/C++ compiler
  and CMake >= 3.16.
- The Liberation Sans font for the appendix figure (package `fonts-liberation` or `liberation-sans-fonts`).
- An NVIDIA GPU for basecalling with dorado. HARICOT reads only the raw signal. The basecalls feed the truth, the
  evaluation and the training labels, and they are the input of Flye and wtdbg2.
- Disk space for the raw reads: 10 GB for D6 below, several TB for all datasets ([All datasets](#all-datasets)).

Compute per set in the manuscript's runs (16 threads, Intel Xeon E5-2650 v4; basecalling, truth and scoring not
included):

| sets | HARICOT | Flye, wtdbg2, Rawsamble |
|---|---|---|
| bacteria (D1-D6) | up to 7 min, 7 GiB | up to 20 min, 10 GiB |
| eukaryotes (D7-D14) | up to 5.9 h, 61 GiB | up to 39 h (Rawsamble on `hg18_both_top`), 46 GiB |

## Installation

Set the work directory `R` in `reproduce/config.sh`. It holds the external tools, the downloads, the read sets, the
basecalls and the runs. `THREADS` (16) sets the threads of each step. Both can also be given on the command line:
`R=/path/to/work bash reproduce/run.sh sd`.

```bash
bash reproduce/0_install.sh
```

This builds HARICOT, builds wtdbg2 and RawHash2 (Rawsamble) at the commits of the manuscript's runs, installs dorado
1.4.0 and 0.9.2 with their models, and creates the conda environments of `envs/`. It downloads the forests and the read
lists of all 86 sets from Zenodo ([10.5281/zenodo.23963126](https://doi.org/10.5281/zenodo.23963126)), checks their
SHA-256 and unpacks them into `models/` and `read_lists/`. The scripts find the environments; none needs to be active.
A rerun skips what is installed.

## One dataset end to end

D6 (*S. dysgalactiae*, key `sd`) has the smallest download, 10.3 GB. Its four sets are `sd_fwd`, `sd_fwd_c12`,
`sd_fwd_c6` and `sd_both_top` ([set names](#all-datasets)).

```bash
bash reproduce/run.sh sd          # steps 1-5
bash reproduce/6_table.sh sd
```

| step | what it does | output |
|---|---|---|
| `1_download.sh` | downloads the reads, writes each set and checks every read's signal against the read list; downloads the reference | `$R/sets/SET.pod5`, `$R/refs/sd.fa` |
| `2_basecall.sh` | basecalls `sd_both_top`, which holds the reads of all four sets, with move tables | `$R/calls/sd.bam` |
| `3_truth.sh` | maps every signal event to its basecalled base and to its place on the reference; writes each set's basecalls | `$R/truth/SET` |
| `4_assemble.sh` | assembles each set with HARICOT, and with Rawsamble, Flye and wtdbg2 ([evaluation/](../evaluation/README.md)) | `$R/runs/SET/ours`, `$R/runs/SET/TOOL` |
| `5_score.sh` | HARICOT's evaluation sequence; QUAST and the connected components of every assembly; the Bandage panels of `sd_fwd` | `$R/runs/SET/score`, `$R/runs/sd_fwd/bandage` |
| `6_table.sh` | the results table and the appendix figure | `$R/table.tsv`, `$R/figure_*.pdf` |

Without a GPU, run step 2 on another machine with the same `R`, or put basecalls made elsewhere at `$R/calls/sd.bam`
(dorado 1.4.0, `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`, `--emit-moves --disable-read-splitting`); the step then skips.

Each step skips the sets it has completed, so a rerun after a failure resumes where it stopped. To redo a step for a
set, delete its output: `$R/runs/SET/ours`, `$R/runs/SET/score`, or `$R/runs/SET/TOOL.done` for Rawsamble, Flye and
wtdbg2, whose failures are kept as their results.

With the published forests, HARICOT's assemblies match the manuscript's runs byte for byte
(`reproduce/expected.sha256`), and step 4 ends with:

```
== 4_assemble.sh sd: done; HARICOT: 4 of 4 sets match the manuscript's runs
```

The manuscript's values for HARICOT and Rawsamble (QUAST against the main contig `CP149485.1`):

| set | tool (`table.tsv` row) | contigs | largest | NG50 | NGA50 | misassemblies | genome fraction (%) |
|---|---|---:|---:|---:|---:|---:|---:|
| `sd_fwd` | HARICOT (`ours_seq`) | 5 | 2,223,364 | 2,223,364 | 666,841 | 6 | 96.585 |
| | Rawsamble (`rawsamble`) | 17 | 459,612 | 213,912 | 44,524 | 59 | 89.215 |
| `sd_fwd_c12` | HARICOT | 2 | 2,208,666 | 2,208,666 | 1,372,513 | 2 | 97.090 |
| | Rawsamble | 9 | 137,731 | - | - | 17 | 20.257 |
| `sd_fwd_c6` | HARICOT | 13 | 584,056 | 354,386 | 210,128 | 6 | 93.740 |
| | Rawsamble | 5 | 73,922 | - | - | 11 | 11.294 |
| `sd_both_top` | HARICOT | 17 | 3,272,293 | 3,272,293 | 548,522 | 16 | 99.733 |
| | Rawsamble | 32 | 537,257 | 146,490 | 42,579 | 83 | 85.219 |

Every tool's QUAST values also depend on the basecalls: they are the input of Flye and wtdbg2, and they give HARICOT's
and Rawsamble's sequences their bases. Flye and wtdbg2 also vary slightly from run to run.

## All datasets

The 56 test sets are the four sets of each of 14 datasets. Set names start with the dataset's key:

| suffix | reads |
|---|---|
| `_fwd` | forward-strand reads, 30x where the data allow |
| `_fwd_c12`, `_fwd_c6` | subsets of `_fwd` at 12x and 6x |
| `_both_top` | `_fwd` plus the reverse-strand reads at the same depth |

```bash
for KEY in cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18; do bash reproduce/run.sh $KEY; done
bash reproduce/6_table.sh
```

`reproduce/datasets.tsv` gives each dataset's reference file, main contig, genome size (the main contig's length, the
`-g` of Flye and wtdbg2), sample rate and QUAST options.

| id | key | species | main contig | sample rate | download |
|---|---|---|---|---|---:|
| D1 | `cj` | *C. jejuni* | `CP149499.1` | 5 kHz | 41 GB |
| D2 | `lm` | *L. monocytogenes* | `CP149495.1` | 5 kHz | 59 GB |
| D3 | `vp1` | *V. parahaemolyticus* chr. 1 | `CP149505.1` | 5 kHz | 85 GB (D3, D4) |
| D4 | `vp2` | *V. parahaemolyticus* chr. 2 | `CP149506.1` | 5 kHz | |
| D5 | `se` | *S. enterica* | `CP149507.1` | 5 kHz | 179 GB |
| D6 | `sd` | *S. dysgalactiae* | `CP149485.1` | 5 kHz | 10 GB |
| D7 | `at2` | *A. thaliana* chr. 2 | `NC_003071.7` | 5 kHz | 164 GB (D7-D9) |
| D8 | `at3` | *A. thaliana* chr. 3 | `NC_003074.8` | 5 kHz | |
| D9 | `at5` | *A. thaliana* chr. 5 | `NC_003076.8` | 5 kHz | |
| D10 | `dm2L` | *D. melanogaster* 2L | `NT_033779.5` | 4 kHz | 228 GB (D10-D13) |
| D11 | `dm2R` | *D. melanogaster* 2R | `NT_033778.4` | 4 kHz | |
| D12 | `dm3L` | *D. melanogaster* 3L | `NT_037436.4` | 4 kHz | |
| D13 | `dmX` | *D. melanogaster* X | `NC_004354.4` | 4 kHz | |
| D14 | `hg18` | human HG002 chr. 18 | `chr18_MATERNAL` | 5 kHz | 1.45 TB |

- Datasets that share a source download it once: files already in `$R/downloads` are kept. Once a dataset's sets are
  written, its downloads can be deleted.
- *A. thaliana* (D7-D9) is one gzip archive, which step 1 downloads whole, once, and checks against its published MD5.
- The *D. melanogaster* files are fast5, which step 1 converts to POD5. These sets are sampled at 4 kHz, and step 2
  basecalls them with dorado 0.9.2 (`dna_r10.4.1_e8.2_400bps_sup@v4.1.0`, `--emit-moves`).
- D14 is scored against the maternal haplotype and, in a second QUAST run, against the paternal chr. 18 (`pat_NGA50`).
  Step 1 downloads both reference files and needs about 4 GB of memory for them.
- Step 2 can run on a GPU machine and the other steps elsewhere, on the same `R`: `run.sh` skips the basecalling once
  `$R/calls/KEY.bam` exists.

## Results table and figure

```bash
bash reproduce/6_table.sh
```

Writes `$R/table.tsv` over the scored sets of the datasets given (default: all 14): one row per set and tool, with
QUAST's values, the connected components of the assembly graph, time and peak memory
([columns](../evaluation/README.md#tabletsv)). It also writes the appendix figure, `$R/figure_bacteria.pdf` (D1-D6)
and `$R/figure_eukaryotes.pdf` (D7-D14), with a PNG preview of each. A dataset not yet scored shows "no panel".

## Retraining the forests

The forests train on 30 sets of five bacteria (`training/pool.tsv`), T1-T5 in the manuscript: *E. coli* (T1,
`ek_fwd`), *S. aureus* (T2, `sa_fwd`), *H. pylori* (T3, `hp2_fwd`), *T. denticola* (T4, `td2_fwd`) and *Nostoc* sp.
PCC 7120 (T5, `an2_fwd`). For each genome there are
its forward-strand reads at 30x, 12x and 6x, plus the same reads cut into 4 kb fragments at the same three depths. The
downloads are 19 GB for `ek_fwd`, `hp2_fwd`, `td2_fwd` and `an2_fwd` together, and 1.65 TB for `sa_fwd`, whose source
files mix 96 barcoded isolates.

```bash
bash reproduce/train.sh
```

For each genome, the script writes the three read sets, basecalls the 30x set (GPU), computes its truth and cuts the
fragments at exact base positions (`truth/fragments.sh`). It writes the table of the 30 training sets
(`$R/train_sets.tsv`) and trains (`training/train.sh`): for each set, HARICOT's candidate detection, a label for each
candidate overlap from the truth, and its features; then each forest five times (seeds 0-4). After basecalling, truth
and training take about 2 hours on 16 cores. `--prepare-only` stops after the table.

The 21 forest files and their checksums go to `$R/train/forests/`, and the script compares the checksums with the
published ones. With basecalls identical to the manuscript's (they are not distributed), the retrained forests equal
the published forests byte for byte. Other basecalls give other labels and other forests. To assemble with the
retrained forests, add `--models $R/train/forests` to `haricot assemble`; it prints `forests=CUSTOM` when they differ
from the published ones.

## Checks

On D6's smallest set, after the steps above:

```bash
. reproduce/config.sh                       # R and THREADS, as the scripts see them
PY=$(conda_env haricot)/bin/python
checks/invariance.sh $R/sets/sd_fwd_c6.pod5 $R/invariance --threads $THREADS --python $PY
target/release/haricot assemble $R/sets/sd_fwd_c6.pod5 -o $R/rerun_sd_fwd_c6 -t $THREADS
checks/cmp_run.sh $R/runs/sd_fwd_c6/ours $R/rerun_sd_fwd_c6
```

`invariance.sh` exits with status 0 when the assemblies of the reads reversed, shuffled and renamed equal the
original's; `cmp_run.sh` when every output file of the two runs is byte-identical. The tests of the code are in
[tests/](../tests/README.md).

After `train.sh`, this command exits with status 0 when training is deterministic and independent of read names and
row order:

```bash
$PY training/train_invariance.py --feat $R/train/work/feat --det $R/train/work/det --bin target/release \
    --exclude-reads $R/train/work/excluded_reads.txt --train hp2_fwd_c6,td2_fwd_c6 --test ek_fwd_c6 --out $R/train_invariance
```

## Software versions

| tool | version |
|---|---|
| Rust | 1.96.0 |
| dorado | 1.4.0 with `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`, `--emit-moves --disable-read-splitting` (5 kHz sets); 0.9.2 with `dna_r10.4.1_e8.2_400bps_sup@v4.1.0`, `--emit-moves` (4 kHz *D. melanogaster* sets) |
| Python | 3.14.4 (`haricot`); 3.12.13 (`haricot-assemblers`) |
| numpy, pandas, scipy, scikit-learn, joblib | 2.4.4, 3.0.2, 1.17.1, 1.9.1, 1.6.0 |
| pod5, pysam, Pillow, matplotlib | 0.3.39, 0.24.0, 12.2.0, 3.10.9 |
| minimap2 | 2.31-r1302 (truth); 2.30 (QUAST) |
| winnowmap, meryl | 2.03, 1.4.1 |
| samtools | 1.23.1 (`haricot`); 1.22.1 (`haricot-assemblers`) |
| seqkit, miniasm, gfatools | 2.13.0, 0.3-r179, 0.5.5 |
| Flye | 2.9.6-b1802 |
| wtdbg2 | 2.5 (`b77c565`) |
| RawHash2 | 2.5 (`92366af`, STORMgroup/RawHash2) |
| QUAST | 5.3.0 |
| BandageNG | 2026.4.1 |
