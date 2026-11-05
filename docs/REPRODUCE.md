# Reproducing the manuscript

This guide rebuilds the manuscript's results from public data. Sections 1-3 cover the requirements, the installation,
and the published forests and read lists. Section 4 runs one test set from download to scores. Section 5 does the same for all
56 test sets, and section 6 builds the results table and the appendix figure. Section 7 retrains the forests, and
section 8 lists the checks.

Each script documents its inputs and outputs in its header and prints its usage with `-h`.

## 1. Requirements

- Linux x86-64, bash, git, curl, GNU time at `/usr/bin/time`, conda, [rustup](https://rustup.rs), a C/C++ compiler
  and CMake >= 3.16.
- The Liberation Sans font for the appendix figure (package `fonts-liberation` or `liberation-sans-fonts`).
- An NVIDIA GPU for basecalling with dorado. HARICOT reads only the raw signal. The basecalls feed the truth, the
  evaluation and the training labels, and they are the input of Flye and wtdbg2.
- Disk space for the raw reads: 10 GB for the example in section 4, and several TB for all sets (section 5 lists
  the download for each dataset).

Compute per set in the manuscript runs (16 threads, Intel Xeon E5-2650 v4; basecalling, truth and scoring not included):

| sets | HARICOT | Flye, wtdbg2, Rawsamble |
|---|---|---|
| bacteria (D1-D6) | up to 7 min, 7 GiB | up to 20 min, 10 GiB |
| eukaryotes (D7-D14) | up to 5.9 h, 61 GiB | up to 39 h (Rawsamble on `hg18_both_top`), 46 GiB |

## 2. Installation

Clone and build HARICOT. `rust-toolchain.toml` pins Rust 1.96.0, which rustup installs on the first build.

```bash
git clone https://github.com/STORMgroup/HARICOT.git
cd HARICOT
cargo build --release --locked
```

The build writes `haricot` and the pipeline programs into `target/release/`. Run every command in this guide from the
root of the checkout.

Create the conda environments. They pin the versions of the manuscript runs ([Software versions](#software-versions)).

```bash
conda env create -f envs/haricot.yml       # haricot: data scripts, truth, training, scoring
conda env create -f envs/assemblers.yml    # haricot-assemblers: Flye and QUAST
conda env create -f envs/bandage.yml       # haricot-bandage: BandageNG, for the appendix figure
```

Choose a work directory `R` for the data, the runs and the external tools. Then build wtdbg2 and RawHash2 (Rawsamble)
at the commits of the manuscript runs. RawHash2's `make cmake` downloads its POD5 library.

```bash
R=/path/to/work
git clone https://github.com/ruanjue/wtdbg2.git $R/tools/wtdbg2
git -C $R/tools/wtdbg2 checkout b77c565
make -C $R/tools/wtdbg2
git clone https://github.com/STORMgroup/RawHash2.git $R/tools/RawHash2
git -C $R/tools/RawHash2 checkout 92366af
git -C $R/tools/RawHash2 submodule update --init --recursive
make -C $R/tools/RawHash2 cmake
```

Install dorado and its basecalling models. The 4 kHz *D. melanogaster* sets need dorado 0.9.2; every other set,
training included, needs dorado 1.4.0.

```bash
curl -L https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz | tar -xz -C $R/tools
curl -L https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.2-linux-x64.tar.gz | tar -xz -C $R/tools
mkdir -p $R/dorado_models                 # dorado download fails without a message if it does not exist
$R/tools/dorado-1.4.0-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
    --models-directory $R/dorado_models
$R/tools/dorado-0.9.2-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.1.0 \
    --models-directory $R/dorado_models
```

### Shell setup

The later sections assume these definitions. Run them in every new shell, from the root of the checkout.

```bash
R=/path/to/work
conda activate haricot
PY=$(command -v python)                  # the haricot environment's Python, for every script
export PATH=$R/tools/wtdbg2:$PATH        # wtdbg2, wtpoa-cns
basecall_5khz() {                        # basecall_5khz READS.pod5 > CALLS.bam
    $R/tools/dorado-1.4.0-linux-x64/bin/dorado basecaller --emit-moves --disable-read-splitting \
        $R/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 "$1"
}
basecall_4khz() {                        # D. melanogaster (dm2L, dm2R, dm3L, dmX)
    $R/tools/dorado-0.9.2-linux-x64/bin/dorado basecaller --emit-moves \
        $R/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.1.0 "$1"
}
mkdir -p $R/sets $R/calls
```

The steps that run Flye, wtdbg2, Rawsamble or QUAST stack `haricot-assemblers` on `haricot`
(`conda activate --stack haricot-assemblers`) and remove it afterwards (`conda deactivate`). Do not run `truth.sh`
while it is stacked: the truth uses minimap2 2.31 from `haricot`, while QUAST uses minimap2 2.30 from
`haricot-assemblers`.

## 3. Forests and read lists

The Zenodo record [10.5281/zenodo.23963126](https://doi.org/10.5281/zenodo.23963126) holds the three trained forests
and the read lists of all 86 sets (56 test, 30 training). Download, verify and unpack both archives in the checkout:

```bash
for f in haricot-forests-v1.tar.gz haricot-read-lists-v1.tar.gz; do
    curl -L -o $f "https://zenodo.org/records/23963126/files/$f?download=1"
done
sha256sum -c - <<'EOF'
2b0c37d21d621539480cf29121184278f70f99308efb206edf1c3bca19e9d4f0  haricot-forests-v1.tar.gz
4c64266e21382a1a47d4e0c75a47b1d27c47b26ec3100a51fc9dd0b8f3b9bd92  haricot-read-lists-v1.tar.gz
EOF
tar -xzf haricot-forests-v1.tar.gz       # -> models/, which haricot assemble reads by default
tar -xzf haricot-read-lists-v1.tar.gz    # -> read_lists/
```

`read_lists/<role>/<set>/` holds a set's read ids in order and, for each read, its public source file and the MD5 of
its raw signal. The training sets of 4 kb fragments (`*_L4000*`) hold the ids and `fragments.tsv` instead, because
`truth/fragments.sh` cuts them from their parent set (section 7). `read_lists/sets.tsv` lists every set. [`datasets/README.md`](../datasets/README.md) describes the
files and the public sources.

Set names start with the dataset's key (`cj`, `at2`, `hg18`, ...):

| suffix | reads |
|---|---|
| `_fwd` | forward-strand reads, 30x where the data allow |
| `_fwd_c12`, `_fwd_c6` | subsets of `_fwd` at 12x and 6x |
| `_both_top` | `_fwd` plus the reverse-strand reads at the same depth |

Each dataset's `_both_top` set therefore contains the reads of its other three sets.

## 4. One test set end to end

This section runs *S. dysgalactiae* at 6x (`sd_fwd_c6`, dataset D6), the test set with the smallest download (10 GB).

```bash
SET=sd_fwd_c6; REF=sd.fa; MAIN=CP149485.1; GSIZE=2204985
```

**Reads and reference.** `download.py` fetches the source files that hold the set's reads. `build_set.py` writes the
set's POD5 and checks every read's signal against the read list; it exits with status 1 on any mismatch.
`references.py` downloads the reference and checks every sequence's length and MD5.

```bash
$PY datasets/download.py $R/downloads --set read_lists/final_test/$SET
$PY datasets/build_set.py read_lists/final_test/$SET $R/downloads -o $R/sets/$SET.pod5
$PY datasets/references.py $R/refs --files $REF
```

**Basecalls** (GPU). The move tables (`--emit-moves`) link each base to its raw signal.

```bash
basecall_5khz $R/sets/$SET.pod5 > $R/calls/$SET.bam
```

**Truth.** `truth.sh` maps every signal event of every read to its basecalled base and to its place on the
reference. HARICOT's evaluation sequence takes its bases from this map. `truth.sh` also writes the set's basecalls
as `reads.fq.gz`, which Flye and wtdbg2 get.

```bash
truth/truth.sh $R/sets/$SET.pod5 $R/calls/$SET.bam $R/refs/$REF $R/truth/$SET --threads 16 --python $PY
```

**HARICOT.**

```bash
target/release/haricot assemble $R/sets/$SET.pod5 -o $R/runs/$SET/ours -t 16
```

**Rawsamble, Flye and wtdbg2.** Each tool gets exactly the set's reads: Rawsamble the raw signal, Flye and wtdbg2 the
basecalls. GNU time records each step in `$R/runs/$SET/<tool>.time.tsv`.

```bash
conda activate --stack haricot-assemblers
evaluation/run_tools.sh $R/sets/$SET.pod5 $R/truth/$SET/reads.fq.gz $R/runs/$SET --genome-size $GSIZE \
    --threads 16 --rawhash2-dir $R/tools/RawHash2 --python $PY
```

**Scores.** `score.sh` writes HARICOT's evaluation sequence (`ours_seq.fa`: each layout step contributes its anchor
read's basecalled bases, with no consensus) and runs QUAST 5.3.0 on every assembly against the main contig, on contigs
of at least 1 kb. It also counts the connected components of each assembly graph. `table.py` collects the results.

```bash
evaluation/score.sh $R/runs/$SET $R/truth/$SET $R/refs/$REF $R/runs/$SET/score --main $MAIN --python $PY
conda deactivate
$PY evaluation/table.py $R/runs $SET --out $R/table_$SET.tsv
```

**Expected result.** With the published forests, HARICOT's assembly matches the manuscript's run byte for byte:

```bash
sha256sum $R/runs/$SET/ours/unitigs.gfa
# 3c24112025d1a20060d3958bfed8b0841713648392e2b025b2427d1477bf2b04
```

The manuscript's values for this set are below. Every tool's QUAST values also depend on the basecalls: they are the
input of Flye and wtdbg2, and they give HARICOT's and Rawsamble's sequences their bases. Flye and wtdbg2 also vary
slightly from run to run.

| tool (`table.py` row) | contigs | largest | NG50 | NGA50 | misassemblies | genome fraction (%) |
|---|---:|---:|---:|---:|---:|---:|
| HARICOT (`ours_seq`) | 13 | 584,056 | 354,386 | 210,128 | 6 | 93.740 |
| Flye (`flye`) | 12 | 707,191 | 281,494 | 281,407 | 1 | 97.890 |
| wtdbg2 (`wtdbg2`) | 12 | 482,994 | 264,338 | 213,352 | 1 | 93.853 |
| wtdbg2's layout without consensus (`wtd_lay`) | 12 | 483,072 | 263,424 | 204,110 | 2 | 87.835 |
| Rawsamble (`rawsamble`) | 5 | 73,922 | - | - | 11 | 11.294 |

## 5. All test sets

The 56 test sets are the four sets of each of 14 datasets. The table gives each dataset's parameters for the commands
below. `MAIN` is the contig QUAST scores against, and `GSIZE` is its length in bases.

| id | `D` | species | `REF` | `MAIN` | `GSIZE` | basecall | `SCORE_OPTS` | download |
|---|---|---|---|---|---:|---|---|---:|
| D1 | `cj` | *C. jejuni* | `cj.fa` | `CP149499.1` | 1768448 | 5 kHz | | 41 GB |
| D2 | `lm` | *L. monocytogenes* | `lm.fa` | `CP149495.1` | 2944530 | 5 kHz | | 59 GB |
| D3 | `vp1` | *V. parahaemolyticus* chr. 1 | `vp.fa` | `CP149505.1` | 3281443 | 5 kHz | | 85 GB (D3, D4) |
| D4 | `vp2` | *V. parahaemolyticus* chr. 2 | `vp.fa` | `CP149506.1` | 1865648 | 5 kHz | | |
| D5 | `se` | *S. enterica* | `se.fa` | `CP149507.1` | 4752118 | 5 kHz | | 179 GB |
| D6 | `sd` | *S. dysgalactiae* | `sd.fa` | `CP149485.1` | 2204985 | 5 kHz | | 10 GB |
| D7 | `at2` | *A. thaliana* chr. 2 | `athal.fa` | `NC_003071.7` | 19698289 | 5 kHz | `--large` | 164 GB (D7-D9) |
| D8 | `at3` | *A. thaliana* chr. 3 | `athal.fa` | `NC_003074.8` | 23459830 | 5 kHz | `--large` | |
| D9 | `at5` | *A. thaliana* chr. 5 | `athal.fa` | `NC_003076.8` | 26975502 | 5 kHz | `--large` | |
| D10 | `dm2L` | *D. melanogaster* 2L | `dmel.fa` | `NT_033779.5` | 23513712 | 4 kHz | `--large` | 228 GB (D10-D13) |
| D11 | `dm2R` | *D. melanogaster* 2R | `dmel.fa` | `NT_033778.4` | 25286936 | 4 kHz | `--large` | |
| D12 | `dm3L` | *D. melanogaster* 3L | `dmel.fa` | `NT_037436.4` | 28110227 | 4 kHz | `--large` | |
| D13 | `dmX` | *D. melanogaster* X | `dmel.fa` | `NC_004354.4` | 23542271 | 4 kHz | `--large` | |
| D14 | `hg18` | human HG002 chr. 18 | `hg002v1.2_matY.fa` | `chr18_MATERNAL` | 78908763 | 5 kHz | `--large --paternal $R/refs/chr18_PAT.fa` | 1.45 TB |

For each dataset, set the variables from its row and run the block below. The block builds the four sets from the
files of `_both_top` and basecalls `_both_top` once, because that set holds the reads of all four. It then runs the
truth and the four assemblers on each set, and scores each set. Use `basecall_4khz` for D10-D13.

```bash
D=cj; REF=cj.fa; MAIN=CP149499.1; GSIZE=1768448; SCORE_OPTS=""
SETS="${D}_fwd ${D}_fwd_c12 ${D}_fwd_c6 ${D}_both_top"

$PY datasets/download.py $R/downloads --set read_lists/final_test/${D}_both_top
for SET in $SETS; do $PY datasets/build_set.py read_lists/final_test/$SET $R/downloads -o $R/sets/$SET.pod5; done
$PY datasets/references.py $R/refs --files $REF
basecall_5khz $R/sets/${D}_both_top.pod5 > $R/calls/$D.bam

for SET in $SETS; do
    truth/truth.sh $R/sets/$SET.pod5 $R/calls/$D.bam $R/refs/$REF $R/truth/$SET --threads 16 --python $PY
    target/release/haricot assemble $R/sets/$SET.pod5 -o $R/runs/$SET/ours -t 16
done

conda activate --stack haricot-assemblers
for SET in $SETS; do
    evaluation/run_tools.sh $R/sets/$SET.pod5 $R/truth/$SET/reads.fq.gz $R/runs/$SET --genome-size $GSIZE \
        --threads 16 --rawhash2-dir $R/tools/RawHash2 --python $PY
    evaluation/score.sh $R/runs/$SET $R/truth/$SET $R/refs/$REF $R/runs/$SET/score --main $MAIN $SCORE_OPTS --python $PY
done
conda deactivate
```

Notes on the downloads:

- Datasets that share a source (D3-D4, D7-D9, D10-D13) download it once: `download.py` keeps the files already
  present.
- *A. thaliana* (D7-D9) is one gzip archive, which `download.py --set` streams in full on every call. Fetch it once
  with `$PY datasets/download.py $R/downloads --source ena_athal` and leave out the `download.py` line of the block.
- The *D. melanogaster* files are fast5; `download.py` converts them to POD5.
- D14 needs both reference files (`--files hg002v1.2_matY.fa chr18_PAT.fa`, about 4 GB of memory).

After scoring, draw the Bandage panels of the dataset's `_fwd` set for the appendix figure. Each graph segment is
coloured by its position on the main contig.

```bash
conda activate --stack haricot-bandage
samtools faidx $R/refs/$REF
evaluation/bandage/panels.sh $R/runs/${D}_fwd $R/refs/$REF $R/runs/${D}_fwd/bandage --main $MAIN --python $PY
conda deactivate
```

For the manuscript, wtdbg2 also ran a second time with `-L 0` (`run_tools.sh --tools flye,wtdbg2,wtdbg2_L0,rawsamble`).
Every test set has a 5 kb read-length floor, so `-L 0` loads the same reads as the `-L 5000` of `-x ont`; the second run
shows wtdbg2's run-to-run variation.

## 6. Results table and appendix figure

```bash
DATASETS="cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18"
SETS=$(for D in $DATASETS; do echo ${D}_fwd ${D}_fwd_c12 ${D}_fwd_c6 ${D}_both_top; done)
$PY evaluation/table.py $R/runs $SETS --out $R/table.tsv
$PY evaluation/bandage/figure.py $R/runs $R/table.tsv $R/figure
```

`table.tsv` has one row per set and tool. Its columns are:

- QUAST's contigs, largest contig, NG50, NGA50, misassemblies and genome fraction (main contig, contigs >= 1 kb), and
  `pat_NGA50` against the paternal chr. 18 for D14;
- the connected components of the assembly graph;
- wall and CPU time summed over the tool's steps, the peak memory of its largest step, and the failed steps;
- `rawsamble_protocol`. When miniasm produces no assembly from RawHash2's overlaps, it is rerun on the overlaps
  filtered by `evaluation/filter_ovl.awk` (the best 25, then 15, then 10 per read). The column records the outcome:
  `PLAIN` (no filter), `DEVIATION` (the filter was needed) or `FAIL` (no assembly).

`figure.py` writes the appendix figure as `$R/figure_bacteria.pdf` (D1-D6) and `$R/figure_eukaryotes.pdf` (D7-D14),
with a PNG preview of each.

## 7. Retraining the forests

The forests train on 30 sets of five bacteria (`training/pool.tsv`): *E. coli* (`ek_fwd`), *S. aureus* (`sa_fwd`),
*H. pylori* (`hp2_fwd`), *T. denticola* (`td2_fwd`) and *Nostoc* sp. PCC 7120 (`an2_fwd`). For each genome there are
its forward-strand reads at 30x, 12x and 6x, plus the same reads cut into 4 kb fragments at the same three depths. The
downloads are 19 GB for `ek_fwd`, `hp2_fwd`, `td2_fwd` and `an2_fwd` together, and 1.65 TB for `sa_fwd`, whose source
files mix 96 barcoded isolates. After download and basecalling, truth and training take about 2 hours on 16 cores.

Build the read sets, the truth and the fragment sets of each genome. `fragments.sh` cuts the reads at exact base
positions into 4 kb fragments, maps the fragments through their parents' truth, and writes `frag.pod5`,
`frag_c12.pod5` and `frag_c6.pod5`.

```bash
GENOMES="ek_fwd sa_fwd hp2_fwd td2_fwd an2_fwd"
for G in $GENOMES; do
    $PY datasets/download.py $R/downloads --set read_lists/training/$G
    for S in $G ${G}_c12 ${G}_c6; do $PY datasets/build_set.py read_lists/training/$S $R/downloads -o $R/sets/$S.pod5; done
    $PY datasets/references.py $R/refs --files $G.fa
    basecall_5khz $R/sets/$G.pod5 > $R/calls/$G.bam
    truth/truth.sh $R/sets/$G.pod5 $R/calls/$G.bam $R/refs/$G.fa $R/truth/$G --threads 16 --python $PY
    truth/fragments.sh $R/sets/$G.pod5 $R/truth/$G $R/frag/$G --base-depth 30 --python $PY
done
```

Write the table of the 30 training sets in the order of `training/pool.tsv`, then train:

```bash
{
    printf 'key\tgenome\tpod5\ttruth\tfragments\n'
    for G in $GENOMES; do
        for S in $G ${G}_c12 ${G}_c6; do
            printf '%s\t%s\t%s\t%s\t-\n' $S $G $R/sets/$S.pod5 $R/truth/$G/truth
        done
        for C in '' _c12 _c6; do
            printf '%s\t%s\t%s\t%s\t%s\n' ${G}_L4000$C $G $R/frag/$G/frag$C.pod5 $R/frag/$G/frag_truth \
                $R/frag/$G/frag.frag.tsv
        done
    done
} > $R/train_sets.tsv
training/train.sh $R/train_sets.tsv $R/train --threads 16 --python $PY
```

For each set, `train.sh` runs HARICOT's candidate detection, labels each candidate overlap from the truth, and computes
its features. It then trains each forest five times (seeds 0-4) and writes the 21 forest files and their checksums to
`$R/train/forests/`. With basecalls identical to the manuscript's (they are not distributed), this reproduces the
published forests byte for byte: `diff $R/train/forests/SHA256SUMS models/SHA256SUMS` prints nothing. Different
basecalls give different labels and different forests. To assemble with the retrained forests, add
`--models $R/train/forests` to `haricot assemble`; it prints `forests=CUSTOM` when they differ from the published ones.

## 8. Checks

| check | command | passes when |
|---|---|---|
| read order and read names | `checks/invariance.sh $R/sets/sd_fwd_c6.pod5 $R/invariance --threads 16 --python $PY` | exit status 0: the assemblies of the reads reversed, shuffled and renamed equal the original's |
| two runs of `haricot assemble` | `checks/cmp_run.sh RUN_A RUN_B` | exit status 0: every output file is byte-identical |
| training | `$PY training/train_invariance.py --feat $R/train/work/feat --det $R/train/work/det --bin target/release --exclude-reads $R/train/work/excluded_reads.txt --train hp2_fwd_c6,td2_fwd_c6 --test ek_fwd_c6 --out $R/train_invariance` | exit status 0: training is deterministic and independent of read names and row order |
| Rust code | `cargo test --release --locked` | all tests pass; the four end-to-end tests report as ignored without their fixture, which is not distributed (`tests/fixtures/README.md` says how to make it) |
| evaluation code | `$PY -m unittest discover -s tests/eval` | all tests pass |

## Software versions

| tool | version |
|---|---|
| Rust | 1.96.0 |
| dorado | 1.4.0 with `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`, `--emit-moves --disable-read-splitting` (5 kHz sets); 0.9.2 with `dna_r10.4.1_e8.2_400bps_sup@v4.1.0`, `--emit-moves` (4 kHz *D. melanogaster* sets) |
| Python | 3.14.4 (`haricot`); 3.12.13 (`haricot-assemblers`) |
| numpy, pandas, scipy, scikit-learn, joblib | 2.4.4, 3.0.2, 1.17.1, 1.9.1, 1.6.0 |
| pod5, pysam, Pillow, matplotlib | 0.3.39, 0.24.0, 12.2.0, 3.10.9 |
| minimap2 | 2.31 (truth); 2.30 (QUAST) |
| winnowmap, meryl | 2.03, 1.4.1 |
| samtools | 1.23.1 (`haricot`); 1.22.1 (`haricot-assemblers`) |
| seqkit, miniasm, gfatools | 2.13.0, 0.3, 0.5.5 |
| Flye | 2.9.6 |
| wtdbg2 | 2.5 (`b77c565`) |
| RawHash2 | `92366af` (STORMgroup/RawHash2) |
| QUAST | 5.3.0 |
| BandageNG | 2026.4.1 |
