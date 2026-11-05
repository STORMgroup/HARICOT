#!/bin/bash
# Score one read set's assemblies as in the manuscript: QUAST on the main contig (contigs >= 1 kb) and the graphs'
# connected components.
#
#   evaluation/score.sh RUN TRUTH REF OUT --main CONTIG [--large] [--paternal PAT.fa] [--threads 4] [--python PY]
#
#   RUN      the set's runs: RUN/ours (haricot assemble -o RUN/ours) and run_tools.sh's RUN/{flye,wtdbg2,rawsamble}
#   TRUTH    truth/truth.sh's output for these reads (truth.evref.bin: the event -> base map; reads.fq.gz)
#   REF      the reference FASTA; --main CONTIG its contig QUAST scores against
#   --large  QUAST --large (the eukaryotes); --paternal PAT.fa a second QUAST run against PAT.fa (HG002)
#   OUT      new directory:
#              ours_seq.fa   HARICOT's sequence: each layout step's anchor read over the step (anchor_seq.py; no consensus)
#              wtd_lay.fa    wtdbg2's layout without its consensus (wtdbg2_lay.py; a control)
#              quast_<label>/report.tsv   per assembly: ours_seq, flye, wtdbg2, wtd_lay, rawsamble (+ wtdbg2_L0)
#              components.tsv  tool, level, graph, segments, components (components.py)
# Tools on PATH: quast (5.3.0, its minimap2 2.30), seqkit; Python with numpy.
set -uo pipefail
usage() { sed -n '5,5p' "$0" | sed 's/^# *//'; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
TH=4; MAIN=""; LARGE=""; PAT=""; PY=python3
ARGS=()
while [ $# -gt 0 ]; do
    case $1 in
        --main) MAIN=$2; shift 2 ;;
        --large) LARGE=--large; shift ;;
        --paternal) PAT=$2; shift 2 ;;
        --threads) TH=$2; shift 2 ;;
        --python) PY=$2; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1"; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} = 4 ] && [ -n "$MAIN" ] || usage
RUN=${ARGS[0]}; TRUTH=${ARGS[1]}; REF=${ARGS[2]}; OUT=${ARGS[3]}
for f in "$RUN/ours/unitig_paths.tsv" "$TRUTH/truth.evref.bin" "$TRUTH/reads.fq.gz" "$REF"; do
    [ -s "$f" ] || { echo "no file $f"; exit 1; }
done
[ ! -e "$OUT" ] || { echo "$OUT exists"; exit 1; }
mkdir -p "$OUT"
seqkit grep -p "$MAIN" "$REF" > "$OUT/ref_main.fa"
[ -s "$OUT/ref_main.fa" ] || { echo "no contig $MAIN in $REF"; exit 1; }

"$PY" -I "$HERE/anchor_seq.py" "$RUN/ours" "$TRUTH/truth.evref.bin" "$TRUTH/reads.fq.gz" "$OUT/ours_seq.fa" 2> "$OUT/ours_seq.log" \
    || echo "anchor_seq.py FAILED ($OUT/ours_seq.log)"
[ -s "$RUN/wtdbg2/asm.ctg.lay.gz" ] && "$PY" -I "$HERE/wtdbg2_lay.py" "$RUN/wtdbg2/asm.ctg.lay.gz" "$OUT/wtd_lay.fa" > /dev/null 2>&1

# QUAST runs with its own installation's directory first on PATH: it aligns with the minimap2 installed next to it
# (QUAST 5.3.0's environment: minimap2 2.30; another minimap2 version changes the reports)
command -v quast > /dev/null || { echo "no quast on PATH"; exit 1; }
QBIN=$(dirname "$(command -v quast)")
quast_one() {  # quast_one LABEL ASM.fa
    [ -s "$2" ] || { echo "$1: no sequence ($2)"; return; }
    PATH="$QBIN:$PATH" quast "$2" -l "$1" -r "$OUT/ref_main.fa" -t "$TH" --min-contig 1000 $LARGE -o "$OUT/quast_$1" \
        > "$OUT/quast_$1.log" 2>&1 || echo "$1: QUAST FAILED ($OUT/quast_$1.log)"
    [ -z "$PAT" ] || PATH="$QBIN:$PATH" quast "$2" -l "$1" -r "$PAT" -t "$TH" --min-contig 1000 $LARGE -o "$OUT/quast_$1_pat" \
        > "$OUT/quast_$1_pat.log" 2>&1 || echo "$1: QUAST (paternal) FAILED"
}
quast_one ours_seq "$OUT/ours_seq.fa"
quast_one flye "$RUN/flye/asm.fa"
quast_one wtdbg2 "$RUN/wtdbg2/asm.fa"
quast_one wtd_lay "$OUT/wtd_lay.fa"
quast_one rawsamble "$RUN/rawsamble/asm.fa"
[ -d "$RUN/wtdbg2_L0" ] && quast_one wtdbg2_L0 "$RUN/wtdbg2_L0/asm.fa"

printf "tool\tlevel\tgraph\tsegments\tcomponents\n" > "$OUT/components.tsv"
for x in ours:layout:ours/unitigs.gfa rawsamble:graph:rawsamble/asm.gfa flye:graph:flye/assembly_graph.gfa \
         wtdbg2:frg_graph:wtdbg2/asm.ctg.dot.gz; do
    IFS=: read -r t l g <<< "$x"
    [ -s "$RUN/$g" ] || continue
    "$PY" -I "$HERE/components.py" "$RUN/$g" | awk -F'\t' -v t=$t -v l=$l -v g=$g '{print t "\t" l "\t" g "\t" $2 "\t" $3}' \
        >> "$OUT/components.tsv"
done
echo "scored: $OUT"
