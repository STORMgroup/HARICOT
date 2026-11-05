#!/bin/bash
# Bandage panels of one read set's assembly graphs, coloured by reference position (the manuscript's appendix figure).
#
#   evaluation/bandage/panels.sh RUN REF OUT --main CONTIG [--threads 16] [--layouts 4] [--python PY]
#
#   RUN   the set's runs and scoring (RUN/ours, RUN/{flye,wtdbg2,rawsamble}, RUN/score/ours_seq.fa: score.sh)
#   REF   the reference FASTA (indexed: REF.fai); --main CONTIG the contig the colours run along
#   OUT   new directory: ours.png, flye.png, wtdbg2.png, rawsamble.png
# Colour (color_gfa.py): a segment's hue is the reference midpoint of its longest primary minimap2 (asm10) block on
# CONTIG, 0 -> 0.85 along it; dark grey: on another contig; light grey: not aligned. HARICOT's unitigs are placed by
# its sequence (score.sh's ours_seq.fa, one record per unitig), wtdbg2's contigs drawn as a graph without links.
# Each panel: BandageNG (2000 x 1500 px) over --layouts random layouts, the least folded kept (pick_best.py).
# Tools on PATH: BandageNG, minimap2; Python with numpy and Pillow.
set -uo pipefail
usage() { sed -n '4,4p' "$0" | sed 's/^# *//'; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
TH=16; NL=4; MAIN=""; PY=python3
ARGS=()
while [ $# -gt 0 ]; do
    case $1 in
        --main) MAIN=$2; shift 2 ;;
        --threads) TH=$2; shift 2 ;;
        --layouts) NL=$2; shift 2 ;;
        --python) PY=$2; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1"; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} = 3 ] && [ -n "$MAIN" ] || usage
RUN=${ARGS[0]}; REF=${ARGS[1]}; OUT=${ARGS[2]}
[ -s "$REF.fai" ] || { echo "no $REF.fai (samtools faidx $REF)"; exit 1; }
mkdir -p "$OUT/work"
export QT_QPA_PLATFORM=offscreen   # BandageNG without a display

render() {  # render NAME GFA COLOURS.csv
    local name=$1 gfa=$2 col=$3 nw i cands=()
    nw=$(awk '/^S\t/{n++} END{w=5*sqrt(n/50); if (w<5) w=5; if (w>80) w=80; printf "%.1f", w}' "$gfa")
    for i in $(seq 1 "$NL"); do
        timeout 1800 BandageNG image "$gfa" "$OUT/work/$name.$i.png" --color "$col" --width 2000 --height 1500 \
            --nodewidth "$nw" > /dev/null 2>&1 &
    done
    wait
    for i in $(seq 1 "$NL"); do [ -s "$OUT/work/$name.$i.png" ] && cands+=("$OUT/work/$name.$i.png"); done
    echo "$name: ${#cands[@]}/$NL layouts, $(grep -c '^S' "$gfa") segments"
    [ ${#cands[@]} -gt 0 ] && "$PY" -I "$HERE/pick_best.py" "$OUT/$name.png" "${cands[@]}"
}
colour() {  # colour NAME SEQS.fa GFA: minimap2 of the sequences -> COLOURS.csv
    minimap2 -t "$TH" -x asm10 --secondary=no "$REF" "$2" 2> /dev/null > "$OUT/work/$1.paf"
    "$PY" -I "$HERE/color_gfa.py" "$OUT/work/$1.paf" --gfa "$3" --fai "$REF.fai" --contig "$MAIN" --out "$OUT/work/$1.csv" | tail -1
}
seq_graph() {  # seq_graph NAME GFA: a graph whose segments carry their sequences
    awk -F'\t' -f "$HERE/gfa_segs.awk" "$2" > "$OUT/work/$1.segs.fa"
    colour "$1" "$OUT/work/$1.segs.fa" "$2"
    render "$1" "$2" "$OUT/work/$1.csv"
}

if [ -s "$RUN/ours/unitigs.gfa" ] && [ -s "$RUN/score/ours_seq.fa" ]; then
    colour ours "$RUN/score/ours_seq.fa" "$RUN/ours/unitigs.gfa"
    render ours "$RUN/ours/unitigs.gfa" "$OUT/work/ours.csv"
fi
[ -s "$RUN/flye/assembly_graph.gfa" ] && seq_graph flye "$RUN/flye/assembly_graph.gfa"
if [ -s "$RUN/wtdbg2/asm.fa" ]; then
    awk -f "$HERE/fa2gfa.awk" "$RUN/wtdbg2/asm.fa" > "$OUT/work/wtdbg2.gfa"
    seq_graph wtdbg2 "$OUT/work/wtdbg2.gfa"
fi
[ -s "$RUN/rawsamble/asm.gfa" ] && seq_graph rawsamble "$RUN/rawsamble/asm.gfa"
rm -f "$OUT"/work/*.png "$OUT"/work/*.segs.fa "$OUT"/work/*.paf
