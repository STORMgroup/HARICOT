#!/bin/bash
# Exact-cut fragments of a read set and their truth, plus depth subsets of the fragments: the training pool's 4 kb sets.
#
#   truth/fragments.sh READS.pod5 TRUTH OUT --base-depth X [--length 4000] [--depths 12,6] [--bin-dir DIR] [--python PY]
#
#   READS.pod5  the read set truth.sh placed (TRUTH = its OUT)
#   OUT         new directory:
#                 frag.pod5, frag.frag.tsv, frag.fq.gz   every read cut at base multiples of --length (fragment.py)
#                 frag.spans.bin, frag_truth.*            their events and placements (event_ref.py --frag-map)
#                 frag_c<D>.pod5                          per depth D: the subset md5("depth:" + id) < D / X (subset.py)
#   --base-depth X   the read set's depth (the subsets' denominator; the training sets: 30)
# Tools: Python with pod5, pysam, numpy.
set -euo pipefail
usage() { sed -n '4,4p' "$0" | sed 's/^# *//'; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
LEN=4000; DEPTHS=12,6; BASE=""; BIN=$HERE/../target/release; PY=python3
ARGS=()
while [ $# -gt 0 ]; do
    case $1 in
        --base-depth) BASE=$2; shift 2 ;;
        --length) LEN=$2; shift 2 ;;
        --depths) DEPTHS=$2; shift 2 ;;
        --bin-dir) BIN=$2; shift 2 ;;
        --python) PY=$2; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1"; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} = 3 ] && [ -n "$BASE" ] || usage
POD=${ARGS[0]}; TRUTH=${ARGS[1]}; OUT=${ARGS[2]}
for f in "$POD" "$TRUTH/moves.bam" "$TRUTH/mm2.bam" "$TRUTH/wm.bam"; do [ -s "$f" ] || { echo "no file $f"; exit 1; }; done
[ ! -e "$OUT" ] || { echo "$OUT exists"; exit 1; }
[ -x "$BIN/haricot-truth-events" ] || { echo "no $BIN/haricot-truth-events: cargo build --release"; exit 1; }
mkdir -p "$OUT"

"$PY" -I "$HERE/fragment.py" "$POD" "$TRUTH/moves.bam" "$OUT/frag" --target-bp "$LEN" 2>&1 | tail -1
"$BIN/haricot-truth-events" "$OUT/frag.pod5" "$OUT/frag.spans.bin" 2>&1 | tail -1
"$PY" -I "$HERE/event_ref.py" "$OUT/frag.spans.bin" "$TRUTH/moves.bam" "$TRUTH/mm2.bam" "$TRUTH/wm.bam" \
    --frag-map "$OUT/frag.frag.tsv" --out "$OUT/frag_truth" 2>&1 | tail -4
for d in ${DEPTHS//,/ }; do
    "$PY" -I "$HERE/subset.py" "$OUT/frag.pod5" "$OUT/frag_c$d.pod5" "$d" "$BASE"
done
