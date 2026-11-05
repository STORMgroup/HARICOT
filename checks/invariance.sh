#!/bin/bash
# Read-order and read-name invariance of haricot assemble on one read set: the assembly must be a function of the reads'
# signals only.
#
#   checks/invariance.sh READS.pod5 OUT [--seed 7] [--rename-seed 11] [--threads N] [--bin-dir DIR] [--python PY]
#
# Runs haricot assemble --keep-work on four inputs:
#   O  READS.pod5 as given
#   R  the same reads in reverse order               (pod5_permute.py --order reverse)
#   S  the same reads in a seeded shuffled order      (pod5_permute.py --order shuffle --seed SEED)
#   N  every read renamed by a seeded random UUID     (pod5_variant.py --mode renr), its ids mapped back (name_back.py)
# and compares R, S and N with O (sym_compare.py): the 8 graph and layout files byte for byte, every candidate cache by
# its canonical dump (haricot-ovl-filter --dump-canonical), every table keyed by read names. PASS: every item the same.
# Python needs pod5.
set -uo pipefail
usage() { sed -n '5,5p' "$0" | sed 's/^# *//'; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
SEED=7; RSEED=11; TH=""; BIN=$HERE/../target/release; PY=python3
ARGS=()
while [ $# -gt 0 ]; do
    case $1 in
        --seed) SEED=$2; shift 2 ;;
        --rename-seed) RSEED=$2; shift 2 ;;
        --threads) TH=$2; shift 2 ;;
        --bin-dir) BIN=$2; shift 2 ;;
        --python) PY=$2; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1"; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} = 2 ] || usage
POD=${ARGS[0]}; OUT=${ARGS[1]}
[ -s "$POD" ] || { echo "no file $POD"; exit 1; }
[ ! -e "$OUT" ] || { echo "$OUT exists"; exit 1; }
for b in haricot haricot-ovl-filter; do [ -x "$BIN/$b" ] || { echo "no $BIN/$b: cargo build --release"; exit 1; }; done
NAME=$(basename "$POD")   # every input keeps this file name: the run's key, so the work files line up
mkdir -p "$OUT/R_in" "$OUT/S_in" "$OUT/N_in"

assemble() {  # assemble INPUT RUN: haricot assemble --keep-work, then the canonical dumps of its candidate caches
    local w=$2/work c
    "$BIN/haricot" assemble "$1" -o "$2" --keep-work ${TH:+-t $TH} > "$2.stdout" 2>&1 || { echo "haricot assemble FAILED ($2.stdout)"; return 1; }
    mkdir -p "$w/canon"
    for c in ovl s2 f ff keb q9c q9c/cand_pool; do
        [ -s "$w/$c.bin" ] || continue
        "$BIN/haricot-ovl-filter" "$w/$c.bin" --dump-canonical "$w/canon/${c//\//_}.tsv" > /dev/null 2>&1 \
            || echo "canonical dump FAILED: $w/$c.bin"
    done
}
result() {  # result LABEL TSV: PASS when every compared item is SAME and the 8 final files are byte-identical
    local n s f
    n=$(($(wc -l < "$2") - 1)); s=$(awk -F'\t' 'NR > 1 && $4 == "SAME"' "$2" | wc -l)
    f=$(awk -F'\t' '$3 == "final" && $4 == "SAME"' "$2" | wc -l)
    echo "$1: $([ "$s" = "$n" ] && [ "$f" = 8 ] && [ "$n" -gt 0 ] && echo PASS || echo FAIL) ($s of $n items the same; final files $f of 8)"
    [ "$s" = "$n" ] && [ "$f" = 8 ] && [ "$n" -gt 0 ]
}

"$PY" -I "$HERE/pod5_permute.py" "$POD" --order reverse --order-out "$OUT/R.order" --dst "$OUT/R_in/$NAME" > "$OUT/inputs.log" 2>&1 \
    && "$PY" -I "$HERE/pod5_permute.py" "$POD" --order shuffle --seed "$SEED" --order-out "$OUT/S.order" --dst "$OUT/S_in/$NAME" \
        >> "$OUT/inputs.log" 2>&1 \
    && "$PY" -I "$HERE/pod5_variant.py" "$POD" --mode renr --seed "$RSEED" --dst "$OUT/N_in/$NAME" --map "$OUT/N.map" \
        >> "$OUT/inputs.log" 2>&1 \
    || { echo "writing the permuted / renamed inputs FAILED ($OUT/inputs.log)"; exit 1; }
assemble "$POD" "$OUT/O" || exit 1
for v in R S N; do assemble "$OUT/${v}_in/$NAME" "$OUT/$v" || exit 1; done
rm -rf "$OUT/R_in" "$OUT/S_in" "$OUT/N_in"
"$PY" -I "$HERE/name_back.py" "$OUT/N" "$OUT/N.map" "$OUT/N.back" > "$OUT/name_back.log" 2>&1 \
    || { echo "name_back.py FAILED ($OUT/name_back.log)"; exit 1; }
ok=0
for x in R:R S:S N:N.back; do
    "$PY" -I "$HERE/sym_compare.py" "$OUT/O" "$OUT/${x#*:}" --label "O-${x%%:*}" --out "$OUT/cmp_${x%%:*}.tsv" > "$OUT/cmp_${x%%:*}.txt" 2>&1
    result "order/name ${x%%:*} vs O" "$OUT/cmp_${x%%:*}.tsv" || ok=1
done
exit $ok
