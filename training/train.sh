#!/bin/bash
# Train the three forests haricot assemble uses (stage1, stage2, rescue) from labelled training sets.
#
#   training/train.sh SETS.tsv OUT [--threads N] [--bin-dir DIR] [--python PY] [--test KEY,...]
#
#   SETS.tsv  one row per training set (header: key genome pod5 truth fragments):
#               key        the set's name (recorded, in table order, in the forests' training list)
#               genome     the genome it belongs to; the row whose key equals its genome is the genome's base set
#               pod5       its reads
#               truth      the prefix of its truth (PREFIX.evref.bin, PREFIX.segs.tsv: truth/truth.sh OUT/truth for a
#                          read set and its subsets, truth/fragments.sh OUT/frag_truth for fragment sets)
#               fragments  a fragment set's frag.tsv (truth/fragments.sh), else -
#   OUT       new directory: OUT/forests/{stage1,stage2,rescue}.* and SHA256SUMS; OUT/work: detections, labels,
#             features, the training logs
#   --test    sets of the table to score with the trained forests (printed only; not used for training)
#
# Steps per set: candidate detection (haricot-detect --production, the first step of haricot assemble), labels (the
# 1 kb signal-exact oracle, haricot-truth-oracle), stage-1 pair features (haricot-pair-features). Then, over all sets:
# the excluded reads (excluded_reads.py), stage 1 (train_stage1.py: 5 forests, seeds 0-4, <= 20,000 labelled rows per
# set), stage-2 graph features on its out-of-bag scores (haricot-graph-features), stage 2 (train_stage2.py), the rescue
# forest (stage 1 without 8 geometry features), each exported for haricot-forest (export_forest.py).
# Python needs numpy, pandas, scikit-learn, joblib, pod5.
set -euo pipefail
usage() { sed -n '4,4p' "$0" | sed 's/^# *//'; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$HERE/../target/release; PY=python3; TEST=""; TH=""
ARGS=()
while [ $# -gt 0 ]; do
    case $1 in
        --threads) TH=$2; shift 2 ;;
        --bin-dir) BIN=$2; shift 2 ;;
        --python) PY=$2; shift 2 ;;
        --test) TEST=$2; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1"; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} = 2 ] || usage
SETS=${ARGS[0]}; OUT=${ARGS[1]}
[ -s "$SETS" ] || { echo "no table $SETS"; exit 1; }
[ ! -e "$OUT" ] || { echo "$OUT exists"; exit 1; }
for b in haricot-detect haricot-truth-oracle haricot-pair-features haricot-graph-features; do
    [ -x "$BIN/$b" ] || { echo "no $BIN/$b: cargo build --release"; exit 1; }
done
[ -z "$TH" ] || export RAYON_NUM_THREADS=$TH
# the recipe of the published forests
SEEDS=0,1,2,3,4; CAP=20000
RESCUE_DROP=z,slope,abs_log_slope,hole_frac,max_hole,resid_med,resid_max,resid_frac1
W=$OUT/work; mkdir -p "$W/det" "$W/labels" "$W/feat" "$W/s2" "$W/model_B" "$W/model_Bq" "$W/logs" "$OUT/forests"
echo "train.sh: $(git -C "$HERE" rev-parse --short HEAD 2> /dev/null || echo '?') sets=$SETS out=$OUT bin=$BIN" \
     "python: $("$PY" -c 'import sys, sklearn, pandas, numpy; print(sys.version.split()[0], "sklearn", sklearn.__version__, "pandas", pandas.__version__, "numpy", numpy.__version__)')"

KEYS=(); declare -A POD TRUTH
while IFS=$'\t' read -r key genome pod truth frag; do
    [ "$key" = key ] && continue
    for f in "$pod" "$truth.evref.bin" "$truth.segs.tsv"; do [ -s "$f" ] || { echo "$key: no file $f"; exit 1; }; done
    KEYS+=("$key"); POD[$key]=$pod; TRUTH[$key]=$truth
done < "$SETS"
for t in ${TEST//,/ }; do [ -n "${POD[$t]:-}" ] || { echo "--test $t is not in $SETS"; exit 1; }; done
TRAIN_ONLY=()
for k in "${KEYS[@]}"; do case ",$TEST," in *",$k,"*) ;; *) TRAIN_ONLY+=("$k") ;; esac; done
TRAIN=$(IFS=,; echo "${TRAIN_ONLY[*]}")
echo "training sets (${#TRAIN_ONLY[@]}): $TRAIN${TEST:+; test: $TEST}"

for k in "${KEYS[@]}"; do
    "$BIN/haricot-detect" "${POD[$k]}" --production --outdir "$W/det/$k.p2f" --save-overlaps "$W/det/$k.ovl.bin" \
        --no-graph > "$W/logs/det.$k.log" 2>&1
    rm -rf "$W/det/$k.p2f"
    "$BIN/haricot-truth-oracle" "$W/det/$k.ovl.bin" "$W/labels/$k.labels.bin" --evref "${TRUTH[$k]}.evref.bin" \
        --segs "${TRUTH[$k]}.segs.tsv" --min-ovl 1000 --min-bins 3 --truth-out "$W/labels/$k.truth1k.tsv" \
        --truth-min 500 > "$W/logs/label.$k.log" 2>&1
    "$BIN/haricot-pair-features" "${POD[$k]}" "$W/det/$k.ovl.bin" --oracle "$W/labels/$k.labels.bin" \
        --truth "$W/labels/$k.truth1k.tsv" --key "$k" --out "$W/feat/$k.tsv" > "$W/logs/feat.$k.log" 2>&1
    echo "[$k] $(grep -E 'prefilter:' "$W/logs/det.$k.log" | tail -1 | cut -c1-100) | $(($(wc -l < "$W/feat/$k.tsv") - 1)) rows"
done

EXCL=$W/excluded_reads.txt
"$PY" "$HERE/excluded_reads.py" "$SETS" "$EXCL"
MB=$W/model_B; MQ=$W/model_Bq
{
    "$PY" "$HERE/train_stage1.py" --train "$TRAIN" --test "$TEST" --tag all5 --features-dir "$W/feat" --out-dir "$MB" \
        --seeds $SEEDS --save-models --cap-per-set $CAP --exclude-reads "$EXCL"
    T1=$("$PY" -c "import json, sys; print(json.load(open(sys.argv[1]))['thresholds']['0.995'])" "$MB/all5.thresholds.json")
    for k in "${TRAIN_ONLY[@]}"; do
        "$BIN/haricot-graph-features" "$W/det/$k.ovl.bin" "$MB/all5__$k.oob.tsv" --t1 "$T1" --out "$W/s2/$k.s2.tsv" > /dev/null 2>&1
    done
    for k in ${TEST//,/ }; do
        "$BIN/haricot-graph-features" "$W/det/$k.ovl.bin" "$MB/all5__$k.scores.tsv" --t1 "$T1" --out "$W/s2/$k.s2.tsv" > /dev/null 2>&1
    done
    "$PY" "$HERE/train_stage2.py" --train "$TRAIN" --test "$TEST" --tag s2all5 --s2-dir "$W/s2" --features-dir "$W/feat" \
        --out-dir "$MB" --seeds $SEEDS --save-models --cap-per-set $CAP --exclude-reads "$EXCL"
    "$PY" "$HERE/export_forest.py" --model-dir "$MB" --tag all5
    "$PY" "$HERE/export_forest.py" --model-dir "$MB" --tag s2all5
    "$PY" "$HERE/train_stage1.py" --train "$TRAIN" --test "$TEST" --tag all5 --features-dir "$W/feat" --out-dir "$MQ" \
        --seeds $SEEDS --save-models --cap-per-set $CAP --drop-features $RESCUE_DROP --exclude-reads "$EXCL"
    "$PY" "$HERE/export_forest.py" --model-dir "$MQ" --tag all5
} > "$W/logs/train.log" 2>&1
grep -E 'out-of-bag AUC' "$W/logs/train.log"

for x in stage1:$MB/all5 stage2:$MB/s2all5 rescue:$MQ/all5; do
    name=${x%%:*}; src=${x#*:}
    for s in ${SEEDS//,/ }; do cp "$src.rf$s.rfb" "$OUT/forests/$name.rf$s.rfb"; done
    cp "$src.rfmeta" "$OUT/forests/$name.rfmeta"
    cp "$src.thresholds.json" "$OUT/forests/$name.thresholds.json"
done
rm -f "$MB"/*.joblib "$MQ"/*.joblib
(cd "$OUT/forests" && sha256sum stage1.* stage2.* rescue.* > SHA256SUMS)
echo "forests: $OUT/forests ($(ls "$OUT/forests" | wc -l) files)"
