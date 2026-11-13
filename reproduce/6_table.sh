#!/bin/bash
# Step 6: the results table and the appendix figure, over the sets step 5 has scored.
#
#   bash reproduce/6_table.sh [KEY ...]
#
#   KEY   datasets of reproduce/datasets.tsv to include (default: all 14)
#
# Writes $R/table.tsv (evaluation/table.py: one row per set and tool) and the appendix figure (evaluation/bandage/
# figure.py: the Bandage panels of each KEY_fwd set): $R/figure_bacteria.pdf (D1-D6) and $R/figure_eukaryotes.pdf
# (D7-D14), with a PNG of each. A dataset without panels shows "no panel".
# Needs: step 5; the haricot environment; the Liberation Sans font.
set -euo pipefail
. "$(dirname "$0")/config.sh"
case ${1:-} in -h | --help) usage; exit 0 ;; esac
need_r
use_haricot
keys=("$@")
if [ $# = 0 ]; then read -ra keys <<< "$KEYS"; fi

sets=()
for k in "${keys[@]}"; do
    dataset "$k"
    for SET in $SETS; do
        if [ -d "$R/runs/$SET/score" ]; then sets+=("$SET"); fi
    done
done
[ ${#sets[@]} -gt 0 ] || die "no scored set of ${keys[*]} in $R/runs: run reproduce/5_score.sh first"
"$PY" -I "$ROOT/evaluation/table.py" "$R/runs" "${sets[@]}" --out "$R/table.tsv"
"$PY" -I "$ROOT/evaluation/bandage/figure.py" "$R/runs" "$R/table.tsv" "$R/figure"
say "== ${0##*/}: ${#sets[@]} sets -> $R/table.tsv, $R/figure_bacteria.pdf, $R/figure_eukaryotes.pdf"
