#!/bin/bash
# Step 3: the signal-exact truth of each of a dataset's four sets.
#
#   bash reproduce/3_truth.sh KEY
#
#   KEY   a dataset of reproduce/datasets.tsv: cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18
#
# For SET in KEY_fwd KEY_fwd_c12 KEY_fwd_c6 KEY_both_top: truth/truth.sh on $R/sets/SET.pod5 (step 1), the basecalls
# $R/calls/KEY.bam (step 2) and the reference (step 1) -> $R/truth/SET: every signal event mapped to its basecalled
# base and to its place on the reference, and the set's basecalls (reads.fq.gz), which Flye and wtdbg2 assemble.
# It runs in the haricot environment alone: the truth aligns with its minimap2 2.31 (QUAST, in step 5, with 2.30).
# A set is written to $R/truth/SET.part and renamed when complete; a complete set is skipped.
# Needs: steps 1 and 2; reproduce/0_install.sh (HARICOT, the haricot environment).
set -euo pipefail
. "$(dirname "$0")/config.sh"
one_key "$@"
use_haricot

calls=$R/calls/$KEY.bam
ref=$R/refs/$REF
[ -s "$calls" ] || die "no $calls: run reproduce/2_basecall.sh $KEY first"
[ -s "$ref" ] || die "no $ref: run reproduce/1_download.sh $KEY first"
for SET in $SETS; do
    out=$R/truth/$SET
    if [ -d "$out" ]; then skip "$SET" "$out"; continue; fi
    pod=$R/sets/$SET.pod5
    [ -s "$pod" ] || die "no $pod: run reproduce/1_download.sh $KEY first"
    rm -rf "$out.part"
    "$ROOT/truth/truth.sh" "$pod" "$calls" "$ref" "$out.part" --threads "$THREADS" --python "$PY"
    mv "$out.part" "$out"
    say "$SET: truth -> $out"
done
say "== ${0##*/} $KEY: done"
