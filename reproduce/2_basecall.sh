#!/bin/bash
# Step 2: the basecalls of one dataset, with move tables (needs an NVIDIA GPU).
#
#   bash reproduce/2_basecall.sh KEY
#
#   KEY   a dataset of reproduce/datasets.tsv: cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18
#
# Basecalls $R/sets/KEY_both_top.pod5 (step 1), which holds the reads of all four sets, into $R/calls/KEY.bam:
# dorado 1.4.0 with dna_r10.4.1_e8.2_400bps_sup@v5.2.0 and --emit-moves --disable-read-splitting for the 5 kHz
# datasets, dorado 0.9.2 with dna_r10.4.1_e8.2_400bps_sup@v4.1.0 and --emit-moves for the 4 kHz D. melanogaster ones.
# Skipped when $R/calls/KEY.bam exists: basecalls made elsewhere with these settings can be put there.
# Needs: reproduce/0_install.sh (dorado and its models), nvidia-smi.
set -euo pipefail
. "$(dirname "$0")/config.sh"
one_key "$@"

out=$R/calls/$KEY.bam
if [ -s "$out" ]; then
    skip "$KEY" "$out"
else
    pod=$R/sets/${KEY}_both_top.pod5
    [ -s "$pod" ] || die "no $pod: run reproduce/1_download.sh $KEY first"
    basecall "$pod" "$RATE" "$out"
fi
say "== ${0##*/} $KEY: done"
