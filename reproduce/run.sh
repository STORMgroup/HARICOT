#!/bin/bash
# Steps 1-5 on one dataset: read sets and reference, basecalls, truth, the four assemblers, scores.
#
#   bash reproduce/run.sh KEY
#
#   KEY   a dataset of reproduce/datasets.tsv: cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18
#
# Each step skips what it has completed, so a rerun after a failure resumes where it stopped. Step 2 needs a GPU unless
# $R/calls/KEY.bam exists. reproduce/6_table.sh then writes the results table and the figure.
set -euo pipefail
. "$(dirname "$0")/config.sh"
one_key "$@"
for step in 1_download 2_basecall 3_truth 4_assemble 5_score; do
    bash "$HERE/$step.sh" "$KEY"
done
