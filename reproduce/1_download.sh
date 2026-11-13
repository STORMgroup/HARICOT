#!/bin/bash
# Step 1: the read sets and the reference of one dataset, from their public sources.
#
#   bash reproduce/1_download.sh KEY
#
#   KEY   a dataset of reproduce/datasets.tsv: cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18
#
# Downloads into $R/downloads the source files that hold the reads of KEY_both_top, which holds the reads of the
# dataset's other three sets (datasets/download.py; files present are kept; the A. thaliana archive is fetched whole,
# once). Writes $R/sets/SET.pod5 for SET in KEY_fwd KEY_fwd_c12 KEY_fwd_c6 KEY_both_top, checking every read's signal
# against the read list (datasets/build_set.py; it stops on any mismatch). Downloads the reference into $R/refs and
# checks every sequence (datasets/references.py). A set or reference present is skipped.
# Needs: reproduce/0_install.sh (the haricot environment, read_lists/).
set -euo pipefail
. "$(dirname "$0")/config.sh"
one_key "$@"
use_haricot

# shellcheck disable=SC2086  # SETS: four set names
get_sets final_test "${KEY}_both_top" $SETS
if [ "$PAT" = - ]; then get_refs "$REF"; else get_refs "$REF" "$PAT"; fi
say "== ${0##*/} $KEY: done"
