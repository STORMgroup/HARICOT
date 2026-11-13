#!/bin/bash
# Step 5: QUAST on every assembly of a dataset's four sets, and the Bandage panels of its KEY_fwd set.
#
#   bash reproduce/5_score.sh KEY
#
#   KEY   a dataset of reproduce/datasets.tsv: cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18
#
# For SET in KEY_fwd KEY_fwd_c12 KEY_fwd_c6 KEY_both_top: evaluation/score.sh -> $R/runs/SET/score: HARICOT's evaluation
# sequence (ours_seq.fa: each layout step contributes its anchor read's basecalled bases, with no consensus), QUAST 5.3.0
# on every assembly against the dataset's main contig, contigs >= 1 kb (--large for the eukaryotes; for HG002 also
# against the paternal chr. 18), and the connected components of each assembly graph. QUAST runs with the
# haricot-assemblers environment stacked on haricot: it aligns with that environment's minimap2 2.30.
# Then evaluation/bandage/panels.sh on KEY_fwd -> $R/runs/KEY_fwd/bandage: one Bandage image per tool, each segment
# coloured by its position on the main contig (haricot-bandage stacked on haricot).
# Outputs are written to NAME.part and renamed when complete; a complete output is skipped.
# Needs: steps 3 and 4; reproduce/0_install.sh (the conda environments).
set -euo pipefail
. "$(dirname "$0")/config.sh"
one_key "$@"
use_haricot
ASM=$(conda_env "$ENV_ASSEMBLERS")
BAN=$(conda_env "$ENV_BANDAGE")
opts=(--main "$MAIN")
if [ "$QUAST" = large ]; then opts+=(--large); fi
if [ "$PAT" != - ]; then opts+=(--paternal "$R/refs/$PAT"); fi

for SET in $SETS; do
    run=$R/runs/$SET
    out=$run/score
    if [ -d "$out" ]; then skip "$SET score" "$out"; continue; fi
    for t in ours flye wtdbg2 wtdbg2_L0 rawsamble; do
        if [ "$t" = ours ]; then f=$run/ours; else f=$run/$t.done; fi
        [ -e "$f" ] || die "no $f: run reproduce/4_assemble.sh $KEY first"
    done
    rm -rf "$out.part"
    PATH="$ASM/bin:$PATH" "$ROOT/evaluation/score.sh" "$run" "$R/truth/$SET" "$R/refs/$REF" "$out.part" "${opts[@]}" \
        --python "$PY"
    mv "$out.part" "$out"
done

run=$R/runs/${KEY}_fwd
out=$run/bandage
if [ -d "$out" ]; then
    skip "${KEY}_fwd Bandage panels" "$out"
else
    [ -s "$R/refs/$REF.fai" ] || samtools faidx "$R/refs/$REF"
    rm -rf "$out.part"
    PATH="$BAN/bin:$PATH" "$ROOT/evaluation/bandage/panels.sh" "$run" "$R/refs/$REF" "$out.part" --main "$MAIN" \
        --threads "$THREADS" --python "$PY"
    mv "$out.part" "$out"
fi
say "== ${0##*/} $KEY: done"
