#!/bin/bash
# Retrain the three forests from the 30 training sets of five bacteria (training/pool.tsv).
#
#   bash reproduce/train.sh [--prepare-only] [GENOME ...]
#
#   GENOME           ek_fwd sa_fwd hp2_fwd td2_fwd an2_fwd (default: all five; training needs all five)
#   --prepare-only   stop after writing the table of training sets
#
# Per genome G: the read sets G, G_c12 and G_c6 ($R/sets; as step 1), the reference $R/refs/G.fa, the basecalls
# $R/calls/G.bam (dorado 1.4.0, GPU; skipped when the file exists, as step 2), the truth $R/truth/G (truth/truth.sh)
# and the 4 kb fragments of G at 30x, 12x and 6x ($R/frag/G; truth/fragments.sh). Then $R/train_sets.tsv: the
# training sets of these genomes in the order of training/pool.tsv, with their POD5, truth and fragment files. Then
# training/train.sh -> $R/train (forests/: the 21 forest files and SHA256SUMS; work/), whose SHA256SUMS is compared
# with the published forests' (models/SHA256SUMS).
# Outputs are written to NAME.part and renamed when complete; a complete output is skipped.
# Needs: reproduce/0_install.sh; a GPU for the basecalls.
set -euo pipefail
. "$(dirname "$0")/config.sh"
prepare_only=0
sel=()
while [ $# -gt 0 ]; do
    case $1 in
        -h | --help) usage; exit 0 ;;
        --prepare-only) prepare_only=1 ;;
        -*) usage >&2; exit 2 ;;
        *)
            case " $GENOMES " in *" $1 "*) sel+=("$1") ;; *) die "unknown genome '$1'; genomes: $GENOMES" ;; esac ;;
    esac
    shift
done
if [ ${#sel[@]} = 0 ]; then read -ra sel <<< "$GENOMES"; fi
need_r
use_haricot
echo "== ${0##*/} ${sel[*]}  R=$R  $(date '+%F %T')"

for G in "${sel[@]}"; do
    get_sets training "$G" "$G" "${G}_c12" "${G}_c6"
    get_refs "$G.fa"
    if [ -s "$R/calls/$G.bam" ]; then
        skip "$G basecalls" "$R/calls/$G.bam"
    else
        basecall "$R/sets/$G.pod5" 5000 "$R/calls/$G.bam"
    fi
    out=$R/truth/$G
    if [ -d "$out" ]; then
        skip "$G truth" "$out"
    else
        rm -rf "$out.part"
        "$ROOT/truth/truth.sh" "$R/sets/$G.pod5" "$R/calls/$G.bam" "$R/refs/$G.fa" "$out.part" --threads "$THREADS" \
            --python "$PY"
        mv "$out.part" "$out"
    fi
    out=$R/frag/$G
    if [ -d "$out" ]; then
        skip "$G fragments" "$out"
    else
        rm -rf "$out.part"
        "$ROOT/truth/fragments.sh" "$R/sets/$G.pod5" "$R/truth/$G" "$out.part" --base-depth 30 --python "$PY"
        mv "$out.part" "$out"
    fi
done

# the training sets: a genome's read sets share its truth; its fragment sets, the fragments' truth (base depth 30)
table=$R/train_sets.tsv
awk -F'\t' -v OFS='\t' -v sel=" ${sel[*]} " -v r="$R" '
    NR == 1 { print "key", "genome", "pod5", "truth", "fragments"; next }
    index(sel, " " $2 " ") == 0 { next }
    $4 == "reads" { print $1, $2, r "/sets/" $1 ".pod5", r "/truth/" $2 "/truth", "-"; next }
    $4 == "fragments" {
        d = ($5 == 30) ? "" : "_c" $5
        print $1, $2, r "/frag/" $2 "/frag" d ".pod5", r "/frag/" $2 "/frag_truth", r "/frag/" $2 "/frag.frag.tsv"
    }' "$ROOT/training/pool.tsv" > "$table.part"
while IFS=$'\t' read -r key genome pod truth frag; do
    if [ "$key" = key ]; then continue; fi
    for f in "$pod" "$truth.evref.bin" "$truth.segs.tsv"; do [ -s "$f" ] || die "$key: no $f"; done
    if [ "$frag" != - ] && [ ! -s "$frag" ]; then die "$key: no $frag"; fi
done < "$table.part"
mv "$table.part" "$table"
say "$table: $(($(wc -l < "$table") - 1)) training sets"
if [ $prepare_only = 1 ]; then
    say "== ${0##*/}: --prepare-only, stopping before the training"
    exit 0
fi

[ ${#sel[@]} = 5 ] || die "the training needs all five genomes ($GENOMES); with fewer, use --prepare-only"
out=$R/train
if [ -d "$out" ]; then
    skip "training" "$out"
else
    rm -rf "$out.part"
    "$ROOT/training/train.sh" "$table" "$out.part" --threads "$THREADS" --python "$PY"
    mv "$out.part" "$out"
fi
if cmp -s "$out/forests/SHA256SUMS" "$ROOT/models/SHA256SUMS"; then
    say "== ${0##*/}: done; the 21 forest files equal the published ones (models/SHA256SUMS)"
else
    say "== ${0##*/}: done; the forests differ from the published ones:"
    diff "$out/forests/SHA256SUMS" "$ROOT/models/SHA256SUMS" || true
fi
