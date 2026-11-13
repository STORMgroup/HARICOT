#!/bin/bash
# Step 4: HARICOT, Rawsamble, Flye and wtdbg2 on each of a dataset's four sets, as in the manuscript's runs.
#
#   bash reproduce/4_assemble.sh KEY
#
#   KEY   a dataset of reproduce/datasets.tsv: cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18
#
# For SET in KEY_fwd KEY_fwd_c12 KEY_fwd_c6 KEY_both_top, into $R/runs/SET:
#   ours/   haricot assemble on $R/sets/SET.pod5 with the published forests. Its unitigs.gfa is compared with the
#           manuscript's run (reproduce/expected.sha256), and the result printed.
#   flye/, wtdbg2/, wtdbg2_L0/, rawsamble/   evaluation/run_tools.sh: Flye (--nano-hq) and wtdbg2 (-x ont;
#           wtdbg2_L0: again with -L 0) on the set's basecalls ($R/truth/SET/reads.fq.gz), Rawsamble (RawHash2's
#           overlaps, miniasm) on the raw signal; per tool TOOL.time.tsv and TOOL.log, and TOOL.done once it has run. A tool's failure is its result,
#           as in the manuscript (table.tsv: failed_steps).
# HARICOT writes to ours.part, renamed when complete. A complete HARICOT run, and a tool with TOOL.done, is skipped;
# delete it to run it again.
# Needs: steps 1 and 3; reproduce/0_install.sh (HARICOT, the forests, wtdbg2, RawHash2, the conda environments).
set -euo pipefail
. "$(dirname "$0")/config.sh"
one_key "$@"
use_haricot
export PATH=$R/tools/wtdbg2:$PATH                # wtdbg2, wtpoa-cns
ASM=$(conda_env "$ENV_ASSEMBLERS")                # Flye, stacked on haricot for run_tools.sh
[ -x "$BIN/haricot" ] || die "no $BIN/haricot: run reproduce/0_install.sh"
for x in flye wtdbg2 wtpoa-cns miniasm gfatools seqkit; do
    PATH="$ASM/bin:$PATH" type -P "$x" > /dev/null || die "no $x: run reproduce/0_install.sh"
done

same=0
for SET in $SETS; do
    pod=$R/sets/$SET.pod5
    run=$R/runs/$SET
    [ -s "$pod" ] || die "no $pod: run reproduce/1_download.sh $KEY first"
    if [ -d "$run/ours" ]; then
        skip "$SET HARICOT" "$run/ours"
    else
        rm -rf "$run/ours.part"
        mkdir -p "$run"
        "$BIN/haricot" assemble "$pod" -o "$run/ours.part" -t "$THREADS"
        mv "$run/ours.part" "$run/ours"
    fi
    want=$(awk -v f="$SET/ours/unitigs.gfa" '$2 == f {print $1}' "$HERE/expected.sha256")
    got=$(sha256sum < "$run/ours/unitigs.gfa" | cut -c1-64)
    if [ "$got" = "$want" ]; then
        say "$SET: unitigs.gfa matches the manuscript's run"
        same=$((same + 1))
    else
        say "$SET: unitigs.gfa DIFFERS from the manuscript's run (SHA-256 $got; manuscript ${want:-none})"
    fi

    fq=$R/truth/$SET/reads.fq.gz
    [ -s "$fq" ] || die "no $fq: run reproduce/3_truth.sh $KEY first"
    for t in flye wtdbg2 wtdbg2_L0 rawsamble; do
        if [ -e "$run/$t.done" ]; then skip "$SET $t" "$run/$t"; continue; fi
        PATH="$ASM/bin:$PATH" "$ROOT/evaluation/run_tools.sh" "$pod" "$fq" "$run" --genome-size "$GSIZE" \
            --threads "$THREADS" --tools "$t" --rawhash2-dir "$R/tools/RawHash2" --python "$PY"
        date '+%F %T' > "$run/$t.done"
    done
done
say "== ${0##*/} $KEY: done; HARICOT: $same of 4 sets match the manuscript's runs"
