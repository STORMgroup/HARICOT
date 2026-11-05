#!/bin/bash
# The comparators on one read set, as in the manuscript: Flye and wtdbg2 on the reads' basecalls, Rawsamble on the raw
# signal. Each step is timed (OUT/<tool>.time.tsv: step, wall s, max RSS KB, user s, sys s, exit; GNU time).
#
#   evaluation/run_tools.sh READS.pod5 CALLS.fq.gz OUT --genome-size BP [--threads 16] [--tools flye,wtdbg2,rawsamble]
#                           [--rawhash2-dir DIR] [--python PY]
#
#   CALLS.fq.gz      basecalls holding (at least) the pod5's reads; the tools get exactly the pod5's reads
#                    (OUT/work/reads.fq.gz). Basecalling is not part of the timed steps.
#   --genome-size    the main contig's length (datasets/reference_sequences.tsv); -g of Flye and wtdbg2, rounded to kb
#   --python         a Python with pod5 (default: python3)
#   --tools          any of flye, wtdbg2, wtdbg2_L0 (wtdbg2 -L 0: also reads < 5 kb), rawsamble
#   --rawhash2-dir   a RawHash2 checkout built with `make cmake` (the manuscript: STORMgroup/RawHash2 92366af,
#                    unmodified): its bin/rawhash2 and the R10.4.1 pore model in extern/local_kmer_models
# Outputs: OUT/flye (Flye's output, asm.fa = assembly.fasta), OUT/wtdbg2 (asm.ctg.lay.gz, asm.ctg.dot.gz, asm.fa =
#   wtpoa-cns's consensus), OUT/rawsamble (ovl.paf.gz, asm.gfa, asm.fa, miniasm_protocol.txt: PLAIN = miniasm on the
#   overlaps as written; DEVIATION = only on overlaps filtered by filter_ovl.awk, the cap given; FAIL = no assembly).
# Tools on PATH: flye, wtdbg2, wtpoa-cns, miniasm, gfatools, seqkit; Python with pod5.
set -uo pipefail
usage() { sed -n '5,6p' "$0" | sed 's/^# *//'; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
TH=16; TOOLS=flye,wtdbg2,rawsamble; GLEN=""; RHD=""; PY=python3
ARGS=()
while [ $# -gt 0 ]; do
    case $1 in
        --genome-size) GLEN=$2; shift 2 ;;
        --threads) TH=$2; shift 2 ;;
        --tools) TOOLS=$2; shift 2 ;;
        --rawhash2-dir) RHD=$2; shift 2 ;;
        --python) PY=$2; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1"; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} = 3 ] && [ -n "$GLEN" ] || usage
POD=${ARGS[0]}; CALLS=${ARGS[1]}; OUT=${ARGS[2]}
for f in "$POD" "$CALLS"; do [ -s "$f" ] || { echo "no file $f"; exit 1; }; done
case ",$TOOLS," in *,rawsamble,*) [ -x "$RHD/bin/rawhash2" ] || { echo "--rawhash2-dir: no $RHD/bin/rawhash2"; exit 1; } ;; esac
RH=$RHD/bin/rawhash2; PORE=$RHD/extern/local_kmer_models/uncalled_r1041_model_only_means.txt
GS=$(( (GLEN + 500) / 1000 ))k
mkdir -p "$OUT/work"

# the pod5's reads: ids, sample rate, their basecalls (not timed)
SR=$("$PY" -I -c "
import pod5, sys
rates = set()
with pod5.Reader(sys.argv[1]) as r, open(sys.argv[2], 'w') as f:
    for x in r.reads():
        f.write(str(x.read_id) + '\n'); rates.add(x.run_info.sample_rate)
print(','.join(map(str, sorted(rates))))" "$POD" "$OUT/work/ids.unsorted") || { echo "reading $POD failed"; exit 1; }
sort -u "$OUT/work/ids.unsorted" > "$OUT/work/ids.txt"; rm -f "$OUT/work/ids.unsorted"
case $SR in *,*|"") echo "sample rates '$SR': want exactly one"; exit 1 ;; esac
FQ=$OUT/work/reads.fq.gz
seqkit grep -f "$OUT/work/ids.txt" "$CALLS" -o "$FQ" 2> /dev/null
nb=$(seqkit stats -T "$FQ" | tail -1 | cut -f4); np=$(wc -l < "$OUT/work/ids.txt")
[ "$nb" = "$np" ] || { echo "$CALLS holds $nb of the pod5's $np reads"; exit 1; }
echo "reads $np, sample rate $SR, -g $GS"

T() {  # T TOOL STEP CMD...: CMD under GNU time, one row of OUT/TOOL.time.tsv
    local t=$1 s=$2; shift 2
    /usr/bin/time -a -o "$OUT/$t.time.tsv" -f "$s\t%e\t%M\t%U\t%S\t%x" "$@"
}
run_flye() {
    T flye flye flye --nano-hq "$FQ" -g $GS -t $TH -o "$OUT/flye" > "$OUT/flye.log" 2>&1 || return 1
    cp "$OUT/flye/assembly.fasta" "$OUT/flye/asm.fa"
}
wtdbg2_run() {  # wtdbg2_run TOOL [wtdbg2 options]
    local t=$1 o=$OUT/$1; shift
    mkdir -p "$o"
    T $t wtdbg2 wtdbg2 -x ont "$@" -g $GS -t $TH -i "$FQ" -fo "$o/asm" > "$OUT/$t.log" 2>&1 || return 1
    T $t wtpoa_cns wtpoa-cns -t $TH -i "$o/asm.ctg.lay.gz" -fo "$o/cns.fa" >> "$OUT/$t.log" 2>&1 || return 1
    mv "$o/cns.fa" "$o/asm.fa"
}
run_wtdbg2() { wtdbg2_run wtdbg2; }
run_wtdbg2_L0() { wtdbg2_run wtdbg2_L0 -L 0; }
run_rawsamble() {
    local o=$OUT/rawsamble n ns kept top
    mkdir -p "$o"
    T rawsamble rh_index "$RH" -x r10ava -t $TH --sample-rate $SR -p "$PORE" -d "$o/idx.ind" "$POD" > "$OUT/rawsamble.log" 2>&1 || return 1
    T rawsamble rh_overlap "$RH" -x r10ava -t $TH --sample-rate $SR "$o/idx.ind" "$POD" > "$o/ovl.paf" 2>> "$OUT/rawsamble.log" || return 1
    rm -f "$o/idx.ind"
    n=$(wc -l < "$o/ovl.paf")
    T rawsamble miniasm miniasm -f "$FQ" "$o/ovl.paf" > "$o/asm.gfa" 2>> "$OUT/rawsamble.log"
    ns=$(grep -c '^S' "$o/asm.gfa")
    if [ "$ns" -gt 0 ]; then
        echo "PLAIN: miniasm ran on the unfiltered PAF ($n records)" > "$o/miniasm_protocol.txt"
    else
        for top in 25 15 10; do
            awk -v MINSPAN=2000 -v MINMATCH=100 -v TOPN=$top -f "$HERE/filter_ovl.awk" "$o/ovl.paf" "$o/ovl.paf" \
                > "$o/ovl.top$top.paf" 2>> "$OUT/rawsamble.log"
            kept=$(wc -l < "$o/ovl.top$top.paf")
            T rawsamble miniasm_top$top miniasm -f "$FQ" "$o/ovl.top$top.paf" > "$o/asm.gfa" 2>> "$OUT/rawsamble.log"
            ns=$(grep -c '^S' "$o/asm.gfa"); rm -f "$o/ovl.top$top.paf"
            if [ "$ns" -gt 0 ]; then
                echo "DEVIATION: miniasm needed the density filter TOPN=$top MINSPAN=2000 MINMATCH=100, $kept of $n records kept" \
                    > "$o/miniasm_protocol.txt"
                break
            fi
        done
    fi
    gzip -1 "$o/ovl.paf"
    [ "$ns" -gt 0 ] || { echo "FAIL: miniasm failed on the plain PAF and at every filter cap" > "$o/miniasm_protocol.txt"; return 1; }
    gfatools gfa2fa "$o/asm.gfa" > "$o/asm.fa"
}

for t in ${TOOLS//,/ }; do
    case $t in flye|wtdbg2|wtdbg2_L0|rawsamble) ;; *) echo "unknown tool $t"; continue ;; esac
    rm -rf "${OUT:?}/$t"; : > "$OUT/$t.time.tsv"
    run_$t || echo "$t: FAILED (log $OUT/$t.log)"
    awk -F'\t' -v t=$t 'NF == 6 {w+=$2; if ($3>m) m=$3} END {printf "%s: wall %.1f s, peak RSS %.2f GB\n", t, w, m/1048576}' "$OUT/$t.time.tsv"
done
rm -f "$FQ"
