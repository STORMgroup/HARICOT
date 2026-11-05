#!/bin/bash
# Signal-exact truth of one read set: every retained event of every read placed on the reference, through the
# basecaller's move table and two aligners. Training labels and ours' evaluation sequence read it.
#
#   truth/truth.sh READS.pod5 CALLS.bam REF.fa OUT [--threads N] [--bin-dir DIR] [--python PY]
#
#   READS.pod5  the read set
#   CALLS.bam   dorado basecalls of (at least) these reads with move tables (unaligned BAM, --emit-moves)
#   REF.fa      the reference
#   OUT         new directory:
#                 spans.bin      raw-signal span of every retained event (haricot-truth-events)
#                 ids.txt        the set's read ids; moves.bam, reads.fq.gz: their basecalls
#                 rep_k15.txt    the reference's repetitive 15-mers (meryl), for winnowmap
#                 mm2.bam, wm.bam  minimap2 and winnowmap alignments (map-ont, primary + supplementary)
#                 truth.evref.bin, truth.segs.tsv, truth.anchored.tsv, truth.summary.txt   (event_ref.py)
#                 truth.junctions.tsv   (junctions.py --min-support 3)
#                 contigs.tsv    reference contig lengths
# Tools on PATH: samtools, minimap2, winnowmap, meryl; Python with pod5, pysam, numpy.
set -euo pipefail
usage() { sed -n '5,5p' "$0" | sed 's/^# *//'; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
TH=4; BIN=$HERE/../target/release; PY=python3
ARGS=()
while [ $# -gt 0 ]; do
    case $1 in
        --threads) TH=$2; shift 2 ;;
        --bin-dir) BIN=$2; shift 2 ;;
        --python) PY=$2; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1"; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} = 4 ] || usage
POD=${ARGS[0]}; CALLS=${ARGS[1]}; REF=${ARGS[2]}; OUT=${ARGS[3]}
for f in "$POD" "$CALLS" "$REF"; do [ -s "$f" ] || { echo "no file $f"; exit 1; }; done
[ ! -e "$OUT" ] || { echo "$OUT exists"; exit 1; }
[ -x "$BIN/haricot-truth-events" ] || { echo "no $BIN/haricot-truth-events: cargo build --release"; exit 1; }
for x in samtools minimap2 winnowmap meryl; do command -v $x > /dev/null || { echo "no $x on PATH"; exit 1; }; done
export LC_ALL=C
mkdir -p "$OUT"

"$BIN/haricot-truth-events" "$POD" "$OUT/spans.bin" 2>&1 | tail -1
"$PY" -I -c "import pod5, sys; [print(x) for x in pod5.Reader(sys.argv[1]).read_ids]" "$POD" | sort -u > "$OUT/ids.txt"
samtools view -@ "$TH" -b -F 0x900 -N "$OUT/ids.txt" "$CALLS" -o "$OUT/moves.bam"
samtools fastq -@ "$TH" "$OUT/moves.bam" 2> /dev/null | gzip -1 > "$OUT/reads.fq.gz"
meryl count k=15 threads="$TH" output "$OUT/ref.merylDB" "$REF" > /dev/null 2>&1
meryl print greater-than distinct=0.9998 "$OUT/ref.merylDB" > "$OUT/rep_k15.txt" 2> /dev/null
rm -rf "$OUT/ref.merylDB"
minimap2 -t "$TH" -ax map-ont --secondary=no "$REF" "$OUT/reads.fq.gz" 2> /dev/null \
    | samtools sort -@ 2 -o "$OUT/mm2.bam" && samtools index "$OUT/mm2.bam"
winnowmap -t "$TH" -W "$OUT/rep_k15.txt" -ax map-ont --secondary=no "$REF" "$OUT/reads.fq.gz" 2> /dev/null \
    | samtools sort -@ 2 -o "$OUT/wm.bam" && samtools index "$OUT/wm.bam"
"$PY" -I "$HERE/event_ref.py" "$OUT/spans.bin" "$OUT/moves.bam" "$OUT/mm2.bam" "$OUT/wm.bam" --out "$OUT/truth" 2>&1 | tail -4
"$PY" -I "$HERE/junctions.py" "$OUT/truth.segs.tsv" --out "$OUT/truth.junctions.tsv" --min-support 3 2>&1 | tail -1
samtools view -H "$OUT/mm2.bam" | awk -F'\t' -f "$HERE/sq_lengths.awk" > "$OUT/contigs.tsv"
echo "$OUT: $(wc -l < "$OUT/ids.txt") reads, $(samtools view -c "$OUT/moves.bam") basecalls"
