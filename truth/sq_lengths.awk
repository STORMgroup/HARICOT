# sq_lengths.awk: a SAM header's @SQ lines -> name<TAB>length (truth.sh's contig table).
# Usage: samtools view -H x.bam | awk -F'\t' -f sq_lengths.awk
$1 == "@SQ" { sub("SN:", "", $2); sub("LN:", "", $3); print $2 "\t" $3 }
