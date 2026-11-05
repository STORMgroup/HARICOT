# GFA S lines with a sequence -> FASTA (segments without a sequence, "*", are skipped).
# Usage: awk -F'\t' -f gfa_segs.awk G.gfa > segs.fa
$1 == "S" && $3 != "*" { print ">" $2; print $3 }
