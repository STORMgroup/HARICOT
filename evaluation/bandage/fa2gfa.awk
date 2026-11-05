# FASTA -> a links-free GFA, one segment per record (for Bandage panels of assemblers that write no GFA).
# Usage: awk -f fa2gfa.awk asm.fa > asm.gfa
BEGIN { print "H\tVN:Z:1.0" }
/^>/ { if (n != "") print "S\t" n "\t" s; n = substr($1, 2); s = ""; next }
{ s = s $0 }
END { if (n != "") print "S\t" n "\t" s }
