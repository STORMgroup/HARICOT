# Density-reducing pre-filter for a Rawsamble r10ava PAF (grouped by query).
# Two passes over ovl.paf:
#   Pass 1: harvest name -> event-length (col2) for the col7 relabel.
#   Pass 2: drop no-hit '*' rows; drop overlaps with qspan<MINSPAN or
#           col10<MINMATCH; relabel col7 (target len) to event space
#           (fix_paf_tlen semantics); then per QUERY read keep only the TOPN
#           highest-scoring overlaps (score = s1:i: chain score, else col10).
# Nothing is silently truncated: END prints all drop counters to stderr.
#   usage: awk -v MINSPAN=.. -v MINMATCH=.. -v TOPN=.. -f filter_ovl.awk ovl.paf ovl.paf > filt.paf
BEGIN {
    FS = OFS = "\t"
    MINSPAN  = MINSPAN  + 0
    MINMATCH = MINMATCH + 0
    TOPN     = TOPN     + 0
}
# Pass 1
NR == FNR {
    if ($1 != "" && $2 ~ /^[0-9]+$/) elen[$1] = $2
    next
}
# Pass 2
{
    if ($3 == "*") { nohit++; next }
    qs = $4 - $3
    nm = $10 + 0
    if (qs < MINSPAN)  { drop_span++;  next }
    if (nm < MINMATCH) { drop_match++; next }
    if (!($6 in elen)) { drop_unknown++; next }
    $7 = elen[$6]
    sc = nm
    for (i = 13; i <= NF; i++) if ($i ~ /^s1:i:/) { split($i, a, ":"); sc = a[3] + 0 }
    q = $1
    if (q != cur) { flush(); cur = q; nb = 0 }
    nb++; buf[nb] = $0; scr[nb] = sc
}
END { flush()
    printf("filter_ovl: kept=%d groups=%d capped=%d | dropped nohit=%d span=%d match=%d unknown=%d\n",
        kept_tot+0, ngroups+0, capped+0, nohit+0, drop_span+0, drop_match+0, drop_unknown+0) > "/dev/stderr"
}
function flush(   i, key, sorted, m, kept, parts, idx) {
    if (nb == 0) return
    ngroups++
    if (TOPN <= 0 || nb <= TOPN) {
        for (i = 1; i <= nb; i++) { print buf[i]; kept_tot++ }
        return
    }
    capped++
    for (i = 1; i <= nb; i++) key[i] = sprintf("%012d|%06d", scr[i], i)
    m = asort(key, sorted)          # ascending by zero-padded score
    kept = 0
    for (i = m; i >= 1 && kept < TOPN; i--) {
        split(sorted[i], parts, "|")
        idx = parts[2] + 0
        print buf[idx]
        kept++; kept_tot++
    }
    delete key; delete sorted
}
