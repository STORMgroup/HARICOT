#!/bin/bash
# cmp_run.sh REF NEW [NREADS]: every file of a recorded haricot assemble run REF against the run NEW, byte for byte (cmp), and the
# files NEW has that REF lacks. Skipped: logs (*.log, *.stdout), the step timings (time.tsv), scorer outputs a recorded
# run may hold (metrics.tsv, pieces.tsv) and the canonical dumps of checks/invariance.sh (work/canon/, sym_compare.py).
# With NREADS (the run's read count), a candidate cache (*.bin) of NEW that is REF's bytes followed by a bincode Vec<u64>
# of NREADS content hashes (8-byte length + 8 bytes per read: pipeline/src/cache.rs) counts
# as "same + hashes" instead of differing.
# Prints each differing / missing / added file and one summary line; exit 0 iff at least one file was compared and
# nothing differs or is missing; exit 2 on a usage error.
usage() { sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
case ${1:-} in -h|--help|'') usage ;; esac
ref=$1 new=${2:-} nreads=${3:-}
[ -d "$ref" ] && [ -d "$new" ] || { echo "cmp_run.sh: REF and NEW must be directories (-h: usage)" >&2; exit 2; }
n=0; same=0; hashed=0; bad=0; miss=0
skip() {
    case $1 in
        *.log|*.stdout|time.tsv|metrics.tsv|pieces.tsv|work/canon/*) return 0 ;;
    esac
    return 1
}
while IFS= read -r f; do
    skip "$f" && continue
    if [ ! -e "$new/$f" ]; then
        miss=$((miss + 1)); echo "  missing: $f"
        continue
    fi
    n=$((n + 1))
    if cmp -s "$ref/$f" "$new/$f"; then
        same=$((same + 1))
        continue
    fi
    if [ -n "$nreads" ] && [ "${f%.bin}" != "$f" ]; then
        rs=$(stat -c %s "$ref/$f"); ns=$(stat -c %s "$new/$f")
        if [ "$ns" -eq $((rs + 8 + 8 * nreads)) ] && cmp -s -n "$rs" "$ref/$f" "$new/$f"; then
            hashed=$((hashed + 1))
            continue
        fi
    fi
    bad=$((bad + 1)); echo "  DIFFERS: $f"
done < <(cd "$ref" && find . -type f | sed 's|^\./||' | sort)
added=$(comm -13 <(cd "$ref" && find . -type f | sed 's|^\./||' | sort) <(cd "$new" && find . -type f | sed 's|^\./||' | sort) \
    | grep -vE '\.log$|\.stdout$|^work/canon/' | tr '\n' ' ')
echo "  added (in NEW only; logs not listed): ${added:-none}"
echo "compared $n files: $same same, $hashed same + content-hash appendix, $bad differ, $miss missing ($new vs $ref)"
[ $bad = 0 ] && [ $miss = 0 ] && [ $n -gt 0 ]
