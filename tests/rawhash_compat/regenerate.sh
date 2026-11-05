#!/usr/bin/env bash
# Regenerate the RawHash2 golden hash fixtures consumed by tests/hash_compat.rs.
#
# Run this ONCE (and re-run only if the signal or the RawHash2 sketch code
# changes). It compiles a tiny C harness against RawHash2's reference sketch
# implementation, generates a deterministic random signal, and dumps the
# resulting hash values per parameter configuration. The committed outputs
# (signal.txt + *.tsv) are what `cargo test` actually checks against, so this
# script is NOT part of the test run and needs no special build deps beyond a
# C compiler.
#
# Usage:  tests/rawhash_compat/regenerate.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
rh="$repo/extern/RawHash2/src"
cc="${CC:-cc}"

if [[ ! -f "$rh/rsketch.c" ]]; then
  echo "error: RawHash2 source not found at $rh (clone RawHash2 into extern/RawHash2)" >&2
  exit 1
fi

bin="$(mktemp -d)/dump_hashes"
echo "building harness -> $bin"
"$cc" -O2 -I "$rh" \
  "$here/dump_hashes.c" "$rh/rsketch.c" "$rh/kalloc.c" \
  -lm -o "$bin"

# 1. Deterministic random signal (2000 values, fixed seed). Spans the low /
#    fine / high quantization regions so every branch of dynamic_quantize and
#    the sig-diff collapse path get exercised.
echo "generating signal.txt"
"$bin" gen 2000 88172645463325252 > "$here/signal.txt"

# 2. One golden per parameter configuration. Columns: e w quant_bit n_buckets
#    k diff fine_min fine_max fine_range. Names match CONFIGS in hash_compat.rs.
#    All use w=0 (dense/reg sketch) so seeds come out in event order, giving an
#    exact, position-by-position comparison. fine_*/quant follow the R10.4.1
#    preset (src/params.rs).
emit() { # <name> <args...>
  local name="$1"; shift
  echo "  -> $name.tsv"
  "$bin" dump "$here/signal.txt" "$@" > "$here/$name.tsv"
}

#     name             e  w  qbit nbkt k  diff fmin fmax frange
emit "e8_q3_diff0"     8  0  3    8    9  0    -2   2    0.75   # R10.4.1 default
emit "e8_q3_nofilter"  8  0  3    8    9  -1   -2   2    0.75   # diff disabled
emit "e4_q3_nofilter"  4  0  3    8    9  -1   -2   2    0.75   # shorter e-mer
emit "e6_q3_diff0"     6  0  3    8    9  0    -2   2    0.75   # mid e + collapse

echo "done. fixtures written under $here"
