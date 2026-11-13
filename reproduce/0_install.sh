#!/bin/bash
# Step 0: HARICOT, the other tools at the manuscript's versions, the published forests and the read lists.
#
#   bash reproduce/0_install.sh
#
# Builds HARICOT (cargo build --release --locked -> target/release/). Into $R/tools: wtdbg2 (b77c565) and RawHash2
# (Rawsamble; 92366af), built from source, and dorado 1.4.0 and 0.9.2; into $R/dorado_models: their basecalling
# models. Creates the conda environments haricot, haricot-assemblers and haricot-bandage from envs/*.yml. Downloads
# the forests and the read lists from Zenodo (10.5281/zenodo.23963126), checks the archives' SHA-256 and unpacks them
# in the checkout: models/ (the forests haricot assemble reads) and read_lists/.
# Each part installed is skipped, so a rerun completes an interrupted installation.
# Needs: git, curl, make, a C/C++ compiler, CMake >= 3.16, GNU time at /usr/bin/time, conda, rustup.
set -euo pipefail
. "$(dirname "$0")/config.sh"
case ${1:-} in -h | --help) usage; exit 0 ;; '') ;; *) usage >&2; exit 2 ;; esac
need_r
for x in git curl make cmake cargo; do
    command -v "$x" > /dev/null || die "no $x on PATH (requirements: reproduce/README.md)"
done
[ -x /usr/bin/time ] || die "no GNU time at /usr/bin/time"
CONDA=${CONDA_EXE:-$(command -v conda || true)}
[ -n "$CONDA" ] || die "no conda: install it, or initialise it in this shell (conda init)"
echo "== ${0##*/}  R=$R  $(date '+%F %T')"
mkdir -p "$R/tools"

say "-- HARICOT"
(cd "$ROOT" && cargo build --release --locked)

at_commit() {  # at_commit DIR COMMIT: DIR is a git checkout at COMMIT
    [ -d "$1/.git" ] && [ "$(git -C "$1" rev-parse HEAD)" = "$(git -C "$1" rev-parse --verify -q "$2^{commit}")" ]
}
d=$R/tools/wtdbg2
if at_commit "$d" "$WTDBG2_COMMIT" && [ -x "$d/wtdbg2" ] && [ -x "$d/wtpoa-cns" ]; then
    skip "wtdbg2 $WTDBG2_COMMIT" "$d"
else
    say "-- wtdbg2 $WTDBG2_COMMIT"
    [ -d "$d/.git" ] || git clone -q "$WTDBG2_URL" "$d"
    git -C "$d" checkout -q "$WTDBG2_COMMIT"
    make -C "$d"
fi
d=$R/tools/RawHash2
if at_commit "$d" "$RAWHASH2_COMMIT" && [ -x "$d/bin/rawhash2" ]; then
    skip "RawHash2 $RAWHASH2_COMMIT" "$d"
else
    say "-- RawHash2 $RAWHASH2_COMMIT"
    [ -d "$d/.git" ] || git clone -q "$RAWHASH2_URL" "$d"
    git -C "$d" checkout -q "$RAWHASH2_COMMIT"
    git -C "$d" submodule update --init --recursive
    make -C "$d" cmake    # also downloads its POD5 library
fi

mkdir -p "$R/dorado_models"   # dorado download fails without a message if the directory does not exist
for rate in 5000 4000; do
    dorado_for "$rate"
    d=$R/tools/dorado-$DV-linux-x64
    if [ -x "$d/bin/dorado" ]; then
        skip "dorado $DV" "$d"
    else
        say "-- dorado $DV"
        rm -rf "$d.part"
        mkdir -p "$d.part"
        curl -fsSL --retry 5 "https://cdn.oxfordnanoportal.com/software/analysis/dorado-$DV-linux-x64.tar.gz" \
            | tar -xz -C "$d.part"
        mv "$d.part/dorado-$DV-linux-x64" "$d"
        rmdir "$d.part"
    fi
    if [ -d "$R/dorado_models/$DM" ]; then
        skip "$DM" "$R/dorado_models/$DM"
    else
        say "-- $DM"
        "$d/bin/dorado" download --model "$DM" --models-directory "$R/dorado_models"
        [ -d "$R/dorado_models/$DM" ] || die "dorado download did not write $R/dorado_models/$DM"
    fi
done

for e in "$ENV_HARICOT haricot" "$ENV_ASSEMBLERS assemblers" "$ENV_BANDAGE bandage"; do
    read -r name yml <<< "$e"
    p=$(env_prefix "$name")
    if [ -n "$p" ]; then
        skip "conda environment $name" "$p"
    else
        say "-- conda environment $name (envs/$yml.yml)"
        "$CONDA" env create -q -n "$name" -f "$ROOT/envs/$yml.yml"
    fi
done

zenodo() {  # zenodo ARCHIVE SHA256 DIR: ARCHIVE from the Zenodo record, checked and unpacked into the checkout's DIR
    local a=$ROOT/$1 sum=$2 d=$ROOT/$3
    if [ -s "$d/SHA256SUMS" ] && (cd "$d" && sha256sum -c --quiet SHA256SUMS > /dev/null 2>&1); then
        skip "$3/ (every file matches its SHA256SUMS)" "$d"
        return
    fi
    if ! echo "$sum  $a" | sha256sum -c --status 2> /dev/null; then
        say "-- $1"
        curl -fsSL --retry 5 -o "$a.part" "$ZENODO/$1?download=1"
        mv "$a.part" "$a"
        echo "$sum  $a" | sha256sum -c --status || die "$a: its SHA-256 differs from the published one"
    fi
    tar -xzf "$a" -C "$ROOT"
    (cd "$d" && sha256sum -c --quiet SHA256SUMS) || die "$d: files differ from its SHA256SUMS"
    say "$3/: unpacked from $1 (SHA-256 checked)"
}
zenodo "$FORESTS_TGZ" "$FORESTS_SHA256" models
zenodo "$READ_LISTS_TGZ" "$READ_LISTS_SHA256" read_lists
say "== ${0##*/}: done"
