# Settings of the reproduction scripts; every step sources this file.
#
# R        the work directory: external tools, downloads, read sets, basecalls, references, truth and runs (several TB
#          for all datasets). Set it here once, or give it in the environment: R=/path/to/work bash reproduce/run.sh sd
# THREADS  threads per step (the manuscript's runs: 16)
R=${R:-}
THREADS=${THREADS:-16}
# the conda environments of envs/*.yml
ENV_HARICOT=haricot
ENV_ASSEMBLERS=haricot-assemblers
ENV_BANDAGE=haricot-bandage

# Nothing to set below this line.
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(dirname "$HERE")
BIN=$ROOT/target/release
case $R in /* | '') ;; *) R=$PWD/$R ;; esac
KEYS="cj lm vp1 vp2 se sd at2 at3 at5 dm2L dm2R dm3L dmX hg18"
GENOMES="ek_fwd sa_fwd hp2_fwd td2_fwd an2_fwd"
WTDBG2_URL=https://github.com/ruanjue/wtdbg2.git
WTDBG2_COMMIT=b77c565
RAWHASH2_URL=https://github.com/STORMgroup/RawHash2.git
RAWHASH2_COMMIT=92366af
ZENODO=https://zenodo.org/records/23963126/files
FORESTS_TGZ=haricot-forests-v1.tar.gz
FORESTS_SHA256=2b0c37d21d621539480cf29121184278f70f99308efb206edf1c3bca19e9d4f0
READ_LISTS_TGZ=haricot-read-lists-v1.tar.gz
READ_LISTS_SHA256=4c64266e21382a1a47d4e0c75a47b1d27c47b26ec3100a51fc9dd0b8f3b9bd92

say() { echo "$*"; }
die() { echo "${0##*/}: $*" >&2; exit 1; }
skip() { echo "$1: done, skipped ($2)"; }
usage() { sed -n '2,/^set -/p' "$0" | sed '$d; s/^# \{0,1\}//'; }

need_r() {  # the work directory must be set; it is created if absent
    [ -n "$R" ] || die "set the work directory R in reproduce/config.sh, or run R=/path/to/work bash $0 ..."
    mkdir -p "$R"
}

dataset() {  # dataset KEY: the dataset's row of datasets.tsv; SETS: its four read sets
    local row
    row=$(awk -F'\t' -v k="$1" 'NR > 1 && $2 == k' "$HERE/datasets.tsv")
    [ -n "$row" ] || die "unknown dataset '$1'; datasets: $KEYS"
    IFS=$'\t' read -r ID KEY SPECIES REF MAIN GSIZE RATE QUAST PAT SOURCE <<< "$row"
    SETS="${KEY}_fwd ${KEY}_fwd_c12 ${KEY}_fwd_c6 ${KEY}_both_top"
}

one_key() {  # one_key "$@": -h prints the script's usage; otherwise exactly one dataset key (dataset)
    case ${1:-} in -h | --help) usage; exit 0 ;; esac
    [ $# = 1 ] || { usage >&2; exit 2; }
    need_r
    dataset "$1"
    echo "== ${0##*/} $KEY ($ID, $SPECIES)  R=$R  $(date '+%F %T')"
}

env_prefix() {  # env_prefix NAME: the prefix of the conda environment NAME; nothing if there is none
    local c=${CONDA_EXE:-$(command -v conda || true)}
    [ -n "$c" ] || die "no conda: install it, or initialise it in this shell (conda init)"
    "$c" env list --json | grep -o '"/[^"]*"' | tr -d '"' | awk -F/ -v n="$1" '$NF == n && !seen++'
}

conda_env() {  # conda_env NAME: as env_prefix, but the environment must exist
    local p
    p=$(env_prefix "$1")
    [ -n "$p" ] && [ -d "$p/bin" ] || die "no conda environment $1: run reproduce/0_install.sh"
    echo "$p"
}

use_haricot() {  # as `conda activate haricot`: its bin first on PATH; PY = its Python, which every script gets.
                 # A step stacks another environment on it for one command: PATH="$(conda_env NAME)/bin:$PATH" CMD
    local p
    p=$(conda_env "$ENV_HARICOT")
    export PATH=$p/bin:$PATH
    PY=$p/bin/python
}

dorado_for() {  # dorado_for RATE: DV, DM, DO = dorado's version, model and options for that sample rate (Hz)
    case $1 in
        5000) DV=1.4.0; DM=dna_r10.4.1_e8.2_400bps_sup@v5.2.0; DO="--emit-moves --disable-read-splitting" ;;
        4000) DV=0.9.2; DM=dna_r10.4.1_e8.2_400bps_sup@v4.1.0; DO="--emit-moves" ;;
        *) die "no basecalling model for the sample rate $1 Hz" ;;
    esac
}

basecall() {  # basecall READS.pod5 RATE OUT.bam: dorado basecalls with move tables (unaligned BAM)
    dorado_for "$2"
    local dorado=$R/tools/dorado-$DV-linux-x64/bin/dorado
    [ -x "$dorado" ] && [ -d "$R/dorado_models/$DM" ] || die "no dorado $DV or its model: run reproduce/0_install.sh"
    nvidia-smi -L > /dev/null 2>&1 \
        || die "basecalling needs an NVIDIA GPU, and nvidia-smi finds none. Run this step on a GPU machine with the same R, or put basecalls made there at $3"
    mkdir -p "$(dirname "$3")"
    say "dorado $DV, $DM: $1 -> $3"
    # shellcheck disable=SC2086  # DO holds several options
    "$dorado" basecaller $DO "$R/dorado_models/$DM" "$1" > "$3.part"
    mv "$3.part" "$3"
}

missing_files() {  # missing_files SET_DIR: how many of the source files the read set needs are not in $R/downloads
    local s f b p n=0
    while IFS=$'\t' read -r s f b; do
        p=$R/downloads/$s/$f
        case $p in *.fast5) p=${p%.fast5}.pod5 ;; esac  # download.py converts fast5 files to POD5
        [ -e "$p" ] || n=$((n + 1))
    done < <("$PY" -I "$ROOT/datasets/download.py" "$R/downloads" --set "$1" --list)
    echo "$n"
}

get_sets() {  # get_sets ROLE BASE SET...: each SET's POD5 in $R/sets, from the source files of BASE, the set of
              # read_lists/ROLE that holds the reads of every SET (download.py, build_set.py); a POD5 present is skipped
    local role=$1 base=$2 s src n todo=()
    shift 2
    for s in "$@"; do
        if [ -s "$R/sets/$s.pod5" ]; then skip "$s" "$R/sets/$s.pod5"; else todo+=("$s"); fi
    done
    [ ${#todo[@]} -gt 0 ] || return 0
    [ -d "$ROOT/read_lists/$role/$base" ] || die "no read list read_lists/$role/$base: run reproduce/0_install.sh"
    n=$(missing_files "$ROOT/read_lists/$role/$base")
    src=$(awk -F'\t' 'NR == 2 {print $2}' "$ROOT/read_lists/$role/$base/index.tsv")
    if [ "$n" = 0 ]; then
        say "$base: every source file is in $R/downloads"
    elif [ "$(awk -F'\t' -v s="$src" '$1 == s {print $2}' "$ROOT/datasets/sources.tsv")" = tgz ]; then
        # a gzip archive is read from its start for any member: fetch it whole, once, and check its MD5
        "$PY" -I "$ROOT/datasets/download.py" "$R/downloads" --source "$src"
    else
        "$PY" -I "$ROOT/datasets/download.py" "$R/downloads" --set "$ROOT/read_lists/$role/$base"
    fi
    mkdir -p "$R/sets"
    for s in "${todo[@]}"; do
        "$PY" -I "$ROOT/datasets/build_set.py" "$ROOT/read_lists/$role/$s" "$R/downloads" -o "$R/sets/$s.pod5"
    done
}

get_refs() {  # get_refs FILE...: the reference files in $R/refs (references.py); a file present is skipped
    local f todo=()
    for f in "$@"; do
        if [ -s "$R/refs/$f" ]; then skip "$f" "$R/refs/$f"; else todo+=("$f"); fi
    done
    [ ${#todo[@]} -gt 0 ] || return 0
    "$PY" -I "$ROOT/datasets/references.py" "$R/refs" --files "${todo[@]}"
}
