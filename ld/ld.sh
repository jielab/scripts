#!/usr/bin/env bash
set -euo pipefail

# 🚩 Help
usage() {
    cat <<'EOF'
Modules:
    run / all  Validate -> chr22 pilot -> automatic D -> autosomes -> BED
    Default: run the complete workflow; stop if validation fails
    --chr LIST  Only compute these autosomes, e.g. 2, 2-21 or 2,3,6

Main options:
    OUT=DIR  Analysis output [/mnt/f/analysis/ld]
    VCF_DIR=DIR  Input VCFs [/mnt/f/gen/1kg/37/vcf]
    PANEL=FILE  Sample panel [/mnt/f/gen/1kg/37/samples.txt]
    PLINK2=FILE  PLINK2 executable [plink2; /mnt/d/software/bin/plink2 fallback]
    THREADS=N  PLINK2 threads [8]
    MEMORY_MB=N  PLINK2 memory [4096 MB]
    D_GRID_KB=LIST  Base D grid [250,500,1000,2000,5000,10000 kb]
    PILOT_MAX_D_KB=N  Extend unfinished chromosomes' grid up to N kb [20000 kb]
    CHRS=LIST  Target autosomes [1-22]; --chr takes priority
    Passed chromosomes keep their validated grid after input/output checks
    Excluding chr22 requires an existing valid pilot; no new chr22 computation
    LD_STORAGE=MODE  LD storage: blocked-zstd (default) or plain; identical computation
    LD_TMP_DIR=DIR  Temporary LD files [OUT/tmp unless configured]; removed after use or failure
    LD_RESUME_DIR=DIR  Reuse a preserved complete LD table; empty disables this checkpoint
    OUT/ld.runtime.env  Optional analysis-specific launch defaults; explicit environment wins
    MIN_FREE_GB=N  Keep this free-space reserve, in decimal GB [5]; LD tables need additional space
    -h, --help  Help; details: /mnt/d/scripts/ld/0.README.md

Example:
    ./ld.sh
    ./ld.sh --chr 2-21
    LD_TMP_DIR=/mnt/f/analysis/ld/tmp ./ld.sh
    PILOT_MAX_D_KB=40000 ./ld.sh
EOF
}

COMMAND=run
COMMAND_SET=0
CHR_ARG=''
while (( $# )); do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --chr|--chr=*)
            if [[ -n "$CHR_ARG" ]]; then
                printf '%s\n' '[LD] Specify --chr only once.' >&2; exit 2
            fi
            if [[ "$1" == --chr ]]; then
                if (( $# < 2 )) || [[ -z "$2" || "$2" == -* ]]; then
                    printf '%s\n' '[LD] --chr requires a list, e.g. --chr 2-21.' >&2; exit 2
                fi
                CHR_ARG="$2"
                shift
            else
                CHR_ARG="${1#--chr=}"
                if [[ -z "$CHR_ARG" ]]; then
                    printf '%s\n' '[LD] --chr requires a nonempty list.' >&2; exit 2
                fi
            fi
            ;;
        run|all)
            if (( COMMAND_SET )); then
                printf '%s\n' '[LD] Specify run or all only once.' >&2; exit 2
            fi
            COMMAND="$1"
            COMMAND_SET=1
            ;;
        *)
            printf '%s\n' '[LD] Use ./ld.sh [run|all] [--chr LIST]. D is selected automatically.' >&2
            exit 2
            ;;
    esac
    shift
done


# 🚩 Paths and common defaults
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
export OUT="${OUT:-/mnt/f/analysis/ld}"
# Keep this analysis's disk/checkpoint choices across new terminal sessions.
# The optional file supplies defaults only; caller-provided values take priority.
if [[ -f "$OUT/ld.runtime.env" ]]; then
    source "$OUT/ld.runtime.env"
fi
export VCF_DIR="${VCF_DIR:-/mnt/f/gen/1kg/37/vcf}"
export PANEL="${PANEL:-/mnt/f/gen/1kg/37/samples.txt}"
export PLINK2="${PLINK2:-plink2}"
export CXX="${CXX:-g++}"
export THREADS="${THREADS:-8}"
export MEMORY_MB="${MEMORY_MB:-4096}"
export D_GRID_KB="${D_GRID_KB:-250,500,1000,2000,5000,10000}"
export PILOT_MAX_D_KB="${PILOT_MAX_D_KB:-20000}"
export STABLE_CONFIRMATIONS="${STABLE_CONFIRMATIONS:-2}"
export REFERENCE_MAX_REL_DEVIATION="${REFERENCE_MAX_REL_DEVIATION:-0.10}"
export COINCIDENT_POLICY="${COINCIDENT_POLICY:-group}"
export LD_STORAGE="${LD_STORAGE:-blocked-zstd}"
export LD_TMP_DIR="${LD_TMP_DIR:-$OUT/tmp}"
export MIN_FREE_GB="${MIN_FREE_GB:-5}"
export CHRS="${CHR_ARG:-${CHRS:-1-22}}"
export LD_TARGETS="${LD_TARGETS:-$SCRIPT_DIR/reference/published_counts.tsv}"


# 🚩 Foreground execution
python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else "[LD] Python 3.11 or newer is required")'
exec python3 "$SCRIPT_DIR/f/ld_pipeline.py" "$COMMAND"
