#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/f/scratch.sh"

#  Help
usage() {
    cat <<'EOF'
Modules:
    run / all  Write chrN.cmd per chromosome and submit with bsub (default)
    --write-cmd  Only write .cmd; do not submit
    --execute  Run analysis inside the allocated job (internal entry)
    --chr LIST  Only compute these autosomes, e.g. 2, 2-21 or 2,3,6
    --check-env  Print Python version, paths and selected chromosomes; do not compute

Main options:
    --date-scratch YYYY-MM-DD  Scratch analysis date [2026-11-03]; retain the same date to resume
    --grch 37|38  Coordinate build [38]; also accepts grch=38 or GRCH=38 environment
    --dir-out DIR --dir-gen DIR --sample-panel FILE  Override output, PFILE directory and panel
    INPUT_FORMAT=pgen|vcf  Input format [pgen]
    OUT=DIR  Analysis output [/scratch/DATE/sph-huangj/analysis/ld/38]
    PFILE_DIR=DIR  chrN.pgen + chrN.pvar[.zst] + chrN.psam [/work/sph-huangj/data/gen/1kg/38/pfile]
    VCF_DIR=DIR  Only used with INPUT_FORMAT=vcf
    PANEL=FILE  Sample panel [/work/sph-huangj/data/gen/1kg/38/samples.txt]; all four-population samples
    PLINK2=FILE  PLINK2 executable [/work/sph-huangj/software/bin/plink2]
    THREADS=N  PLINK2 threads [4]
    MEMORY_MB=N  PLINK2 memory [8192 MB]
    D_GRID_KB=LIST  Base D grid [250,500,1000,2000,5000,10000 kb]
    MAX_D_KB=N  Initial D budget [40000 kb]; unstable chromosomes extend automatically
    AUTO_EXTEND_D=TRUE|FALSE  Double D up to observed chromosome span if needed [TRUE]
    REFERENCE_COUNT_POLICY=MODE  off for GRCh38; warn/error allowed for original GRCh37 cohort
    REFERENCE_MAX_REL_DEVIATION=N  Count comparison threshold [0.10; not a statistical test]
    CHRS=LIST  Target autosomes [1-22]; --chr takes priority
    Passed chromosomes keep their validated grid after input/output checks
    LD_STORAGE=MODE  LD storage: blocked-zstd (default) or plain; identical computation
    LD_TMP_DIR=DIR  Temporary LD files [OUT/tmp unless configured]; removed after use or failure
    LD_RESUME_DIR=DIR  Reuse a preserved complete LD table; empty disables this checkpoint
    OUT/ld.runtime.env  Optional analysis-specific launch defaults; explicit environment wins
    OUT/rebuild_v2/logs/  Run status/STOP, partition errors and native PLINK logs
    MIN_FREE_GB=N  Filesystem/group-quota reserve, decimal GB [100]
    LD_QUEUE=NAME  Submission queue [HPC_CPU_QUEUE, normally 33]
    -h, --help  Show options and examples below beside this script

Usage examples:
    cd /work/sph-huangj/scripts/ld

    # Fresh GRCh38: all autosomes plus automatic dependent merge.
    ./ld.sh --grch 38 --chr 1-22

    # Alternative queue only; do not submit both copies.
    LD_QUEUE=38 ./ld.sh --grch 38 --chr 1-22

    # Fresh QiMing run; inputs and new outputs follow 0f/0os.sh.
    ./ld.sh --grch 38 --chr 1-22 --date-scratch 2026-11-03

    # Inspect paths without submitting or computing.
    ./ld.sh --check-env --grch 38 --chr 1-22

    # Write command files only; no submission.
    ./ld.sh --grch 38 --chr 1-22 --write-cmd

    # Default output: /scratch/2026-11-03/sph-huangj/analysis/ld/38
    # Each chromosome uses chrN/ and compressed temporary files in chrN/tmp/.
    # Inputs: /work/sph-huangj/data/gen/1kg/38/pfile (never removed).
    # Use --date-scratch YYYY-MM-DD only to select a different scratch date.
    # Merge waits for done("ld38.step1.chr*") and requires chr1-22 completion.
    # Empty bjobs is not proof of success: check completion records and logs.

EOF
}

COMMAND=run
COMMAND_SET=0
CHR_ARG=''
GRCH_ARG=''
CHECK_ENV=0
EXECUTE=0
WRITE_ONLY=0
while (( $# )); do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --dir-out|--dir-gen|--sample-panel)
            [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "[LD] $1 requires a path" >&2; exit 2; }
            case $1 in
                --dir-out) export OUT=$2 ;;
                --dir-gen) export PFILE_DIR=$2 ;;
                --sample-panel) export PANEL=$2 ;;
            esac
            shift ;;
        --check-env) CHECK_ENV=1 ;;
        --write-cmd) WRITE_ONLY=1 ;;
        --execute) EXECUTE=1 ;;
        --grch|--grch=*|grch=*)
            if [[ -n "$GRCH_ARG" ]]; then
                printf '%s\n' '[LD] Specify --grch only once.' >&2; exit 2
            fi
            if [[ "$1" == --grch ]]; then
                if (( $# < 2 )); then
                    printf '%s\n' '[LD] --grch requires 37 or 38.' >&2; exit 2
                fi
                GRCH_ARG="$2"
                shift
            else
                GRCH_ARG="${1#*=}"
            fi
            if [[ "$GRCH_ARG" != 37 && "$GRCH_ARG" != 38 ]]; then
                printf '%s\n' '[LD] --grch must be 37 or 38.' >&2; exit 2
            fi
            ;;
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
            printf '%s\n' '[LD] Use ./ld.sh [run|all] [--grch 37|38] [--chr LIST]. D is selected automatically.' >&2
            exit 2
            ;;
    esac
    shift
done


#  Paths and common defaults
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
export GENOME_BUILD="${GRCH_ARG:-${GRCH:-${GENOME_BUILD:-38}}}"
if [[ "$GENOME_BUILD" != 37 && "$GENOME_BUILD" != 38 ]]; then
    printf '%s\n' '[LD] GRCH must be 37 or 38.' >&2; exit 2
fi
export GRCH="$GENOME_BUILD"
source "$SCRIPT_DIR/f/env.sh"
if [[ -f "$OUT/ld.runtime.env" ]]; then
    # Source defaults, then restore explicitly resolved caller settings.
    declare -A ld_saved=()
    for ld_key in OUT PFILE_DIR PANEL VCF_DIR D_GRID_KB MAX_D_KB AUTO_EXTEND_D STABLE_CONFIRMATIONS LD_TMP_DIR; do
        [[ ! -v $ld_key ]] || ld_saved[$ld_key]=${!ld_key}
    done
    source "$OUT/ld.runtime.env"
    for ld_key in "${!ld_saved[@]}"; do export "$ld_key=${ld_saved[$ld_key]}"; done
    unset ld_saved ld_key
fi
export D_GRID_KB="${D_GRID_KB:-250,500,1000,2000,5000,10000}"
export MAX_D_KB="${MAX_D_KB:-${PILOT_MAX_D_KB:-40000}}"
export STABLE_CONFIRMATIONS="${STABLE_CONFIRMATIONS:-2}"
export REFERENCE_MAX_REL_DEVIATION="${REFERENCE_MAX_REL_DEVIATION:-0.10}"
export REFERENCE_COUNT_POLICY="${REFERENCE_COUNT_POLICY:-off}"
export COINCIDENT_POLICY="${COINCIDENT_POLICY:-group}"
export LD_STORAGE="${LD_STORAGE:-blocked-zstd}"
export LD_TMP_DIR="${LD_TMP_DIR:-$OUT/tmp}"
export MIN_FREE_GB="${MIN_FREE_GB:-100}"
export CHRS="${CHR_ARG:-${CHRS:-1-22}}"
export LD_TARGETS="${LD_TARGETS:-}"


#  Foreground execution
if (( CHECK_ENV )); then
    "$PYTHON" -c 'import json,os,sys; print(json.dumps(dict(python=sys.executable, version=sys.version.split()[0], **{k:os.environ.get(k) for k in ("GRCH","GENOME_BUILD","OUT","PFILE_DIR","PANEL","CHRS")}), indent=2))'
    exit 0
fi
"$PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else "[LD] Python 3.11 or newer is required")'
if (( ! EXECUTE )); then
    exec bash "$SCRIPT_DIR/f/jobs.sh" "$WRITE_ONLY"
fi
exec "$PYTHON" "$SCRIPT_DIR/f/ld_pipeline.py" "$COMMAND"
