#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Usage: ./2.herv_census.sh [OPTIONS]
    Variant-conditioned sequence reconstruction and RepeatMasker census for one sample.
    HG00096 represents one sample; its results do not establish population allele patterns.

Main options:
    --dat NAME --build 37|38 --sample ID [1kg/38/HG00096]
    --out DIR --pfile-dir DIR | --vcf-dir DIR --fasta FILE --library FILE
    --replace TRUE|FALSE               Recompute or reuse complete outputs [TRUE]
    --jobs N --cores N                 Parallel chunks / CPU budget [8 / 16]
    --submit auto|local|hpc --queue NAME --job-cores N --job-memory-gb N --job-hours N
    --job-dir DIR --write-cmd --wait TRUE|FALSE --plan-only
    --help-full                        Show all engine options

Can run alongside scan, PCA and GWAS; not a prerequisite for herv_enrichment.
Population census requires a stratified cohort plan and aggregation of all selected samples.

Examples:
    ./2.herv_census.sh --dat 1kg --build 38 --sample HG00096 \
      --out /work/sph-huangj/analysis/allele3/herv_census/1kg/38/HG00096 \
      --replace TRUE --jobs 8 --cores 16 \
      --submit hpc --queue 33 --job-cores 40 --job-memory-gb 160 --job-hours 120
PUBLIC_ENTRY_HELP
    exit 0
fi
# Independent individual-genome HERV census; use --help-full for parameters.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=herv_census
# shellcheck source=f/0.run.sh
source "$ROOT/f/0.run.sh"

usage() {
    bash "${BASH_SOURCE[0]}" --help
}

main() {
    local arg
    for arg in "$@"; do
        case "$arg" in -h|--help|help) usage; return 0 ;; esac
    done
    local config=${HERV_CONFIG:-"$ROOT/herv.env"}
    if [[ -f $config ]]; then
        set -a
        # shellcheck disable=SC1090
        source "$config"
        set +a
    fi
    local analysis_root=${HERV_ANALYSIS_ROOT:-${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}}
    local out=${HERV_OUT:-} gen=${HERV_GEN:-1kg} sample=${HERV_SAMPLE:-HG00096}
    local build=${HERV_BUILD:-38} pfile=${HERV_PFILE_DIR:-}
    local vcf=${HERV_VCF_DIR:-}
    local fasta=${HERV_FASTA:-}
    local jobs=${HERV_JOBS:-8} cores=${HERV_CORES:-16}
    local python=${HERV_PYTHON:-${ALLELE_PYTHON:-}}
    local memory=${HERV_MEMORY_LIMIT_GB:-${ALLELE_MEMORY_LIMIT_GB:-16}}
    local interval=${HERV_SCREEN_LOG_INTERVAL:-${ALLELE_SCREEN_LOG_INTERVAL:-300}}
    local verbose=0 dry_run=0 full_help=0 replace=${ALLELE_REPLACE:-TRUE} explicit_out=0
    local -a extra=() args=()
    [[ -z ${HERV_LIBRARY:-} ]] || extra+=(--library "$HERV_LIBRARY")
    [[ -z ${HERV_SCRATCH:-} ]] || extra+=(--scratch "$HERV_SCRATCH")
    while (($#)); do
        if [[ $1 == --*=* ]]; then
            arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"
        fi
        case "$1" in
            --dry-run) dry_run=1; shift ;;
            --help-full) full_help=1; shift ;;
            --verbose-screen) verbose=1; shift ;;
            --dir-analysis-root) need_value "$@"; analysis_root=$2; shift 2 ;;
            --out) need_value "$@"; out=$2; explicit_out=1; shift 2 ;;
            --dat) need_value "$@"; gen=$2; shift 2 ;;
            --build) need_value "$@"; build=$2; shift 2 ;;
            --replace) need_value "$@"; replace=${2^^}; shift 2 ;;
            --sample) need_value "$@"; sample=$2; shift 2 ;;
            --vcf-dir) need_value "$@"; vcf=$2; shift 2 ;;
            --pfile-dir) need_value "$@"; pfile=$2; shift 2 ;;
            --fasta) need_value "$@"; fasta=$2; shift 2 ;;
            --jobs) need_value "$@"; jobs=$2; shift 2 ;;
            --cores) need_value "$@"; cores=$2; shift 2 ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --input-format|--pgen-python|--library|--repeatmasker|--chunk-bp|--halo|--scratch)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            *) fail "Unknown option: $1 (see --help)" ;;
        esac
    done
    [[ $replace == TRUE || $replace == FALSE ]] || fail '--replace must be TRUE or FALSE'
    [[ $gen =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && $sample =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || fail '--gen and --sample must be single names (letters, digits, dots, underscores or hyphens)'
    [[ $build == 37 || $build == 38 ]] || fail '--build must be 37 or 38'
    export ALLELE_ACTIVE_HERV_BUILD=$build
    local sample_out=${DIR_DATA:-${DIR_BIG}}/gen/$gen/$build/herv/$sample
    if [[ ${submit:-local} == hpc ]]; then sample_out=$DIR_SCRATCH/allele3/herv/$gen/$build/$sample; fi
    if [[ -z $out ]]; then
        out=$analysis_root/herv_census/$gen/$build/$sample
    fi
    vcf=${vcf:-${DIR_DATA:-${DIR_BIG}}/gen/$gen/$build/vcf}
    pfile=${pfile:-${DIR_DATA:-${DIR_BIG}}/gen/$gen/$build/pfile}
    if [[ -z $fasta ]]; then
        if [[ $build == 37 ]]; then fasta=${ALLELE_FASTA37:-${DIR_DATA:-${DIR_BIG}}/gen/1hgp/37/GRCH37.fasta.gz}; else fasta=${ALLELE_FASTA:-${DIR_DATA:-${DIR_BIG}}/gen/1hgp/38/GRCH38.fasta.gz}; fi
    fi
    extra+=(--pfile-dir "$pfile")
    positive_integer --jobs "$jobs"
    positive_integer --cores "$cores"
    workflow_cpu_check "$jobs" 8
    if [[ ${submit:-local} == hpc ]]; then workflow_cpu_check "$cores" 16; fi
    local conda_root=${HERV_CONDA_ROOT:-${ALLELE_CONDA_ROOT:-${METAG_CONDA_ROOT:-$HOME/anaconda3}}}
    python=$(workflow_python "$python" "$conda_root/envs/herv-pgen/bin/python" \
        "$conda_root/envs/allele3-herv-pgen/bin/python" "$conda_root/envs/allele3-evolution/bin/python")
    args=(--out "$out" --replace "$replace" --gen "$gen" --sample "$sample" --fasta "$fasta"
          --vcf-dir "$vcf" --jobs "$jobs" --cores "$cores" "${extra[@]}")
    if ((full_help)); then args=(--help); fi
    printf '[HERV dataset] dat=%s; build=%s; sample=%s; PGEN=%s\n' "$gen" "$build" "$sample" "$pfile"
    local mode="FRESH RECALCULATION"
    [[ $replace == TRUE ]] || mode="REUSE COMPLETE OUTPUT OR REBUILD"
    printf '[mode] %s; replace=%s; output=%s\n' "$mode" "$replace" "$out"
    local entry=$ROOT/f/workflow_data_check.py
    local -a entry_args=(census)
    if ((full_help)); then entry=$ROOT/f/census.py; entry_args=(); fi
    export ALLELE_LOG_STYLE=herv WORKFLOW_CONSOLE=$ROOT/f/console.py HERV_SCREEN_OUT=$out
    workflow_launch herv "$python" "$entry" "$out" \
        "$memory" "$interval" "$verbose" "$dry_run" "$((1-full_help))" "${entry_args[@]}" "${args[@]}"
}

main "$@"

