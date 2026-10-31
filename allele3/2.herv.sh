#!/usr/bin/env bash
# Independent individual-genome HERV census. See 0.README.md for dependencies.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=herv
# shellcheck source=f/0.run.sh
source "$ROOT/f/0.run.sh"

usage() {
    cat <<'HELP'
Usage: ./2.herv.sh [OPTIONS]

Reconstruct phased genomes and annotate HERVs with RepeatMasker; no subcommands.
Default: recalculate into /mnt/f/gen/DAT/BUILD/herv/SAMPLE.

Options:
  --out DIR              Override the dataset/build/sample output directory
  --dat NAME, --gen NAME Dataset [1kg]
  --build 37|38          Reference build for inputs and output [38]
  --sample ID            Individual [HG00096]
  --replace TRUE|FALSE   Recalculate or reuse checkpoints [TRUE]
  --reuse-from DIR       Reuse a completed cache with --replace FALSE
  --input-format auto|vcf|pgen  Input format [auto]
  --dir-vcf DIR / --dir-pfile DIR  Genotype input directories
  --fasta FILE / --library FILE   Reference genome / repeat library
  --jobs N / --cores N   Concurrent chunks / CPU limit [8 / 16]
  --screen-log-interval N  Progress interval in seconds [300]; full logs in OUT/logs
  --verbose-screen      Show the full workload output on screen
  --dry-run              Print the command only
  --help-full            Show detailed options
  -h, --help             Show this help

Examples:
  ./2.herv.sh [OPTIONS]
  ./2.herv.sh --out /mnt/f/gen/1kg/38/herv/HG00096
HELP
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
    # Legacy environment names remain aliases; the main workflow configuration is never loaded here.
    local analysis_root=${HERV_ANALYSIS_ROOT:-}
    local out=${HERV_OUT:-${ALLELE_HERV_OUT:-}} gen=${HERV_GEN:-1kg} sample=${HERV_SAMPLE:-HG00096}
    local build=${HERV_BUILD:-38} pfile=${HERV_PFILE_DIR:-}
    local vcf=${HERV_VCF_DIR:-${ALLELE_HERV_VCF_DIR:-}}
    local fasta=${HERV_FASTA:-${ALLELE_FASTA:-}}
    local jobs=${HERV_JOBS:-${ALLELE_HERV_JOBS:-8}} cores=${HERV_CORES:-${ALLELE_HERV_CORES:-16}}
    local python=${HERV_PYTHON:-${ALLELE_PYTHON:-}}
    local memory=${HERV_MEMORY_LIMIT_GB:-${ALLELE_MEMORY_LIMIT_GB:-16}}
    local interval=${HERV_SCREEN_LOG_INTERVAL:-${ALLELE_SCREEN_LOG_INTERVAL:-300}}
    local verbose=0 dry_run=0 full_help=0 replace=${HERV_REPLACE:-TRUE} explicit_out=0 explicit_reuse=0 no_reuse=0
    local -a extra=() args=()
    local reuse=${HERV_REUSE_DIR:-}
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
            --dir-out|--out) need_value "$@"; out=$2; explicit_out=1; shift 2 ;;
            --gen|--dat) need_value "$@"; gen=$2; shift 2 ;;
            --build) need_value "$@"; build=$2; shift 2 ;;
            --reuse-from) need_value "$@"; reuse=$2; explicit_reuse=1; shift 2 ;;
            --replace) need_value "$@"; replace=${2^^}; shift 2 ;;
            --no-reuse) reuse=''; no_reuse=1; shift ;;
            --sample) need_value "$@"; sample=$2; shift 2 ;;
            --dir-vcf|--vcf-dir) need_value "$@"; vcf=$2; shift 2 ;;
            --dir-pfile|--pfile-dir) need_value "$@"; pfile=$2; shift 2 ;;
            --fasta|--fasta38) need_value "$@"; fasta=$2; shift 2 ;;
            --jobs) need_value "$@"; jobs=$2; shift 2 ;;
            --cores) need_value "$@"; cores=$2; shift 2 ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --input-format|--pgen-python|--library|--repeatmasker|--chunk-bp|--halo|--scratch|--dir-scratch)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            *) fail "Unknown option: $1 (see --help)" ;;
        esac
    done
    [[ $replace == TRUE || $replace == FALSE ]] || fail '--replace must be TRUE or FALSE'
    if [[ $replace == TRUE && $full_help == 0 ]]; then
        ((explicit_reuse == 0)) || fail 'Use --replace FALSE with --reuse-from; TRUE requests recalculation'
        reuse=''
    fi
    [[ $gen =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && $sample =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || fail '--gen and --sample must be single names (letters, digits, dots, underscores or hyphens)'
    [[ $build == 37 || $build == 38 ]] || fail '--build must be 37 or 38'
    local sample_out=/mnt/f/gen/$gen/$build/herv/$sample
    if [[ -z $out ]]; then
        if [[ -n $analysis_root ]]; then out=$analysis_root/herv/$gen/$sample; else out=$sample_out; fi
    fi
    vcf=${vcf:-/mnt/f/gen/$gen/$build/vcf}
    pfile=${pfile:-/mnt/f/gen/$gen/$build/pfile}
    fasta=${fasta:-/mnt/f/gen/1hgp/$build/GRCH$build.fasta.gz}
    extra+=(--pfile-dir "$pfile")
    # The canonical sample directory can be recalculated or resumed in place.
    # Keep other cache sources protected; do not attach a directory to itself.
    local protected_cache=${HERV_REUSE_DIR:-$sample_out}
    if [[ $(readlink -m -- "$protected_cache") == $(readlink -m -- "$out") ]]; then protected_cache=''; fi
    export HERV_REUSE_DIR=$protected_cache
    if [[ $replace == FALSE && $explicit_reuse == 0 && $no_reuse == 0 && -z $reuse ]]; then reuse=$protected_cache; fi
    if [[ -n $reuse && $(readlink -m -- "$reuse") == $(readlink -m -- "$out") ]]; then reuse=''; fi
    positive_integer --jobs "$jobs"
    positive_integer --cores "$cores"
    ((jobs <= 8)) || fail 'Use at most 8 workers on this shared host'
    local conda_root=${HERV_CONDA_ROOT:-${ALLELE_CONDA_ROOT:-${METAG_CONDA_ROOT:-$HOME/anaconda3}}}
    python=$(workflow_python "$python" "$conda_root/envs/herv-pgen/bin/python" \
        "$conda_root/envs/allele3-herv-pgen/bin/python" "$conda_root/envs/allele3-evolution/bin/python")
    args=(--out "$out" --replace "$replace" --gen "$gen" --sample "$sample" --fasta "$fasta"
          --vcf-dir "$vcf" --jobs "$jobs" --cores "$cores" "${extra[@]}")
    [[ -z $reuse ]] || args+=(--reuse-from "$reuse")
    if ((full_help)); then args=(--help); fi
    printf '[HERV dataset] dat=%s; build=%s; sample=%s; PGEN=%s\n' "$gen" "$build" "$sample" "$pfile"
    local mode="FRESH RECALCULATION"
    [[ $replace == TRUE ]] || mode="CHECKPOINT RESUME / CACHE REUSE"
    printf '[mode] %s; replace=%s; output=%s\n' "$mode" "$replace" "$out"
    local entry=$ROOT/runtime/herv_runtime.py
    if ((full_help)); then entry=$ROOT/f/1.herv.py; fi
    export WORKFLOW_CONSOLE=$ROOT/runtime/herv_console.py HERV_SCREEN_OUT=$out
    workflow_launch herv "$python" "$entry" "$out" \
        "$memory" "$interval" "$verbose" "$dry_run" "$((1-full_help))" "${args[@]}"
}

main "$@"
