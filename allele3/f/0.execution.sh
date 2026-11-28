#!/usr/bin/env bash
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case ${1:-} in
        -h|--help|help)
            echo 'Internal runtime setup for the allele3 workflows.'
            exit 0 ;;
    esac
fi
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    echo 'This file is a sourced library. Use a root workflow .sh; see -h.' >&2
    exit 2
fi
# Scheduler/path bootstrap. No genomic file is opened on a submission host.

allele_execution_init() {
    local os_config=${ALLELE_OS_CONFIG:-$ROOT/../0f/0os.sh} arg
    local help=0
    for arg in "$@"; do
        case "$arg" in -h|--help|help|--help-full|--version) help=1 ;; esac
    done
    if [[ ! -f $os_config ]]; then
        if ((help)); then ALLELE_EXECUTION_ARGS=("$@"); return 0; fi
        fail "Shared OS configuration missing: $os_config. Install beside scripts/0f/0os.sh, or set ALLELE_OS_CONFIG."
    fi
    local -a expanded=()
    for arg in "$@"; do
        if [[ $arg == --*=* ]]; then expanded+=("${arg%%=*}" "${arg#*=}"); else expanded+=("$arg"); fi
    done
    # Resolve the machine before sourcing the shared path table.
    local i chosen_mode=${submit:-${SUBMIT:-auto}}
    for ((i=0;i<${#expanded[@]};i++)); do
        case "${expanded[i]}" in
            --submit) [[ $((i+1)) -lt ${#expanded[@]} ]] || fail '--submit needs auto, local or hpc'; chosen_mode=${expanded[i+1]} ;;
        esac
    done
    case "$chosen_mode" in auto|local|hpc) ;; *) fail '--submit needs auto, local or hpc' ;; esac
    # 0os.sh recomputes tracked defaults while preserving explicit root overrides.
    export submit=$chosen_mode SUBMIT=$chosen_mode
    # The shared file owns the OS roots. No separate OS table is maintained.
    source "$os_config"
    workflow_args "${expanded[@]}"
    set -- "${WORKFLOW_ARGS[@]}"
    export ALLELE_REPLACE=${ALLELE_REPLACE:-FALSE}
    ALLELE_EXECUTION_ARGS=()
    while (($#)); do
        case "$1" in
            --dir-out-allele3_scan|--dir-out-gnomad_detailed|--dir-out-avi_scan|--dir-out-herv_scan)
                [[ $# -ge 2 && -n $2 && $2 != --* ]] || fail "$1 requires a directory"
                local product_dir=$2
                [[ $product_dir == /* ]] || product_dir=$PWD/$product_dir
                case "$1" in
                    --dir-out-allele3_scan) export ALLELE_SUMMARY_STORE=$product_dir ;;
                    --dir-out-gnomad_detailed) export ALLELE_SCAN_STORE=$product_dir ;;
                    --dir-out-avi_scan) export ALLELE_AVI_STORE=$product_dir ;;
                    --dir-out-herv_scan) export ALLELE_HERV_ARCHIVE=$product_dir ;;
                esac
                shift 2 ;;
            --cores)
                case "${0##*/}" in
                    allele3.sh|allele3.paper.sh|3.allele3_scan.sh) shift 2 ;;
                    *) ALLELE_EXECUTION_ARGS+=("$1" "$2"); shift 2 ;;
                esac ;;
            --grch)
                [[ $# -ge 2 && $2 == 38 ]] || fail '--grch requires 38 for this workflow'
                export ALLELE_GRCH=38; shift 2 ;;
            --replace)
                [[ $# -ge 2 && ( $2 == TRUE || $2 == FALSE ) ]] || fail '--replace requires TRUE or FALSE'
                export ALLELE_REPLACE=$2; shift 2 ;;
            --write-cmd) export ALLELE_WRITE_CMD=1; shift ;;
            --retry-failed|--dry-run|--preflight)
                fail "$1 was removed; use --plan-only to inspect or rerun with --replace FALSE" ;;
            --queue|--job-cores|--job-memory-gb|--job-hours|--max-running|--job-dir|--lsf-memory-scope)
                [[ $# -ge 2 && -n $2 && $2 != --* ]] || fail "$1 requires a value"
                case "$1" in
                    --queue) export ALLELE_JOB_QUEUE=$2 ;;
                    --job-cores) export ALLELE_JOB_CORES=$2 ;;
                    --job-memory-gb) export ALLELE_JOB_MEMORY_GIB=$2 ;;
                    --job-hours) export ALLELE_JOB_HOURS=$2 ;;
                    --max-running) export ALLELE_MAX_RUNNING=$2 ;;
                    --job-dir) export ALLELE_JOB_DIR=$2 ;;
                    --lsf-memory-scope) export ALLELE_LSF_MEMORY_SCOPE=$2 ;;
                esac
                shift 2 ;;
            *) ALLELE_EXECUTION_ARGS+=("$1"); shift ;;
        esac
    done
    export ALLELE_GRCH=${ALLELE_GRCH:-38}
    [[ $ALLELE_GRCH == 38 ]] || fail 'ALLELE_GRCH must be 38'
    if [[ ${ALLELE_BATCH_WORKER:-0} != 1 ]]; then
        export ALLELE_ANALYSIS_ROOT=${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}
        case "$ALLELE_ANALYSIS_ROOT" in
            */37|*/37/*) fail 'Analysis root conflicts with --grch 38' ;;
            *) ;; # Build is recorded in manifests; do not add another directory level.
        esac
    fi
    export ALLELE_EXECUTION_MODE=$submit
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 && $submit == hpc && -z ${LSB_JOBID:-} ]]; then
        fail 'HPC worker requires an LSF allocation (LSB_JOBID); use bsub to submit the generated .cmd'
    fi
    if (( ! help )); then
        printf '[execution] mode=%s worker=%s jobs_N=%s max_running=%s\n' "$submit" "${ALLELE_BATCH_WORKER:-0}" "$jobs_N" "${ALLELE_MAX_RUNNING:-32}"
        printf '[paths] DIR_DATA=%s ANALYSIS_ROOT=%s\n' "$DIR_DATA" "$ANALYSIS_ROOT"
    fi
}

allele_execution_dispatch() {
    local arg help=0 batch=0 planner dry=0
    for arg in "$@"; do
        case "$arg" in -h|--help|help|--help-full|--version) help=1 ;; --plan-only) dry=1 ;; esac
    done
    [[ ${1:-} != batch ]] || batch=1
    if ((help && ! batch)); then return 0; fi
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 ]]; then return 0; fi
    # Both backends first freeze the same scientific command as a .cmd file.
    # Only the receipt-aware dispatcher decides bash versus bsub execution.
    planner=${ALLELE_PLANNER_PYTHON:-${SHARED_PYTHON:-python3}}
    command -v "$planner" >/dev/null 2>&1 || planner=python3
    command -v "$planner" >/dev/null 2>&1 || fail 'Python 3 is required to generate the job plan'
    if [[ ${BASH_SOURCE[2]##*/} == 2.herv_scan.sh ]] && (( ! dry && ! help )); then
        if "$planner" "$ROOT/f/herv_archive.py" "$@"; then exit 0; fi
    fi
    local -a planned=("$planner" "$ROOT/f/batch_workflows.py" --entry "${BASH_SOURCE[2]}" --mode "${submit:-local}" -- "$@")
    exec "${planned[@]}"

}

workflow_memory_check() {
    local memory=$1 allocated=${ALLELE_ALLOCATED_MEMORY_GIB:-}
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 && ${submit:-local} == hpc ]]; then
        [[ -n ${LSB_JOBID:-} ]] || fail 'Missing LSF job allocation'
    else
        local budget=${ALLELE_ALLOCATED_MEMORY_GIB:-${ALLELE_JOB_MEMORY_GIB:-32}}
        ((memory <= budget)) || fail "Requested memory $memory GiB exceeds local budget $budget GiB; set --job-memory-gb"
    fi
}

workflow_cpu_check() {
    local count=$1 local_cap=${ALLELE_ALLOCATED_CORES:-$(getconf _NPROCESSORS_ONLN)} allocated=${ALLELE_ALLOCATED_CORES:-0}
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 && ${submit:-local} == hpc ]]; then
        ((count <= allocated)) || fail "Requested $count CPU slots exceeds allocated $allocated"
    else
        ((count <= local_cap)) || fail "Use at most $local_cap workers/threads on this local shared host"
    fi
}
