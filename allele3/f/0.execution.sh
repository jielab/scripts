#!/usr/bin/env bash
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
    # Resolve these two selectors before the first source. In 0os.sh itself,
    # switching mode unsets DIR_SCRATCH; a following --date-scratch would then
    # read it under set -u. Preselection avoids that order-sensitive case
    # without changing the shared file.
    local i chosen_mode=${submit:-${SUBMIT:-auto}} chosen_date=${DATE_SCRATCH:-2026-11-03}
    for ((i=0;i<${#expanded[@]};i++)); do
        case "${expanded[i]}" in
            --submit) [[ $((i+1)) -lt ${#expanded[@]} ]] || fail '--submit needs local or hpc'; chosen_mode=${expanded[i+1]} ;;
            --date-scratch) [[ $((i+1)) -lt ${#expanded[@]} ]] || fail '--date-scratch needs a date'; chosen_date=${expanded[i+1]} ;;
        esac
    done
    case "$chosen_mode" in auto|local|hpc) ;; *) fail '--submit needs local or hpc' ;; esac
    if [[ ${submit:-${SUBMIT:-auto}} != auto && $chosen_mode != "${submit:-$SUBMIT}" ]]; then
        unset DIR_WORK DIR_DATA DIR_SCRATCH DIR_SOFTWARE SHARED_PYTHON
    elif [[ ${DIR_SCRATCH:-} == /scratch/${DATE_SCRATCH:-2026-11-03}/sph-huangj && $chosen_date != ${DATE_SCRATCH:-2026-11-03} ]]; then
        unset DIR_SCRATCH
    fi
    export submit=$chosen_mode SUBMIT=$chosen_mode DATE_SCRATCH=$chosen_date
    # The shared file owns the OS roots. No separate OS table is maintained.
    source "$os_config"
    workflow_args "${expanded[@]}"
    set -- "${WORKFLOW_ARGS[@]}"
    ALLELE_EXECUTION_ARGS=()
    while (($#)); do
        case "$1" in
            --write-cmd) export ALLELE_WRITE_CMD=1; shift ;;
            --retry-failed) export ALLELE_RETRY_FAILED=1; shift ;;
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
    export ALLELE_EXECUTION_MODE=$submit
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 && $submit == hpc && -z ${LSB_JOBID:-} ]]; then
        fail 'HPC worker requires an LSF allocation (LSB_JOBID); use bsub to submit the generated .cmd'
    fi
    if (( ! help )); then
        printf '[execution] mode=%s worker=%s date_scratch=%s jobs_N=%s max_running=%s\n' "$submit" "${ALLELE_BATCH_WORKER:-0}" "$DATE_SCRATCH" "$jobs_N" "${ALLELE_MAX_RUNNING:-32}"
        printf '[paths] DIR_WORK=%s DIR_DATA=%s DIR_SCRATCH=%s\n' "$DIR_WORK" "$DIR_DATA" "$DIR_SCRATCH"
    fi
}

allele_execution_dispatch() {
    local arg help=0 batch=0 planner dry=0
    for arg in "$@"; do
        case "$arg" in -h|--help|help|--help-full|--version) help=1 ;; --dry-run|--preflight) dry=1 ;; esac
    done
    [[ ${1:-} != batch ]] || batch=1
    if ((help && ! batch)); then return 0; fi
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 ]]; then return 0; fi
    if [[ ${submit:-local} == local && $batch == 0 && ${ALLELE_WRITE_CMD:-0} != 1 ]]; then
        if [[ -z ${LSB_JOBID:-} ]] && command -v bsub >/dev/null 2>&1; then
            fail 'Local computation on an LSF submission host requires an allocation; choose --submit hpc'
        fi
        return 0
    fi
    planner=${ALLELE_PLANNER_PYTHON:-${SHARED_PYTHON:-python3}}
    command -v "$planner" >/dev/null 2>&1 || planner=python3
    command -v "$planner" >/dev/null 2>&1 || fail 'Python 3 is required to generate the job plan'
    local -a planned=("$planner" "$ROOT/f/batch_workflows.py" --entry "${BASH_SOURCE[2]}" --mode "${submit:-local}" -- "$@")
    if [[ ${submit:-local} == local && ${ALLELE_WRITE_CMD:-0} != 1 ]] && (( ! help && ! dry )); then
        export XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-/run/user/$(id -u)}
        export DBUS_SESSION_BUS_ADDRESS=${DBUS_SESSION_BUS_ADDRESS:-unix:path=$XDG_RUNTIME_DIR/bus}
        command -v systemd-run >/dev/null 2>&1 || fail 'Local batches require the process-tree memory cap (systemd-run)'
        [[ -r /sys/fs/cgroup/cgroup.controllers ]] || fail 'Local batches require cgroup v2'
        timeout -k 1s 5s systemctl --user show-environment >/dev/null 2>&1 || fail 'User systemd is unavailable; local batch was not started without a memory cap'
        local bytes=$((32 * 1024 * 1024 * 1024))
        exec systemd-run --user --scope --quiet --unit="allele3-batch-$(date +%s)-$$-$RANDOM" \
            -p "MemoryMax=$bytes" -p MemorySwapMax=0 -p OOMPolicy=kill -p TasksMax=512 \
            bash -c "$(declare -f verify_limits)"$'\nset -euo pipefail\nverify_limits "$@"' \
            allele3-batch-guard "$bytes" "${planned[@]}"
    fi
    exec "${planned[@]}"
}

workflow_memory_check() {
    local memory=$1 allocated=${ALLELE_ALLOCATED_MEMORY_GIB:-}
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 && ${submit:-local} == hpc ]]; then
        [[ -n ${LSB_JOBID:-} && $allocated =~ ^[1-9][0-9]*$ ]] || fail 'Missing verified HPC memory allocation'
        ((memory <= allocated)) || fail "Requested memory $memory GiB exceeds allocated $allocated GiB"
    else
        ((memory <= 32)) || fail 'Memory cap must be at most 32 GiB on this local shared host'
    fi
}

workflow_cpu_check() {
    local count=$1 local_cap=$2 allocated=${ALLELE_ALLOCATED_CORES:-0}
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 && ${submit:-local} == hpc ]]; then
        ((count <= allocated)) || fail "Requested $count CPU slots exceeds allocated $allocated"
    else
        ((count <= local_cap)) || fail "Use at most $local_cap workers/threads on this local shared host"
    fi
}
