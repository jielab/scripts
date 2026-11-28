#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
source ./phenotype_jobs.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
# Shared YGB launcher: resolve paths, freeze one .cmd, then run it or use LSF.
ygb_job_init() {
    local directory argument
    directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
    while [[ ! -f "$directory/0f/0os.sh" ]]; do
        [[ $directory != / ]] || { echo 'ERROR: cannot locate 0f/0os.sh' >&2; return 2; }
        directory=$(dirname -- "$directory")
    done
    shared_dir=$directory/0f
    # Preserve the explicit local-copy convention; scripts/ygb uses 0os defaults.
    case ${BASH_SOURCE[1]} in */0LOCAL/*) export submit=local SUBMIT=local ;; esac
    source "$shared_dir/0os.sh" || return
    workflow_args "$@" || return
    set -- "${WORKFLOW_ARGS[@]}"
    YGB_ARGS=()
    YGB_EXECUTE=0
    YGB_WRITE_ONLY=0
    while (($#)); do
        case "$1" in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
        case "$1" in
            --cores) shift 2 ;;
            --execute) YGB_EXECUTE=1; shift ;;
            --write-cmd|--dry-run) YGB_WRITE_ONLY=1; shift ;;
            --background)
                export submit=hpc SUBMIT=hpc
                source "$shared_dir/0os.sh" || return
                shift ;;
            --run-cmd)
                case ${2:-} in
                    TRUE|true) YGB_WRITE_ONLY=0 ;;
                    FALSE|false) YGB_WRITE_ONLY=1 ;;
                    *) echo 'ERROR: --run-cmd requires TRUE or FALSE' >&2; return 2 ;;
                esac
                shift 2 ;;
            -q|--queue|-n|--slots|--memory-gb|--memory-gib|--walltime|--threads|--dir-jobs)
                [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }
                case "$1" in
                    -q|--queue) YGB_QUEUE=$2 ;;
                    -n|--slots) YGB_SLOTS=$2 ;;
                    --memory-gb|--memory-gib) YGB_MEMORY_GIB=$2 ;;
                    --walltime) YGB_WALLTIME=$2 ;;
                    --threads) YGB_THREADS=$2 ;;
                    --dir-jobs) YGB_JOB_DIR=$2 ;;
                esac
                shift 2 ;;
            *) YGB_ARGS+=("$1"); shift ;;
        esac
    done
    YGB_THREADS=${YGB_THREADS:-$WORKFLOW_CORES}
    (( YGB_THREADS <= WORKFLOW_CORES )) || { echo "--threads exceeds --cores" >&2; return 2; }
    YGB_SLOTS=${YGB_SLOTS:-$WORKFLOW_CORES}
    (( YGB_SLOTS >= YGB_THREADS )) || { echo "slots must cover threads" >&2; return 2; }
    YGB_MEMORY_GIB=${YGB_MEMORY_GIB:-$MEMORY_DEFAULT_GIB}
    YGB_QUEUE=${YGB_QUEUE:-${HPC_CPU_QUEUE:-33}}
    YGB_WALLTIME=${YGB_WALLTIME:-${HPC_WALLTIME:-24:00}}
    YGB_PYTHON=${YGB_PYTHON:-$SHARED_PYTHON}
    if [[ $submit == hpc ]]; then
        YGB_RSCRIPT=${YGB_RSCRIPT:-Rscript}
    else
        YGB_RSCRIPT=${YGB_RSCRIPT:-$(command -v Rscript || printf Rscript)}
    fi
    for argument in YGB_THREADS YGB_SLOTS YGB_MEMORY_GIB; do
        [[ ${!argument} =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: $argument must be a positive integer" >&2; return 2; }
    done
    [[ $YGB_QUEUE =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo 'ERROR: invalid queue' >&2; return 2; }
    [[ $YGB_WALLTIME =~ ^[0-9]+:[0-5][0-9]$ && $YGB_WALLTIME != 0:00 ]] || { echo 'ERROR: walltime must be HHH:MM and positive' >&2; return 2; }
    if (( YGB_EXECUTE )) && [[ $submit == hpc && -z ${LSB_JOBID:-} ]]; then
        echo 'ERROR: HPC workers require an LSF allocation; submit the .cmd with bsub' >&2
        return 2
    fi
}

ygb_dispatch() {
    local job_name=$1 root=$2 entry=$3 resources name cmd run_dir
    shift 3
    (( YGB_EXECUTE == 0 )) || return 0
    (( YGB_SLOTS >= YGB_THREADS )) || YGB_SLOTS=$YGB_THREADS
    if [[ $submit == hpc ]]; then
        resources=$("$YGB_PYTHON" "$shared_dir/qiming_resources.py" --queue "$YGB_QUEUE" \
            --cores "$YGB_SLOTS" --memory-gib "$YGB_MEMORY_GIB" --walltime "$YGB_WALLTIME" --format shell) || return
        read -r YGB_QUEUE YGB_SLOTS <<< "$resources"
    fi
    root=${YGB_JOB_DIR:-$root/logs}
    [[ $root != *$'\n'* && $root != *'"'* && $root != *'\'* ]] || { echo 'ERROR: job directory contains an unsupported character' >&2; return 2; }
    mkdir -p -- "$root"
    root=$(cd -- "$root" && pwd -P)
    job_name=${job_name//[^a-zA-Z0-9_.-]/_}
    run_dir=$("$YGB_PYTHON" "$shared_dir/workflow_layout.py" "$root") || return
    cmd=$run_dir/$job_name.cmd
    {
        printf '#!/bin/bash -l\n'
        if [[ $submit == hpc ]]; then
            printf '#BSUB -q %s\n#BSUB -n %s\n' "$YGB_QUEUE" "$YGB_SLOTS"
            printf '#BSUB -R "span[hosts=1] rusage[mem=%sGB/host]"\n' "$YGB_MEMORY_GIB"
:
            printf '#BSUB -J ygb.%s\n#BSUB -o "%s/worker.out"\n#BSUB -e "%s/worker.err"\n' "$job_name" "$run_dir" "$run_dir"
        fi
        printf 'set -euo pipefail\ncd %q\n' "$PWD"
        # Freeze only workflow settings, never credentials or the entire environment.
        for name in "${OS_PATH_KEYS[@]}" submit SUBMIT jobs_N SUBMIT_WAIT \
            YGB_THREADS YGB_SLOTS YGB_MEMORY_GIB YGB_QUEUE YGB_WALLTIME YGB_PYTHON YGB_RSCRIPT \
            DIR0 GEN_DIR PGS_DIR LIFTOVER_BIN PGS_PLINK_MEMORY_MB YGB_PLINK_MEMORY_MB \
            GP_REBUILD_FULL_FST SRD_YEAR_CAP ADULT_ENABLE ADULT_TRAITS ADULT_CIP_FROM_YGB \
            ADULT_CIP_FILE ADULT_TRAIT_FILE ADULT_CASE_MODE ADULT_STANDARDIZE ADULT_MIN_CIP \
            PHEWAS_ROSTER PHEWAS_COVERAGE_DEFINITION PHEWAS_REBUILD_OUTDIR PHEWAS_REPLACE; do
            case $name in DIR_WINDOWS*|OS_RSCRIPT_WINDOWS|OS_RTOOLS_WINDOWS) [[ $submit != hpc ]] || continue ;; esac
            [[ ! -v $name ]] || printf 'export %s=%q\n' "$name" "${!name}"
        done
        printf 'export WORKFLOW_WORKER=1 SCRIPT_CONSOLE_ACTIVE=1\n'
        if [[ $submit == hpc ]]; then
            printf '[[ -n ${LSB_JOBID:-} ]] || { echo "Submit this .cmd with bsub" >&2; exit 2; }\n'
            printf '(( ${LSB_DJOB_NUMPROC:-0} >= %s )) || { echo "Insufficient allocated CPU slots" >&2; exit 2; }\n' "$YGB_SLOTS"
        fi
        if [[ $submit == hpc ]]; then
            printf 'source %q\n' "$SCRIPTS_ROOT/0f/0os.sh"
        fi
        printf 'exec bash %q --execute'  "$entry"
        printf ' %q' "$@"
        printf '\n'
    } > "$cmd"
    chmod +x "$cmd"
    bash -n "$cmd"
    printf '[YGB] Written %s\n' "$cmd"
    if (( YGB_WRITE_ONLY )); then
        printf '[YGB] Command file only; submit=%s\n' "$submit"
        exit 0
    fi
    if [[ $submit == local ]]; then
        exec bash "$cmd"
    fi
    local runner=("$YGB_PYTHON" "$shared_dir/submit_jobs.py" --submit hpc --jobs-N "$jobs_N")
    [[ $SUBMIT_WAIT != TRUE ]] || runner+=(--wait)
    exec "${runner[@]}" "$cmd"
}

ygb_worker_runtime() {
    # Absolute runtimes avoid whichever unrelated conda environment a login shell uses.
    command -v "$YGB_PYTHON" >/dev/null || { echo "ERROR: Python not found: $YGB_PYTHON" >&2; return 2; }
    local python_executable
    python_executable=$(command -v "$YGB_PYTHON")
    PATH="$DIR_SOFTWARE/bin:$(dirname -- "$python_executable"):$PATH"
    export PATH
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=$YGB_THREADS R_DATATABLE_NUM_THREADS=$YGB_THREADS
    export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    export TMPDIR=${TMPDIR:-$DIR_RUNTIME_TMP}
    YGB_PLINK_MEMORY_MB=${YGB_PLINK_MEMORY_MB:-$((YGB_MEMORY_GIB * 1024 * 3 / 4))}
    export YGB_PLINK_MEMORY_MB
    [[ $YGB_PLINK_MEMORY_MB =~ ^[1-9][0-9]*$ ]] && (( YGB_PLINK_MEMORY_MB <= YGB_MEMORY_GIB * 1024 )) || {
        echo 'ERROR: YGB_PLINK_MEMORY_MB exceeds the job RAM budget or is invalid' >&2; return 2;
    }
}
