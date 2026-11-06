#!/usr/bin/env bash
# Shared UKB launcher: resolve paths, freeze one .cmd, then run it or use LSF.
ukb_workflow_init() {
    local directory argument
    directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
    while [[ ! -f "$directory/0f/0os.sh" ]]; do
        [[ $directory != / ]] || { echo 'ERROR: cannot locate 0f/0os.sh' >&2; return 2; }
        directory=$(dirname -- "$directory")
    done
    shared_dir=$directory/0f
    # Preserve the explicit local-copy convention; scripts/ukb uses 0os defaults.
    case ${BASH_SOURCE[1]} in */0LOCAL/*) export submit=local SUBMIT=local ;; esac
    source "$shared_dir/0os.sh" || return
    workflow_args "$@" || return
    set -- "${WORKFLOW_ARGS[@]}"
    UKB_ARGS=()
    UKB_EXECUTE=0
    UKB_WRITE_ONLY=0
    while (($#)); do
        case "$1" in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
        case "$1" in
            --execute) UKB_EXECUTE=1; shift ;;
            --write-cmd|--dry-run) UKB_WRITE_ONLY=1; shift ;;
            --background)
                export submit=hpc SUBMIT=hpc
                source "$shared_dir/0os.sh" || return
                shift ;;
            --run-cmd)
                case ${2:-} in
                    TRUE|true) UKB_WRITE_ONLY=0 ;;
                    FALSE|false) UKB_WRITE_ONLY=1 ;;
                    *) echo 'ERROR: --run-cmd requires TRUE or FALSE' >&2; return 2 ;;
                esac
                shift 2 ;;
            -q|--queue|-n|--slots|--memory-gb|--memory-gib|--walltime|--threads|--dir-jobs)
                [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }
                case "$1" in
                    -q|--queue) UKB_QUEUE=$2 ;;
                    -n|--slots) UKB_SLOTS=$2 ;;
                    --memory-gb|--memory-gib) UKB_MEMORY_GIB=$2 ;;
                    --walltime) UKB_WALLTIME=$2 ;;
                    --threads) UKB_THREADS=$2 ;;
                    --dir-jobs) UKB_JOB_DIR=$2 ;;
                esac
                shift 2 ;;
            *) UKB_ARGS+=("$1"); shift ;;
        esac
    done
    UKB_THREADS=${UKB_THREADS:-4}
    UKB_SLOTS=${UKB_SLOTS:-$UKB_THREADS}
    UKB_MEMORY_GIB=${UKB_MEMORY_GIB:-$MEMORY_DEFAULT_GIB}
    UKB_QUEUE=${UKB_QUEUE:-${HPC_CPU_QUEUE:-33}}
    UKB_WALLTIME=${UKB_WALLTIME:-${HPC_WALLTIME:-24:00}}
    UKB_PYTHON=${UKB_PYTHON:-$SHARED_PYTHON}
    if [[ $submit == hpc ]]; then
        UKB_RSCRIPT=${UKB_RSCRIPT:-${LE8_RSCRIPT:-$DIR_WORK/.conda/envs/le8/bin/Rscript}}
    else
        UKB_RSCRIPT=${UKB_RSCRIPT:-$(command -v Rscript || printf Rscript)}
    fi
    for argument in UKB_THREADS UKB_SLOTS UKB_MEMORY_GIB; do
        [[ ${!argument} =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: $argument must be a positive integer" >&2; return 2; }
    done
    [[ $UKB_QUEUE =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo 'ERROR: invalid queue' >&2; return 2; }
    [[ $UKB_WALLTIME =~ ^[0-9]+:[0-5][0-9]$ && $UKB_WALLTIME != 0:00 ]] || { echo 'ERROR: walltime must be HHH:MM and positive' >&2; return 2; }
    if (( UKB_EXECUTE )) && [[ $submit == hpc && -z ${LSB_JOBID:-} ]]; then
        echo 'ERROR: HPC workers require an LSF allocation; submit the .cmd with bsub' >&2
        return 2
    fi
}

ukb_dispatch() {
    local job_name=$1 root=$2 entry=$3 resources name cmd run_dir
    shift 3
    (( UKB_EXECUTE == 0 )) || return 0
    (( UKB_SLOTS >= UKB_THREADS )) || UKB_SLOTS=$UKB_THREADS
    if [[ $submit == hpc ]]; then
        resources=$("$UKB_PYTHON" "$shared_dir/qiming_resources.py" --queue "$UKB_QUEUE" \
            --cores "$UKB_SLOTS" --memory-gib "$UKB_MEMORY_GIB" --walltime "$UKB_WALLTIME" --format shell) || return
        read -r UKB_QUEUE UKB_SLOTS <<< "$resources"
    fi
    root=${UKB_JOB_DIR:-$root/jobs}
    [[ $root != *$'\n'* && $root != *'"'* && $root != *'\'* ]] || { echo 'ERROR: job directory contains an unsupported character' >&2; return 2; }
    mkdir -p -- "$root"
    root=$(cd -- "$root" && pwd -P)
    job_name=${job_name//[^a-zA-Z0-9_.-]/_}
    run_dir=$(mktemp -d "$root/$job_name.$(date +%Y%m%dT%H%M%S).XXXXXX")
    cmd=$run_dir/$job_name.cmd
    {
        printf '#!/bin/bash -l\n'
        if [[ $submit == hpc ]]; then
            printf '#BSUB -q %s\n#BSUB -n %s\n#BSUB -W %s\n' "$UKB_QUEUE" "$UKB_SLOTS" "$UKB_WALLTIME"
            printf '#BSUB -M %sGB\n#BSUB -R "span[hosts=1] rusage[mem=%sGB/host]"\n' "$UKB_MEMORY_GIB" "$UKB_MEMORY_GIB"
            printf '#BSUB -J ukb.%s\n#BSUB -o "%s/%%J.out"\n#BSUB -e "%s/%%J.err"\n' "$job_name" "$run_dir" "$run_dir"
        fi
        printf 'set -euo pipefail\ncd %q\n' "$PWD"
        # Freeze only workflow settings, never credentials or the entire environment.
        for name in "${OS_PATH_KEYS[@]}" submit SUBMIT jobs_N SUBMIT_WAIT \
            UKB_THREADS UKB_SLOTS UKB_MEMORY_GIB UKB_QUEUE UKB_WALLTIME UKB_PYTHON UKB_RSCRIPT \
            DIR0 GEN_DIR PGS_DIR LIFTOVER_BIN PGS_PLINK_MEMORY_MB UKB_PLINK_MEMORY_MB \
            GP_REBUILD_FULL_FST SRD_YEAR_CAP ADULT_ENABLE ADULT_TRAITS ADULT_CIP_FROM_UKB \
            ADULT_CIP_FILE ADULT_TRAIT_FILE ADULT_CASE_MODE ADULT_STANDARDIZE ADULT_MIN_CIP \
            PHEWAS_ROSTER PHEWAS_COVERAGE_DEFINITION PHEWAS_REBUILD_OUTDIR PHEWAS_REPLACE; do
            case $name in DIR_WINDOWS*|OS_RSCRIPT_WINDOWS|OS_RTOOLS_WINDOWS) [[ $submit != hpc ]] || continue ;; esac
            [[ ! -v $name ]] || printf 'export %s=%q\n' "$name" "${!name}"
        done
        printf 'export WORKFLOW_WORKER=1 SCRIPT_CONSOLE_ACTIVE=1\n'
        if [[ $submit == hpc ]]; then
            printf '[[ -n ${LSB_JOBID:-} ]] || { echo "Submit this .cmd with bsub" >&2; exit 2; }\n'
            printf '(( ${LSB_DJOB_NUMPROC:-0} >= %s )) || { echo "Insufficient allocated CPU slots" >&2; exit 2; }\n' "$UKB_SLOTS"
        fi
        printf 'exec bash %q --execute' "$entry"
        printf ' %q' "$@"
        printf '\n'
    } > "$cmd"
    chmod +x "$cmd"
    bash -n "$cmd"
    printf '[UKB] Written %s\n' "$cmd"
    if (( UKB_WRITE_ONLY )); then
        printf '[UKB] Command file only; submit=%s\n' "$submit"
        exit 0
    fi
    if [[ $submit == local ]]; then
        exec bash "$cmd"
    fi
    local runner=("$UKB_PYTHON" "$shared_dir/submit_jobs.py" --submit hpc --jobs-N "$jobs_N")
    [[ $SUBMIT_WAIT != TRUE ]] || runner+=(--wait)
    exec "${runner[@]}" "$cmd"
}

ukb_worker_runtime() {
    # Absolute runtimes avoid whichever unrelated conda environment a login shell uses.
    command -v "$UKB_PYTHON" >/dev/null || { echo "ERROR: Python not found: $UKB_PYTHON" >&2; return 2; }
    local python_executable
    python_executable=$(command -v "$UKB_PYTHON")
    PATH="$DIR_SOFTWARE/bin:$(dirname -- "$python_executable"):$PATH"
    export PATH
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=$UKB_THREADS R_DATATABLE_NUM_THREADS=$UKB_THREADS
    export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    export TMPDIR=${TMPDIR:-$DIR_RUNTIME_TMP}
    UKB_PLINK_MEMORY_MB=${UKB_PLINK_MEMORY_MB:-$((UKB_MEMORY_GIB * 1024 * 3 / 4))}
    export UKB_PLINK_MEMORY_MB
    [[ $UKB_PLINK_MEMORY_MB =~ ^[1-9][0-9]*$ ]] && (( UKB_PLINK_MEMORY_MB <= UKB_MEMORY_GIB * 1024 )) || {
        echo 'ERROR: UKB_PLINK_MEMORY_MB exceeds the job RAM budget or is invalid' >&2; return 2;
    }
}
