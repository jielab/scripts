#!/usr/bin/env bash
# Reuse the repository's path table; a standalone unpacked demo also works.
flux_runtime() {
    local requested=hpc previous='' arg
    for arg in "$@"; do
        if [[ $previous == --machine ]]; then requested=$arg; fi
        case "$arg" in --machine=*) requested=${arg#*=} ;; esac
        previous=$arg
    done
    if [[ -f "$FLUX_ROOT/../0f/0os.sh" ]]; then
        submit=$requested
        source "$FLUX_ROOT/../0f/0os.sh"
    else
        case "$requested" in
            auto) if [[ -n ${LSB_JOBID:-}${LSF_ENVDIR:-}${LSF_SERVERDIR:-} ]] || command -v bsub >/dev/null 2>&1; then submit=hpc; else submit=local; fi ;;
            local|hpc) submit=$requested ;;
            *) echo 'ERROR: --machine must be auto, local or hpc' >&2; return 2 ;;
        esac
        if [[ $submit == hpc ]]; then
            export DIR_WORK=${DIR_WORK:-/work/sph-huangj}
        else
            export DIR_WORK=${DIR_WORK:-/mnt/d}
        fi
        export DIR_UKB=${DIR_UKB:-$DIR_WORK/phe/ukb}
        export ANALYSIS_ROOT=${ANALYSIS_ROOT:-$DIR_WORK/analysis}
    fi
    export FLUX_MACHINE=$submit
    export FLUX_ENV=${FLUX_ENV:-$DIR_WORK/software/flux}
    if [[ -z ${FLUX_PYTHON:-} ]]; then
        if [[ -x "$FLUX_ENV/bin/python" ]]; then FLUX_PYTHON=$FLUX_ENV/bin/python;
        else FLUX_PYTHON=$(command -v python3); fi
    fi
    export FLUX_PYTHON
}
