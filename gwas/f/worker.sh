#!/usr/bin/env bash
set -euo pipefail
if [[ ${submit:-local} == hpc && -z ${LSB_JOBID:-} ]]; then
    echo 'Submit this worker with bsub; it cannot run on a login node.' >&2
    exit 2
fi
if [[ -n ${GWAS_RUNTIME_ENV:-} ]]; then
    [[ -r $GWAS_RUNTIME_ENV ]] || { echo "Missing runtime: $GWAS_RUNTIME_ENV" >&2; exit 2; }
    PS1=${PS1-}
    set +u
    source "$GWAS_RUNTIME_ENV"
    set -euo pipefail
fi
export GWAS_WORKER=1 GWAS_PLAN_ONLY=FALSE
export OMP_NUM_THREADS=${GWAS_PROCESS_THREADS:-1} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=${GWAS_PROCESS_THREADS:-1}
exec bash "$@"
