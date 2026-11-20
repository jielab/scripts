#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
bash ./worker.sh -h  # 内部工具/单独阶段；调用要求见下文
# 从所属项目主入口选择模块；本文件不代表整个项目
# 必填输入沿用现有接口，不自动补造数据或调度参数

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
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
