#!/bin/bash -l

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
bash ./install-hpc.sh -h  # 内部工具/单独阶段；调用要求见下文
# 从所属项目主入口选择模块；本文件不代表整个项目
# 必填输入沿用现有接口，不自动补造数据或调度参数

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
set -eo pipefail
: "${LSB_JOBID:?Submit with bsub -q short}"
export PS1="${PS1-}"
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
source f/hpc.runtime.env
export CONDA_ENVS_PATH=$GU_CONDA_BASE/envs CONDA_PKGS_DIRS=$GU_SOFT/pkgs
export CONDARC=$GU_SOFT/condarc
export GU_INSTALL_ENV_READY=1 PIP_NO_INDEX=1 PYTHONDONTWRITEBYTECODE=1
export GU_BUILD_JOBS=8 CMAKE_BUILD_PARALLEL_LEVEL=8 MAKEFLAGS=-j8
trap 's=$?; printf "%s\n" "$s" > "$GU_SOFT/logs/install.exit"' EXIT
bash install.sh --ibdmix-run cpu
bash install.sh --check
bash f/install-arg-hpc.sh
export GU_TEST_OUTPUT_ROOT=$GU_SOFT/validation
mkdir -p "$GU_TEST_OUTPUT_ROOT"
export IBDMIX_UPSTREAM=$GU_SOFT/ibdmix IBDMIX_PREPARE_BIN=$GU_SOFT/ibdmix-prepare/prepare-v3
printf 'GU core software installed and checked\n' > "$GU_SOFT/logs/core.complete"
