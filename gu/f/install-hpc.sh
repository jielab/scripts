#!/bin/bash -l
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
bash install.sh --GPU FALSE
bash install.sh --check
bash f/install-arg-hpc.sh
export GU_TEST_OUTPUT_ROOT=$GU_SOFT/validation
mkdir -p "$GU_TEST_OUTPUT_ROOT"
export IBDMIX_UPSTREAM=$GU_SOFT/ibdmix IBDMIX_PREPARE_BIN=$GU_SOFT/ibdmix-prepare/prepare-v3
python -m pytest -q -p no:cacheprovider f/tests -rs
printf 'GU core software installed and checked\n' > "$GU_SOFT/logs/core.complete"
