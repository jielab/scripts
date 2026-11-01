#!/usr/bin/env bash
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)/0f/0os.sh"
if [[ $submit == hpc ]]; then
    if ! type module >/dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
    module load gcc/11.2.0
    export METAG_PYTHON=${METAG_PYTHON:-$DIR_WORK/.conda/envs/py311/bin/python3}
    export CC=/share/base/gcc/11.2.0/bin/gcc CXX=/share/base/gcc/11.2.0/bin/g++
else
    export METAG_PYTHON=${METAG_PYTHON:-python3}
    export CC=${CC:-gcc} CXX=${CXX:-g++}
fi
export OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset OMP_NUM_THREADS
export PYTHONDONTWRITEBYTECODE=1
