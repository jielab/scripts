#!/usr/bin/env bash
# Standalone on either cluster. Select modules in the login shell or an explicit
# METAG_ENV_FILE; no TaiYi-specific compiler/Python path is assumed.
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)/0f/0os.sh"
if [[ -n ${METAG_ENV_FILE:-} ]]; then
    source "$METAG_ENV_FILE"
fi
if command -v bsub >/dev/null 2>&1; then
    if [[ -z ${METAG_PYTHON:-} && -x $OS_GU_CONDA/bin/python ]]; then
        export METAG_PYTHON=$OS_GU_CONDA/bin/python
    fi
    # QiMing's system GCC 8 needs a newer C++17 filesystem runtime.
    if [[ $(g++ -dumpversion 2>/dev/null) == 8* ]]; then
        if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
            source /etc/profile.d/modules.sh
        fi
        if type module >/dev/null 2>&1; then module load gcc/12.2.0; fi
    fi
fi
export METAG_PYTHON=${METAG_PYTHON:-python3}
export CC=${CC:-gcc} CXX=${CXX:-g++}
export OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset OMP_NUM_THREADS
export PYTHONDONTWRITEBYTECODE=1
