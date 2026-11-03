#!/usr/bin/env bash
# Standalone on either cluster. Select modules in the login shell or an explicit
# METAG_ENV_FILE; no TaiYi-specific compiler/Python path is assumed.
if [[ -n ${METAG_ENV_FILE:-} ]]; then
    source "$METAG_ENV_FILE"
fi
export METAG_PYTHON=${METAG_PYTHON:-python3}
export CC=${CC:-gcc} CXX=${CXX:-g++}
export OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset OMP_NUM_THREADS
export PYTHONDONTWRITEBYTECODE=1
