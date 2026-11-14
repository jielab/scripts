#!/usr/bin/env bash
# Source from the account shell before starting Python/R or submitting jobs.
# Workflow-specific TMPDIR/cache paths can still override these defaults.
# Keep this separate from frozen workflow files so queued plans stay valid.
case "${BASH_SOURCE[0]}" in
    /work/sph-huangj/*) _runtime_cache_default=/work/sph-huangj/tmp ;;
    /mnt/d/*) _runtime_cache_default=/mnt/d/tmp ;;
    *) echo 'runtime_cache.sh: unsupported checkout location' >&2; return 2 2>/dev/null || exit 2 ;;
esac
export RUNTIME_CACHE_ROOT=${RUNTIME_CACHE_ROOT:-$_runtime_cache_default}
unset _runtime_cache_default
export DIR_RUNTIME_TMP=$RUNTIME_CACHE_ROOT
export TMPDIR=$RUNTIME_CACHE_ROOT/runtime
export TMP=$TMPDIR TEMP=$TMPDIR
export XDG_CACHE_HOME=$RUNTIME_CACHE_ROOT/cache
export PYTHONPYCACHEPREFIX=$XDG_CACHE_HOME/python
export PIP_CACHE_DIR=$XDG_CACHE_HOME/pip
export MPLCONFIGDIR=$XDG_CACHE_HOME/matplotlib
export NUMBA_CACHE_DIR=$XDG_CACHE_HOME/numba
export TORCH_HOME=$XDG_CACHE_HOME/torch
export TRITON_CACHE_DIR=$XDG_CACHE_HOME/triton
export HF_HOME=$XDG_CACHE_HOME/huggingface
export CUDA_CACHE_PATH=$XDG_CACHE_HOME/cuda
mkdir -p -- "$TMPDIR" "$PYTHONPYCACHEPREFIX" "$PIP_CACHE_DIR" \
    "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$TORCH_HOME" \
    "$TRITON_CACHE_DIR" "$HF_HOME" "$CUDA_CACHE_PATH"
