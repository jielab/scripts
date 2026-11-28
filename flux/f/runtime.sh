#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
source ./runtime.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
# Resolve paths without borrowing a Python environment.
flux_paths() {
    local requested=auto previous='' arg
    for arg in "$@"; do
        if [[ $previous == --machine ]]; then requested=$arg; fi
        case "$arg" in --machine=*) requested=${arg#*=} ;; esac
        previous=$arg
    done
    submit=$requested
    source "$FLUX_ROOT/../0f/0os.sh" || return
    export FLUX_MACHINE=$submit
    if [[ -z ${FLUX_ENV:-} ]]; then
        local candidate
        for candidate in "${CONDA_PREFIX:-/nonexistent}" "$HOME/anaconda3/envs/ai" "$HOME/miniforge3/envs/ai" "$HOME/.conda/envs/ai" "$OS_GU_CONDA/envs/ai" "$DIR_WORK/.conda/envs/ai"; do
            if [[ ${candidate##*/} == ai && -x $candidate/bin/python && -d $candidate/conda-meta ]]; then
                export FLUX_ENV=$candidate
                break
            fi
        done
        export FLUX_ENV=${FLUX_ENV:-$HOME/.conda/envs/ai}
    fi
}

flux_runtime() {
    flux_paths "$@" || return
    if [[ -n ${FLUX_PYTHON:-} && $FLUX_PYTHON != "$FLUX_ENV/bin/python" ]]; then
        echo "ERROR: FLUX_PYTHON must be $FLUX_ENV/bin/python; unset the old override." >&2
        return 2
    fi
    export FLUX_PYTHON=$FLUX_ENV/bin/python
    if [[ ! -x $FLUX_PYTHON || ! -d $FLUX_ENV/conda-meta ]]; then
        echo "ERROR: Existing ai environment missing: $FLUX_ENV. Set FLUX_ENV to its absolute path." >&2
        return 2
    fi
    export PYTHONNOUSERSITE=1
    export LD_LIBRARY_PATH="$FLUX_ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    "$FLUX_PYTHON" -I "$FLUX_ROOT/f/check_install.py" --preflight --prefix "$FLUX_ENV" || return
}
