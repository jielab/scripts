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
# Reuse the repository's path table; a standalone unpacked demo also works.
flux_runtime() {
    local requested=auto previous='' arg
    for arg in "$@"; do
        if [[ $previous == --machine ]]; then requested=$arg; fi
        case "$arg" in --machine=*) requested=${arg#*=} ;; esac
        previous=$arg
    done
    submit=$requested
    source "$FLUX_ROOT/../0f/0os.sh" || return
    export FLUX_MACHINE=$submit
    export FLUX_ENV=${FLUX_ENV:-$DIR_WORK/software/flux}
    if [[ -z ${FLUX_PYTHON:-} ]]; then
        if [[ -x "$FLUX_ENV/bin/python" ]]; then FLUX_PYTHON=$FLUX_ENV/bin/python;
        else FLUX_PYTHON=$(command -v python3); fi
    fi
    export FLUX_PYTHON
}
