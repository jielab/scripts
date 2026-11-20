#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
source ./cxx.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
set -euo pipefail
args=()
if [[ -n ${LD_ZSTD_PREFIX:-} ]]; then
    args+=("-I$LD_ZSTD_PREFIX/include" "-L$LD_ZSTD_PREFIX/lib" "-Wl,-rpath,$LD_ZSTD_PREFIX/lib")
fi
exec "${LD_CXX:-g++}" "${args[@]}" "$@" -pthread -fopenmp

