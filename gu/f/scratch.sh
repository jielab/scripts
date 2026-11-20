#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
source ./scratch.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
# Shared OS defaults and common submission options.
submit=${submit:-${SUBMIT:-auto}}
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)/0f/0os.sh"
workflow_args "$@"
_gu_args=()
for ((_gu_i=0;_gu_i<${#WORKFLOW_ARGS[@]};_gu_i++)); do
    if [[ ${WORKFLOW_ARGS[_gu_i]} == --cores ]]; then ((_gu_i+=1)); continue; fi
    _gu_args+=("${WORKFLOW_ARGS[_gu_i]}")
done
set -- "${_gu_args[@]}"
unset _gu_args _gu_i
