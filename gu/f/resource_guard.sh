#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
source ./resource_guard.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
# The scheduler and OS enforce resource policy. No extra process memory cap.
gu_resource_guard() {
    export GU_RESOURCE_GUARD=auto
}
