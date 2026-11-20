#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
bash ./qc_entry.sh -h  # 内部工具/单独阶段；调用要求见下文
# 从所属项目主入口选择模块；本文件不代表整个项目
# 必填输入沿用现有接口，不自动补造数据或调度参数

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
module=$1
shift
if [[ ${1:-} == --help-full ]]; then
    shift
    source "$HERE/f/env.sh"
    exec "$METAG_PYTHON" "$HERE/f/${module}.py" "$@" --help
fi
for arg in "$@"; do
    if [[ $arg == -h || $arg == --help ]]; then
        bash "$HERE/f/usage.sh" "$module"
        exit 0
    fi
done
source "$HERE/f/env.sh"
exec "$METAG_PYTHON" "$HERE/f/qc_submit.py" "$module" "$@"

