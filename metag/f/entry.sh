#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
bash ./entry.sh -h  # 内部工具/单独阶段；调用要求见下文
# 从所属项目主入口选择模块；本文件不代表整个项目
# 必填输入沿用现有接口，不自动补造数据或调度参数

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
set -euo pipefail
METAG_CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "$METAG_CODE_ROOT/f/env.sh"
workflow_args "$@"
set -- "${WORKFLOW_ARGS[@]}"
"$METAG_PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else "metag HPC needs Python >= 3.10; load a site Python module or set METAG_PYTHON in env.sh")'
"$METAG_PYTHON" -c 'import yaml' || { echo 'PyYAML is required in METAG_PYTHON (python -m pip install PyYAML).' >&2; exit 1; }
exec "$METAG_PYTHON" "$METAG_CODE_ROOT/f/parallel_cli.py" "$@"

