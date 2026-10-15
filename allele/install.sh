#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=f/common.sh
source "$ROOT/f/common.sh"
allele_load_config "$ROOT"

usage(){ cat <<'HELP'
用法
  ./install.sh          创建/更新 conda 环境
  ./install.sh --check  检查环境、依赖和数据路径
HELP
}

check(){
  allele_require_env
  allele_python "$ROOT/f/validate.py" --root "$ROOT" --out "$ALLELE_OUT"
}

case "${1:-install}" in
  -h|--help) usage; exit 0;;
  --check) [[ $# -eq 1 ]] || { usage; exit 2; }; check; exit $?;;
  install) [[ $# -le 1 ]] || { usage; exit 2; };;
  *) usage; exit 2;;
esac

allele_find_conda
if allele_env_exists; then
  "$ALLELE_CONDA" env update -n "$ALLELE_ENV_NAME" -f "$ROOT/environment.yml"
else
  "$ALLELE_CONDA" env create -n "$ALLELE_ENV_NAME" -f "$ROOT/environment.yml"
fi
check
