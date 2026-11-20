#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
case ${1:-create} in
    -h|--help) echo 'Usage: bash install.sh [create|update]. Creates/updates conda env grid; Optional deep/PRSformer GPU dependencies are separate.';exit 0 ;;
    create|update) action=${1:-create} ;;
    *) echo 'Use create or update' >&2;exit 2 ;;
esac
if command -v mamba >/dev/null 2>&1; then manager=mamba
elif command -v conda >/dev/null 2>&1; then manager=conda
else echo 'Load conda or mamba, then rerun install.sh.' >&2;exit 2
fi
"$manager" env "$action" -f "$ROOT/env.yml"
echo 'Activate grid, then run: bash grid.sh check'
echo 'Set GRID_PYTHON to this environment python if it is outside the standard paths.'
