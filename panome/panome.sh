#!/usr/bin/env bash
set -euo pipefail
PANOME_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${PANOME_PYTHON:-}" ]]; then
  python_bin="$PANOME_PYTHON"
elif [[ -x "$PANOME_HOME/.venv/bin/python" ]]; then
  python_bin="$PANOME_HOME/.venv/bin/python"
else
  python_bin=python3
fi
export OMP_NUM_THREADS="${PANOME_THREADS:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMEXPR_NUM_THREADS="$OMP_NUM_THREADS"
export CUBLAS_WORKSPACE_CONFIG="${CUBLAS_WORKSPACE_CONFIG:-:4096:8}"
exec "$python_bin" "$PANOME_HOME/f/panome.py" "$@"
