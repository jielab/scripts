#!/usr/bin/env bash
set -euo pipefail
PANOME_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PANOME_PYTHON:-$PANOME_HOME/.venv/bin/python}"
if [[ ! -x "$python_bin" ]]; then python_bin=python3; fi
export OMP_NUM_THREADS="${PANOME_THREADS:-2}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
exec "$python_bin" "$PANOME_HOME/f/panome.py" "$@"
