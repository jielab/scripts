#!/usr/bin/env bash
# Public entry for GRID twin analysis, benchmark dispatch, and the viewer.
set -euo pipefail
umask 077
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
# Match the existing noninteractive GRID environment selection without importing
# upstream pipeline configuration or changing the caller's R libraries.
runtime=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
if [[ -x "$runtime/bin/python3" ]]; then
  default_python="$runtime/bin/python3"
else
  default_python=python3
fi
PYTHON=${GRID_PYTHON:-${ABM_PYTHON:-$default_python}}
export PYTHONDONTWRITEBYTECODE=1
exec "$PYTHON" "$ROOT/f/twin.py" "$@"
