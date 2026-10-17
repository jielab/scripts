#!/usr/bin/env bash
# Public entry for GRID twin analysis, benchmark dispatch, and the viewer.
# Usage and copy-ready examples: ./grid.sh --help (defined in f/twin.py).
set -euo pipefail
# Preserve the public command's raw output while supervising every worker.
SCRIPT_CONSOLE_PASSTHROUGH=1 source /mnt/d/scripts/0f/console.sh
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
