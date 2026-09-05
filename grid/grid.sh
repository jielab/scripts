#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

# Make the supported environment usable for non-interactive GRID calls without
# requiring `conda activate` first.
GRID_CONDA_ENV=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid-argneedle}
if [[ -x "$GRID_CONDA_ENV/bin/python3" && ${CONDA_PREFIX:-} != "$GRID_CONDA_ENV" ]]; then
  export PATH="$GRID_CONDA_ENV/bin:$PATH"
fi
if [[ -d "$GRID_CONDA_ENV/lib/R/library" ]]; then
  export R_ENVIRON_USER=/dev/null
  export R_LIBS_USER="$GRID_CONDA_ENV/lib/R/library"
fi
export PYTHONNOUSERSITE=1
if [[ -z ${ARG_NEEDLE_HOME:-} && -x "$GRID_CONDA_ENV/bin/python3" ]]; then
  ARG_NEEDLE_HOME=$(
    "$GRID_CONDA_ENV/bin/python3" - <<'PY' 2>/dev/null || true
import importlib.util
from pathlib import Path
spec = importlib.util.find_spec("arg_needle.scripts.infer_args_advanced")
if spec is not None and spec.origin is not None:
    print(Path(spec.origin).resolve().parents[2])
PY
  )
  export ARG_NEEDLE_HOME
fi
module=${1:-help}
case "$module" in
  -h|--help) module=help ;;
  pca|ancestry|csx|disco|arg|grid|eval|all|help) ;;
  *) echo "ERROR: unknown GRID module '$module'" >&2; bash "$ROOT/f/help.sh" >&2; exit 2 ;;
esac
shift || true
exec bash "$ROOT/f/${module}.sh" "$@"
