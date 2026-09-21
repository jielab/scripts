#!/usr/bin/env bash
# Make the supported environment usable for non-interactive GRID calls without
# requiring `conda activate` first.
GRID_CONDA_ENV=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
if [[ -x "$GRID_CONDA_ENV/bin/python3" && ${CONDA_PREFIX:-} != "$GRID_CONDA_ENV" ]]; then
  grid_conda_init="$(dirname -- "$(dirname -- "$GRID_CONDA_ENV")")/etc/profile.d/conda.sh"
  if [[ -f "$grid_conda_init" ]]; then
    # Third-party activation/deactivation hooks may read unset variables (e.g.
    # Anaconda's GeoTIFF hook). Relax nounset only while Conda runs its hooks,
    # then restore the caller's setting on both success and failure.
    grid_conda_nounset=FALSE
    [[ $- != *u* ]] || grid_conda_nounset=TRUE
    set +u
    if source "$grid_conda_init" && conda activate "$GRID_CONDA_ENV"; then
      grid_conda_status=0
    else
      grid_conda_status=$?
    fi
    if [[ $grid_conda_nounset == TRUE ]]; then set -u; fi
    unset grid_conda_nounset
    if (( grid_conda_status != 0 )); then
      printf 'ERROR: Could not activate Conda environment: %s\n' "$GRID_CONDA_ENV" >&2
      unset grid_conda_init
      return "$grid_conda_status"
    fi
    unset grid_conda_status
  else
    export PATH="$GRID_CONDA_ENV/bin:$PATH"
  fi
  unset grid_conda_init
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
