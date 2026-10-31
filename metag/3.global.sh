#!/usr/bin/env bash
# Independent global-index route. Does not call or modify the MetaCache route.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1
export PYTHONPATH="$HERE/f${PYTHONPATH:+:$PYTHONPATH}"

METAG_CONDA=""
for candidate in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" "$(command -v conda 2>/dev/null || true)"; do
  if [[ -n "$candidate" && -x "$candidate" ]]; then METAG_CONDA="$candidate"; break; fi
done
if [[ -n "$METAG_CONDA" ]]; then
  METAG_PREFIX="$(dirname "$(dirname "$METAG_CONDA")")/envs/metag"
  if [[ -f "$METAG_PREFIX/conda-meta/history" ]]; then
    set +u
    source "$(dirname "$(dirname "$METAG_CONDA")")/etc/profile.d/conda.sh"
    conda activate "$METAG_PREFIX"
    set -u
  fi
fi
if [[ "${1:-}" == "install" ]]; then
  shift
  exec bash "$HERE/f/install_global_kraken.sh" "$@"
fi
exec python3 "$HERE/f/global_pipeline.py" "$@"
