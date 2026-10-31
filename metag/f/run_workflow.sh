#!/usr/bin/env bash
# Shared environment setup for the public entry scripts.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1
export PYTHONPATH="$HERE/f${PYTHONPATH:+:$PYTHONPATH}"

CONDA=""
for c in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" "$(command -v conda 2>/dev/null || true)"; do
  if [[ -n "$c" && -x "$c" ]]; then CONDA="$c"; break; fi
done
if [[ -n "$CONDA" ]]; then
  ENV_PREFIX="$(dirname "$(dirname "$CONDA")")/envs/metag"
  if [[ -f "$ENV_PREFIX/conda-meta/history" ]]; then
    set +u
    source "$(dirname "$(dirname "$CONDA")")/etc/profile.d/conda.sh"
    conda activate "$ENV_PREFIX"
    set -u
  fi
fi
source "$HERE/f/cuda_env.sh"
metag_use_cuda || true
if [[ "${1:-}" == "composition" ]]; then
  shift
  exec python3 "$HERE/f/sample_report.py" "$@"
fi
exec python3 "$HERE/f/workflow.py" "$@"
