#!/usr/bin/env bash
set -euo pipefail


ENV_NAME="${PFI_CONDA_ENV_NAME:-pfi}"
if command -v mamba >/dev/null 2>&1; then
  mamba env create -f environment.yml || mamba env update -f environment.yml
elif command -v conda >/dev/null 2>&1; then
  conda env create -f environment.yml || conda env update -f environment.yml
else
  echo "ERROR: conda or mamba is required." >&2
  exit 1
fi
conda run -n "$ENV_NAME" python - <<'PY'
import torch, transformers, duckdb, pyarrow, pandas, sklearn
print('OK: pfi Python environment is ready')
PY
