#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
ENV_NAME=${ENV_NAME:-grid-argneedle}
PREFIX=${ARG_NEEDLE_INSTALL_DIR:-$HOME/software/arg-needle}
[[ -x "$HOME/miniforge3/bin/conda" ]] && export PATH="$HOME/miniforge3/bin:$PATH"
if command -v mamba >/dev/null 2>&1; then solver=mamba; elif command -v conda >/dev/null 2>&1; then solver=conda; else echo 'ERROR: install Miniforge/Conda first.' >&2; exit 2; fi
if "$solver" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then "$solver" env update -n "$ENV_NAME" -f "$ROOT/install/environment.grid.yml" --prune
else "$solver" env create -n "$ENV_NAME" -f "$ROOT/install/environment.grid.yml"; fi
mkdir -p "$PREFIX"; cd "$PREFIX"
[[ -d arg-needle-scripts/.git ]] || git clone --recursive https://github.com/PalamaraLab/arg-needle-scripts.git
[[ -d arg-needle/.git ]] || git clone --recursive https://github.com/PalamaraLab/arg-needle.git arg-needle
git -C arg-needle submodule update --init --recursive
# ARG-Needle 1.0.3 is the newest official release compatible with the NumPy
# <2 / tskit 0.6 ABI required by the bundled PRS-CSx workflow. Install only
# runtime wheels; plotting/notebook extras declared by ASMC are not used.
"$solver" run -n "$ENV_NAME" python -m pip install --no-deps \
  'tskit==0.6.4' 'msprime==1.3.4' 'fastcluster==1.2.6' \
  'arg-needle-lib==1.1.3' 'asmc-preparedecoding==2.2.5' \
  'asmc-asmc==1.4.0' 'arg-needle==1.0.3'
"$solver" run -n "$ENV_NAME" python -m pip install click newick demes psutil

# The current official arg-needle package contains the advanced CLI.  Do not put
# the historical scripts checkout first on PYTHONPATH: its empty arg_needle
# package would shadow the installed implementation.
ARG_HOME=$("$solver" run -n "$ENV_NAME" python - <<'PY'
import importlib.util
from pathlib import Path
spec = importlib.util.find_spec("arg_needle.scripts.infer_args_advanced")
if spec is None or spec.origin is None:
    raise SystemExit(1)
print(Path(spec.origin).resolve().parents[2])
PY
)
if ! "$solver" run -n "$ENV_NAME" python -c 'import arg_needle,arg_needle_lib,tskit; from arg_needle import build_arg'; then
  cat >&2 <<MSG
ARG-Needle compiled modules are still unavailable. Follow the upstream installation instructions in:
  $PREFIX/arg-needle/README.md
Then validate with:
  conda run -n $ENV_NAME python -c 'import arg_needle,arg_needle_lib'
MSG
  exit 1
fi
cat <<MSG
Environment and ARG-Needle imports validated.
Activate the environment before running GRID (grid.sh also auto-detects it):
  conda activate $ENV_NAME
  export ARG_NEEDLE_HOME=$ARG_HOME
MSG
