#!/usr/bin/env bash
set -euo pipefail


# 🚩 Command-line help
usage() {
	cat <<'HELP'
cd /mnt/d/scripts/grid

./install_grid.sh

Installs the shared CSx / Disco / GRID CPU environment and the pinned RDS reader.
PRSformer uses a separate GPU environment; see README.prsformer.md.
Optional ARG tools are managed by the gu workflow and are not required by GRID.
HELP
}
case "${1:-}" in -h | --help | help)
	usage
	exit 0
	;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
ENV_NAME=${ENV_NAME:-grid}
[[ -x "$HOME/miniforge3/bin/conda" ]] && export PATH="$HOME/miniforge3/bin:$PATH"
if command -v mamba >/dev/null 2>&1; then solver=mamba; elif command -v conda >/dev/null 2>&1; then solver=conda; else
	echo 'ERROR: install Miniforge/Conda first.' >&2
	exit 2
fi
if "$solver" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
	"$solver" env update -n "$ENV_NAME" -f "$ROOT/environment.yml"
else "$solver" env create -n "$ENV_NAME" -f "$ROOT/environment.yml"; fi

# GRID uses PLINK for genotypes and a pure Python RDS reader; no ARG build is needed.
"$solver" run -n "$ENV_NAME" python -m pip install 'rdata==1.1.0'
"$solver" run -n "$ENV_NAME" python -c 'import numpy,pandas,scipy,sklearn,matplotlib,openpyxl,rdata,joblib'
"$solver" run -n "$ENV_NAME" plink2 --version
cat <<MSG
Environment imports and PLINK validated.
grid.sh finds the default grid environment and places its bin directory on PATH.
Show workflow examples:
  ./grid.sh -h
MSG
