#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Default  Create/update the PFI Conda environment

Main options:
    PFI_CONDA_ENV_NAME  Conda environment [pfi]
    -h, --help  Show help

Example:
    ./01_create_conda_env.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 01_create_conda_env
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
