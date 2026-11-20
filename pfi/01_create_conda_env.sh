#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
PFI: Create/update the PFI Conda environment.

Usage:
  ./01_create_conda_env.sh run [options]

Modules:
  run                              Create/update the PFI Conda environment.

Options:
  --submit auto|local|hpc          Execution mode [auto].
  --plan-only                      Write a local/LSF .cmd without running
  --job-hours N --queue NAME       
  -h, --help                       Show this help.

Examples:
  ./01_create_conda_env.sh run --plan-only
  ./01_create_conda_env.sh run
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
submit=${submit:-${SUBMIT:-auto}}
source "$_paths_directory/0f/0os.sh"
unset _paths_directory


# Resolve --dir-dat/--dir-big/--dir-tmp through the shared path table.
source "$SCRIPTS_ROOT/0f/project_entry.sh"
project_workflow_init "$@" || exit $?
set -- "${WORKFLOW_ARGS[@]}"
if [[ ${1:-} == run ]]; then shift; fi

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "${SCRIPTS_ROOT}/0f/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			bash "${BASH_SOURCE[0]}" --help
			exit 0
			;;
	esac
fi


# 01_create_conda_env
set -euo pipefail

ENV_NAME="${PFI_CONDA_ENV_NAME:-pfi}"
if command -v mamba >/dev/null 2>&1; then
	mamba env create -f env.yml || mamba env update -f env.yml
elif command -v conda >/dev/null 2>&1; then
	conda env create -f env.yml || conda env update -f env.yml
else
	echo "ERROR: conda or mamba is required." >&2
	exit 1
fi
conda run -n "$ENV_NAME" python - <<'PY'
import torch, transformers, duckdb, pyarrow, pandas, sklearn
print('OK: pfi Python environment is ready')
PY
