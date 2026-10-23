#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# Install one shared environment and validate both GRID and PRSformer.
set -euo pipefail
umask 077

usage() {
	cat <<'HELP'
Modules:
    Default  Install and validate the shared GRID environment

Main options:
    GRID_ENV_NAME  Conda environment [grid]
    PRSFORMER_UPSTREAM_DIR  Source [/mnt/f/software/PRSformer]
    --dry-run  Print commands without changes
    -h, --help  Help; dependency details: 0.README.md

Example:
    # After installation, validate the environment and run with prepared inputs:
    cd /mnt/d/scripts/grid
    bash grid.sh abm --check
    bash grid.sh prsformer --standalone --check-runtime

    With PCA/CSx/Disco inputs ready (default traits: height,ldl,t2dm):
    bash grid.sh prsformer
    bash grid.sh abm
    bash grid.sh final
    Skip PRSformer training when the desired completed benchmark already exists.
    See bash grid.sh --help for upstream commands.

    ./install.sh --dry-run
    ./install.sh
HELP
}

DRY_RUN=FALSE
while (($#)); do
	case $1 in
		-h|--help|help) usage; exit 0 ;;
		--dry-run) DRY_RUN=TRUE; shift ;;
		*) printf 'ERROR: unknown argument %s; see --help\n' "$1" >&2; exit 2 ;;
	esac
done

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
ENV_NAME=${GRID_ENV_NAME:-${ENV_NAME:-grid}}
UPSTREAM=${PRSFORMER_UPSTREAM_DIR:-/mnt/f/software/PRSformer}
REVISION=7dea4e1bb27975885c937f2be82f9243bc82bda7
[[ -x "$HOME/miniforge3/bin/conda" ]] && export PATH="$HOME/miniforge3/bin:$PATH"
if command -v mamba >/dev/null 2>&1; then solver=mamba
elif command -v conda >/dev/null 2>&1; then solver=conda
else
	printf 'ERROR: install Miniforge/Conda first.\n' >&2
	exit 2
fi

run() {
	printf 'RUN '; printf '%q ' "$@"; printf '\n'
	[[ $DRY_RUN == FALSE ]] || return 0
	"$@"
}

# Check an existing checkout before changing the environment. Never reset it.
if [[ -e $UPSTREAM ]]; then
	if [[ $DRY_RUN == FALSE ]]; then
		command -v git >/dev/null || { printf 'ERROR: git is required.\n' >&2; exit 2; }
		[[ $(git -C "$UPSTREAM" rev-parse HEAD) == "$REVISION" ]] || {
			printf 'ERROR: existing PRSformer source must be at %s.\n' "$REVISION" >&2; exit 2;
		}
		[[ -z $(git -C "$UPSTREAM" status --porcelain) ]] || {
			printf 'ERROR: PRSformer checkout has local changes; inspect them first.\n' >&2; exit 2;
		}
	fi
fi

if "$solver" env list | awk -v name="$ENV_NAME" '$1 == name {found=1} END {exit !found}'; then
	run "$solver" env update -n "$ENV_NAME" -f "$ROOT/environment.yml"
else
	run "$solver" env create -n "$ENV_NAME" -f "$ROOT/environment.yml"
fi

if [[ ! -e $UPSTREAM ]]; then
	run mkdir -p -- "$(dirname -- "$UPSTREAM")"
	run "$solver" run -n "$ENV_NAME" git clone https://github.com/23andMe/PRSformer.git "$UPSTREAM"
	run "$solver" run -n "$ENV_NAME" git -C "$UPSTREAM" checkout --detach "$REVISION"
fi

run "$solver" run -n "$ENV_NAME" python -m pip check
run "$solver" run -n "$ENV_NAME" python -c 'import numpy,pandas,scipy,sklearn,matplotlib,openpyxl,rdata,joblib,torch,natten,pgenlib,pyreadr'
run "$solver" run -n "$ENV_NAME" Rscript --vanilla -e 'p <- c("shiny","DT","DBI","RSQLite","ggplot2","jsonlite","testthat"); stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))'
run "$solver" run -n "$ENV_NAME" plink2 --version
run "$solver" run -n "$ENV_NAME" python "$ROOT/f/main.py" --check-device
run "$solver" run -n "$ENV_NAME" bash "$ROOT/grid.sh" prsformer --standalone \
	--python python --upstream-dir "$UPSTREAM" --check-runtime

if [[ $DRY_RUN == TRUE ]]; then
	printf 'PLAN ONLY: no environment, source checkout or data was changed.\n'
else
	printf 'Shared GRID / PRSformer environment %s and both GPU checks passed.\n' "$ENV_NAME"
fi
