#!/usr/bin/env bash
# Install one shared environment and validate both GRID and PRSformer.
set -euo pipefail
umask 077

usage() {
	cat <<'HELP'
Usage: bash install_grid.sh [--dry-run]

Optional first-time setup or dependency repair. Working installations do not
need this command after code updates. No analysis command installs software.

Installs all PCA / CSx / Disco / GRID / PRSformer dependencies from environment.yml
into one Conda environment, then checks PLINK and both GPU implementations.
The supported stack is Python 3.11, PyTorch 2.12.0 / CUDA 13.0 / NATTEN 0.21.7.
PRSformer's official source is pinned and installed separately from its runtime.
No real cohort is read and no analysis is trained by the installation checks.

Overrides:
  GRID_ENV_NAME          Shared Conda environment name (default grid).
  PRSFORMER_UPSTREAM_DIR  Official source directory (default /mnt/f/software/PRSformer).
  --dry-run              Print installation/check commands without making changes.

After installation:
  cd /mnt/d/scripts/grid
  bash grid.sh abm --check-device
  bash grid.sh prsformer --standalone --check-runtime

With PCA/CSx/Disco inputs ready (default traits: height,ldl,t2dm):
  bash grid.sh prsformer
  bash grid.sh abm
  bash grid.sh final
Skip PRSformer training when the desired completed benchmark already exists.
See bash grid.sh --help for upstream commands.
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
run "$solver" run -n "$ENV_NAME" plink2 --version
run "$solver" run -n "$ENV_NAME" python "$ROOT/f/main.py" --check-device
run "$solver" run -n "$ENV_NAME" bash "$ROOT/grid.sh" prsformer --standalone \
	--python python --upstream-dir "$UPSTREAM" --check-runtime

if [[ $DRY_RUN == TRUE ]]; then
	printf 'PLAN ONLY: no environment, source checkout or data was changed.\n'
else
	printf 'Shared GRID / PRSformer environment %s and both GPU checks passed.\n' "$ENV_NAME"
fi
