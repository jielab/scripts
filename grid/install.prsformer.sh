#!/usr/bin/env bash
# Install the independent PRSformer runtime, then exercise the official model.
set -euo pipefail
umask 077
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
VENV=${GRID_PRSFORMER_VENV:-${HOME}/.venvs/grid-prsformer}
UPSTREAM=${PRSFORMER_UPSTREAM_DIR:-/mnt/f/software/PRSformer}
BOOTSTRAP=${GRID_PRSFORMER_BOOTSTRAP_PYTHON:-}
REVISION=7dea4e1bb27975885c937f2be82f9243bc82bda7
if [[ ${1:-} == --help ]]; then
	cat <<'HELP'
Usage: bash install.prsformer.sh
Installs Python dependencies and the pinned official PRSformer checkout, then
runs a synthetic CUDA forward/backward check (does not read UKB data).
Overrides: GRID_PRSFORMER_VENV, GRID_PRSFORMER_BOOTSTRAP_PYTHON,
PRSFORMER_UPSTREAM_DIR. Supports Python 3.11 or 3.12 on Linux/WSL x86_64.
The modern stack is torch 2.12.0 / CUDA 13.0 / NATTEN 0.21.7.
HELP
	exit 0
fi
[[ $# == 0 ]] || { printf 'ERROR: unknown argument; see --help\n' >&2; exit 2; }
command -v git >/dev/null || { printf 'ERROR: git is required\n' >&2; exit 2; }
if [[ ! -x $VENV/bin/python ]]; then
	if [[ -z $BOOTSTRAP ]]; then
		for candidate in python3.12 python3.11 "$HOME/anaconda3/envs/ai/bin/python" "$HOME/miniforge3/envs/grid/bin/python"; do
			if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys; assert sys.version_info[:2] in ((3,11),(3,12))' 2>/dev/null; then
				BOOTSTRAP=$candidate
				break
			fi
		done
	fi
	[[ -n $BOOTSTRAP ]] || { printf 'ERROR: set GRID_PRSFORMER_BOOTSTRAP_PYTHON to Python 3.11/3.12\n' >&2; exit 2; }
	"$BOOTSTRAP" -m venv "$VENV"
fi
PYTHON="$VENV/bin/python"
"$PYTHON" -c 'import sys; assert sys.version_info[:2] in ((3,11),(3,12)), "Use Python 3.11 or 3.12"'
if [[ ! -e $UPSTREAM ]]; then
	mkdir -p -- "$(dirname -- "$UPSTREAM")"
	git clone https://github.com/23andMe/PRSformer.git "$UPSTREAM"
	git -C "$UPSTREAM" checkout --detach "$REVISION"
else
	[[ $(git -C "$UPSTREAM" rev-parse HEAD) == "$REVISION" ]] || { printf 'ERROR: existing upstream checkout must be at %s; inspect it before changing revisions\n' "$REVISION" >&2; exit 2; }
	[[ -z $(git -C "$UPSTREAM" status --porcelain) ]] || { printf 'ERROR: existing upstream checkout has local changes; inspect them before installation\n' >&2; exit 2; }
fi
"$PYTHON" -m pip install --upgrade pip
"$PYTHON" -m pip install --only-binary=:all: -r "$ROOT/requirements.prsformer.txt" -f https://whl.natten.org
"$PYTHON" -m pip check
bash "$ROOT/3.prsformer.sh" --python "$PYTHON" --upstream-dir "$UPSTREAM" --check-runtime
printf '\nPRSformer installation and GPU check completed.\n'
