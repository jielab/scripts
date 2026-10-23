#!/usr/bin/env bash

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Library  Shared functions/configuration; source from the owning entry script

Main options:
    No standalone analysis options. Read the variables/functions below.
    -h, --help  Show this library help

Example:
    source ./conda_env.sh
HELP
			exit 0
			;;
	esac
fi
# Shared discovery; sourcing this file does not activate or install anything.
metag_find_conda() {
  local candidate executable
  if [[ -n "${METAG_CONDA_ROOT:-}" ]]; then
    realpath -m -- "$METAG_CONDA_ROOT"
    return
  fi
  for candidate in "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/mambaforge"; do
    if [[ -r "$candidate/etc/profile.d/conda.sh" && -x "$candidate/bin/conda" ]]; then
      printf '%s\n' "$candidate"
      return
    fi
  done
  executable="${CONDA_EXE:-$(type -P conda || true)}"
  if [[ -n "$executable" ]]; then
    candidate="$(dirname -- "$(dirname -- "$(realpath -- "$executable")")")"
    if [[ -r "$candidate/etc/profile.d/conda.sh" && -x "$candidate/bin/conda" ]]; then
      printf '%s\n' "$candidate"
      return
    fi
  fi
  return 1
}

metag_activate() {
  local conda_root conda_sh restore_nounset=0 status=0
  conda_root="$(metag_find_conda)" || {
    echo 'Conda not found. Run ./install.sh first.' >&2
    return 1
  }
  conda_sh="$conda_root/etc/profile.d/conda.sh"
  [[ -r "$conda_sh" ]] || { printf 'Conda initialization missing: %s\n' "$conda_sh" >&2; return 1; }
  [[ "$-" == *u* ]] && restore_nounset=1
  # Some Conda deactivation hooks read variables absent in a noninteractive shell.
  set +u
  # shellcheck source=/dev/null
  if source "$conda_sh" && conda activate "$conda_root/envs/metag" >&2; then
    export METAG_CONDA_ROOT="$conda_root"
    export METAG_PYTHON="${METAG_PYTHON:-${PYTHON_BIN:-$CONDA_PREFIX/bin/python3}}"
  else
    status=1
  fi
  if ((restore_nounset)); then set -u; fi
  return "$status"
}
