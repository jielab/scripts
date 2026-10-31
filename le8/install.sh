#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# Dependency entry point for LE8. This installer is never called by le8.sh.
# CIGMA uses a separate environment because its pandas/scipy requirements differ.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

# Standard Conda YAML, printed on demand; pip/cran comments select light installs.
environment_recipe() {
    case "$1" in
    le8) cat "$ROOT/environment.yml"
        ;;
    cigma|le8-cigma) cat <<'CIGMA_ENV'
# Full upstream CIGMA environment for new installations. The prepared-matrix
# adapter has also been numerically tested in the local fit-only runtime.
name: le8-cigma
channels:
  - conda-forge
dependencies:
  - python=3.12
  - pip
  - r-base
  - r-optparse
  - r-numderiv
  - r-matrix
  - pip:
      - cigma==1.1.0
CIGMA_ENV
        ;;
    *) echo "ERROR: environment must be le8 or cigma" >&2; return 2 ;;
    esac
}

profile=report
mode=pip
environment_name=le8
dry_run=FALSE
python_bin=${LE8_REPORT_PYTHON:-${PYTHON_BIN:-${HOME}/anaconda3/envs/le8/bin/python3}}
r_bin=${R_BIN:-Rscript}
install_r=TRUE
while (($#)); do
    case "$1" in
    --report) profile=report; shift ;;
    --abm) profile=abm; shift ;;
    --all) profile=abm; shift ;;
    --validation) profile=validation; shift ;;
    --env) mode=conda; environment_name=${2:?--env requires le8 or cigma}; shift 2 ;;
    --print-environment) mode=print; environment_name=${2:?--print-environment requires le8 or cigma}; shift 2 ;;
    --python) python_bin=${2:?--python requires FILE}; shift 2 ;;
    --r-bin) r_bin=${2:?--r-bin requires FILE}; shift 2 ;;
    --no-r) install_r=FALSE; shift ;;
    --dry-run) dry_run=TRUE; shift ;;
    -h|--help)
	cat <<'HELP'
Modules:
    --report / --abm / --all / --validation  Dependency profiles
    --env le8|cigma  Create/update a complete Conda environment

Main options:
    --python FILE --r-bin FILE  Interpreters
    --no-r  Skip R packages
    --print-environment le8|cigma  Print environment recipe
    --dry-run  Print commands without installing
    -h, --help  Help; isolated CIGMA setup: 0.README.md

Example:
    # Dependency groups (choose the group needed for this installation):
    ./install.sh --report
    ./install.sh --abm
    ./install.sh --all
    ./install.sh --validation
    ./install.sh --print-environment le8
    ./install.sh --print-environment cigma

    ./install.sh --env le8
    ./install.sh --env cigma --dry-run
HELP
        exit 0 ;;
    *) echo "ERROR: Unknown installation option: $1" >&2; exit 2 ;;
    esac
done
if [[ $mode == print ]]; then environment_recipe "$environment_name"; exit 0; fi
requirement_file=$(mktemp)
trap 'rm -f -- "$requirement_file"' EXIT
run() {
    if [[ $dry_run == TRUE ]]; then printf '%q ' "$@"; printf '\n'; else "$@"; fi
}
if [[ $mode == conda ]]; then
    environment_recipe "$environment_name" > "$requirement_file"
    [[ $environment_name != cigma ]] || environment_name=le8-cigma
    if [[ $dry_run == TRUE ]]; then cat "$requirement_file"; fi
    conda_bin=${CONDA_EXE:-conda}
    run "$conda_bin" env update --name "$environment_name" --file "$requirement_file"
    exit 0
fi
if ! command -v "$python_bin" >/dev/null 2>&1; then python_bin=python3; fi
environment_recipe le8 | awk -v profile="$profile" '
  / # pip:/ {
    split($0, fields, " # pip:"); group=fields[2]; package=fields[1]
    sub(/^[[:space:]]*-[[:space:]]*/, "", package)
    if (package ~ /^[A-Za-z0-9_.-]+=[^=]/) sub(/=/, "==", package)
    if (group=="report" || (profile!="report" && group=="abm") ||
        (profile=="validation" && group=="validation")) print package
  }
' > "$requirement_file"
if [[ $dry_run == TRUE ]]; then cat "$requirement_file"; fi
run "$python_bin" -m pip install -r "$requirement_file"
if [[ $install_r == TRUE ]]; then
    mapfile -t r_packages < <(environment_recipe le8 | sed -n 's/.* # cran://p')
    run "$r_bin" - "${r_packages[@]}" <<'R_INSTALL'
# Explicit dependency installation; app.R never installs software.
p <- commandArgs(trailingOnly=TRUE)
m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
if (length(m)) install.packages(m, repos="https://cloud.r-project.org")
R_INSTALL
fi
