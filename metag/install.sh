#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Metagenomics installation
Install the pinned CPU tools through the local / LSF setup workflow.

Main options:
  --software-only     Install/verify tools without reference or cohort data (put first)
  --bootstrap-python  Create a local .venv and install Python dependencies from env.yml (put first)
  --qc-tools          Install pinned minimap2 for optional species refinement (put first)
  --config FILE       Runtime JSON paths/resources (auto local/HPC profile by default)
  --submit hpc|local|auto  Execution machine [auto]
  --plan-only         Generate commands without running or submitting them
  submit=auto|local|hpc  Optional environment override; automatic by default
  jobs_N=N            Environment variable controlling local concurrency

Examples:
    ./install.sh --plan-only
    ./install.sh
    # Check library commands after installation.
    ./1.lib.sh -h
    # Optional local Python bootstrap.
    ./install.sh --bootstrap-python
PUBLIC_ENTRY_HELP
    exit 0
fi
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ ${1:-} == --qc-tools ]]; then
    shift
    source "$HERE/f/env.sh"
    workflow_args "$@"
    set -- "${WORKFLOW_ARGS[@]}"
    exec "$METAG_PYTHON" "$HERE/f/software_install.py" --qc-only "$@"
fi
if [[ ${1:-} == --bootstrap-python ]]; then
    python3 -m venv "$HERE/.venv"
    mapfile -t python_packages < <(awk '/# pip:bootstrap/ {print $2}' "$HERE/env.yml")
    "$HERE/.venv/bin/python" -m pip install "${python_packages[@]}"
    echo "Python ready: $HERE/.venv/bin/python"
    exit 0
fi
if [[ ${1:-} == --software-only ]]; then
    shift
    source "$HERE/f/env.sh"
    workflow_args "$@"
    set -- "${WORKFLOW_ARGS[@]}"
    exec "$METAG_PYTHON" "$HERE/f/software_install.py" "$@"
fi
case "${1:-}" in
    -h|--help|help)
        bash "${BASH_SOURCE[0]}" --help
        exit 0 ;;
esac
if [[ $# == 0 || ${1:0:1} == - ]]; then
    set -- setup "$@"
fi
if [[ $1 != setup ]]; then
    echo "install.sh prepares software/metadata only; builds use 1.build-lib.sh and analyses use metag.sh." >&2
    exit 2
fi
exec bash "$HERE/f/entry.sh" "$@"
