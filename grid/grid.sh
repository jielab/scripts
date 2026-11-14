#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
GRID ABM analysis (CPU)

Modules:
    all           Full ABM workflow [default]
    prepare       Prepare inputs and freeze reference geometry
    twins         Predict using eligible reference neighbours
    evaluate      Evaluate held-out predictions
    architecture  Summarize the fixed representation
    final         Aggregate completed results
    cleanup       Preview scratch cleanup; --apply deletes files

Main options:
    --traits CSV              Traits [height,ldl,t2dm]
    --analysis-root DIR --config FILE  Results and settings [config/abm.yml]
    --replace TRUE|FALSE      Recompute or reuse verified outputs [TRUE]
    --submit hpc|local|auto    Machine selection [auto]
    --queue NAME --job-cores N --job-memory-gb N --job-hours N --threads N
    --plan-only --wait TRUE|FALSE  Preview or wait for completion

Each trait runs twins -> evaluate; architecture -> final follows all traits.
Use queue auto for the default 256 GiB RAM request. ABM runs on CPU.
Fitted benchmarks use 2.csx.sh, 3.disco.sh and 4.prsformer.sh.

Examples:
    ./grid.sh all --traits height,ldl,t2dm --queue auto --job-cores 40 --plan-only
    ./grid.sh all --traits height,ldl,t2dm --queue auto --job-cores 40
PUBLIC_ENTRY_HELP
    exit 0
fi
# GRID transparent agent/reference counting; fitted benchmarks have separate entries.
set -euo pipefail
umask 077
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

source "$ROOT/f/workflow_entry.sh"
grid_workflow_entry grid "$@"
set -- "${GRID_ENTRY_ARGS[@]}"

usage() {
bash "${BASH_SOURCE[0]}" --help
}

case "${1:-all}" in -h|--help|help) usage; exit 0 ;; esac
runtime=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
if [[ -x "$runtime/bin/python3" ]]; then default_python="$runtime/bin/python3"; else default_python=python3; fi
PYTHON=${GRID_PYTHON:-$default_python}
PYTHON_PATH=$(command -v -- "$PYTHON")
export GRID_PYTHON="$PYTHON_PATH"
export PATH="$(dirname -- "$PYTHON_PATH"):$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$ROOT/f${PYTHONPATH:+:$PYTHONPATH}"
module=all
if [[ $# -gt 0 && $1 != -* ]]; then module=$1; shift; fi
case "$module" in
  prepare|geometry|radius|radius_twins|twins|evaluate|architecture|final|all)
    exec "$PYTHON_PATH" -m abm.count "$module" "$@"
    ;;
  cleanup)
    exec "$PYTHON_PATH" "$ROOT/f/twin.py" "$module" "$@"
    ;;
  *) echo "ERROR: unknown GRID module: $module" >&2; usage >&2; exit 2 ;;
esac
