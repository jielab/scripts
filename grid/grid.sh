#!/usr/bin/env bash
# GRID transparent agent/reference counting; fitted benchmarks have separate entries.
set -euo pipefail
umask 077
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

source "$ROOT/f/workflow_entry.sh"
grid_workflow_entry grid "$@"
set -- "${GRID_ENTRY_ARGS[@]}"

usage() {
cat <<'HELP'
GRID transparent count ABM:
  prepare       Verify/reuse shared cohort and raw CSx LD-block inputs; freeze
                training-only geometry, donor bank and one reference radius.
  twins         Identify ALL eligible in-radius training donors. Predict by
                their arithmetic mean (binary: cases / all donors). No fitting.
  evaluate      Attach outcomes AFTER prediction freeze; assess full population,
                supported people, uncertainty, actual coverage and composition.
  architecture  Summarize this frozen geometry and available CSx regional weights.
  final         Read only current count_abm completion receipts; aggregate report.
  all           prepare -> twins -> evaluate -> architecture -> final (default).
  self-check    Repository checks.
  cleanup       Existing cache maintenance utility (preview by default).

Examples:
  ./grid.sh all --submit local --traits height,ldl,t2dm --plan-only
  ./grid.sh all --submit hpc --traits height,ldl,t2dm --plan-only
  ./grid.sh all --traits height,ldl,t2dm --device auto
  ./grid.sh prepare --traits height --ld-block-files ALL.37.bed
  ./grid.sh twins --traits height --device cpu
  ./grid.sh evaluate --traits height
  ./grid.sh final --traits height,ldl,t2dm
  ./grid.sh all --input-manifest '/path/{trait}.input.json' --dry-run

Common options:
  --submit local|hpc --jobs-N N --wait TRUE|FALSE --plan-only
  --hpc-config FILE --queue NAME --gpu-queue NAME --job-cores N --job-memory-gb N
  --traits CSV          Default: height,ldl,t2dm
  --analysis-root DIR   Default: /mnt/d/analysis/grid
  --config FILE_OR_DIR  Default: config/count_abm.json
  --block-config FILE   Producer block-input configuration (one or four LD maps)
  --prepared FILE       Existing trusted GRID prepared cohort
  --device auto|cpu|cuda; --threads N; --dry-run
  <module> --help       Complete options

Deployment:
  Set DIR_WORK / DIR_DATA / DIR_SCRATCH for actual mounts and GRID_CONDA_ENV
  (or GRID_PYTHON / GRID_RSCRIPT) for the installed scientific environment.
  Generate plans on the execution machine. After fixing a failed LSF job,
  create a new plan; resume does not automatically retry failed jobs.
  Final report: count_abm/final/index.html; comparisons: count_vs_benchmarks.tsv.
  Comparisons require identical test IDs/outcomes and compatible source receipts.

Individual identities, complete donor memberships and counterexamples are
written only under count_abm/<trait>/_LOCAL_ONLY. Unsupported queries abstain;
a clearly labelled global donor-mean fallback preserves population evaluation.
Support and donor dispersion do not guarantee an individual's prediction error.

Fitted methods remain separate: 2.csx.sh blockprepare/blockridge/blocknet for
block regression; 1.pca.sh, 3.disco.sh and 4.prsformer.sh for their benchmarks.
This entry does not fit blocknet, PRSformer, CountHeads or an error regression.
Historical abm/ and final/ outputs are not used by the count_abm final report.
HELP
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
  prepare|twins|evaluate|architecture|final|all)
    exec "$PYTHON_PATH" "$ROOT/f/count_abm.py" "$module" "$@"
    ;;
  assess|abm)
    echo "GRID assess/abm now refers only to evaluation of frozen transparent counts; use 'all' to build them." >&2
    exec "$PYTHON_PATH" "$ROOT/f/count_abm.py" evaluate "$@"
    ;;
  blocknet)
    echo "Block regression moved to ./2.csx.sh blockridge or blocknet (inputs: blockprepare)." >&2
    echo "The GRID innovation entry predicts directly from observed training-donor counts." >&2
    exit 2
    ;;
  cleanup|self-check)
    exec "$PYTHON_PATH" "$ROOT/f/twin.py" "$module" "$@"
    ;;
  *) echo "ERROR: unknown GRID module: $module" >&2; usage >&2; exit 2 ;;
esac
