#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
PRSformer analysis

Modules:
    all           Prepare inputs, train both methods, and generate reports
    prepare       Shared cohort and genotype preparation (CPU)
    prsformer     Train PRSformer and generate its report
    prsformer2    Train PRSformer2 and generate its report
    report        Regenerate both reports from completed models (CPU)
    Default: all

Main options:
    --submit hpc|local|auto    Machine selection [auto]
    --traits CSV              Traits [height,ldl,t2dm]
    --replace TRUE|FALSE      Overwrite requested results [TRUE]; FALSE reuses complete outputs
    --analysis-root DIR       Analysis output directory
    --config FILE             Method parameters [config/benchmarks.yml]
    --queue NAME              CPU queue
    --job-cores N             CPU slots per CPU job
    --job-memory-gb N         Host RAM per CPU job (GiB)
    --job-hours N             CPU allocation time limit (hours)
    --threads N               CPU worker threads
    --gpu-queue NAME          GPU queue [4v100-16-e5]
    --gpu-cores N             CPU slots per GPU job
    --gpu-memory-gb N         GPU-node host RAM [48 GiB]; not GPU VRAM
    --gpu-hours N             GPU allocation time limit for both models combined
    --gpu-threads N           CPU worker threads in GPU jobs
    --jobs-N N                Concurrent local jobs / active submissions [8]
    --wait TRUE|FALSE         Wait for the complete submitted DAG [FALSE]
    --plan-only               Write reviewable commands without executing/submitting
    --with-grid               Also run independent CPU GRID analysis (all only)

Training requires one GPU. Scheduling is automatic: one shared CPU preparation
job, then two jobs per trait (GPU: both models in sequence; CPU: both reports).
PRSformer2 is an exploratory extension. Job names use prsformer., without grid.

Low-memory defaults: microbatch 1, accumulation 64 (effective batch 64),
0 loader workers, FP16 and gradient checkpointing. Override via
PRSFORMER_BATCH_SIZE / PRSFORMER_ACCUMULATION_STEPS / PRSFORMER_WORKERS.

Examples:
    ./4.prsformer.sh --traits height,ldl,t2dm --submit hpc --gpu-queue 4v100-16-e5 --gpu-memory-gb 48
    ./4.prsformer.sh --traits height,ldl,t2dm
    ./4.prsformer.sh --traits height,ldl,t2dm --plan-only
    ./4.prsformer.sh report --traits height,ldl,t2dm
PUBLIC_ENTRY_HELP
    exit 0
fi
# PRSformer + PRSformer2 exploratory.
#
# This is deliberately separate from grid.sh:
#   PRSformer  = independent-trait adaptation
#   PRSformer2 = local exploratory extension, NOT part of the GRID method claim
#   prepare    = shared cohort/genotype preparation that GRID may read
set -euo pipefail
umask 077

usage() {
bash "${BASH_SOURCE[0]}" --help
}

case "${1:-all}" in
  -h|--help|help) usage; exit 0 ;;
esac

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/workflow_entry.sh"
grid_workflow_entry prsformer "$@"
set -- "${GRID_ENTRY_ARGS[@]}"
export SCRIPT_FOREGROUND=1
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-${GRID_OUTPUT_ROOT:-${ANALYSIS_ROOT}/grid}/logs}
source "$ROOT/../0f/console.sh"

runtime=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
if [[ -x "$runtime/bin/python3" ]]; then
  default_python="$runtime/bin/python3"
else
  default_python=python3
fi
PYTHON=${GRID_PYTHON:-$default_python}
PYTHON_PATH=$(command -v -- "$PYTHON")
export GRID_PYTHON="$PYTHON_PATH"
export PATH="$(dirname -- "$PYTHON_PATH"):$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$ROOT/f${PYTHONPATH:+:$PYTHONPATH}"

# Retain the large shared preparation until GRID has had an opportunity to
# build/verify its own producer-bound block cache. GRID may later release it
# through its normal verified cleanup.
export PRSFORMER_KEEP_GENOTYPES=1
export GRID_KEEP_SHARED_INPUTS=1

module=all
if [[ $# -gt 0 && $1 != -* ]]; then module=$1; shift; fi
args=("$@")

case "$module" in
  prepare)
    "$PYTHON_PATH" "$ROOT/f/shared_cohort.py" prepare "${args[@]}"
    ;;
  prsformer)
    exec "$PYTHON_PATH" "$ROOT/f/prsformer_workflow.py" prsformer "${args[@]}"
    ;;
  prsformer2)
    exec "$PYTHON_PATH" "$ROOT/f/prsformer_workflow.py" prsformer2 "${args[@]}"
    ;;
  all)
    "$0" prepare "${args[@]}"
    "$0" prsformer "${args[@]}"
    "$0" prsformer2 "${args[@]}"
    ;;
  *)
    echo "ERROR: unknown 4.prsformer.sh module: $module" >&2
    usage >&2
    exit 2
    ;;
esac
