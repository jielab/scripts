#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
GRID: evaluate completed predictions and aggregate reports.

Usage:
  ./grid.paper.sh run|evaluate|final [options]

Modules:
  run                             Per-trait evaluate jobs -> one final aggregate job.
  evaluate                        Evaluate completed predictions for selected traits.
  final                           Aggregate completed evaluations and reports.

Options:
  --grch 38                       Target build [38 only].
  --traits CSV                    Traits [config: height,ldl,t2dm].
  --config FILE                   Same YAML/JSON config as analysis [config/grid.yml].
  --analysis-root DIR             Same root as analysis [ANALYSIS_ROOT/grid/38].
  --bootstrap N                   Evaluation bootstrap count [config: 200].
  --with-deep | --without-deep    Match the enabled deep prediction selection.
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    CPU queue [auto].
  --cores N                       Per-job CPU budget from 0f/0os.sh.
  --threads N                     Scientific CPU budget [--cores].
  --job-memory-gb N               Host RAM [HPC 128 / local 32 GiB].
  --job-hours N                   Job duration [120 hours].
  --jobs N                        Concurrent local jobs [1].
  --replace TRUE|FALSE            Recompute / reuse matching outputs [FALSE].
  --wait TRUE|FALSE               Wait for HPC workflow [FALSE]; local waits.
  --plan-only                     Write plan without submitting; skips input validation.
  -h, --help                      Show this help.

Examples:
  ./grid.paper.sh run --grch 38 --traits height,ldl,t2dm --plan-only
  ./grid.paper.sh run --grch 38 --traits height --submit hpc --cores 40
  ./grid.paper.sh evaluate --grch 38 --traits height --without-deep --bootstrap 200
  # Wait for evaluation completion before final when submitted separately.
  ./grid.paper.sh final --grch 38 --traits height
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
case "$1" in
    evaluate|final|run) ;;
    *) printf 'Unknown module: %s\n' "$1" >&2; usage >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ $1 == run ]]; then shift; set -- paper-run "$@"; fi
exec bash "$ROOT/f/run_entry.sh" grid "$@"
