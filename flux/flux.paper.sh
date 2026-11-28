#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    ""|help|-h|--help) cat <<'HELP'
HealthFlux | summarize a completed analysis and generate publication figures
Usage: ./flux.paper.sh run [--dir-out DIR] [--cores N] [--config FILE] [--plan-only]
Run after the analysis jobs finish. One CPU job (flux.paper, normally 2 cores)
reads the fitted outputs, builds the main and available ablation summaries,
then renders the publication tables and figures.
--dir-out is the same analysis root used by flux.sh.
Outputs
  main/reports/report.html           Cohort and model summaries
  paper/main/reports/index.html      Publication report and PNG figures
  paper/main/results/                Publication tables
Examples
  ./flux.paper.sh run
  ./flux.paper.sh run --plan-only
HELP
        exit 0 ;;
    run) shift ;;
    *) echo 'Use ./flux.paper.sh -h' >&2; exit 2 ;;
esac
for arg in "$@"; do case "$arg" in -h|--help) exec "$0" -h;; esac; done
FLUX_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export FLUX_ROOT PYTHONDONTWRITEBYTECODE=1
source "$FLUX_ROOT/f/runtime.sh"
flux_runtime "$@"
exec "$FLUX_PYTHON" "$FLUX_ROOT/f/flux.py" paper "$@"
