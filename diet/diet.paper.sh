#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
Diet publication: render figures, assemble tables and build the offline report.

Usage:
  ./diet.paper.sh run|final|web [options]

Modules:
  run                             Run final, then web.
  final                           Render publication figures and numbered .xlsx tables.
  web                             Refresh the viewer from existing analysis output.

Options:
  --study dose|maha|all            Study selection [all].
  --cohort ukb|nhanes|chns|charls|all
                                  Cohort selection [all].
  --config FILE                   Pipeline configuration [config/diet.json].
  --dir-output DIR                Results root containing completed analyses.
  --dir-files DIR                 Dictionary/resource directory.
  --collect-only                  Collect existing figures without rerendering composites.
  --replace TRUE|FALSE            Recompute / reuse verified outputs [TRUE].
  --rscript FILE                  R interpreter; alias --r-bin.
  --python FILE                   Python interpreter; alias --python-bin.
  --submit auto|local|hpc         Execution environment [auto: detect local/LSF].
  --cores N                       CPU budget [local: 16; HPC: 40].
  --queue NAME                    HPC CPU queue [33].
  --wait TRUE|FALSE               Wait for workflow completion [FALSE].
  --job-memory-gb N               Job memory limit in GiB.
  --job-hours N                   Job duration limit in hours.
  --dir-jobs DIR                  Directory for generated commands and plans.
  --plan-only                     Write a plan without executing or submitting.
  --dry-run                       Print scientific commands without executing.
  -h, --help                      Show this help.

Examples:
  # Finish diet.sh analyses before building the paper.
  ./diet.paper.sh run --study all --cohort all
  ./diet.paper.sh run --study all --cohort all --plan-only
  ./diet.paper.sh final --study dose --cohort ukb --submit local --dry-run
  ./diet.paper.sh run --study all --cohort all --submit hpc --wait TRUE
  ./diet.paper.sh web --study all --cohort all
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
if [[ $1 == run ]]; then shift; set -- final,web "$@"; fi
case "$1" in
    final|web|final,web) ;;
    *) printf 'Unknown module: %s\n' "$1" >&2; usage >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
diet_entry="$ROOT/diet.paper.sh"
source "$ROOT/f/entry.sh"
