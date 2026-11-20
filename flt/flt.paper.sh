#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    ""|-h|--help|help)
        cat <<'HELP'
FLT publication: render manuscript figures and export source tables.

Usage:
  ./flt.paper.sh run|fig1|fig2|fig3|fig4|fig5|tables [options]

Modules:
  run                             Render fig1 through fig5, then export tables.
  fig1 | fig2 | fig3 | fig4 | fig5 Render the selected manuscript figure.
  tables                          Export publication source tables.

Options:
  --dir-output DIR                Results root containing completed FLT analyses.
  --submit auto|local|hpc         Execution environment [auto: detect local/LSF].
  --cores N                       CPU budget [local: 16; HPC: 40].
  --queue NAME                    HPC CPU queue [33].
  --wait TRUE|FALSE               Wait for workflow completion [FALSE].
  --job-memory-gb N               Job memory limit in GiB.
  --job-hours N                   Job duration limit in hours.
  --dir-jobs DIR                  Directory for generated commands and plans.
  --plan-only                     Write a plan without executing or submitting.
  --dry-run                       Print scientific commands without executing.
  --check-inputs                  Check external files only; do not run or submit.
  -h, --help                      Show this help.

Environment:
  FLT_EXTERNAL_ROOT              External resources [DIR_FILES/flt].
  FLT_ATLAS_DIR                  GWAS Atlas rsID.csv directory.
  FLT_METASCAPE_RESULTS           Native Metascape Enrichment_GO/_FINAL_GO.csv.
  FLT_RSCRIPT, FLT_PYTHON         Override scientific interpreters.
  FLT_FUMA_DIR                    FUMA return directory for genetic figures.
  Shared scripts/0f selects local/HPC paths; keep diet/ beside flt/ and 0f/.
  Run the required flt.sh analyses and obtain external FUMA results first.

Examples:
  ./flt.paper.sh run
  ./flt.paper.sh run --check-inputs
  ./flt.paper.sh run --plan-only
  ./flt.paper.sh fig1 --submit local --dry-run
  ./flt.paper.sh run --submit hpc --queue 33 --cores 40 --wait TRUE
  ./flt.paper.sh tables
HELP
        exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
source "$(dirname -- "${BASH_SOURCE[0]}")/f/entry.sh"
