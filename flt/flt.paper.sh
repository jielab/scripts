#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    ""|help|-h|--help)
        cat <<'HELP'
Usage: ./flt.paper.sh run [fig1|fig2|fig3|fig4|fig5|genetics|tables] [OPTIONS]

Read completed analyses and write manuscript PNGs, source tables and a
coverage report under paper/main. Individual figure jobs report missing inputs.
  --dir-output DIR               Analysis root [ANALYSIS_ROOT/flt]
  --panel primary|extended      Match the analysis run [primary]
  --cores N                      Maximum CPUs per job
  --queue NAME                   LSF queue [33]
  --wait                         Wait for figure jobs
  --plan-only | --dry-run         Inspect the workflow
  --check-inputs                 Report availability of manuscript inputs
HELP
        exit 0 ;;
esac
source "$(dirname -- "${BASH_SOURCE[0]}")/f/entry.sh"
