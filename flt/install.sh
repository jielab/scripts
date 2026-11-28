#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    ""|-h|--help|help)
        cat <<'HELP'
FLT analysis dependencies (same local/HPC runtime as flt.sh).
  ./install.sh run                 Install missing R packages and BWMR; check inputs.
  ./install.sh --check             Check only; no downloads or installation.
  ./install.sh run --queue 33 --cores 1
Options: --cores N, --queue NAME,
         --plan-only, --dry-run, --dir-output DIR, --wait.
Installation waits for completion by default, including on HPC via bsub.
Uses FLT_RSCRIPT, FLT_PYTHON, FLT_R_LIBRARY and FLT_BWMR overrides.
Pre-existing software is reused. Raw data, LD panels and credentials must exist.
This prepares flt.sh run; publication-only MRI rendering is outside this installer.
HELP
        exit 0 ;;
    run) shift; set -- run --wait "$@" ;;
    --check|check) [[ $# == 0 ]] || shift; set -- check --wait "$@" ;;
    *) echo 'Usage: ./install.sh run|--check [options]' >&2; exit 2 ;;
esac
source "$(dirname -- "${BASH_SOURCE[0]}")/f/entry.sh"
