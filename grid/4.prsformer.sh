#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
GRID: PRSformer preparation, training and prediction.

Usage:
  ./4.prsformer.sh run|MODULE [options]

Modules:
  run                              Prepare -> train -> predict.
  plan | prepare | train | predict Run one PRSformer stage.

Options:
  --grch 38                        Target genome build [38].
  --config FILE                    YAML/JSON configuration [config/grid.yml].
  --analysis-root DIR              Analysis output root.
  --traits CSV                     height,ldl,t2dm [configured traits].
  --genotype-dir DIR               Imputed chrN pfiles; alias: --dir-imp.
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME --cores N           Queue and CPU budget from shared 0f settings.
  --job-hours N                    Maximum job duration.
  --replace TRUE|FALSE             Recompute or reuse completed outputs [FALSE].
  --plan-only                      Write commands without running or submitting.
  --threads N                      Scientific CPU budget [--cores].
  --jobs N                         Concurrent local workers [1].
  --wait TRUE|FALSE                Wait for HPC completion [FALSE].
  -h, --help                       Show this help.

Examples:
  ./4.prsformer.sh run --grch 38 --plan-only
  ./4.prsformer.sh run --grch 38
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- all "$@"; fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
stage=all
case ${1:-} in all|plan|prepare|train|predict) stage=$1;shift ;; esac
exec bash "$ROOT/f/run_entry.sh" prsformer all --prsformer-stage "$stage" "$@"
