#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
PRS-CSx / PRS-CS: chromosome inference, scoring and fitting.

Usage:
  ./2.csx.sh run|MODULE [options]

Modules:
  run                              Cohort + GWAS inputs -> chromosome inference -> merge -> score -> fit.
  infer                            GWAS inputs -> chromosome inference -> merge.
  prepare | inputs                 Prepare cohort or GWAS manifest.
  score | fit                      Use completed upstream outputs.

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
  --chr, --chrs LIST               Autosomes 1-22 or a subset [1-22]; X unsupported.
  --infer-jobs N                   Concurrent model chains per chromosome [4].
  --blas-threads N                 BLAS threads per chain [1].
  --infer-memory-gb N              RAM per chromosome/input/merge worker [16 GiB].
  --split-file FILE|none           Frozen split roster [none; configured seed].
  --gwas-manifest FILE             Explicit GWAS manifest [auto].
  --phi-grid CSV                   Positive fixed phi values [0.01]; also fits auto.
  -h, --help                       Show this help.

Examples:
  ./2.csx.sh run --grch 38 --plan-only
  ./2.csx.sh run --grch 38 --chr 1-22 --traits t2dm
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- all "$@"; fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/run_entry.sh" csx "$@"
