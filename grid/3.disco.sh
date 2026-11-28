#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
GRID: DiscoDivas prediction from completed CSX outputs.

Usage:
  ./3.disco.sh run [options]

Modules:
  run | disco                      Fit DiscoDivas using completed cohort, score and fit outputs.

Options:
  --grch 38                        Target genome build [38].
  --config FILE                    YAML/JSON configuration [config/grid.yml].
  --analysis-root DIR              Analysis output root.
  --traits CSV                     height,ldl,t2dm [default: all three].
  --genotype-dir DIR               Imputed chrN pfiles; alias: --dir-imp.
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME --cores N           Queue and CPU budget from shared 0f settings.
  --replace TRUE|FALSE             Recompute or reuse completed outputs [FALSE].
  --plan-only                      Write commands without running or submitting.
  --threads N                      Scientific CPU budget [--cores].
  --jobs N                         Concurrent local workers [1].
  --wait TRUE|FALSE                Wait for HPC completion [FALSE].
  -h, --help                       Show this help.

Outputs:
  analysis/grid/3.disco/main/{height,ldl,t2dm}; existing 3.disco/<trait> trees are reused in place.

Examples:
  ./3.disco.sh run --grch 38 --traits height,ldl,t2dm --plan-only
  ./3.disco.sh run --grch 38 --traits height,ldl,t2dm
PUBLIC_ENTRY_HELP
            exit 0 ;;
        run|disco) ;;
        *) echo "Unknown module: $1; use --help" >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- all "$@"; fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/run_entry.sh" disco "$@"
