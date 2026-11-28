#!/usr/bin/env bash
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd)/output_dirs.sh"
script_output_directory_args "${BASH_SOURCE[0]}" "$@" || { return 2 2>/dev/null || exit 2; }
set -- "${SCRIPT_OUTPUT_DIRECTORY_ARGS[@]}"

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
GRID: project genotypes onto the reference PCA.

Usage:
  ./1.pca.sh run [options]

Modules:
  run | projection | pca           Project PCA and assign ancestry.

Permanent output directories (legacy options remain accepted):
  --dir-out-pca_projection DIR
                                  YGB PCA/ancestry tables [${DIR_GEN}/ygb/38/pca_proj].

Options:
  --grch 38                        Target genome build [38].
  --config FILE                    YAML/JSON configuration [config/grid.yml].
  --analysis-root DIR              Analysis output root.
  --traits CSV                     height,ldl,t2dm [configured traits].
  --genotype-dir DIR               Imputed chrN pfiles; alias: --dir-imp.
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME --cores N           Queue and CPU budget from shared 0f settings.
  --replace TRUE|FALSE             Recompute or reuse completed outputs [FALSE].
  --plan-only                      Write commands without running or submitting.
  --threads N                      Scientific CPU budget [--cores].
  --jobs N                         Concurrent local workers [1].
  --wait TRUE|FALSE                Wait for HPC completion [FALSE].
  -h, --help                       Show this help.

Examples:
  ./1.pca.sh run --grch 38 --plan-only
  ./1.pca.sh run --grch 38
PUBLIC_ENTRY_HELP
            exit 0 ;;
        run|projection|pca) ;;
        *) echo "Unknown module: $1; use --help" >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- all "$@"; fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/run_entry.sh" pca "$@"
