#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
GRID: regional models, ancestry-aware prediction and diagnostics.

Usage:
  ./grid.sh run|MODULE [options]

Modules:
  run                             regions -> architecture / abm / enabled block models;
                                  twins follows abm; includes configured abm_deep.
  regions                         Build regional block scores from existing CSx weights.
  architecture                    Summarize regional architecture.
  block_models                    Fit regional glmnet and LightGBM models (CPU).
  abm                             Fit ancestry-aware neighbour model (CPU).
  abm_deep                        Fit block-token deep model (CPU or GPU).
  twins                           Export ABM neighbour diagnostics.
  robustness                      Geometry sensitivity analysis (opt-in).
  check                           Check required inputs without submitting analysis.
  Publication                     Use ./grid.paper.sh for evaluate and final.

Options:
  --grch 38                       Target genome build [38 only]
  --traits CSV                    height,ldl,t2dm [all three]
  --config FILE                   YAML/JSON configuration [config/grid.yml]
  --analysis-root DIR             Output/input root [${ANALYSIS_ROOT}/grid/38]
  --genotype-dir DIR              Imputed chrN pfiles; alias: --dir-imp
  --consensus-dir DIR             Existing GRCh38 chr*.bed catalog; no LD rebuilding
  --geometry MODE                 consensus [default], regional, global_prs, legacy_unit_block
  --bootstrap N                   Evaluation bootstrap replicates [200]
  --with-deep                     Include abm_deep [enabled in default config]
  --without-deep                  Omit abm_deep from run
  --abm-deep-run MODE             cpu|gpu|auto; auto means GPU [config: cuda]
                                  cpu keeps deep enabled but runs it on CPU
  --submit MODE                   auto|hpc|local [auto]
  --queue NAME                    CPU queue [auto]
  --cores N                       CPU slots per job [shared 0f/0os.sh WORKFLOW_CORES]
  --gpu-queue NAME                Deep GPU queue [config: 2v100-32-e5]
  --gpu-cores N                   CPU slots per deep GPU job [8]
  --gpu-memory-limit-gb N Deep host RAM cap [HPC_GPU_MEMORY_GIB or 160 GiB];
                                  actual reservation is min(cap, configured job memory)
  --job-hours N                   Maximum duration per job [120 h]
  --threads N                     Scientific CPU budget [--cores]
  --jobs N                        Concurrent local jobs [1]; not an HPC concurrency limit
  --replace TRUE|FALSE            Recompute stages / reuse matching completed outputs [FALSE]
  --wait TRUE|FALSE               Wait for submitted HPC workflow [FALSE]; local waits
  --plan-only                     Write commands/plan without submitting; skips input validation

  --job-memory-gb N               Per-job host RAM [HPC 128 / local 32 GiB].
  -h, --help                      Show this help.

Examples:
  ./grid.sh run --grch 38 --traits height,ldl,t2dm --plan-only
  ./grid.sh run --grch 38 --traits height,ldl,t2dm --submit hpc
  ./grid.sh run --grch 38 --traits height --without-deep
  ./grid.sh run --grch 38 --traits height --abm-deep-run cpu
  ./grid.sh run --grch 38 --traits height --gpu-queue 2v100-32-e5 --gpu-cores 8
  ./grid.sh abm_deep --grch 38 --traits height --abm-deep-run gpu
  # After selected model jobs finish:
  ./grid.paper.sh run --grch 38 --traits height --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
case "$1" in
    regions|block_models|abm|twins|abm_deep|architecture|robustness|check|run) ;;
    *) printf 'Unknown module: %s\n' "$1" >&2; usage >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/run_entry.sh" grid "$@"
