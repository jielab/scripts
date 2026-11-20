#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
METAG: group screening and selected-species splitting.

Usage:
  ./metag.sh run|split|MODULE[,MODULE...] [options]

Modules:
  run | split                     split-group -> split-species; no evaluation/report.
  split-group                     Coarse group split only.
  split-species                   Refine selected coarse groups; requires group outputs.
  Publication                     Use ./metag.paper.sh for evaluate and report.

Options:
  --index-use optimized|legacy|full  Index profile [optimized].
                                  optimized: compact lib/group/optimized screening;
                                  legacy: preserved old group index;
                                  full: complete union, expensive random I/O.
  --dat-list CSV                  Dataset names.
  --sample DATASET/SAMPLE         Select one sample; repeatable; excludes --batch/--batch-size.
  --batch SELECTION               Per-dataset batches: 1-2, 1,3-5, 3-, all [all].
  --batch-size N                  Samples per job [configured; default 10].
  --config FILE                   Project configuration override.
  --root DIR                      Project data root override.
  --replace TRUE|FALSE            Recompute / reuse compatible owned outputs [FALSE].
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    LSF queue [38].
  --cores N                       CPU budget [HPC 40 / local 16].
  --threads N                     Classifier threads [max(1, cores-2)], <= cores.
  --mem-gib N                     Host RAM reservation [448 GiB].
  --hours N                       Job duration [120 hours].
  --plan-only                     Generate commands without executing/submitting.
  --allow-local                   Allow local execution where guarded by configuration.
  --groups CSV                    Detailed groups [human,plant]; all supported.
                                  Requires a workflow containing split-species.
  --include-uncertain             Refine uncertain bins with full union (costly);
                                  requires split-species.
  -h, --help                      Show this help.

Examples:
  ./metag.sh run --dat-list Ndesert,Weurope,Wglobal,Wxiamen,changj,emp,neon,ocean,scotland,sub --batch 1 --batch-size 10 --index-use optimized --plan-only
  ./metag.sh run --dat-list Ndesert --batch 1 --index-use legacy
  ./metag.sh split-group --dat-list Ndesert --batch 1 --plan-only
  ./metag.sh split-species --sample Ndesert/SRR23862437 --groups human,plant --index-use optimized --plan-only
  # Wait for splitting to finish before evaluation/report.
  ./metag.paper.sh run --dat-list Ndesert --batch 1 --index-use optimized --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|--help-full|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
if [[ $1 == run || $1 == split ]]; then shift; set -- split-group,split-species "$@"; fi
if [[ $1 == *, ]]; then usage >&2; exit 2; fi
IFS=, read -ra selected_modules <<< "$1"
for module in "${selected_modules[@]}"; do
    case "$module" in
        split-group|split-species) ;;
        *) printf 'Unknown module: %s\n' "$module" >&2; usage >&2; exit 2 ;;
    esac
done
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/env.sh"
exec "$METAG_PYTHON" "$ROOT/f/metag_submit.py" "$@"
