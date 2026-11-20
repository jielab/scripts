#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
METAG: evaluate existing splits and refresh completion/count reports.

Usage:
  ./metag.paper.sh run|MODULE[,MODULE...] [options]

Modules:
  run                             evaluate -> report for selected samples/batches.
  evaluate                        Evaluate existing group/species splits.
  report                          Refresh completion/count reports without splitting.

Options:
  --index-use optimized|legacy|fullIndex profile [optimized].
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
  --evaluation basic|full         basic: QC/reference controls [default]; full also runs
                                  full-library agreement and competitive alignment.
  --stage group|species|both      Evaluation stage selection [both].
  --truth FILE                    Optional truth TSV: ordinal,true_taxid/true_group.
  --evaluation-pairs N            Reservoir limit per group/status stratum (>0).
  --out PATH                      Cohort report destination; requires report in workflow.
  -h, --help                      Show this help.

Examples:
  ./metag.paper.sh run --dat-list Ndesert --batch 1 --index-use optimized --plan-only
  ./metag.paper.sh evaluate --sample Ndesert/SRR23862437 --index-use optimized --evaluation full --plan-only
  ./metag.paper.sh report --dat-list Ndesert --index-use optimized --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|--help-full|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
if [[ $1 == run ]]; then shift; set -- evaluate,report "$@"; fi
if [[ $1 == *, ]]; then usage >&2; exit 2; fi
IFS=, read -ra selected_modules <<< "$1"
for module in "${selected_modules[@]}"; do
    case "$module" in
        evaluate|report) ;;
        *) printf 'Unknown module: %s\n' "$module" >&2; usage >&2; exit 2 ;;
    esac
done
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/env.sh"
exec "$METAG_PYTHON" "$ROOT/f/metag_submit.py" "$@"
