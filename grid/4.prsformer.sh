#!/usr/bin/env bash
set -euo pipefail
usage() {
    cat <<'HELP'
PRSformer, using the fixed GRID cohort and the official 23andMe architecture.

  ./4.prsformer.sh run

  run       For each trait prepare dosages on CPU, then train and export predictions.
  plan      Estimate resources using the target SNP manifests and roster.
  prepare   Build or resume TRAIN-QC'd PGEN dosages.
  train     Train, resume a compatible checkpoint, or recover the export.
  predict   Score from the completed checkpoint and export for GRID.

  --config FILE                 Default config/grid.yml
  --analysis-root DIR           Default ANALYSIS_ROOT/grid
  --genotype-dir DIR            Imputed GRCh38 chromosome PGENs (alias --dir-imp)
  --traits CSV                  height,ldl,t2dm
  --grch 38                     Genotype build
  --cores N --threads N         CPU slots and scientific threads
  --queue NAME                  Preparation queue
  --gpu-queue NAME              Training queue [8v100-32-sc]
  --gpu-cores N                 CPU slots for a GPU job [8]
  --jobs N                      Local concurrency [1]
  --plan-only                   Write the task graph without dispatch
  --wait TRUE|FALSE             Wait for LSF jobs [FALSE]. Local runs always wait.
  --replace TRUE|FALSE          Rebuild selected stage [FALSE]
  --submit auto|local|hpc       Auto-detected by 0f/0os.sh
  -h, --help                    Show this help

Fresh 4.prsformer outputs:
  main/<trait>/                 Model, training history and frozen predictions
  tmp/<trait>/                  Dosage matrix and preparation progress
  qc/                           Per-trait export and completion records
  jobs/ and logs/               Submission plans and stdout/stderr
A clean run builds its own cohort and SNP maps from phenotype/PCA/family inputs
and PGENs, then trains into main/<trait>/. No 2.csx outputs are required.
Full-genome neighborhood attention uses a GPU for training and inference.
By default the worker benchmarks full-length microbatches and checkpointing
on the allocated GPU. runtime_tuning.json records timings; TRAIN reports real throughput.
Effective batch size, full variant set, sample splits and TUNE stopping are preserved.
A job interrupted in this version can continue from its saved epoch with the train subcommand.
Set prsformer.upstream_dir to the official checkout before training.
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
    run) stage=all ;;
    plan|prepare|train|predict) stage=$1 ;;
    *) echo "Unknown PRSformer module: $1" >&2; usage >&2; exit 2 ;;
esac
shift
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/run_entry.sh" prsformer all --prsformer-stage "$stage" "$@"
