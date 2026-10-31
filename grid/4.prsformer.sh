#!/usr/bin/env bash
# PRSformer + PRSformer2 exploratory.
#
# This is deliberately separate from grid.sh:
#   PRSformer  = independent-trait adaptation
#   PRSformer2 = local exploratory extension, NOT part of the GRID method claim
#   prepare    = shared cohort/genotype preparation that GRID may read
set -euo pipefail
umask 077

usage() {
cat <<'HELP'
4.prsformer.sh modules:
  prepare     Prepare and RETAIN the shared cohort/genotype cache used by
              PRSformer and, when required, GRID LD-block construction.
              Does NOT train PRSformer.
  prsformer   Run/resume PRSformer (independent traits).
  prsformer2  Run/resume PRSformer2 exploratory with separate calibration.
  all         prepare -> prsformer -> prsformer2.

Default module: all (prepare -> prsformer -> prsformer2)

Examples:
  ./4.prsformer.sh --traits height,ldl,t2dm
  ./4.prsformer.sh prepare --traits height,ldl,t2dm
  ./4.prsformer.sh prsformer --traits height,ldl,t2dm --device cuda
  ./4.prsformer.sh prsformer2 --traits height,ldl,t2dm --device cuda
  ./4.prsformer.sh all --traits height,ldl,t2dm --device cuda

Training runtime: GPU batch 16, gradient accumulation 4 (64 subjects/update),
2 loader workers and 2 CPU threads. Existing runs resume from their last saved
epoch, retaining the optimizer and 30-epoch learning-rate schedule.
Override runtime with PRSFORMER_BATCH_SIZE, PRSFORMER_ACCUMULATION_STEPS,
PRSFORMER_WORKERS and --threads; resume requires the same effective batch.
30 epochs is a maximum; validation early stopping remains enabled (patience 3).

Default permanent inputs (outside the analysis output directory):
  READ genotypes: /mnt/f/gen/ukb/37/hap
  READ ancestry: /mnt/f/gen/ukb/pca_proj/ukb.ancestry.auto.tsv.gz
  READ SNP reference: /mnt/f/LD/csx/snpinfo_mult_1kg_hm3
  READ upstream software: /mnt/f/software/PRSformer (PRSFORMER_UPSTREAM_DIR override)
  There is currently NO separate permanent F-drive model/genotype-cache output.
  Prepared data, checkpoints and models remain under --analysis-root.
  KEEP_GENOTYPES retains them there; it does not copy them to F:.
  Deleting the analysis root therefore does not preserve these caches/models.

Scientific boundary:
  GRID blocknet/twin-count does NOT consume PRSformer or PRSformer2 predictions.
  It may consume the shared preparation files retained by "prepare".
  An active local/kinship/ACTIVE.json applies the verified family grouping
  and excludes development relatives without copying the genotype matrix.
HELP
}

case "${1:-all}" in
  -h|--help|help) usage; exit 0 ;;
esac

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export SCRIPT_FOREGROUND=1
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-/mnt/d/analysis/grid/logs}
source "$ROOT/../0f/console.sh"

runtime=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
if [[ -x "$runtime/bin/python3" ]]; then
  default_python="$runtime/bin/python3"
else
  default_python=python3
fi
PYTHON=${GRID_PYTHON:-$default_python}
PYTHON_PATH=$(command -v -- "$PYTHON")
export GRID_PYTHON="$PYTHON_PATH"
export PATH="$(dirname -- "$PYTHON_PATH"):$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$ROOT/f${PYTHONPATH:+:$PYTHONPATH}"

# Retain the large shared preparation until GRID has had an opportunity to
# build/verify its own producer-bound block cache. GRID may later release it
# through its normal verified cleanup.
export PRSFORMER_KEEP_GENOTYPES=1
export GRID_KEEP_SHARED_INPUTS=1

module=all
if [[ $# -gt 0 && $1 != -* ]]; then module=$1; shift; fi
args=("$@")

case "$module" in
  prepare)
    exec "$PYTHON_PATH" "$ROOT/f/prsformer_prepare_shared.py" "${args[@]}"
    ;;
  prsformer)
    exec "$PYTHON_PATH" "$ROOT/f/twin.py" prsformer "${args[@]}"
    ;;
  prsformer2)
    exec "$PYTHON_PATH" "$ROOT/f/twin.py" prsformer2 "${args[@]}"
    ;;
  all)
    "$0" prepare "${args[@]}"
    "$0" prsformer "${args[@]}"
    "$0" prsformer2 "${args[@]}"
    ;;
  *)
    echo "ERROR: unknown 4.prsformer.sh module: $module" >&2
    usage >&2
    exit 2
    ;;
esac
