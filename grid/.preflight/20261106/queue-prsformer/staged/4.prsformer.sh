#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# CPU planning flags for the analysis command below:
#   --submit hpc --queue 33 --job-cores 40 --job-memory-gb 160 --plan-only
# New CPU plans auto-round slots. GPU stages keep their separate --gpu-queue resources.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

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
  all         One job per trait: prepare -> prsformer -> prsformer2.

Default module: all (three trait jobs; each runs both methods sequentially).
Scheduler: --submit local|hpc --jobs-N N --wait TRUE|FALSE --plan-only
Resources: --hpc-config FILE --queue NAME --gpu-queue NAME
           --job-cores N --job-memory-gb N --job-hours N

Examples:
  # QiMing: PCA from 1.pca.sh must be complete. CSx (2.csx.sh) and Disco
  # (3.disco.sh) do not need to finish: this workflow uses raw genotypes and
  # phenotypes, and creates/reuses the shared cohort cache itself.
  cd /work/sph-huangj/scripts/grid
  queueinfo
  bash 4.prsformer.sh all --submit hpc --hpc-config config/hpc.qiming.json \
    --config config/abm.qiming-v100.json \
    --queue v3-6t --gpu-queue 8v100-32-sc --gpu-count 1 \
    --device cuda --threads 2 --traits height,ldl,t2dm --wait FALSE
  # Three GPU jobs: height.all, ldl.all, t2dm.all. Each reserves the maximum
  # CPU/RAM required by preparation and either training method:
  #   -q 8v100-32-sc -n 32 -W 120:00 -M 256GB
  #   -R "span[hosts=1] rusage[mem=256GB/host]"
  #   -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"
  # 2026-11-06: 8v100-32-sc is Open (V100 32-GiB GPUs, 384-GiB host RAM),
  # while 4a100-40 is Closed. All 16 V100s were allocated at this check;
  # Open permits submission but does not promise an immediate start.
  # Recheck queueinfo; queues reserved for other groups are not alternatives.
  # The V100 config changes PRSformer2 AMP from BF16 to FP16; V100 cannot
  # run BF16. PRSformer already uses FP16. The installed CUDA 12.6/NATTEN
  # stack passed the official-model forward/backward check on this queue.
  # Shared preparation uses a blocking lock; completed common inputs are reused.
  # 120 hours covers the entire trait job; --job-hours overrides it.
  # Add --plan-only to preview without submitting or preparing inputs.
  # The entry calls bsub once per trait; do not add an outer bsub.
  # The separate prsformer/prsformer2 modules may still run independently;
  # either alone uses the GPU profile's 8 cores / 160 GiB after inputs exist.
  # --queue selects CPU-only modules; --gpu-queue selects these GPU bundles.

  # Other module options:
  ./4.prsformer.sh --traits height,ldl,t2dm
  ./4.prsformer.sh prepare --traits height,ldl,t2dm
  ./4.prsformer.sh prsformer --traits height,ldl,t2dm --device cuda
  ./4.prsformer.sh prsformer2 --traits height,ldl,t2dm --device cuda
  ./4.prsformer.sh all --traits height,ldl,t2dm --device cuda

PRSformer runtime: GPU batch 16, accumulation 4 (64 subjects/update),
2 loader workers; the HPC example explicitly selects 2 CPU threads.
Existing runs resume from their last saved epoch, retaining the optimizer
and 30-epoch learning-rate schedule. Override PRSformer runtime with
PRSFORMER_BATCH_SIZE, PRSFORMER_ACCUMULATION_STEPS, PRSFORMER_WORKERS and
--threads; resume requires the same effective batch. Validation early
stopping remains enabled (patience 3; 30 epochs maximum).
PRSformer2 settings come from the config's "prsformer" section. The V100
profile uses FP16, batch 1, accumulation 64, and retains its existing
learning rate, seed, 30-epoch maximum and early-stopping patience 5.

Default permanent inputs (outside the analysis output directory):
  READ genotypes: ${DIR_GEN}/ukb/37/typ (hap fallback if typ chr1 is absent;
                  override with GRID_TARGET_DIR)
  READ phenotypes: ${DIR_UKB}/Rdata/all.rds
  READ PCA: ${DIR_GEN}/ukb/pca_proj/ukb.discodivas.pca.tsv.gz
  READ family relationships: ${DIR_GEN}/ukb/37/typ/ukb_rel.dat
  READ ancestry: ${DIR_GEN}/ukb/pca_proj/ukb.ancestry.auto.tsv.gz
  READ SNP reference: ${DIR_LD}/csx/snpinfo_mult_1kg_hm3
  READ upstream software: ${DIR_SOFTWARE}/PRSformer (PRSFORMER_UPSTREAM_DIR override)
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
source "$ROOT/f/workflow_entry.sh"
grid_workflow_entry prsformer "$@"
set -- "${GRID_ENTRY_ARGS[@]}"
export SCRIPT_FOREGROUND=1
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-${GRID_OUTPUT_ROOT:-${ANALYSIS_ROOT}/grid}/logs}
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
    exec "$PYTHON_PATH" "$ROOT/f/shared_cohort.py" prepare "${args[@]}"
    ;;
  prsformer)
    exec "$PYTHON_PATH" "$ROOT/f/prsformer_dispatch.py" prsformer "${args[@]}"
    ;;
  prsformer2)
    exec "$PYTHON_PATH" "$ROOT/f/prsformer_dispatch.py" prsformer2 "${args[@]}"
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
