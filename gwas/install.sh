#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# This installer/checker runs directly; for a CPU installation worker, for example:
#   bsub -q 33 -n 10 -W 4:00 -M 32GB -R "span[hosts=1] rusage[mem=32GB/host]" bash install.sh
# GPU builds/checks require a suitable GPU allocation instead.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$ROOT/../0f/0os.sh"
workflow_args "$@"; set -- "${WORKFLOW_ARGS[@]}"
case ${1:---check} in
  -h|--help) echo 'Usage: bash install.sh [--check|--install] [--submit hpc|local]'; exit 0 ;;
  --check) ;;
  --install) "$OS_GU_CONDA/bin/conda" env update -p "${SAIGE_ENV_DIR:-$OS_GU_CONDA/envs/saige}" -f "$ROOT/environment.yml" ;;
  *) echo "Unknown option: $1" >&2; exit 2 ;;
esac
source "$ROOT/f/runtime.env"
for command in plink2 regenie bcftools samtools Rscript; do
    command -v "$command" >/dev/null || { echo "missing executable $command in PATH=$PATH" >&2; exit 2; }
done
"${GWAS_LD_PYTHON:-python3}" -c 'import numpy,pandas,scipy,yaml; print("GWAS Python dependencies PASS")'
saige=${SAIGE_ENV_DIR:-$OS_GU_CONDA/envs/saige}
[[ -x $saige/bin/Rscript ]] || { echo "missing $saige/bin/Rscript" >&2; exit 2; }
env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE "$saige/bin/Rscript" --vanilla -e 'library(SAIGE); cat("SAIGE", as.character(packageVersion("SAIGE")), "PASS\n")'
Rscript --vanilla -e 'p<-c("data.table","ggplot2","jsonlite"); for(x in p) stopifnot(requireNamespace(x,quietly=TRUE)); cat("GWAS R dependencies PASS\n")'
plink2 --version
bcftools --version | head -2
samtools --version | head -2
echo 'GWAS software check PASS (no cohort data needed)'
