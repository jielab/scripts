#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
GWAS dependencies

Main options:
    --check    Check installed tools and dependencies
    --install  Install the configured environment

Default: --check. --install updates the configured SAIGE environment.

Examples:
    ./install.sh --check
    ./install.sh --install
PUBLIC_ENTRY_HELP
    exit 0
fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$ROOT/../0f/0os.sh"
workflow_args "$@"; set -- "${WORKFLOW_ARGS[@]}"
case ${1:---check} in
  -h|--help) echo 'Usage: bash install.sh [--check|--install] [--submit auto|local|hpc] (default: auto)'; exit 0 ;;
  --check) ;;
  --install) "$OS_GU_CONDA/bin/conda" env update -p "${SAIGE_ENV_DIR:-$OS_GU_CONDA/envs/saige}" -f "$ROOT/env.yml" ;;
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
