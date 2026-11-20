#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
GWAS: check or install the analysis environment.

Usage:
  ./install.sh run [options]

Modules:
  run                              Check configured tools and dependencies.

Options:
  --check                          Check installed tools and dependencies
  --install                        Install the configured environment
  -h, --help                       Show this help.

Examples:
  ./install.sh run
  ./install.sh --install
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- --check "$@"; fi
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
saige_rscript=$saige/bin/Rscript
if [[ $submit == hpc ]]; then saige_rscript=$WORKFLOW_RSCRIPT; fi
[[ -x $saige_rscript ]] || { echo "missing $saige_rscript" >&2; exit 2; }
env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE "$saige_rscript" --vanilla -e 'library(SAIGE); cat("SAIGE", as.character(packageVersion("SAIGE")), "PASS\n")'
"${WORKFLOW_RSCRIPT:-Rscript}" --vanilla -e 'p<-c("data.table","ggplot2","jsonlite"); for(x in p) stopifnot(requireNamespace(x,quietly=TRUE)); cat("GWAS R dependencies PASS\n")'
plink2 --version
bcftools --version | head -2
samtools --version | head -2
echo 'GWAS software check PASS (no cohort data needed)'
