#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    ""|-h|--help|help)
        cat <<'HELP'
FLT: prepare phenotypes and run observational and genetic analyses.

Usage:
  ./flt.sh run|MODULE [options]

Modules (run executes in this order):
  prepare                         Prepare exposure, outcome and lifestyle data.
  cox                             Cox association analysis.
  intake                          Dietary intake correlations.
  behavior                        Other phenotype associations.
  protein                         Protein associations.
  mri                             MRI associations.
  uvmr | mvmr                     Univariable / multivariable Mendelian randomization.
  ldsc                            LDSC / HDL genetic correlations.
  placo                           PLACO analysis.
  smr                             Prepare and run SMR.
  phewas                          PheWAS using returned FUMA and GWAS Atlas exports.
  run                             Execute all modules above sequentially.

Options:
  --dir-output DIR                Results root [ANALYSIS_ROOT/flt].
  --submit auto|local|hpc         Execution environment [auto: detect local/LSF].
  --cores N                       CPU budget [local: 16; HPC: 40].
  --queue NAME                    HPC CPU queue [33].
  --wait TRUE|FALSE               Wait for workflow completion [FALSE].
  --job-memory-gb N               Job memory limit in GiB.
  --job-hours N                   Job duration limit in hours.
  --dir-jobs DIR                  Directory for generated commands and plans.
  --plan-only                     Write a plan without executing or submitting.
  --dry-run                       Print scientific commands without executing.
  --check-inputs                  Check external files only; do not run or submit.
  -h, --help                      Show this help.

Environment:
  FLT_EXTERNAL_ROOT              External resources [DIR_FILES/flt].
  FLT_ATLAS_DIR                  GWAS Atlas rsID.csv directory.
  FLT_METASCAPE_RESULTS           Native Metascape Enrichment_GO/_FINAL_GO.csv.
  FLT_UKB, FLT_GWAS               Override input data directories.
  FLT_DICTIONARY                  Override the FLT dictionary workbook.
  FLT_FUMA_DIR                    FUMA return directory containing leadSNPs.txt.
  FLT_RSCRIPT, FLT_PYTHON         Override scientific interpreters.
  Shared scripts/0f selects local/HPC paths; keep diet/ beside flt/ and 0f/.
  phewas requires external FUMA returns and GWAS Atlas exports; resume it once ready.

Examples:
  ./flt.sh run
  ./flt.sh run --check-inputs
  ./flt.sh run --plan-only
  ./flt.sh run --submit local --cores 16 --dry-run
  ./flt.sh run --submit hpc --queue 33 --cores 40 --wait TRUE
  ./flt.sh phewas
  # After the required analyses finish:
  ./flt.paper.sh run
HELP
        exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
source "$(dirname -- "${BASH_SOURCE[0]}")/f/entry.sh"
