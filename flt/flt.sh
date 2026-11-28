#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    ""|help|-h|--help)
        cat <<'HELP'
Usage: ./flt.sh run [MODULE] [OPTIONS]

Modules: prepare, cox, cox-sensitivity, intake, behavior, protein, mri,
         uvmr, mvmr, ldsc, placo, coloc, smr, phewas.
Without MODULE, run the analysis graph. FUMA-dependent stages run when their
inputs are available. Use ./flt.paper.sh run to generate the manuscript outputs.

  --dir-output DIR              Analysis root [ANALYSIS_ROOT/flt]
  --dir-out-phewas-cache DIR    Persistent cache [DIR_FILES/flt/phewas_cache]
  --panel primary|extended      GWAS panel [primary: five diseases plus edu]
  --cores N                     Maximum analysis CPUs per job [HPC 20 / local 16]
  --whole-node                  Reserve the full queue-33 node for heavy jobs
  --queue NAME                  LSF queue [33]
  --wait                        Wait for LSF jobs and return their status
  --plan-only                   Write the job graph and commands
  --dry-run                     Print the graph without writing files
  --check-inputs                Check the starting inputs for the selected graph

GWAS inputs: FLT_GWAS/{flt,brain,covs}/common/ID/gwas/ID.gz.
Use FLT_YGB, FLT_DICTIONARY and FLT_EXTERNAL_ROOT to override inputs.
HELP
        exit 0 ;;
esac
source "$(dirname -- "${BASH_SOURCE[0]}")/f/entry.sh"
