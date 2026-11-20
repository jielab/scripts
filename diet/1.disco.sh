#!/usr/bin/env bash
case "${1:-}" in
    ""|-h|--help|help)
        cat <<'HELP'
Diet scoring: generate outcome-free diet scores before diet.sh analyses.

Usage:
  ./1.disco.sh run|COHORT[,COHORT...] [options]

Modules:
  run                             Score all configured cohorts.
  ukb | nhanes | chns | charls    Score one cohort or a comma-separated list.

Options:
  --score all|NAME,NAME            Score totals: diet.XXX.sum names or XXX [all].
  --study dose|maha|all            Study selection [all].
  --config FILE                   Pipeline configuration [config/diet.json].
  --dir-output DIR                Results root [configured output directory].
  --dir-files DIR                 Dictionary/resource directory.
  --dir-ukb DIR                   UKB input directory.
  --dir-nhanes DIR                NHANES input directory.
  --dir-nhanes-dose DIR           NHANES DOSE root containing cdc/ and FPED/MPED.
  --dir-chns DIR                  CHNS input directory.
  --dir-charls DIR                CHARLS input directory.
  --charls-input FILE             CHARLS input override.
  --charls-phe FILE               CHARLS phenotype override.
  --charls-config FILE            CHARLS configuration override.
  --foods FILE                    Food input override.
  --replace TRUE|FALSE            Recompute / reuse verified outputs [TRUE].
  --rscript FILE                  R interpreter; alias --r-bin.
  --python FILE                   Python interpreter; alias --python-bin.
  --submit auto|local|hpc         Execution environment [auto: detect local/LSF].
  --cores N                       CPU budget [local: 16; HPC: 40].
  --queue NAME                    HPC CPU queue [33].
  --wait TRUE|FALSE               Wait for workflow completion [FALSE].
  --job-memory-gb N               Job memory limit in GiB.
  --job-hours N                   Job duration limit in hours.
  --dir-jobs DIR                  Directory for generated commands and plans.
  --plan-only                     Write a plan without executing or submitting.
  --dry-run                       Print scientific commands without executing.
  -h, --help                      Show this help.

Examples:
  ./1.disco.sh ukb,nhanes,chns
  ./1.disco.sh run --plan-only
  ./1.disco.sh ukb --submit local --cores 16 --dry-run
  ./1.disco.sh ukb,nhanes,chns --submit hpc --queue 33 --cores 40 --wait TRUE
  # After scoring finishes successfully:
  ./diet.sh run --study all --cohort all
HELP
        exit 0 ;;
    run) shift; set -- all "$@" ;;
    --plan-only) set -- all "$@" ;;
    all) if [[ ${PROJECT_JOB_ACTIVE:-0} != 1 && ${WORKFLOW_WORKER:-0} != 1 ]]; then echo 'Use run for all cohorts.' >&2; exit 2; fi ;;
esac
# Outcome-free diet scoring. Local/HPC roots and submission use scripts/0f.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export DIET_ROOT="$ROOT"
HELPERS=${DIET_HELPERS:-$(dirname -- "$ROOT")/0f}
[[ -f "$HELPERS/0os.sh" ]] || { echo "Place diet/ beside scripts/0f or set DIET_HELPERS." >&2; exit 2; }
source "$HELPERS/0os.sh"
workflow_args "$@"
set -- "${WORKFLOW_ARGS[@]}"
source "$ROOT/f/runtime_env.sh" "$@"
# Preserve original scientific environment settings in reproducible diet job plans.
if [[ -z ${DIET_ENGINE_ENV:-} ]]; then
    DIET_ENGINE_ENV=$("${DIET_PYTHON:-python3}" - <<'PYENV'
import json, os
print(json.dumps({k:v for k,v in os.environ.items() if k.startswith(("DOSE_","MAHA_","CHNS_","NHANES_"))},separators=(",",":")))
PYENV
)
    export DIET_ENGINE_ENV
fi
preview=0
for arg in "$@"; do case "$arg" in -h|--help|help|--dry-run|--list-steps) preview=1;; esac; done
if ((preview==0)) && (($#>0)); then
    plan_only=0
    for arg in "$@"; do case "$arg" in --plan-only|--write-cmd) plan_only=1;; esac; done
    if ((plan_only)) || [[ ${WORKFLOW_WORKER:-0} != 1 && ${PROJECT_JOB_ACTIVE:-0} != 1 ]]; then
        exec "${DIET_PYTHON:-python3}" "$HELPERS/project_job.py" --project diet --entry "$ROOT/1.disco.sh" -- "$@"
    fi
    [[ $submit != hpc || -n ${LSB_JOBID:-} ]] || {
        echo 'HPC workers require an LSF allocation.' >&2; exit 2;
    }
fi
exec "${DIET_PYTHON:-python3}" "$ROOT/f/pipeline.py" --entry disco "$@"
