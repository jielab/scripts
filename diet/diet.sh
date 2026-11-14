#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Diet analysis

Modules:
    all                    Full analysis including FoodWAS, ML and molecular follow-up
    all --core             Core only; omits inputs required by DOSE Fig2, Fig3a and Fig5
    ukb|nhanes|chns|charls   Cohort phenotypes; scores come from 1.disco.sh
    assoc|bayes|final|web    Individual downstream stages
    foodwas|ml|mechanism|signature  UKB DOSE extensions

Main options:
    --cohort ukb --study dose|maha|all  Cohort and study
    --replace TRUE|FALSE               Recompute or reuse verified outputs [TRUE]
    --jobs N --job-cores N --job-memory-gb N --queue NAME
    --dir-output DIR --dir-files DIR   Results and dictionary directories
    --plan-only                        Preview commands without submitting

Local/HPC paths and R/Python follow shared 0f settings. No arguments show help.
Scoring job: diet.disco; downstream job: diet. --core excludes extensions.
Complete DOSE Fig1 and Fig3f also require NHANES analysis.
MAHA figures require --study maha or --study all and the corresponding cohorts.

Examples:
    ./diet.sh all --cohort ukb --study dose --core --plan-only
    # Generate scores first; wait for successful completion.
    ./1.disco.sh ukb --study dose --replace FALSE
    # Run after scoring is DONE; skip if an existing job already includes this step.
    ./diet.sh all --cohort ukb --study dose --core --replace FALSE
    # Refresh the viewer from existing results.
    ./diet.sh web
PUBLIC_ENTRY_HELP
    exit 0
fi
# Cohorts -> association -> Bayes -> manuscript figures -> offline web report.
set -euo pipefail
usage() {
    bash "${BASH_SOURCE[0]}" --help
}
# Help is pure text: no runtime activation, planning, or analysis.
if (($# == 0)); then usage; exit 0; fi
for arg in "$@"; do
    case "$arg" in -h|--help|help) usage; exit 0 ;; esac
done
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
        exec "${DIET_PYTHON:-python3}" "$ROOT/f/submit.py" --entry diet -- "$@"
    fi
    [[ $submit != hpc || -n ${LSB_JOBID:-} ]] || {
        echo 'HPC workers require an LSF allocation.' >&2; exit 2;
    }
fi
exec "${DIET_PYTHON:-python3}" "$ROOT/f/pipeline.py" --entry diet "$@"
