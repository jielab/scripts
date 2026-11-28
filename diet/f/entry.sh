# Shared execution environment for diet analysis and publication stages.
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
        exec "${DIET_PYTHON:-python3}" "$HELPERS/submit.py" --project diet --entry "$diet_entry" -- "$@"
    fi
    [[ $submit != hpc || -n ${LSB_JOBID:-} ]] || {
        echo 'HPC workers require an LSF allocation.' >&2; exit 2;
    }
fi
exec "${DIET_PYTHON:-python3}" "$ROOT/f/pipeline.py" --entry diet "$@"
