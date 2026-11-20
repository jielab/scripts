#!/usr/bin/env bash

# Omitted action uses this entry's configured workflow; explicit actions are unchanged.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    if (($# == 0)); then
        set -- all
    elif [[ $1 == --* && $1 != --help && $1 != --help-full && $1 != --help-all && $1 != --execute && $1 != --worker ]]; then
        set -- all "$@"
    fi
fi

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
./1.disco.sh  # 默认评分流程，全部队列和评分
./1.disco.sh ukb,nhanes,chns,charls  # 可删除不需要的队列
./1.disco.sh --score all  # 也可 --score medito,dash,mind；完整名称见下文

Diet scoring
UKB MAHA diet*.rds are published to ${DIR_UKB}/Rdata after successful scoring.
DOSE keeps its separate analysis cache.

Modules:
    all or comma-separated ukb,nhanes,chns,charls

Main options:
    --score all|NAME,NAME    Select diet.XXX.sum totals by XXX (default: all)
    --replace TRUE|FALSE     Recompute or reuse verified outputs
    --cores N --queue NAME
                            CPU default follows 0f/0os.sh (WORKFLOW_CORES)
    --plan-only              Preview commands without submitting

Local/HPC paths and runtimes follow shared 0f settings. Job name: diet.disco.
CHNS provenance: DIET_CHNS/chns.manifest.tsv beside cohort data.
A missing provenance file is logged; scoring does not regenerate it.
Scoring diagnostics stay in .cache; no maha/COHORT directory is created.
Each scoring engine requires its own configured inputs. Wait for successful scoring
before running ./diet.sh all --cohort ukb --study dose --core.

Examples:
    ./1.disco.sh all
    ./1.disco.sh ukb,nhanes,chns,charls # same as "all"
    ./1.disco.sh all --score medito,medi24,dash,mind,hpdi,ahei,phdi,modern,digm,maha,hlcd,hlfd,redii,fds,maha_bal,maha_strict,maha_nodairy,maha_noprotein,upf,dose,dose11,dose9,dose9half,maha_processed,maha_nomeat,medi,amed,chdi,maha_density,maha_residual,dash_fung,dash_density --replace FALSE

The example lists all 32 supported score names across cohorts.
Unsupported scores are skipped and logged. Shared engines may compute additional
scores internally. *.disco.rds contains IDs and selected totals; *.raw.rds keeps
full scoring data for downstream analyses. Same-name totals from different
profiles remain separate (for example, diet.dash.sum and dose.diet.dash.sum).
PUBLIC_ENTRY_HELP
    exit 0
fi
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
