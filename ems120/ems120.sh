#!/usr/bin/env bash
# EMS120: default QiMing LSF submission; --submit local runs the same worker.
# This entry point performs only lightweight planning on the login node.
set -euo pipefail
umask 077

export EMS120_SCRIPT_DIR
EMS120_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export EMS120_SHARED_DIR=${EMS120_SHARED_DIR:-$EMS120_SCRIPT_DIR/../0f}
if [[ ! -f $EMS120_SHARED_DIR/0os.sh ]]; then
    echo "Missing $EMS120_SHARED_DIR/0os.sh. Install this project as scripts/ems120 beside the existing scripts/0f directory." >&2
    exit 2
fi
EMS120_SHARED_DIR=$(cd -- "$EMS120_SHARED_DIR" && pwd -P)

# Set the project default before reading the shared table. An explicit local
# CLI option is resolved by workflow_args and recomputes all derived roots.
export submit=${EMS120_SUBMIT:-${submit:-${SUBMIT:-hpc}}}
export SUBMIT=$submit
case $submit in hpc|local) ;; *) echo "EMS120 submit must be hpc or local" >&2; exit 2 ;; esac
source "$EMS120_SHARED_DIR/0os.sh"
args=()
while (($#)); do
    case $1 in
        --submit)
            case ${2:-} in hpc|local) args+=("$1" "$2"); shift 2 ;;
                *) echo "--submit must be hpc or local; automatic local fallback is disabled." >&2; exit 2 ;; esac ;;
        --submit=*)
            case ${1#*=} in hpc|local) args+=("$1"); shift ;;
                *) echo "--submit must be hpc or local; automatic local fallback is disabled." >&2; exit 2 ;; esac ;;
        --wait)
            case ${2:-} in TRUE|FALSE|true|false) args+=("$1" "$2"); shift 2 ;;
                *) args+=(--wait TRUE); shift ;; esac ;;
        *) args+=("$1"); shift ;;
    esac
done
workflow_args "${args[@]}"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

planner=${EMS120_PLANNER_PYTHON:-${WORKFLOW_PYTHON:-python3}}
if ! command -v "$planner" >/dev/null 2>&1; then
    if [[ -z ${EMS120_PLANNER_PYTHON:-}${WORKFLOW_PYTHON:-} && -x ${SHARED_PYTHON:-} ]]; then
        planner=$SHARED_PYTHON
    else
        echo "Planning requires Python 3.9+; set EMS120_PLANNER_PYTHON to its executable." >&2
        exit 2
    fi
fi
exec "$planner" "$EMS120_SCRIPT_DIR/f/dispatch.py" "${WORKFLOW_ARGS[@]}"
