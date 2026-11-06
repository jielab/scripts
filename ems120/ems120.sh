#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# EMS120 CPU stages obey the CPU-quarter rule; only the ML stage requests a GPU.
# Review: bash ems120.sh --submit hpc --queue 4a100-40 --cores 8 --memory-gb 96 --gpus 1 --plan-only
# GPU request also includes: -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

# EMS120: CPU prepare -> GPU ML -> CPU report; --submit local runs the same stages.
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
if [[ $submit == hpc && -f $EMS120_SCRIPT_DIR/f/hpc.env ]]; then
    source "$EMS120_SCRIPT_DIR/f/hpc.env"
fi
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

planner=${EMS120_PLANNER_PYTHON:-${WORKFLOW_PYTHON:-${SHARED_PYTHON:-python3}}}
if ! command -v "$planner" >/dev/null 2>&1; then
    if [[ -z ${EMS120_PLANNER_PYTHON:-}${WORKFLOW_PYTHON:-} && -x ${SHARED_PYTHON:-} ]]; then
        planner=$SHARED_PYTHON
    else
        echo "Planning requires Python 3.9+; set EMS120_PLANNER_PYTHON to its executable." >&2
        exit 2
    fi
fi
exec "$planner" "$EMS120_SCRIPT_DIR/f/dispatch.py" "${WORKFLOW_ARGS[@]}"
