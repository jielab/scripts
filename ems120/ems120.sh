#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Usage: ./ems120.sh [OPTIONS]
    Prepare -> phone -> geo -> dx -> report. Every stage defaults to CPU.

Main options:
    --submit hpc|local|auto         Execution environment [auto]
    --phone-run cpu|gpu|auto        Phone XGBoost processor [cpu]
    --geo-run cpu|gpu|auto          Address Transformer processor [cpu]
    --dx-run cpu|gpu|auto           Disease Transformer/CV processor [cpu]
    --replace TRUE|FALSE           Overwrite requested outputs [TRUE]
    --report-only                  CPU report using complete data, no fitting
    --preflight                    Input/runtime checks with synthetic training
    --cpu-queue NAME --cpu-cores N --cpu-threads N
    --cpu-memory-gb N --cpu-hours N
    --gpu-queue NAME --gpu-cores N --gpu-memory-gb N --gpu-hours N
    --gpu-model NAME               Optional LSF model constraint
    --dir-dat DIR --dir-big DIR    Shared roots
    --dir-data DIR --dir-analysis DIR --dir-base-model DIR
    --dir-model DIR --dir-files DIR --dir-jobs DIR
    --conda-env NAME --conda-exe FILE --rscript FILE --python-lib DIR
    --start-step all|fig1|fig2|fig3|fig4|fig5|fig6|figs
    --end-step fig1|fig2|fig3|fig4|fig5|fig6|figs
    --run-cv TRUE|FALSE             Five-fold disease validation [TRUE]
    --wait TRUE|FALSE               Wait for submitted work [FALSE]
    --plan-only                     Write commands without executing
    --help                          Show this help

'auto' resolves to CPU for all three methods. GPU selection allocates one GPU
only to that method. Phone selection still compares its original CPU models.

Examples:
    ./ems120.sh --plan-only
    ./ems120.sh --geo-run gpu --dx-run gpu
    ./ems120.sh --report-only --start-step fig2 --end-step fig2
PUBLIC_ENTRY_HELP
    exit 0
fi
# EMS120: sequential prepare -> phone -> geo -> dx -> report, CPU by default.
# This entry point performs only lightweight planning on the login node.
set -euo pipefail

usage() {
    bash "${BASH_SOURCE[0]}" --help
}

case ${1:-} in
    -h|--help) usage; exit 0 ;;
esac
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
export submit=${EMS120_SUBMIT:-${submit:-${SUBMIT:-auto}}}
export SUBMIT=$submit
case $submit in auto|hpc|local) ;; *) echo "EMS120 submit must be auto, hpc or local" >&2; exit 2 ;; esac
source "$EMS120_SHARED_DIR/0os.sh"
args=()
while (($#)); do
    case $1 in
        --submit)
            case ${2:-} in auto|hpc|local) args+=("$1" "$2"); shift 2 ;;
                *) echo "--submit must be auto, hpc or local." >&2; exit 2 ;; esac ;;
        --submit=*)
            case ${1#*=} in auto|hpc|local) args+=("$1"); shift ;;
                *) echo "--submit must be auto, hpc or local." >&2; exit 2 ;; esac ;;
        --wait)
            case ${2:-} in TRUE|FALSE|true|false) args+=("$1" "$2"); shift 2 ;;
                *) echo "--wait requires TRUE or FALSE." >&2; exit 2 ;; esac ;;
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
exec "$planner" "$EMS120_SCRIPT_DIR/f/submit.py" "${WORKFLOW_ARGS[@]}"
