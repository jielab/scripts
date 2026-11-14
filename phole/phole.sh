#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
PHOLE analysis

Modules:
    cross-check|manifest|specificity-check|fingerprint-check  Input checks
    prepare|edge-check|fit|discover|validate|plot             Analysis and figures
    fingerprint-plan|rebuild-plan                           Review plans

Main options:
    --queue NAME --cores N --memory-gb N  Resources [33 / 40 / 32]
    --dir-out DIR --protein-map FILE     Results and protein mapping
    --replace TRUE|FALSE                 Recompute or reuse outputs [FALSE]
    --plan-only                          Preview commands without submitting

Local/HPC is detected automatically. Comma-separated modules run in one job.
Choose either the full workflow or the resume command; do not submit both.

Examples:
    ./phole.sh cross-check --plan-only
    # Full workflow using existing summary data.
    ./phole.sh cross-check,manifest,specificity-check,fingerprint-check,prepare,edge-check,fit,discover,validate,plot
    # Alternative: resume after the input checks have succeeded.
    ./phole.sh prepare,edge-check,fit,discover,validate,plot
PUBLIC_ENTRY_HELP
    exit 0
fi
# Public PHOLE workflow and submission interface.
usage() {
    bash "${BASH_SOURCE[0]}" --help
}
set -euo pipefail
# Print one concise help page and stop before the Python parser.
for _phole_arg in "$@"; do
    case "$_phole_arg" in -h|--help) usage; exit 0 ;; esac
done
unset _phole_arg
PHOLE_ENTRY_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
phole_runtime_override=${PHOLE_PYTHON:-${WORKFLOW_PYTHON:-}}
source "$PHOLE_ENTRY_DIR/../0f/workflow_env.sh"
workflow_init phole "$@"
# Use the tested project runtime on HPC without activating or changing shared environments.
# Explicit --python and PHOLE_PYTHON remain authoritative; other machines retain shared defaults.
if [[ $submit == hpc && -z ${PHOLE_PYTHON:-} && -x "$DIR_SOFTWARE/venvs/phole-v5/bin/python" ]]; then
    export PHOLE_PYTHON="$DIR_SOFTWARE/venvs/phole-v5/bin/python"
elif [[ $submit == local && -z $phole_runtime_override ]]; then
    source "$PHOLE_ENTRY_DIR/../0f/local_runtime.sh"
    export PHOLE_PYTHON=$(workflow_local_runtime phole python "$(workflow_local_runtime '' python python3)")
fi
unset phole_runtime_override
exec "${PHOLE_PLAN_PYTHON:-python3}" "$PHOLE_ENTRY_DIR/f/pipeline.py" "${WORKFLOW_ARGS[@]}"
