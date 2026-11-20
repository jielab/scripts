#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
PHOLE: protein interaction input checks, analysis and validation.

Usage:
  ./phole.sh run|MODULE [options]

Modules:
  run                              Input checks -> prepare -> edge-check -> fit -> discover -> validate -> plot.
  MODULE[,MODULE...]               Select cross-check, manifest, specificity-check, fingerprint-check,
                                   prepare, edge-check, fit, discover, validate or plot.
  fingerprint-plan | rebuild-plan  Write reviewable plans.

Options:
  --queue NAME --cores N --memory-gb N Resources [33 / 40 / 32]
  --dir-out DIR --protein-map FILE Results and protein mapping
  --replace TRUE|FALSE             Recompute or reuse outputs [FALSE]
  --plan-only                      Preview commands without submitting
  -h, --help                       Show this help.

Examples:
  ./phole.sh run --plan-only
  ./phole.sh cross-check --plan-only
  ./phole.sh prepare,edge-check,fit,discover,validate,plot
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- cross-check,manifest,specificity-check,fingerprint-check,prepare,edge-check,fit,discover,validate,plot "$@"; fi
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
