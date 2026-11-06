#!/usr/bin/env bash
# Called before Conda/model imports/memory scopes from every public GRID entry.
grid_workflow_entry() {
    local entry=$1
    shift
    local helper_root planner_python token help=false
    helper_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
    source "$helper_root/../0f/workflow_env.sh"
    workflow_init grid "$@" || return
    GRID_ENTRY_ARGS=("${WORKFLOW_ARGS[@]}")
    # Workers execute one already scheduled scientific stage, never submit again.
    if [[ ${WORKFLOW_WORKER:-0} == 1 ]]; then
        if [[ $submit == hpc && ( -z ${LSB_JOBID:-} || ${WORKFLOW_RESOURCE_MANAGER:-} != lsf ) ]]; then
            echo 'HPC scientific workers require an allocated workflow .cmd job; no calculation was started.' >&2
            return 2
        fi
        if [[ -n ${GRID_TRAIT_STEPS:-} ]]; then
            "${WORKFLOW_PYTHON:-python3}" "$helper_root/f/trait_workflow.py"
            exit $?
        fi
        return 0
    fi
    for token in "${GRID_ENTRY_ARGS[@]}"; do
        case "$token" in -h|--help|help) help=true ;; esac
    done
    [[ $help == false ]] || return 0
    # A local preview may name a worker interpreter that exists only on HPC.
    # Do not rewrite that intended worker path to this computer's interpreter.
    planner_python=${GRID_PLANNER_PYTHON:-python3}
    command -v "$planner_python" >/dev/null 2>&1 || { echo "Missing planning Python: $planner_python" >&2; return 2; }
    "$planner_python" "$helper_root/f/hpc_workflow.py" --entry "$entry" "${GRID_ENTRY_ARGS[@]}"
    exit $?
}
