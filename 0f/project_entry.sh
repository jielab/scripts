#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Shared project command and submission helpers

Main options:
    Source this library from a project entry script.
    project_workflow_init "$@"  Initialize the caller workflow.

This is a sourced helper library, not an analysis entry point.

Examples:
    source ./0os.sh
    source ./project_entry.sh
PUBLIC_ENTRY_HELP
    exit 0
fi
# Small adapter for sequential project entries. The scientific script is the
# same in a local .cmd and an LSF .cmd; nested stages reuse the allocation.
project_workflow_init() {
    local helper entry project arg planning_python plan_only=0 help_requested=0 preview=
    helper=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
    entry=$(realpath -- "${BASH_SOURCE[1]}")
    project=$(basename -- "$(dirname -- "$entry")")
    source "$helper/0os.sh" || return
    [[ $(dirname -- "$entry") != "$SCRIPTS_ROOT" ]] || project=$(basename -- "$entry" .sh)
    workflow_args "$@" || return
    if [[ $project == pfi ]]; then
        export PFI_ANALYSIS_ROOT=${PFI_ANALYSIS_ROOT:-$ANALYSIS_ROOT/pfi}
        if [[ $submit == local && -z ${PFI_PYTHON_BIN:-} ]]; then
            source "$helper/local_runtime.sh"
            export PFI_PYTHON_BIN=$(workflow_local_runtime pfi python python3)
        fi
    fi
    [[ $(basename -- "$entry") != pfi.sh || ${#WORKFLOW_ARGS[@]} -gt 0 ]] || return 0
    # Some legacy stage scripts only recognize help as their first argument.
    # Normalize it before worker re-entry and never let a planning flag fall
    # through to a stage which does not parse arbitrary options itself.
    [[ ${WORKFLOW_ARGS[0]:-} != help ]] || help_requested=1
    for arg in "${WORKFLOW_ARGS[@]}"; do
        [[ $arg != -- ]] || break
        case "$arg" in
            -h|--help) help_requested=1 ;;
            --plan-only|--write-cmd) plan_only=1 ;;
            --list-steps|--dry-run) preview=$arg ;;
        esac
    done
    if ((help_requested)); then
        WORKFLOW_ARGS=(--help)
        return 0
    fi
    if ((plan_only == 0)); then
        if [[ -n $preview ]]; then
            # DOSE implements these previews itself, before loading its runtime.
            # Other adapted legacy stages may silently ignore unknown options.
            [[ $project == dose ]] || {
                echo "$preview is not supported by this entry; use --plan-only." >&2
                return 2
            }
            return 0
        fi
        if [[ ${WORKFLOW_WORKER:-0} == 1 || ${PROJECT_JOB_ACTIVE:-0} == 1 ]]; then
            [[ $submit != hpc || -n ${LSB_JOBID:-} ]] || {
                echo 'HPC workers require an LSF allocation.' >&2; return 2;
            }
            return 0
        fi
    fi
    planning_python=${WORKFLOW_PYTHON:-python3}
    command -v "$planning_python" >/dev/null 2>&1 || {
        echo "Planning requires Python 3: $planning_python" >&2; return 2;
    }
    "$planning_python" "$helper/project_job.py" --project "$project" --entry "$entry" -- "${WORKFLOW_ARGS[@]}"
    local code=$?
    exit "$code"
}
