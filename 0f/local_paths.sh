#!/usr/bin/env bash
# Shared root initialization for scripts/0LOCAL entry points; never submit HPC work.
local_workflow_init() {
    local helper
    helper=$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")
    export submit=local SUBMIT=local
    source "$helper/0os.sh"
    workflow_args "$@"
    [[ $submit == local ]] || { echo 'scripts/0LOCAL supports --submit local only' >&2; return 2; }
    export SCRIPTS_LOCAL_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[1]}")/.." && pwd -P)"
}
