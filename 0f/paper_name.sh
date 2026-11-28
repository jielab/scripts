#!/usr/bin/env bash
# All paper jobs use the entry filename plus its module, without option values.
workflow_paper_name() {
    local entry=${1##*/}; shift
    local action=
    [[ ${1:-} == -* || -z ${1:-} ]] || action=$1
    local name=${entry%.sh}
    [[ -z $action ]] || name+=.$action
    name=${name//[^a-zA-Z0-9_.-]/_}
    export WORKFLOW_JOB_NAME=${name:0:120}
}
