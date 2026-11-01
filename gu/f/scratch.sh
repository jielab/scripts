#!/usr/bin/env bash
# Shared OS defaults and common submission options.
submit=${submit:-${SUBMIT:-auto}}
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)/0f/0os.sh"
workflow_args "$@"
set -- "${WORKFLOW_ARGS[@]}"
