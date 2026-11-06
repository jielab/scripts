#!/usr/bin/env bash
# Cheap transfer inventory. Does not submit jobs, open genotype contents or create output directories.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$ROOT/0os.sh"
workflow_args "$@"; set -- "${WORKFLOW_ARGS[@]}"
exec "${SHARED_PYTHON:-python3}" "$ROOT/check_paths.py" "$@"
