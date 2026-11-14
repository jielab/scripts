#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Check default input paths

Main options:
    --project NAME  all|allele3|grid|gu|gwas|ld|le8|metag [all]
    --chr LIST      Autosomes to inspect [1-22]
    --json          Print machine-readable results

Checks file presence without reading genotype contents or submitting jobs.
Module-specific input validation is still required.

Examples:
    ./check_paths.sh --project all
    ./check_paths.sh --project grid --json
PUBLIC_ENTRY_HELP
    exit 0
fi
# Cheap transfer inventory. Does not submit jobs, open genotype contents or create output directories.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$ROOT/0os.sh"
workflow_args "$@"; set -- "${WORKFLOW_ARGS[@]}"
exec "${SHARED_PYTHON:-python3}" "$ROOT/check_paths.py" "$@"
