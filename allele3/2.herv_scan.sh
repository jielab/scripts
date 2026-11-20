#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
ALLELE3: HERV sequence reconstruction and summary.

Usage:
  ./2.herv_scan.sh run|MODULE [options]

Modules:
  run | cohort                     Run the configured HERV cohort study.
  plan | sample | summarize        Plan, process one sample, or summarize.

Options:
  --max-running N                  Concurrent sample jobs [4].
  --sample-cores N                 Cores per sample [16].
  --grch 37|38                     Genome build [38].
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME --cores N           Queue and CPU budget from shared 0f settings.
  --job-hours N                    Maximum job duration.
  --replace TRUE|FALSE             Recompute or reuse completed outputs [FALSE].
  --plan-only                      Write commands without running or submitting.
  --help-full                      Show all engine options.
  -h, --help                       Show this help.

Examples:
  ./2.herv_scan.sh run --grch 38 --plan-only
  ./2.herv_scan.sh summarize --grch 38
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
set -euo pipefail
if [[ ${1:-} == highlights ]]; then
    echo "This module belongs to ./3.allele3_scan.sh highlights" >&2; exit 2
fi
# Public descriptive alias; retain the existing engine action for compatibility.
if [[ ${1:-} == cohort ]]; then shift; set -- run "$@"; fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=herv-cohort
export ALLELE_REPLACE=${ALLELE_REPLACE:-FALSE}
export HERV_COHORT_ONLY=1
source "$ROOT/f/0.run.sh"

main() {
    local action=run python=${HERV_PYTHON:-${ALLELE_PYTHON:-}}
    local out=${HERV_OUT:-${ALLELE_HERV_SCAN_STORE:-${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}/herv_scan}}
    local memory=${HERV_MEMORY_LIMIT_GB:-${ALLELE_MEMORY_LIMIT_GB:-32}}
    local interval=${HERV_SCREEN_LOG_INTERVAL:-${ALLELE_SCREEN_LOG_INTERVAL:-300}}
    local verbose=0 full_help=0 cores=${HERV_CORES:-${ALLELE_ALLOCATED_CORES:-32}} arg
    local -a args=()
    if (($#)) && [[ $1 != -* ]]; then action=$1; shift; fi
    case "$action" in plan|run|sample|summarize) ;; help) bash "${BASH_SOURCE[0]}" --help; return ;; *) fail "Unknown command: $action" ;; esac
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            -h|--help) bash "${BASH_SOURCE[0]}" --help; return ;;
            --help-full) full_help=1; shift ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --out) need_value "$@"; out=$2; shift 2 ;;
            --cores) need_value "$@"; cores=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --verbose-screen) verbose=1; shift ;;
            --replace) need_value "$@"; export ALLELE_REPLACE=$2; shift 2 ;;
            *) args+=("$1"); shift ;;
        esac
    done
    python=$(workflow_python "$python")
    if ((full_help)); then exec "$python" "$ROOT/f/herv_cohort.py" --help; fi
    positive_integer --cores "$cores"
    workflow_cpu_check "$cores"
    export ALLELE_ACTIVE_HERV_BUILD=38 HERV_BUILD=38
    export WORKFLOW_LOG_DIR="$out/logs" WORKFLOW_NAME=herv-cohort
    workflow_launch herv-cohort "$python" "$ROOT/f/herv_cohort.py" "$out" \
        "$memory" "$interval" "$verbose" 0 1 "$action" --out "$out" \
        --cores "$cores" --replace "$ALLELE_REPLACE" "${args[@]}"
}
main "$@"
