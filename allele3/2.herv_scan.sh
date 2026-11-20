#!/usr/bin/env bash
# Stratified sequence reconstruction + RepeatMasker; completed work is reused.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
--grch 38  Target build [38]; new analysis output: ANALYSIS_ROOT/allele3/38/
./2.herv_scan.sh --grch 38 # 52 人 HERV 分析及汇总
./2.herv_scan.sh MODULE --grch 38 # cohort|plan|sample|summarize
./2.herv_scan.sh --max-running 4 --sample-cores 16 --grch 38

Usage: ./2.herv_scan.sh [cohort|plan|summarize] [OPTIONS]
No action: 52-person HERV sequence study (26 populations x 2).
  --max-running N       Concurrent sample jobs [4]
  --sample-cores N      Cores per sample [16]
  --replace FALSE|TRUE  Reuse completed matching work [FALSE]
  --plan-only           Inspect without submitting
  --queue NAME --job-hours N
  --help-full           All engine options
HPC: 52 sample jobs + dependent summary. Job name: allele3.herv_scan
Enrichment is handled by 3.allele3_scan.sh.
PUBLIC_ENTRY_HELP
    exit 0
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
