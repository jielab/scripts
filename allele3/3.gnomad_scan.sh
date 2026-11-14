#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Usage: ./3.gnomad_scan.sh [scan|verify|scan,verify] [OPTIONS]
    scan: scan gnomAD PASS SNVs. verify: validate completed scan outputs.
    scan,verify: scan then verify in the same allocation; fail if either fails.
    Default output: $ALLELE_ANALYSIS_ROOT/gnomad_scan (or $ALLELE_SCAN_STORE).

Options:
    --out DIR --vcf FILE --vcf-dir DIR --chr autosomes|all|1,2,X
    --jobs N --replace TRUE|FALSE --wait TRUE|FALSE
    --queue NAME --job-cores N --job-dir DIR --write-cmd
    -h, --help

--wait TRUE waits for successful job completion, not just submission.

Examples:
    ./3.gnomad_scan.sh scan,verify --wait TRUE && ./allele3.sh
    ./3.gnomad_scan.sh scan
    ./3.gnomad_scan.sh verify --wait TRUE
PUBLIC_ENTRY_HELP
    exit 0
fi
# Expensive gnomAD preparation is explicit and independent of allele3.sh.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=allele3-scan
source "$ROOT/f/0.run.sh"

usage() {
    bash "${BASH_SOURCE[0]}" --help
}

main() {
    local action=scan explicit_action=0 arg python=${ALLELE_PYTHON:-} native=0 full_help=0 verbose=0
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local store=${ALLELE_SCAN_STORE:-${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}/gnomad_scan}
    local -a args=()
    local explicit_store=0 replace=${ALLELE_REPLACE:-TRUE}
    if (($#)) && [[ $1 != -* ]]; then action=$1; explicit_action=1; shift; fi
    case "$action" in
        help) usage; return 0 ;;
        verify|scan|scan,verify) ;;
        *) fail "Unknown command: $action (see --help)" ;;
    esac
    args+=("$action")
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            -h|--help) usage; return 0 ;;
            --help-full) full_help=1; native=1; shift ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --verbose-screen) verbose=1; shift ;;
            --out) need_value "$@"; store=$2; explicit_store=1; args+=(--out "$2"); shift 2 ;;
            --replace) need_value "$@"; replace=${2^^}; shift 2 ;;
            --dry-run|--plan-only) native=1; args+=("$1"); shift ;;
            *) args+=("$1"); shift ;;
        esac
    done
    [[ $replace == TRUE || $replace == FALSE ]] || fail '--replace must be TRUE or FALSE'
    ((explicit_store)) || args+=(--out "$store")
    args+=(--replace "$replace")
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    local entry=$ROOT/f/scan_store.py
    if ((full_help)); then entry=$ROOT/f/scan.py; args=(--help); fi
    printf '[configured scan store] %s\n' "$store"
    printf '[mode] scan-%s; replace=%s; selected store/output=%s\n' "$action" "$replace" "$store"
    launch_workflow "scan-${action//,/-}" "$python" "$entry" "${ALLELE_SCAN_LOG_OUT:-$store}" "$memory" "$interval" \
        "$verbose" 0 "$((1-native))" "$native" "${args[@]}"
}

main "$@"

