#!/usr/bin/env bash
# Expensive gnomAD preparation is explicit and independent of allele3.sh.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=allele3-scan
source "$ROOT/f/0.run.sh"

usage() {
    cat <<'HELP'
Usage: ./3.scan.sh [MODULE] [OPTIONS]

Modules:
  run      Run a new scan (default; requires --out DIR)
  plan     Inventory existing results and annotation cache
  adopt    Copy and verify completed results
  verify   Check the permanent store

Options:
  --out DIR, --store DIR  Output/store [/mnt/f/annot/gnomad/allele3; run requires --out]
  --source DIR           Adoption source [/mnt/d/analysis/allele3/final/gnomad_local]
  --replace TRUE|FALSE   Recalculate or resume [TRUE]
  --dry-run, --preflight Inspect metadata only
  --help-full            Show detailed scanner options
  -h, --help             Show this help

Examples:
  ./3.scan.sh [run|plan|adopt|verify]
  ./3.scan.sh --out /mnt/f/annot/gnomad/allele3
HELP
}

main() {
    local action=run explicit_action=0 arg python=${ALLELE_PYTHON:-} native=0 full_help=0 verbose=0
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local store=${ALLELE_SCAN_STORE:-/mnt/f/annot/gnomad/allele3}
    local -a args=()
    local explicit_store=0 replace=${ALLELE_SCAN_REPLACE:-TRUE}
    if (($#)) && [[ $1 != -* ]]; then action=$1; explicit_action=1; shift; fi
    case "$action" in
        help) usage; return 0 ;;
        plan) native=1 ;;
        adopt|verify|run) ;;
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
            --store|--dir-out|--out) need_value "$@"; store=$2; explicit_store=1; args+=(--store "$2"); shift 2 ;;
            --replace) need_value "$@"; replace=${2^^}; shift 2 ;;
            --source) need_value "$@"; args+=(--source "$2"); shift 2 ;;
            --dry-run|--preflight) native=1; args+=("$1"); shift ;;
            *) args+=("$1"); shift ;;
        esac
    done
    [[ $replace == TRUE || $replace == FALSE ]] || fail '--replace must be TRUE or FALSE'
    if [[ $action == run && $full_help == 0 ]]; then
        ((explicit_store)) || fail 'Fresh scan requires --out DIR (or --store DIR). Use --replace FALSE to resume, or verify to check completed results.'
    fi
    args+=(--replace "$replace")
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    ((memory <= 32)) || fail 'Memory cap must be at most 32 GiB on this shared host'
    python=$(workflow_python "$python" /home/huangj/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    local entry=$ROOT/f/scan_store.py
    if ((full_help)); then entry=$ROOT/f/scan.py; args=(--help); fi
    printf '[configured scan store] %s\n' "${ALLELE_SCAN_STORE:-/mnt/f/annot/gnomad/allele3}"
    printf '[mode] scan-%s; replace=%s; selected store/output=%s\n' "$action" "$replace" "$store"
    launch_workflow "scan-$action" "$python" "$entry" "${ALLELE_SCAN_LOG_OUT:-${ALLELE_ANALYSIS_ROOT:-/mnt/d/analysis/allele3}/scan_operations}" "$memory" "$interval" \
        "$verbose" 0 "$((1-native))" "$native" "${args[@]}"
}

main "$@"
