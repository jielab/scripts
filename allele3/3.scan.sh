#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# Add to the analysis command below for a 160-GiB job:
#   --submit hpc --queue 33 --job-cores 40 --job-memory-gb 160 --job-hours 24 --write-cmd
# New HPC plans auto-round slots; remove --write-cmd to submit the reviewed scope.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

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
  --out DIR, --store DIR  Output/store [${ANALYSIS_ROOT}/allele3/stores/gnomad_scan; run requires --out]
  --source DIR           Adoption source [${ANALYSIS_ROOT}/allele3/final/gnomad_local]
  --replace TRUE|FALSE   Recalculate or resume [TRUE]
  --dry-run, --preflight Inspect metadata only
  --help-full            Show detailed scanner options
  -h, --help             Show this help

Execution (shared with ld/metag 0f/0os.sh):
  --submit local|hpc    Local process cap or .cmd/bsub DAG submission
  --date-scratch DATE   Dated HPC scratch root; specify the active allocation date
  --write-cmd           Generate reviewable .cmd files without submitting
  --max-running N       HPC concurrency ceiling [32]; independent of --jobs-N
  --jobs-N N            Local parallel jobs [1; maximum 8]
  --wait TRUE|FALSE     Wait for submitted jobs [FALSE]
  --queue NAME --job-cores N --job-memory-gb N --job-hours N
  --job-dir DIR         Stable submission receipts; new plan/input -> new directory
  --retry-failed        Explicit EXIT retry with pending dependency repair

Examples:
    source "${SCRIPTS_ROOT}/0f/0os.sh"
  ./3.scan.sh [run|plan|adopt|verify]
  ./3.scan.sh --out ${ANALYSIS_ROOT}/allele3/stores/gnomad_scan
HELP
}

main() {
    local action=run explicit_action=0 arg python=${ALLELE_PYTHON:-} native=0 full_help=0 verbose=0
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local store=${ALLELE_SCAN_STORE:-${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}/stores/gnomad_scan}
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
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    local entry=$ROOT/f/scan_store.py
    if ((full_help)); then entry=$ROOT/f/scan.py; args=(--help); fi
    printf '[configured scan store] %s\n' "$store"
    printf '[mode] scan-%s; replace=%s; selected store/output=%s\n' "$action" "$replace" "$store"
    launch_workflow "scan-$action" "$python" "$entry" "${ALLELE_SCAN_LOG_OUT:-${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}/scan_operations}" "$memory" "$interval" \
        "$verbose" 0 "$((1-native))" "$native" "${args[@]}"
}

main "$@"

