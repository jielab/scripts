#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
ALLELE3: summary/detailed scans and HERV highlights.

Usage:
  ./3.allele3_scan.sh run|MODULE [options]

Modules:
  run                              Three scans, then HERV highlights; requires completed HERV cohort.
  scan | merge | verify            Process or validate one selected source.
  highlights                       Run HERV enrichment from completed scans.

Options:
  --target gnomad|allofus|custom   Source for a single scan.
  --mode summary|detailed          Scan mode; detailed supports gnomad.
  --out DIR                        Scan output directory.
  --jobs N                         Worker count.
  --grch 37|38                     Genome build [38].
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME --cores N           Queue and CPU budget from shared 0f settings.
  --job-hours N                    Maximum job duration.
  --replace TRUE|FALSE             Recompute or reuse completed outputs [FALSE].
  --plan-only                      Write commands without running or submitting.
  --help-full                      Show all engine options.
  -h, --help                       Show this help.

Examples:
  ./3.allele3_scan.sh run --grch 38 --plan-only
  ./3.allele3_scan.sh scan --target gnomad --mode detailed --grch 38
  ./3.allele3_scan.sh highlights --grch 38
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- all "$@"; fi
set -euo pipefail
# Unqualified invocation selects all prepared stores; explicit source/mode selects one.
if (($# == 0)) || [[ $1 == -* ]]; then
    _single=0
    for _arg in "$@"; do
        case "$_arg" in --target|--target=*|--mode|--mode=*|--help-full) _single=1 ;; esac
    done
    if ((_single)); then set -- scan "$@"; else set -- all "$@"; fi
    unset _single _arg
fi
# Reject removed action spellings before loading machine configuration or
# creating a scheduler plan. Dataset names belong only to --target.
case "${1:-scan}" in
    all|scan|merge|verify|highlights|help|-*) ;;
    gnomad|allofus|custom)
        printf 'Use scan --target %s; data sources are not action names.\n' "$1" >&2; exit 2 ;;
    deep)
        printf 'Use scan --target gnomad --mode detailed. Actions are scan, merge, verify.\n' >&2; exit 2 ;;
    *) printf 'Unknown action: %s. Use scan, merge or verify; select data with --target.\n' "$1" >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=allele3-scan
source "$ROOT/f/0.run.sh"

usage() { bash "${BASH_SOURCE[0]}" --help; }

main() {
    local action=scan target=gnomad scan_mode=summary arg python=${ALLELE_PYTHON:-}
    local native=0 full_help=0 verbose=0 explicit_store=0 jobs=
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local replace=${ALLELE_REPLACE:-FALSE} store=
    local -a args=()
    if (($#)) && [[ $1 != -* ]]; then
        case "$1" in
            scan|merge|verify|highlights) action=$1; shift ;;
            help) usage; return 0 ;;
            *) fail "Unknown command: $1 (see --help)" ;;
        esac
    fi
    case "$action" in scan|merge|verify|highlights) ;; *) fail "Unknown action: $action" ;; esac
    if [[ $action == highlights ]]; then target=gnomad,allofus; fi
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            -h|--help) usage; return 0 ;;
            --help-full) full_help=1; native=1; shift ;;
            --target) need_value "$@"; target=$2; shift 2 ;;
            --mode) need_value "$@"; scan_mode=$2; shift 2 ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --verbose-screen) verbose=1; shift ;;
            --out) need_value "$@"; store=$2; explicit_store=1; shift 2 ;;
            --replace) need_value "$@"; replace=${2^^}; shift 2 ;;
            --jobs) need_value "$@"; jobs=$2;
                if [[ $action == highlights ]]; then args+=(--enrichment-jobs "$2"); else args+=(--jobs "$2"); fi
                shift 2 ;;
            --dry-run|--plan-only) native=1; args+=(--plan-only); shift ;;
            *) args+=("$1"); shift ;;
        esac
    done
    [[ $action == highlights ]] || [[ $target == gnomad || $target == allofus || $target == custom ]] || fail '--target must be gnomad, allofus or custom'
    [[ $scan_mode == summary || $scan_mode == detailed ]] || fail '--mode must be summary or detailed'
    [[ $replace == TRUE || $replace == FALSE ]] || fail '--replace must be TRUE or FALSE'
    if [[ -n $jobs ]]; then positive_integer --jobs "$jobs"; workflow_cpu_check "$jobs"; fi
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" "${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3" \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    local entry=$ROOT/f/summary_scan.py label=summary-scan
    if [[ $action == highlights ]]; then
        ((explicit_store)) || store=${HERV_OUT:-${ALLELE_HERV_SCAN_STORE:-${ALLELE_ANALYSIS_ROOT}/herv_scan}}
        export WORKFLOW_LOG_DIR=$store/logs
        entry=$ROOT/f/herv_cohort.py; label=herv-highlights
        args=(highlights --target "$target" --out "$store" --replace "$replace" "${args[@]}")
    elif [[ $scan_mode == detailed ]]; then
        [[ $target == gnomad ]] || fail '--mode detailed supports --target gnomad only'
        [[ $action != merge ]] || fail '--mode detailed supports scan and verify; scan performs its own completion checks'
        ((explicit_store)) || store=${ALLELE_SCAN_STORE:-${ALLELE_ANALYSIS_ROOT}/gnomad_scan}
        export WORKFLOW_LOG_DIR=$store/logs
        entry=$ROOT/f/scan_store.py; label=scan-${action//,/-}
        args=("$action" --out "$store" --replace "$replace" "${args[@]}")
        if ((full_help)); then entry=$ROOT/f/scan.py; args=(--help); fi
    else
        ((explicit_store)) || store=${ALLELE_SUMMARY_STORE:-${ALLELE_ANALYSIS_ROOT}/allele3_scan}/$target
        export WORKFLOW_LOG_DIR=$store/logs
        args=("$action" --target "$target" --out "$store" --replace "$replace" "${args[@]}")
        if ((full_help)); then args=(--help); fi
    fi
    printf '[configured scan store] %s\n' "$store"
    printf '[mode] %s/%s; target=%s; replace=%s\n' "$label" "$action" "$target" "$replace"
    launch_workflow "$label" "$python" "$entry" "$store" "$memory" "$interval" \
        "$verbose" 0 "$((1-native))" "$native" "${args[@]}"
}

main "$@"
