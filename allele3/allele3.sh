#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Usage: ./allele3.sh [all|hotspots|density|herv_enrichment|base_simulation|base_occupancy MODULE] [OPTIONS]

Modules:
    all              One CPU job: downstream analyses using completed scan and PCA
    hotspots         Genome-wide hotspot discovery and matched comparisons [default]
    density          Continuous multi-ALT density from a completed scan
    herv_enrichment  HERV enrichment after scan and PCA
    base_occupancy   Observed nucleotide states, functional and ancient strata
    base_simulation  Four-base mutation, selection and drift simulation

Main options:
    --submit hpc|local|auto
    --out DIR --replace TRUE|FALSE [TRUE]
    --queue NAME --job-cores N --job-memory-gb N --job-hours N
    --job-dir DIR --write-cmd --wait TRUE|FALSE --plan-only
    all: --config JSON [config/workflow.json]
    all checks required PCA files and scan manifest exist; no size/hash audit.
    Missing inputs: print the independent script to run and exit.
    hotspots: --scan-dir DIR --with-atlas auto|TRUE|FALSE --atlas-resources TSV
    -h, --help; --help-full for scientific options

Examples:
    ./allele3.sh all --plan-only
    ./allele3.sh all
    # Alternative: scan and run hotspot analysis only.
    ./3.gnomad_scan.sh scan,verify --wait TRUE && ./allele3.sh
PUBLIC_ENTRY_HELP
    exit 0
fi
# Genome-first discovery; the permanent completed scan is a read-only input.
set -euo pipefail

usage() {
case "${1:-allele3}" in
allele3) bash "${BASH_SOURCE[0]}" --help
;;
herv_enrichment) cat <<'HELP'
Usage: ./allele3.sh herv_enrichment [STAGE] [OPTIONS]

Stages:
    primary -> summary -> spatial -> windows -> alleles -> time
            -> selection -> candidates -> genotypes -> ld -> report
    IBDmix is disabled by default. Inputs include a completed scan and shared PCA.

Options:
    --out DIR --scan-dir DIR --stages CSV --ibdmix TRUE|FALSE
    --metadata TSV --coordinates FILE --basis-dir DIR
    --replace TRUE|FALSE [TRUE]
    --submit hpc|local|auto --queue NAME --job-cores N --job-memory-gb N --job-hours N
    --job-dir DIR --write-cmd --wait TRUE|FALSE --plan-only --help-full

Run scan, PCA and census separately with their numbered entry scripts.
Use allele3.sh all for downstream analyses using the completed inputs.
Use this entry after its inputs exist to run only the HERV analysis.
HELP
;;
*) printf 'Unknown help topic: %s\n' "$1" >&2; exit 2 ;;
esac

}
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Summary help needs only Bash, even without the scientific environment.
case "${1:-}" in
    -h|--help|help) usage allele3; exit 0 ;;
    all)
        case "${2:-}" in -h|--help|help) exec python3 "$ROOT/f/all_workflow.py" --help ;; esac ;;
    herv_enrichment)
        case "${2:-}" in -h|--help|help) usage herv_enrichment; exit 0 ;; esac ;;
esac
# Three primary workflows plus explicit base_occupancy modules; other worker routes are internal.
case "${1:-}" in
    all|hotspots|density|herv_enrichment|base_occupancy|base_simulation|burden|batch|""|-*|help) ;;
    pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report)
        [[ ${ALLELE_BATCH_WORKER:-0} == 1 ]] || { printf 'Use hotspots, density or herv_enrichment; see --help.\n' >&2; exit 2; } ;;
    *) printf 'Unknown module: %s. Use hotspots, density, herv_enrichment or base_occupancy; GWAS uses ../gwas/gwas.sh.\n' "$1" >&2; exit 2 ;;
esac
if [[ ${1:-} == burden ]]; then
    shift
    source "$ROOT/../0f/workflow_env.sh"
    workflow_init allele3-burden "$@"
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    if [[ ${WORKFLOW_ARGS[0]:-} == demo ]]; then
        exec "${ALLELE_PYTHON:-/work/sph-huangj/.conda/envs/allele3-evolution/bin/python}" "$ROOT/f/burden.py" demo "${WORKFLOW_ARGS[@]:1}"
    fi
    exec "${ALLELE_PLANNER_PYTHON:-python3}" "$ROOT/f/burden_plan.py" "${WORKFLOW_ARGS[@]}"
fi
WORKFLOW_NAME=allele3
# shellcheck source=f/0.run.sh
source "$ROOT/f/0.run.sh"

usage_herv() { usage herv_enrichment; }


base_occupancy_main() {
    local module=${1:-help} arg python=${ALLELE_PYTHON:-}
    (($# == 0)) || shift
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local analysis_root=${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}
    local out=$analysis_root/base_occupancy/${module//-/_}
    case "$module" in hotspots|density|base_simulation) out=$analysis_root/$module ;; esac
    if [[ $module == cohort && ${1:-} =~ ^(plan|summarize)$ ]]; then out=$out/$1; fi
    local dry_run=0 native=0 verbose=0 out_explicit=0
    local -a args=("$module")
    case "$module" in
        help|-h|--help) native=1 ;;
        hotspots|genome-windows|density|hotspot-characterize|window-selection|data|allofus|allofus-cache|mfi|mapped-targets|panel|atlas|atlas-allofus|balance|space|allele-time|domains|cohort|archaic-alleles|decompose|opportunity|opportunity-model|functional|strata|joint-counts|regions|base_simulation|review) ;;
        *) fail "Unknown base_occupancy module: $module (use base_occupancy --help)" ;;
    esac
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --verbose-screen) verbose=1; shift ;;
            --dry-run) dry_run=1; shift ;;
            --out) need_value "$@"; out=$2; out_explicit=1; args+=("$1" "$2"); shift 2 ;;
            -h|--help|--help-full) native=1; args+=(--help); shift ;;
            *) args+=("$1"); shift ;;
        esac
    done
    if (( ! native && ! out_explicit )); then args+=(--out "$out"); fi
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    export ALLELE_LOG_STYLE=pipeline WORKFLOW_CONSOLE="$ROOT/f/console.py"
    # launch_workflow separates module logs from atomic stage outputs.
    launch_workflow "$module" "$python" "$ROOT/f/base_occupancy.py" "$out" "$memory" "$interval" \
        "$verbose" "$dry_run" "$((1-native))" "$native" "${args[@]}"
}

herv_main() {
    if [[ ${1:-} == base_occupancy ]]; then shift; base_occupancy_main "$@"; return; fi
    case "${1:-}" in
        -h|--help|help) usage_herv; return 0 ;;
        --version) printf '8.0.0\n'; return 0 ;;
        pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report|"") ;;
        -*) ;;
        *) fail "Unknown analysis stage: $1. PCA, individual HERV and scan have separate 1.pca3D.sh / 2.herv_census.sh / 3.gnomad_scan.sh entry points." ;;
    esac
    local command=pipeline arg python=${ALLELE_PYTHON:-}
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local analysis_root=${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}
    local out=${ALLELE_HERV_FOCUS_OUT:-} verbose=0 native=0
    local -a args=()
    if (($#)) && [[ $1 != -* ]]; then command=$1; shift; fi
    args+=("$command")
    [[ -z ${ALLELE_METADATA:-} ]] || args+=(--metadata "$ALLELE_METADATA")
    [[ -z ${ALLELE_COORDINATES:-} ]] || args+=(--coordinates "$ALLELE_COORDINATES")
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            -h|--help) usage_herv; return 0 ;;
            --help-full) args=(--help); native=1; shift ;;
            --version) printf '8.0.0\n'; return 0 ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --verbose-screen) verbose=1; shift ;;
            --preflight|--dry-run) args+=("$1"); native=1; shift ;;
            --dir-analysis-root) need_value "$@"; analysis_root=$2; args+=("$1" "$2"); shift 2 ;;
            --dir-out|--out) need_value "$@"; out=$2; shift 2 ;;
            *) args+=("$1"); shift ;;
        esac
    done
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    out=$("$python" -c 'import sys; sys.path.insert(0, sys.argv[1]); from herv import resolve_output_dir; print(resolve_output_dir(sys.argv[2], sys.argv[3]))' \
        "$ROOT/f" "$analysis_root" "$out")
    args+=(--dir-out "$out")
    export ALLELE_LOG_STYLE=pipeline WORKFLOW_CONSOLE="$ROOT/f/console.py"
    launch_workflow "$command" "$python" "$ROOT/f/herv.py" "$out" "$memory" "$interval" \
        "$verbose" 0 "$((1-native))" "$native" pipeline "${args[@]}"
}

main() {
    case "${1:-}" in
        -h|--help|help) usage ;;
        --version) printf '8.1.0\n' ;;
        all) fail "The all workflow must be launched through its generated all.cmd" ;;
        base_occupancy) shift; base_occupancy_main "$@" ;;
        hotspots) shift; base_occupancy_main hotspots "$@" ;;
        density) shift; base_occupancy_main density "$@" ;;
        base_simulation) shift; base_occupancy_main base_simulation "$@" ;;
        herv_enrichment) shift; herv_main "$@" ;;
        pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report)
            herv_main "$@" ;;
        ""|-*) base_occupancy_main hotspots "$@" ;;
        *) fail "Unknown module: $1 (use --help; the default is hotspots)" ;;
    esac
}

main "$@"

