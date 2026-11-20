#!/usr/bin/env bash

# Omitted action uses this entry's configured workflow; explicit actions are unchanged.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    if (($# == 0)); then
        set -- all
    elif [[ $1 == --* && $1 != --help && $1 != --help-full && $1 != --help-all && $1 != --execute && $1 != --worker ]]; then
        set -- all "$@"
    fi
fi

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
--grch 38  Target build [38]; new analysis output: ANALYSIS_ROOT/allele3/38/
./allele3.sh --grch 38 # 完整下游流程；先完成 PCA 和 gnomAD detailed
./allele3.sh MODULE --grch 38 # all|enrich_test|allele_test|hotspots|density|herv_enrichment|base_simulation|base_occupancy
./allele3.sh --replace FALSE --grch 38

# 完整下游流程（复制运行）：
bash allele3.sh --grch 38
# 可选：--cores N；默认读取 0f/0os.sh（HPC 40 / local 16）

Usage: ./allele3.sh [all|enrich_test|allele_test|hotspots|density|...] [OPTIONS]
  all          Default downstream workflow; wait for PCA + gnomAD detailed scan
  enrich_test  Region enrichment; needs target summary + reference
  allele_test  Within-site tests; needs target summary + existing Atlas VCFs
  --out DIR --replace FALSE|TRUE [FALSE]
  --plan-only --queue NAME --cores N --job-hours N
  --help-full  All options for the selected module
HPC is detected automatically. Job name: allele3
The all workflow does not require HERV cohort or AVI scan outputs.
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
    --replace TRUE|FALSE [FALSE]
    --submit hpc|local|auto --queue NAME --cores N --job-hours N
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
        case "${2:-}" in -h|--help|help)
            source "$ROOT/../0f/0os.sh"
            exec "$ALLELE_PLANNER_PYTHON" "$ROOT/f/all_workflow.py" --help ;; esac ;;
    herv_enrichment)
        case "${2:-}" in -h|--help|help) usage herv_enrichment; exit 0 ;; esac ;;
esac
# Three primary workflows plus explicit base_occupancy modules; other worker routes are internal.
case "${1:-}" in
    enrich_test|allele_test|all|hotspots|density|herv_enrichment|base_occupancy|base_simulation|burden|batch|""|-*|help) ;;
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
        exec "${ALLELE_PYTHON:-${SHARED_PYTHON:-python3}}" "$ROOT/f/burden.py" demo "${WORKFLOW_ARGS[@]:1}"
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

research_main() {
    local module=$1; shift
    local target=gnomad reference=alphagenome out= python=${ALLELE_PYTHON:-}
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local -a args=()
    while (($#)); do
        case "$1" in
            --target) need_value "$@";target=$2;args+=("$1" "$2");shift 2;;
            --ref) need_value "$@";reference=$2;args+=("$1" "$2");shift 2;;
            --out) need_value "$@";out=$2;shift 2;;
            --python) need_value "$@";python=$2;shift 2;;
            --memory-limit-gb) need_value "$@";memory=$2;shift 2;;
            --screen-log-interval) need_value "$@";interval=$2;shift 2;;
            --help-full|-h|--help) exec "${python:-python3}" "$ROOT/f/research_tests.py" --help;;
            *) args+=("$1");shift;;
        esac
    done
    if [[ -z $out ]]; then
        out="$ALLELE_ANALYSIS_ROOT/$module/$target"
        [[ $module != enrich_test ]] || out="$ALLELE_ANALYSIS_ROOT/$module/${reference}_${target}"
    fi
    python=$(workflow_python "$python" "$HOME/anaconda3/envs/allele3-evolution/bin/python")
    export WORKFLOW_LOG_DIR="$out/logs"
    workflow_launch "$module" "$python" "$ROOT/f/research_tests.py" "$out" "$memory" "$interval" 0 0 1 \
        "$module" --out "$out" --replace "${ALLELE_REPLACE:-FALSE}" "${args[@]}"
}

herv_main() {
    if [[ ${1:-} == base_occupancy ]]; then shift; base_occupancy_main "$@"; return; fi
    case "${1:-}" in
        -h|--help|help) usage_herv; return 0 ;;
        --version) printf '2026.11.15.3\n'; return 0 ;;
        pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report|"") ;;
        -*) ;;
        *) fail "Unknown analysis stage: $1. PCA, individual HERV and scan have separate 1.ancient_pca.sh / 2.herv_scan.sh / 3.allele3_scan.sh entry points." ;;
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
            --version) printf '2026.11.15.3\n'; return 0 ;;
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
        "$verbose" 0 "$((1-native))" "$native" "${args[@]}"
}

main() {
    case "${1:-}" in
        -h|--help|help) usage ;;
        --version) printf '2026.11.15.3\n' ;;
        enrich_test|allele_test) local module=$1;shift;research_main "$module" "$@" ;;
        all) fail "The all workflow must be launched through its generated all.cmd" ;;
        base_occupancy) shift; base_occupancy_main "$@" ;;
        hotspots) shift; base_occupancy_main hotspots "$@" ;;
        density) shift; base_occupancy_main density "$@" ;;
        base_simulation) shift; base_occupancy_main base_simulation "$@" ;;
        herv_enrichment) shift; herv_main "$@" ;;
        pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report)
            herv_main "$@" ;;
        ""|-*) research_main enrich_test "$@" ;;
        *) fail "Unknown module: $1 (use --help; the default is enrich_test)" ;;
    esac
}

main "$@"

