#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
allele3: genome-first downstream analyses on prepared scan/PCA inputs.

Usage:
  ./allele3.sh run|MODULE [options]

Modules:
  run                             Ordered downstream steps from config/workflow.json:
                                  environment check, hotspots, density, HERV enrichment,
                                  allele_time, space and base_simulation.
  hotspots                        Discover/characterize multi-allelic hotspots.
  density                         Genome-window density analysis.
  herv_enrichment                 HERV analysis; optional stage token after module.
  base_simulation                 Base-occupancy simulation workflow.
  base_occupancy SUBCOMMAND       Specialized occupancy/space/time analyses.
  burden                          Burden workflow; use burden --help for its options.
  batch                           Batch workflows; use batch --help for its options.
  Publication                     Use ./allele3.paper.sh for enrich_test and allele_test.

Options:
  --grch 38                       Target build [38]; uses the shared build-specific root.
  --out DIR                       Output root/directory for the selected workflow.
  --replace TRUE|FALSE            Recompute / reuse completed results [FALSE].
  --python FILE                   Scientific interpreter [ALLELE_PYTHON].
  --job-memory-gb N               Allocation RAM [run: HPC 160 / local 32 GiB].
  --memory-limit-gb N             Scientific module RAM (individual modules only).
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    LSF queue.
  --cores N                       CPU budget from shared 0f/0os.sh.
  --job-hours N                   Job duration limit.
  --job-dir DIR                   Generated command/job directory.
  --write-cmd                     Write a reviewable command file.
  --wait TRUE|FALSE               Wait for submitted work.
  --plan-only                     Write plan without executing/submitting.
  --help-full                     Native full help for the selected scientific module.
  --config FILE                   run only: ordered steps [config/workflow.json].
  --scan-dir DIR                  Module input: completed detailed scan.
  --metadata TSV                  HERV/module sample metadata.
  --coordinates FILE              HERV/module shared PCA coordinates.
  --basis-dir DIR                 HERV/module PCA basis directory.
  --stages CSV                    HERV stage selection.
  --ibdmix TRUE|FALSE             HERV IBDmix branch [FALSE].
  --screen-log-interval N         Scientific progress interval.
  --verbose-screen                Detailed scientific progress.
  -h, --help                      Show this help.

Examples:
  ./allele3.sh run --grch 38 --replace FALSE --plan-only
  ./allele3.sh run --grch 38 --submit hpc --queue 33 --cores 40
  ./allele3.sh hotspots --grch 38 --scan-dir /work/sph-huangj/analysis/allele3/gnomad_scan --plan-only
  ./allele3.sh herv_enrichment --help
  ./allele3.sh base_occupancy --help
  ./allele3.paper.sh run --target gnomad --grch 38 --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
# Genome-first discovery; the permanent completed scan is a read-only input.
set -euo pipefail

module_usage() {
case "${1:-}" in
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

Use allele3.sh run for downstream analyses using the completed inputs.
Use this entry after its inputs exist to run only the HERV analysis.
HELP
;;
*) printf 'Unknown help topic: %s\n' "$1" >&2; exit 2 ;;
esac

}
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Summary help needs only Bash, even without the scientific environment.
case "${1:-}" in
    run)
        case "${2:-}" in -h|--help|help)
            source "$ROOT/../0f/0os.sh"
            exec "$ALLELE_PLANNER_PYTHON" "$ROOT/f/all_workflow.py" --help ;; esac ;;
    herv_enrichment)
        case "${2:-}" in -h|--help|help) module_usage herv_enrichment; exit 0 ;; esac ;;
esac
# Three primary workflows plus explicit base_occupancy modules; other worker routes are internal.
case "${1:-}" in
    run|hotspots|density|herv_enrichment|base_occupancy|base_simulation|burden|batch|--version|help) ;;
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

usage_herv() { module_usage herv_enrichment; }


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
        run) fail "The workflow must be launched through its generated command file" ;;
        base_occupancy) shift; base_occupancy_main "$@" ;;
        hotspots) shift; base_occupancy_main hotspots "$@" ;;
        density) shift; base_occupancy_main density "$@" ;;
        base_simulation) shift; base_occupancy_main base_simulation "$@" ;;
        herv_enrichment) shift; herv_main "$@" ;;
        pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report)
            herv_main "$@" ;;
        *) fail "Unknown module: $1 (use --help)" ;;
    esac
}

main "$@"

