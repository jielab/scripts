#!/usr/bin/env bash
# HERV analyses; the permanent completed scan is a read-only input.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=allele3
# shellcheck source=f/0.run.sh
source "$ROOT/f/0.run.sh"

usage() {
    cat <<'HELP'
Usage: ./allele3.sh [MODULE] [OPTIONS]

Modules:
  pipeline         Run all analysis stages (default)
  primary          Compute primary enrichment from current scan counters
  archaic          Prepare GU/IBDmix cohorts and calculate tract enrichment
  summary / spatial / windows   Enrichment and spatial analyses
  alleles / time / selection     Allele catalog and ancient analyses
  candidates / genotypes / ld    Candidates, genotype validation and LD
  batch MODULE       allofus / atlas / mfi / herv / pipeline (LSF or local DAG)
  research MODULE    data / allofus / allofus-cache / mfi / mapped-targets / panel / atlas / atlas-allofus / balance / space / domains / cohort / archaic-alleles
  report           Assemble the results

Options:
  --out DIR              Output [/mnt/d/analysis/allele3]
  --scan-dir DIR         Scan input [/mnt/f/annot/gnomad/allele3]
  --archaic-root DIR     GU raw cohorts [/mnt/f/analysis/gu]
  --archaic-datasets LIST  Names/globs [1kg]; auto selects all complete cohorts
  --archaic-references S  Altai,Chagyr,Denisova,Denisova25,Vindija by default
  --archaic-chromosomes S  auto includes available X/Y separately; or autosomes,X,Y / X,Y
  --archaic-cache DIR    Per-cohort cache [SCAN-DIR-archaic]
  --stages LIST          Comma-separated stages; only requested stages run
  --metadata TSV         Sample metadata [allele3.env: /mnt/f/gen/ancient/allele3/pca3d/metadata/sample_metadata.tsv]
  --coordinates FILE    Shared PCA coordinates [allele3.env: /mnt/f/gen/ancient/allele3/pca3d/pca/all.sscore]
  --preflight, --dry-run Check paths and dependencies without running analysis
  --screen-log-interval S  Progress interval in seconds [300]
  --verbose-screen       Show full output on screen (always saved in logs)
  --help-full            Show detailed options
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
  ./allele3.sh [pipeline|primary|archaic|summary|spatial|windows|alleles|time|selection|candidates|genotypes|ld|report]
  ./allele3.sh --out /mnt/d/analysis/allele3
  ./allele3.sh --archaic-datasets 1kg
  ./allele3.sh --archaic-datasets '1kg,ukb00[1-9],ukb01[0-3]'

  ./allele3.sh batch atlas --submit hpc --date-scratch 2026-11-02 \
    --targets /scratch/2026-11-02/sph-huangj/analysis/allele3/final/herv_focus_20261028/candidates/targets.tsv \
    --out /scratch/2026-11-02/sph-huangj/analysis/allele3/research/atlas \
    --job-dir /scratch/2026-11-02/sph-huangj/analysis/allele3/jobs/atlas

  ./allele3.sh batch atlas --submit local --jobs-N 1 \
    --targets /mnt/d/analysis/allele3/final/herv_focus_20261028/candidates/targets.tsv \
    --out /mnt/d/analysis/allele3/research/atlas \
    --job-dir /mnt/d/analysis/allele3/jobs/atlas

  Atlas requires extracted alphagenome_variant_impact_score_snvs.tsv.gz
  and a matching .tbi or .csi under DIR_DATA/annot/allele3/alphagenome/hg38/.
  Missing required inputs report full paths and exit nonzero before submission.
  Add --write-cmd to generate scripts only; generate HPC scripts on the HPC.
HELP
}

research_main() {
    local module=${1:-help} arg python=${ALLELE_PYTHON:-}
    (($# == 0)) || shift
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local out=${ALLELE_ANALYSIS_ROOT:-/mnt/d/analysis/allele3}/research/$module
    if [[ $module == cohort && ${1:-} =~ ^(plan|summarize)$ ]]; then out=$out/$1; fi
    local dry_run=0 native=0 verbose=0 out_explicit=0
    local -a args=("$module")
    case "$module" in
        help|-h|--help) native=1 ;;
        data|allofus|allofus-cache|mfi|mapped-targets|panel|atlas|atlas-allofus|balance|space|domains|cohort|archaic-alleles) ;;
        *) fail "Unknown research module: $module (use research --help)" ;;
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
            -h|--help) native=1; args+=("$1"); shift ;;
            *) args+=("$1"); shift ;;
        esac
    done
    if (( ! native && ! out_explicit )); then args+=(--out "$out"); fi
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" /home/huangj/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    WORKFLOW_CONSOLE="$ROOT/runtime/pipeline_console.py"
    # Derived modules atomically replace their result directory. Keep a running
    # log outside that directory so its published path remains valid.
    WORKFLOW_LOG_DIR="${out%/}.logs"
    launch_workflow "research-$module" "$python" "$ROOT/f/research.py" "$out" "$memory" "$interval" \
        "$verbose" "$dry_run" "$((1-native))" "$native" "${args[@]}"
}

main() {
    if [[ ${1:-} == research ]]; then shift; research_main "$@"; return; fi
    case "${1:-}" in
        -h|--help|help) usage; return 0 ;;
        --version) printf '6.0.1\n'; return 0 ;;
        pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report|"") ;;
        -*) ;;
        *) fail "Unknown analysis stage: $1. PCA, individual HERV and scan have separate 1.pca3D.sh / 2.herv.sh / 3.scan.sh entry points." ;;
    esac
    local command=pipeline arg python=${ALLELE_PYTHON:-}
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local analysis_root=${ALLELE_ANALYSIS_ROOT:-/mnt/d/analysis/allele3}
    local out=${ALLELE_HERV_FOCUS_OUT:-} verbose=0 native=0
    local -a args=()
    if (($#)) && [[ $1 != -* ]]; then command=$1; shift; fi
    args+=("$command")
    [[ -z ${ALLELE_METADATA:-} ]] || args+=(--metadata "$ALLELE_METADATA")
    [[ -z ${ALLELE_COORDINATES:-} ]] || args+=(--coordinates "$ALLELE_COORDINATES")
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            -h|--help) usage; return 0 ;;
            --help-full) args=(--help); native=1; shift ;;
            --version) printf '6.0.1\n'; return 0 ;;
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
    python=$(workflow_python "$python" /home/huangj/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    out=$("$python" -c 'import sys; sys.path.insert(0, sys.argv[1]); from herv_pipeline import resolve_output_dir; print(resolve_output_dir(sys.argv[2], sys.argv[3]))' \
        "$ROOT/f" "$analysis_root" "$out")
    args+=(--dir-out "$out")
    WORKFLOW_CONSOLE="$ROOT/runtime/pipeline_console.py"
    launch_workflow "$command" "$python" "$ROOT/f/herv_pipeline.py" "$out" "$memory" "$interval" \
        "$verbose" 0 "$((1-native))" "$native" "${args[@]}"
}

main "$@"

