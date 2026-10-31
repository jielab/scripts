#!/usr/bin/env bash
# Standalone fixed-basis PCA and PC1/PC2/age visualization.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=pca3d
# shellcheck source=f/0.run.sh
source "$ROOT/f/0.run.sh"

usage() {
    cat <<'HELP'
Usage: ./1.pca3D.sh [MODULE] [OPTIONS]

Modules:
  run / plan       Fit/project a new PCA (default), or generate the plan only
  render           Render existing PCA coordinates
  atlas            Build an atlas from existing HTML or coordinates
  prepare-basis    Prepare a projection basis without refitting
  project / project-plan        Run external projection, or inspect its plan
  trajectory / trajectory-plan  Analyze temporal structure, or audit inputs

Options:
  --out DIR              Output [/mnt/f/gen/ancient/allele3/pca3d/]
  --replace TRUE|FALSE   Recalculate or resume [TRUE; run only]
  --pgen-prefix PREFIX   Input PGEN/PVAR/PSAM prefix
  --metadata TSV         Sample metadata
  --reference-metadata TSV  Projection reference metadata (optional override)
  --target-metadata TSV  Optional projection target metadata snapshot
  --eigenvec FILE        Existing PCA scores for rendering
  --basis-dir DIR        Fixed PCA basis for projection
  --coordinate-receipt JSON  Verified score receipt for temporal analysis
  --dry-run              Print the command only
  --help-full            Show detailed options
  -h, --help             Show this help

Examples:
  ./1.pca3D.sh [render|run|plan|atlas|prepare-basis|project|project-plan|trajectory|trajectory-plan]
  ./1.pca3D.sh --out /mnt/f/gen/ancient/allele3/pca3d/
HELP
}

main() {
    local arg
    for arg in "$@"; do
        case "$arg" in -h|--help|help) usage; return 0 ;; --help-full) "$(workflow_python "${PCA3D_PYTHON:-}" /home/huangj/anaconda3/envs/allele3-evolution/bin/python3)" "$ROOT/f/0.pca3D.py" --help; return 0 ;; esac
    done
    if [[ ${PCA3D_CONFIG+x} && ! -f $PCA3D_CONFIG ]]; then fail "Explicit PCA3D_CONFIG does not exist: $PCA3D_CONFIG"; fi
    local config=${PCA3D_CONFIG:-"$ROOT/pca3d.env"}
    if [[ -f $config ]]; then
        set -a
        # shellcheck disable=SC1090
        source "$config"
        set +a
    fi
    local command=run
    if (($#)) && [[ $1 != -* ]]; then command=$1; shift; fi
    case "$command" in
        pca-plan) command=plan ;; pca-time) command=render ;;
        run|plan|render|atlas|prepare-basis|project-plan|project|trajectory-plan|trajectory) ;; *) fail "Unknown command: $command" ;;
    esac
    local analysis_root=${PCA3D_ANALYSIS_ROOT:-/mnt/d/analysis/allele3} out=${PCA3D_OUT:-}
    local ancient=${PCA3D_ANCIENT_DIR:-/mnt/f/gen/ancient} cache=${PCA3D_PGEN_CACHE_DIR:-}
    local python=${PCA3D_PYTHON:-} memory=${PCA3D_MEMORY_LIMIT_GB:-16}
    local threads=${PCA3D_THREADS:-4} pca_memory=${PCA3D_MEMORY_MB:-8192}
    local interval=${PCA3D_SCREEN_LOG_INTERVAL:-300} verbose=0 dry_run=0 full_help=0
    local replace=${PCA3D_REPLACE:-TRUE}
    local -a extra=() args=()
    local -A supplied=()
    for arg in "$@"; do [[ $arg != --* ]] || supplied[${arg%%=*}]=1; done
    if [[ $command == run && -z ${supplied[--out]:-${supplied[--dir-out]:-}} && -z $out ]]; then
        fail 'PCA calculation requires --out DIR; see --help for the default output path'
    fi
    while (($#)); do
        if [[ $1 == --*=* ]]; then
            arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"
        fi
        case "$1" in
            --dry-run) dry_run=1; shift ;;
            --replace) need_value "$@"; replace=${2^^}; shift 2 ;;
            --help-full) full_help=1; shift ;;
            --verbose-screen) verbose=1; shift ;;
            --dir-analysis-root) need_value "$@"; analysis_root=$2; shift 2 ;;
            --dir-out|--out) need_value "$@"; out=$2; shift 2 ;;
            --dir-ancient|--ancient-dir) need_value "$@"; ancient=$2; shift 2 ;;
            --dir-pgen-cache) need_value "$@"; cache=$2; shift 2 ;;
            --no-pgen-cache) cache=''; shift ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --pca-threads|--threads|--cores) need_value "$@"; threads=$2; shift 2 ;;
            --pca-memory-mb) need_value "$@"; pca_memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --simulated|--skip-atlas|--allow-unknown-info|--allow-palindromic|--no-plots|--date-interval-uniform) extra+=("$1"); shift ;;
            --reference-metadata|--target-metadata|--date-draws|--coordinate-receipt|--provenance-manifest|--training-ids|--supplement-metadata|--date-distributions|--bin-years|--max-age-bp|--min-samples|--pcs|--bootstrap|--display-limit|--pgen-prefix|--metadata|--hq-variants|--eigenvec|--html|--source-html|--supplement-xlsx|--qpadm|--ancestry-model|--selected-variants|--target-variants|--pca-basis|--pca-basis-label|--pca-basis-min-age-bp|--pca-basis-max-age-bp|--pca-maf|--pca-geno|--pca-max-per-group|--pca-components|--pca-target-flank-bp|--pca-time-bin-years|--pca-exclude-bed|--pca-max-variants|--pca-layer-years|--pca-time-view|--pca-color-by|--pca-time-scale|--pca-display-max-age-bp|--basis-dir|--reuse-basis-dir|--target-glob|--project-out|--target-info|--min-info|--min-overlap|--min-call-rate|--count-kind|--plink2|--seed|--title)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            *) fail "Unknown option: $1 (see --help)" ;;
        esac
    done
    [[ $replace == TRUE || $replace == FALSE ]] || fail '--replace must be TRUE or FALSE'
    positive_integer --memory-limit-gb "$memory"
    positive_integer --pca-threads "$threads"
    positive_integer --pca-memory-mb "$pca_memory"
    ((threads <= 8)) || fail 'Use at most 8 threads on this shared host'
    if ((pca_memory > memory * 512)); then pca_memory=$((memory * 512)); fi
    ((pca_memory >= 640)) || fail 'PLINK needs at least 640 MiB; increase --memory-limit-gb / --pca-memory-mb'
    local conda_root=${PCA3D_CONDA_ROOT:-$HOME/anaconda3}
    python=$(workflow_python "$python" "$conda_root/envs/pca3d/bin/python" \
        "$conda_root/envs/allele3-evolution/bin/python" /home/huangj/anaconda3/envs/allele3-evolution/bin/python3)
    # Add defaults only after recording explicit input modes. CLI values always win.
    local temporal=0 key value
    [[ $command != trajectory && $command != trajectory-plan ]] || temporal=1
    local -a defaults=()
    if [[ -z ${supplied[--source-html]:-} ]]; then
        if [[ $command == render || $command == atlas || $temporal == 1 ]]; then
            if [[ -z ${supplied[--metadata]:-} && -n ${PCA3D_METADATA:-} ]]; then defaults+=(--metadata "$PCA3D_METADATA"); fi
            if [[ -z ${supplied[--eigenvec]:-} && -n ${PCA3D_EIGENVEC:-} ]]; then defaults+=(--eigenvec "$PCA3D_EIGENVEC"); fi
        fi
    fi
    if [[ $command == render && -z ${supplied[--eigenvec]:-} && -z ${supplied[--reuse-basis-dir]:-} && -n ${PCA3D_REUSE_BASIS:-} ]]; then
        defaults+=(--reuse-basis-dir "$PCA3D_REUSE_BASIS")
    fi
    if ((temporal)); then
        [[ -n ${supplied[--basis-dir]:-} || -n ${supplied[--eigenvec]:-} || -z ${PCA3D_REUSE_BASIS:-} ]] || defaults+=(--basis-dir "$PCA3D_REUSE_BASIS")
        out=${out:-$analysis_root/pca3d/temporal}
        args=("$command" --out "$out" "${defaults[@]}" "${extra[@]}")
        ((full_help == 0)) || args=(--help)
        workflow_launch pca-trajectory "$python" "$ROOT/f/pca_trajectory.py" "$out" "$memory" "$interval" "$verbose" "$dry_run" "$((1-full_help))" "${args[@]}"
        return
    fi
    for key in --pgen-prefix --hq-variants --target-variants --plink2; do
        case "$key" in
            --pgen-prefix) value=${PCA3D_PGEN_PREFIX:-} ;;
            --hq-variants) value=${PCA3D_HQ_VARIANTS:-} ;;
            --target-variants) value=${PCA3D_TARGET_VARIANTS:-} ;;
            --plink2) value=${PCA3D_PLINK2:-} ;;
        esac
        [[ -n ${supplied[$key]:-} || -z $value ]] || defaults+=("$key" "$value")
    done
    if [[ $command == project-plan || $command == project ]]; then
        [[ -n ${supplied[--basis-dir]:-} || -z ${PCA3D_PROJECT_BASIS:-} ]] || defaults+=(--basis-dir "$PCA3D_PROJECT_BASIS")
        [[ -n ${supplied[--target-info]:-} || -z ${PCA3D_PROJECT_INFO:-} ]] || defaults+=(--target-info "$PCA3D_PROJECT_INFO")
        [[ -n ${supplied[--target-glob]:-} || -z ${PCA3D_PROJECT_TARGET_GLOB:-} ]] || defaults+=(--target-glob "$PCA3D_PROJECT_TARGET_GLOB")
    fi
    extra=("${defaults[@]}" "${extra[@]}")
    out=${out:-$ancient/allele3/pca3d}
    args=("$command" --out "$out" --replace "$replace" --ancient-dir "$ancient" --pca-threads "$threads"
        --pca-memory-mb "$pca_memory" "${extra[@]}")
    [[ -z $cache ]] || args+=(--dir-pgen-cache "$cache")
    if ((full_help)); then args=(--help); fi
    workflow_launch pca3d "$python" "$ROOT/f/0.pca3D.py" "$out" "$memory" \
        "$interval" "$verbose" "$dry_run" "$((1-full_help))" "${args[@]}"
}

main "$@"
