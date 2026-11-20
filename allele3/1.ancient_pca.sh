#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
--grch 38  Target build [38]; new analysis output: ANALYSIS_ROOT/allele3/38/
./1.ancient_pca.sh --grch 38 # 计算默认 PCA
./1.ancient_pca.sh MODULE --grch 38 # pca|study|render|stability|drivers
./1.ancient_pca.sh --replace FALSE --pca-threads 4 --grch 38

Usage: ./1.ancient_pca.sh [pca|study|render|stability|drivers] [OPTIONS]
No action: compute PCA from the prepared ancient/38/ancient.p* panel.
  study                 PCA sensitivity experiments after input preparation
  --pca-threads N --pca-memory-mb N --out DIR
  --replace FALSE|TRUE   Reuse completed matching work [FALSE]
  --plan-only           Inspect without submitting
  --queue NAME --job-hours N
  --help-full           All engine options
HPC is detected automatically. Job name: allele3.ancient_pca
PUBLIC_ENTRY_HELP
    exit 0
fi
# Standalone fixed-basis PCA and PC1/PC2/age visualization.
set -euo pipefail
# Public descriptive alias; retain the existing engine action for compatibility.
if [[ ${1:-} == pca ]]; then shift; set -- run "$@"; fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=pca3d
# shellcheck source=f/0.run.sh
source "$ROOT/f/0.run.sh"

usage() {
    bash "${BASH_SOURCE[0]}" --help
}

main() {
    local arg
    local help_entry="$ROOT/f/pca.py"
    [[ ${1:-} != space ]] || help_entry="$ROOT/f/allele_space.py"
    case ${1:-} in study|study-plan|stability|drivers) help_entry="$ROOT/f/pca_study.py" ;; esac
    for arg in "$@"; do
        case "$arg" in -h|--help|help) usage; return 0 ;; --help-full) "$(workflow_python "${PCA3D_PYTHON:-}" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3)" "$help_entry" --help; return 0 ;; esac
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
        run|plan|render|atlas|prepare-basis|project-plan|project|trajectory-plan|trajectory|space|study|study-plan|stability|drivers) ;; *) fail "Unknown command: $command" ;;
    esac
    local analysis_root=${PCA3D_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3} out=${PCA3D_OUT:-}
    local ancient=${PCA3D_ANCIENT_DIR:-${DIR_GEN}/ancient} cache=${PCA3D_PGEN_CACHE_DIR:-}
    local python=${PCA3D_PYTHON:-${ALLELE_PYTHON:-}} memory=${PCA3D_MEMORY_LIMIT_GB:-${ALLELE_MEMORY_LIMIT_GB:-16}}
    local threads=${PCA3D_THREADS:-4} pca_memory=${PCA3D_MEMORY_MB:-16384}
    local interval=${PCA3D_SCREEN_LOG_INTERVAL:-300} verbose=0 dry_run=0 full_help=0
    local replace=${ALLELE_REPLACE:-FALSE}
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
            --simulated|--skip-atlas|--atlas-no-static|--allow-unknown-info|--allow-palindromic|--no-plots|--date-interval-uniform|--joint-sites|--allow-unverified-coordinates|--keep-intermediates) extra+=("$1"); shift ;;
            --scenarios|--scenario-file|--study-dir|--runtime-root|--window-bp|--neighbors|--max-neighbor-samples|--components|--exclude-bed)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            --sample-features|--feature-schema|--atlas-theme)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            --observations|--coordinates|--time-dir|--alleles-dir|--genome-build|--evidence-kind|--pc-bins|--min-earliest-carriers|--chunk-size|--plot-alleles|--date-tolerance-years|--date-mismatch-policy|--min-joint-ac)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            --reference-metadata|--target-metadata|--date-draws|--coordinate-receipt|--provenance-manifest|--training-ids|--supplement-metadata|--date-distributions|--bin-years|--max-age-bp|--min-samples|--pcs|--bootstrap|--display-limit|--pgen-prefix|--metadata|--hq-variants|--eigenvec|--html|--source-html|--supplement-xlsx|--qpadm|--ancestry-model|--selected-variants|--target-variants|--pca-basis|--pca-basis-label|--pca-basis-min-age-bp|--pca-basis-max-age-bp|--pca-maf|--pca-geno|--pca-max-per-group|--pca-components|--pca-target-flank-bp|--pca-time-bin-years|--pca-exclude-bed|--pca-max-variants|--pca-layer-years|--pca-time-view|--pca-color-by|--pca-time-scale|--pca-display-max-age-bp|--basis-dir|--reuse-basis-dir|--target-glob|--project-out|--target-info|--min-info|--min-overlap|--min-call-rate|--count-kind|--plink2|--seed|--title)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            *) fail "Unknown option: $1 (see --help)" ;;
        esac
    done
    [[ $replace == TRUE || $replace == FALSE ]] || fail '--replace must be TRUE or FALSE'
    positive_integer --memory-limit-gb "$memory"
    positive_integer --pca-threads "$threads"
    positive_integer --pca-memory-mb "$pca_memory"
    workflow_cpu_check "$threads" 8
    if ((pca_memory > memory * 512)); then pca_memory=$((memory * 512)); fi
    ((pca_memory >= 640)) || fail 'PLINK needs at least 640 MiB; increase --memory-limit-gb / --pca-memory-mb'
    local conda_root=${PCA3D_CONDA_ROOT:-$HOME/anaconda3}
    python=$(workflow_python "$python" "$conda_root/envs/pca3d/bin/python" \
        "$conda_root/envs/allele3-evolution/bin/python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3)
    export WORKFLOW_LOG_DIR=${WORKFLOW_LOG_DIR:-${PCA3D_RUNTIME_ROOT:-$analysis_root/pca3d}/logs}
    case "$command" in
        study|study-plan|stability|drivers)
            local study_runtime="${PCA3D_RUNTIME_ROOT:-$analysis_root/pca3d}/study" i
            for ((i=0; i<${#extra[@]}; i++)); do
                [[ ${extra[i]} != --runtime-root ]] || study_runtime=${extra[i+1]}
            done
            if [[ -z ${supplied[--out]:-${supplied[--dir-out]:-}} ]]; then
                out=${PCA3D_STUDY_OUT:-${PCA3D_OUT:-$ancient/38/pca3d}/study}
            fi
            # Study resumes complete scenarios unless a replacement was explicitly requested.
            [[ -n ${supplied[--replace]:-} ]] || replace=FALSE
            args=("$command" --out "$out" --ancient-dir "$ancient" --replace "$replace"
                --threads "$threads" --memory-mb "$pca_memory")
            [[ -n ${supplied[--runtime-root]:-} ]] || args+=(--runtime-root "$study_runtime")
            if [[ $command == study || $command == study-plan ]]; then
                [[ -n ${supplied[--pgen-prefix]:-} || -z ${PCA3D_PGEN_PREFIX:-} ]] || args+=(--pgen-prefix "$PCA3D_PGEN_PREFIX")
                [[ -n ${supplied[--hq-variants]:-} || -z ${PCA3D_HQ_VARIANTS:-} ]] || args+=(--hq-variants "$PCA3D_HQ_VARIANTS")
                [[ -n ${supplied[--metadata]:-} || -n ${supplied[--supplement-xlsx]:-} || ! -f ${PCA3D_METADATA:-} ]] || args+=(--metadata "$PCA3D_METADATA")
            fi
            [[ -n ${supplied[--plink2]:-} || -z ${PCA3D_PLINK2:-} ]] || args+=(--plink2 "$PCA3D_PLINK2")
            args+=("${extra[@]}")
            export WORKFLOW_LOG_DIR="$study_runtime/logs"
            workflow_launch pca-study "$python" "$ROOT/f/pca_study.py" "$out" "$memory" "$interval" \
                "$verbose" "$dry_run" 1 "${args[@]}"
            return ;;
    esac
    if [[ $command == space ]]; then
        [[ -n ${supplied[--metadata]:-} || -z ${PCA3D_METADATA:-} ]] || extra+=(--metadata "$PCA3D_METADATA")
        [[ -n ${supplied[--coordinates]:-} || -z ${PCA3D_EIGENVEC:-} ]] || extra+=(--coordinates "$PCA3D_EIGENVEC")
        out=${out:-$analysis_root/base_occupancy/space}
        workflow_launch base_occupancy-space "$python" "$ROOT/f/allele_space.py" "$out" "$memory" "$interval" \
            "$verbose" "$dry_run" 1 --out "$out" "${extra[@]}"
        return
    fi
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
        out=${out:-$ancient/38/pca3d/temporal}
        args=(trajectory "$command" --out "$out" "${defaults[@]}" "${extra[@]}")
        ((full_help == 0)) || args=(trajectory --help)
        workflow_launch pca-trajectory "$python" "$ROOT/f/pca_plot.py" "$out" "$memory" "$interval" "$verbose" "$dry_run" "$((1-full_help))" "${args[@]}"
        return
    fi
    if [[ $command == run || $command == plan || $command == prepare-basis ]]; then
        if [[ -z ${supplied[--metadata]:-} && -z ${supplied[--supplement-xlsx]:-} ]]; then
            if [[ -f ${PCA3D_METADATA:-} ]]; then defaults+=(--metadata "$PCA3D_METADATA")
            else defaults+=(--supplement-xlsx "$ancient/files/Nature.xlsx"); fi
        fi
    fi
    export WORKFLOW_LOG_DIR=${WORKFLOW_LOG_DIR:-${PCA3D_RUNTIME_ROOT:-$analysis_root/pca3d}/logs}
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
    out=${out:-$ancient/38/pca3d}
    args=("$command" --out "$out" --replace "$replace" --ancient-dir "$ancient" --pca-threads "$threads"
        --pca-memory-mb "$pca_memory" "${extra[@]}")
    [[ -z $cache ]] || args+=(--dir-pgen-cache "$cache")
    if ((full_help)); then args=(--help); fi
    workflow_launch pca3d "$python" "$ROOT/f/pca.py" "$out" "$memory" \
        "$interval" "$verbose" "$dry_run" "$((1-full_help))" core "${args[@]}"
}

main "$@"

