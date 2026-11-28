#!/usr/bin/env bash
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case ${1:-} in
        -h|--help|help)
            cat <<'ALLELE_HELP'
source ./0.run.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

内部共享函数库，不是独立分析命令。
由顶层入口 source，定义运行、日志和资源限制。
运行分析：bash /work/sph-huangj/scripts/allele3/allele3.sh -h
ALLELE_HELP
            exit 0 ;;
    esac
fi
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    echo 'This file is a sourced library. Use a root workflow .sh; see -h.' >&2
    exit 2
fi
# Shared process-tree limits and logging for the three root entry points.

fail() { printf '%s: %s\n' "${WORKFLOW_NAME:-allele3}" "$*" >&2; exit 2; }
source "$ROOT/f/0.execution.sh"
allele_execution_init "$@"
set -- "${ALLELE_EXECUTION_ARGS[@]}"
# Help remains available even when an explicitly configured path is missing.
_allele_help=0
for _allele_arg in "$@"; do
    case "$_allele_arg" in -h|--help|help|--help-full|--version) _allele_help=1 ;; esac
done
if ((_allele_help == 0)); then
    if [[ ${ALLELE_CONFIG+x} && ! -f $ALLELE_CONFIG ]]; then fail "Explicit ALLELE_CONFIG does not exist: $ALLELE_CONFIG"; fi
    if [[ -n ${ROOT:-} && -f ${ALLELE_CONFIG:-"$ROOT/allele3.env"} ]]; then
        set -a
        source "${ALLELE_CONFIG:-"$ROOT/allele3.env"}"
        set +a
    fi
fi
unset _allele_arg _allele_help

need_value() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || fail "$1 requires a value"; }
positive_integer() { [[ $2 =~ ^[1-9][0-9]*$ && ${#2} -le 6 ]] || fail "$1 must be a positive integer"; }

verify_limits() {
    local expected=${1:?memory bytes required}
    shift
    local cgroup
    cgroup=$(awk -F: '$1 == "0" { print $3 }' /proc/self/cgroup)
    local base=/sys/fs/cgroup$cgroup
    if [[ -z $cgroup || ! -r $base/memory.max || $(<"$base/memory.max") != "$expected" ||
          $(<"$base/memory.swap.max") != 0 || $(<"$base/memory.oom.group") != 1 ]]; then
        printf '%s: cgroup RAM/swap/OOM-group verification failed; workload was not started\n' "${WORKFLOW_NAME:-allele3}" >&2
        exit 2
    fi
    printf '[%s] verified cgroup=%s memory.max=%s swap.max=0 oom.group=1\n' "${WORKFLOW_NAME:-allele3}" "$cgroup" "$expected" >&2
    exec "$@"
}

launch_workflow() {
    local command=$1 python=$2 entry=$3 out=$4 memory=$5 screen_interval=$6
    local verbose_screen=$7 dry_run=$8 logged=$9 native_help=${10}
    local log_out=${out%/} log_default
    local layout_root=${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT:-$HOME/analysis}/allele3}
    if [[ $log_out == "$layout_root/"* ]]; then
        local relative=${log_out#"$layout_root/"} module
        module=${relative%%/*}
        log_default="$layout_root/$module/logs"
        [[ $relative == "$module" ]] || log_default+=/${relative#*/}
    else
        log_default="$log_out/logs"
    fi
    local log_dir=${WORKFLOW_LOG_DIR:-$log_default}
    local product_log
    product_log=$("$python" "$ROOT/f/output_layout.py" "$out" --kind logs --permanent-only)
    [[ -z $product_log ]] || log_dir=$product_log
    shift 10
    local -a args=("$@")
    local _path_index _path_option
    for ((_path_index=0; _path_index<${#args[@]}-1; _path_index++)); do
        _path_option=${args[_path_index]}
        case "$_path_option" in
            --metadata|--coordinates|--eigenvec|--basis-dir|--reuse-basis-dir|--target-glob|--target-info|--reuse-from|--fasta|--vcf-dir|--dir-vcf|--pfile-dir|--library|--ancient-dir|--dir-ancient|--source|--store|--scan-dir|--dir-modern)
                printf '[resolved path] %s = %s\n' "$_path_option" "${args[_path_index+1]}" ;;
        esac
    done
    if ((dry_run)); then
        printf 'Working-set budget=%s GiB; no extra OS memory cap; output=%s\n' "$memory" "$out"
        if ((logged)); then
            printf 'Screen progress interval=%ss; verbose=%s; full logs=%s/\n' "$screen_interval" "$verbose_screen" "$log_dir"
        fi
        printf '%q ' "$python" "$entry" "${args[@]}"; printf '\n'
        return 0
    fi
    command -v "$python" >/dev/null 2>&1 || fail "Python not found: $python. Install env.yml or supply --python. Help: ./allele3.sh -h"
    # Use matching native tools (bcftools/PLINK) from the selected environment.
    python=$(command -v "$python")
    PATH="$(dirname -- "$python"):$PATH"
    export PATH
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export MPLCONFIGDIR=${MPLCONFIGDIR:-"$HOME/.cache/allele3-evolution/matplotlib"}
    if ((native_help)); then exec "$python" "$entry" "${args[@]}"; fi
    "$python" "$ROOT/f/output_layout.py" "$out"
    local -a launch=("$python" "$ROOT/f/workflow_data_check.py" --workflow "$command" --module-out "$out" --
        "$python" "$entry" "${args[@]}")
    # These modules already own their shallow output layout and perform the
    # source-specific checks themselves. Do not relocate their reusable shards
    # or 52-sample census into the historical publication/tmp layout.
    case "$command" in
        pca-study|herv-cohort|herv-highlights|summary-scan|avi-scan|scan-scan|scan-verify|scan-scan-verify|enrich_test|allele_test|avi_window_enrichment|herv_region_enrichment|archaic_tract_enrichment|alt_impact_similarity)
            launch=("$python" "$entry" "${args[@]}") ;;
    esac
    if ((logged)); then
        local -a console=("$python" "${WORKFLOW_CONSOLE:-$ROOT/f/console.py}" --log-dir "$log_dir"
            --label "$command" --interval "$screen_interval")
        ((verbose_screen == 0)) || console+=(--verbose)
        launch=("${console[@]}" -- "${launch[@]}")
    fi
    if [[ ${ALLELE_BATCH_WORKER:-0} == 1 && ${submit:-local} == hpc ]]; then
        workflow_memory_check "$memory"
        printf '[allocation] LSF job=%s cores=%s; LSF controls memory and runtime\n' "$LSB_JOBID" "$ALLELE_ALLOCATED_CORES" >&2
        exec "${launch[@]}"
    fi
    exec "${launch[@]}"

}

workflow_python() {
    local requested=$1 candidate
    shift
    if [[ -n $requested ]]; then printf '%s\n' "$requested"; return; fi
    if [[ ${submit:-local} == hpc ]]; then
        printf '%s\n' "${DIR_WORK}/.conda/envs/allele3-evolution/bin/python"
        return
    fi
    if [[ -n ${CONDA_PREFIX:-} && ${CONDA_DEFAULT_ENV:-} != base && -x $CONDA_PREFIX/bin/python ]]; then
        printf '%s\n' "$CONDA_PREFIX/bin/python"; return
    fi
    for candidate in "$@"; do
        if [[ -x $candidate ]]; then printf '%s\n' "$candidate"; return; fi
    done
    if [[ -n ${CONDA_PREFIX:-} && -x $CONDA_PREFIX/bin/python ]]; then
        printf '%s\n' "$CONDA_PREFIX/bin/python"; return
    fi
    printf 'python3\n'
}


workflow_launch() {
    local label=$1 python=$2 entry=$3 out=$4 memory=$5 interval=$6
    local verbose=$7 dry_run=$8 logged=$9
    shift 9
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    export WORKFLOW_NAME=$label
    launch_workflow "$label" "$python" "$entry" "$out" "$memory" "$interval" \
        "$verbose" "$dry_run" "$logged" "$((1-logged))" "$@"
}

allele_execution_dispatch "$@"

