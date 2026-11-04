#!/usr/bin/env bash
# Portable entry-point setup. Planning must precede genotype/GWAS discovery.

gwas_runtime_init() {
    local entry=${1:?entry required} root token value help=0 legacy_root=
    shift
    root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
    export GWAS_CODE_ROOT=$root
    source "$root/../0f/workflow_env.sh"
    # The historical flag now uses the same resource-aware dispatcher.
    local -a requested=()
    while (($#)); do
        case "$1" in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
        if [[ $1 == --submit-bsub ]]; then
            case ${2:-} in
                TRUE|true) requested+=(--submit hpc) ;;
                FALSE|false) requested+=(--submit local) ;;
                *) echo '--submit-bsub requires TRUE or FALSE' >&2; return 2 ;;
            esac
            shift 2
        else requested+=("$1"); shift; fi
    done
    workflow_init gwas "${requested[@]}" || return
    source "$root/../0f/cli_dirs.sh"
    script_directory_args "$root/$entry" "${WORKFLOW_ARGS[@]}"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    GWAS_ARGS=()
    export GWAS_PLAN_ONLY=${GWAS_PLAN_ONLY:-FALSE}
    while (($#)); do
        token=$1
        case "$token" in
            --plan-only) export GWAS_PLAN_ONLY=TRUE; shift ;;
            --queue|--cores|--memory-gb|--hours|--runtime-env|--plan-dir|--memory-scope|--tool-threads|--tool-memory-mb)
                [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "$1 requires a value" >&2; return 2; }
                value=$2
                case "$token" in
                    --queue) export GWAS_QUEUE=$value ;;
                    --cores) export GWAS_CORES=$value ;;
                    --memory-gb) export GWAS_MEMORY_GIB=$value ;;
                    --hours) export GWAS_HOURS=$value ;;
                    --runtime-env) export GWAS_RUNTIME_ENV=$(realpath -m -- "$value") ;;
                    --plan-dir) export GWAS_PLAN_DIR=$(realpath -m -- "$value") ;;
                    --memory-scope) export GWAS_MEMORY_SCOPE=$value ;;
                    --tool-threads) export GWAS_TOOL_THREADS=$value ;;
                    --tool-memory-mb) export GWAS_PLINK_MEMORY_MB=$value ;;
                esac
                shift 2 ;;
            --data-root)
                if [[ $entry == format.sh ]]; then
                    [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "$1 requires a value" >&2; return 2; }
                    legacy_root=$(realpath -m -- "$2")
                    GWAS_ARGS+=("$1" "$legacy_root"); shift 2
                else GWAS_ARGS+=("$1"); shift; fi ;;
            *) GWAS_ARGS+=("$1"); shift ;;
        esac
    done
    export GWAS_DATA_ROOT=${GWAS_DATA_ROOT:-$DIR_DATA}
    export GWAS_FILES_ROOT=${GWAS_FILES_ROOT:-${legacy_root:+$legacy_root/files}}
    export GWAS_FILES_ROOT=${GWAS_FILES_ROOT:-$OS_FILES_ROOT}
    export GWAS_SOFTWARE_ROOT=${GWAS_SOFTWARE_ROOT:-${legacy_root:+$legacy_root/software}}
    export GWAS_SOFTWARE_ROOT=${GWAS_SOFTWARE_ROOT:-$DIR_SOFTWARE}
    export GWAS_ANALYSIS_ROOT=${GWAS_ANALYSIS_ROOT:-$ANALYSIS_ROOT/gwas}
    export GWAS_PHENOTYPE_FILE=${GWAS_PHENOTYPE_FILE:-${OS_SNP_LIST%/*}/ukb.phe}
    export CHECK_GRCH_SNP_LIST=${CHECK_GRCH_SNP_LIST:-${legacy_root:+$legacy_root/ukb/common/snp.lst}}
    export CHECK_GRCH_SNP_LIST=${CHECK_GRCH_SNP_LIST:-$OS_SNP_LIST}
    export GWAS_TOOL_THREADS=${GWAS_TOOL_THREADS:-2} GWAS_PLINK_MEMORY_MB=${GWAS_PLINK_MEMORY_MB:-8192}
    export PGS_PLINK_MEMORY_MB=${PGS_PLINK_MEMORY_MB:-8192}
    export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-$GWAS_ANALYSIS_ROOT/logs}
    export PATH="$GWAS_SOFTWARE_ROOT/bin:$PATH"
    for token in "${GWAS_ARGS[@]}"; do
        case "$token" in -h|--help|--help-all|help) help=1 ;; esac
    done
    ((${#GWAS_ARGS[@]})) || help=1
    if ((help == 0)) && [[ ${GWAS_WORKER:-0} != 1 ]] &&
        { [[ $submit == hpc ]] || [[ $GWAS_PLAN_ONLY == TRUE ]]; }; then
        # --submit hpc inside a running job uses that allocation, never nests bsub.
        if [[ -z ${LSB_JOBID:-} || $GWAS_PLAN_ONLY == TRUE ]]; then
            exec "${GWAS_PLAN_PYTHON:-python3}" "$root/f/dispatch.py" "$entry" "${GWAS_ARGS[@]}"
        else
            exec "${GWAS_PLAN_PYTHON:-python3}" "$root/f/dispatch.py" --in-allocation "$entry" "${GWAS_ARGS[@]}"
        fi
    fi
    if [[ $submit == hpc && ${GWAS_WORKER:-0} == 1 && -z ${LSB_JOBID:-} ]]; then
        echo 'HPC analysis requires an LSF allocation; submit the generated .cmd.' >&2
        return 2
    fi
    if ((help == 0)) && [[ -n ${GWAS_RUNTIME_ENV:-} && ${GWAS_WORKER:-0} != 1 ]]; then
        exec bash "$root/f/worker.sh" "$root/$entry" "${GWAS_ARGS[@]}"
    fi
}

gwas_runtime_usage() {
    cat <<'HELP'

Local/HPC execution (all three entry points):
    --submit auto|local|hpc       Default hpc; explicit auto detects LSF
    --plan-only                  Write .cmd/plan.json without scanning data or running
    --run-cmd TRUE|FALSE          Execute or only prepare (format/gwas)
    --queue NAME|auto --cores N   LSF queue/CPU allocation [auto / derived from workers]
                                 QiMing auto fits summed RAM and cores on one node
    --memory-gb N --hours N       Total job RAM GiB / walltime hours [derived / 120]
    --tool-threads N --tool-memory-mb N  format PLINK/GCTA threads and PLINK RAM MiB [2/8192]
    --runtime-env FILE            Optional compute-node environment setup script
    --plan-dir DIR --wait TRUE|FALSE  Plan location / wait for successful completion
    --date-scratch YYYY-MM-DD     Shared scratch date; reuse the same date when resuming
    --memory-scope site_host|host|task  LSF memory reservation convention [site_host]
    DIR_WORK, DIR_DATA, DIR_SCRATCH and GWAS_* exports override default paths.
    --jobs is within-job parallelism; --threads is per GWAS process.
    HPC keeps each requested pipeline in one allocation and preserves its internal
    dependencies/checkpoints. Explicitly activate the scientific environment or
    supply --runtime-env; no software is installed automatically.
HELP
}
