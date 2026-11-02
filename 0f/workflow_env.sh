#!/usr/bin/env bash
# Lightweight path/resource initialization. Never imports numerical libraries,
# reads cohort matrices, activates Conda, or executes an analysis on a login node.
# Call workflow_init PROJECT "$@"; consume the returned WORKFLOW_ARGS array.

workflow_init() {
    local project=${1:?project required}
    shift
    local shared_root token next requested_mode= requested_date= i
    shared_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
    local -a original=("$@")
    # Select the OS table BEFORE sourcing it, preserving explicit DIR_* overrides.
    for ((i=0; i<${#original[@]}; i++)); do
        token=${original[$i]}
        case "$token" in
            --submit)
                next=${original[$((i+1))]:-}
                [[ $next == local || $next == hpc || $next == auto ]] || { echo '--submit requires local|hpc|auto' >&2; return 2; }
                requested_mode=$next ;;
            --submit=*) requested_mode=${token#*=} ;;
            --date-scratch) requested_date=${original[$((i+1))]:-} ;;
            --date-scratch=*) requested_date=${token#*=} ;;
        esac
    done
    [[ -z $requested_mode ]] || export submit=$requested_mode SUBMIT=$requested_mode
    if [[ -n $requested_date ]]; then
        # A date change moves an inherited default scratch path; an explicitly
        # customized directory remains the caller's authoritative choice.
        if [[ ${submit:-${SUBMIT:-auto}} == hpc && ${DIR_SCRATCH:-} == /scratch/${DATE_SCRATCH:-2026-10-31}/sph-huangj ]]; then
            unset DIR_SCRATCH
            [[ ${GRID_OUTPUT_ROOT:-} != /scratch/${DATE_SCRATCH:-2026-10-31}/sph-huangj/analysis/grid ]] || unset GRID_OUTPUT_ROOT
            [[ ${GRID_ANALYSIS_ROOT:-} != /scratch/${DATE_SCRATCH:-2026-10-31}/sph-huangj/analysis/grid ]] || unset GRID_ANALYSIS_ROOT
            [[ ${SCRIPT_LOG_DIR:-} != /scratch/${DATE_SCRATCH:-2026-10-31}/sph-huangj/analysis/grid/logs ]] || unset SCRIPT_LOG_DIR
            [[ ${LE8_ANALYSIS_ROOT:-} != /scratch/${DATE_SCRATCH:-2026-10-31}/sph-huangj/analysis/le8 ]] || unset LE8_ANALYSIS_ROOT
        fi
        export DATE_SCRATCH=$requested_date
    fi
    source "$shared_root/0os.sh" || return
    WORKFLOW_ARGS=()
    while (($#)); do
        case "$1" in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
        case "$1" in
            --submit) shift 2 ;;
            --jobs-N|--jobs_N)
                [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; return 2; }
                export jobs_N=$2; shift 2 ;;
            --wait)
                if [[ ${2:-} == TRUE || ${2:-} == FALSE || ${2:-} == true || ${2:-} == false ]]; then
                    export SUBMIT_WAIT=${2^^}; shift 2
                else export SUBMIT_WAIT=TRUE; shift; fi ;;
            --date-scratch) shift 2 ;;
            *) WORKFLOW_ARGS+=("$1"); shift ;;
        esac
    done
    [[ $jobs_N =~ ^[1-9][0-9]*$ ]] && ((jobs_N <= 10000)) || { echo '--jobs-N requires 1..10000' >&2; return 2; }
    [[ $SUBMIT_WAIT == TRUE || $SUBMIT_WAIT == FALSE ]] || { echo '--wait requires TRUE or FALSE' >&2; return 2; }
    [[ $DATE_SCRATCH =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] && [[ $(date -d "$DATE_SCRATCH" +%F 2>/dev/null) == "$DATE_SCRATCH" ]] || { echo 'Invalid scratch date' >&2; return 2; }
    export WORKFLOW_PROJECT=$project
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    if [[ $project == grid ]]; then
        export GRID_OUTPUT_ROOT=${GRID_OUTPUT_ROOT:-$ANALYSIS_ROOT/grid}
        export GRID_ANALYSIS_ROOT=${GRID_ANALYSIS_ROOT:-$GRID_OUTPUT_ROOT}
        export GRID_DATA_ROOT=${GRID_DATA_ROOT:-$DIR_WORK}
        export GRID_GWAS_DIR=${GRID_GWAS_DIR:-$DIR_DATA/gwas/4grid/common}
        if [[ -z ${GRID_TARGET_DIR:-} ]]; then
            GRID_TARGET_DIR=$DIR_DATA/gen/ukb/37/typ
            if [[ ! -s $GRID_TARGET_DIR/chr1.pgen && ! -s $GRID_TARGET_DIR/chr1.bed && -s $DIR_DATA/gen/ukb/37/hap/chr1.pgen ]]; then
                GRID_TARGET_DIR=$DIR_DATA/gen/ukb/37/hap
            fi
            export GRID_TARGET_DIR
        fi
        export GRID_IMP_DIR=${GRID_IMP_DIR:-$DIR_DATA/gen/ukb/37/imp}
        export GRID_PHE_FILE=${GRID_PHE_FILE:-$DIR_WORK/data/ukb/Rdata/phe.rds}
        export GRID_SCORE_DIR=${GRID_SCORE_DIR:-$DIR_DATA/gen/ukb/pgs}
        export GRID_CSX_REF_DIR=${GRID_CSX_REF_DIR:-$DIR_DATA/LD/csx}
        export GRID_REMOVE=${GRID_REMOVE:-$OS_FILES_ROOT/ukb.exclude.id}
        export GRID_REL_FILE=${GRID_REL_FILE:-$DIR_DATA/gen/ukb/37/typ/ukb_rel.dat}
        export GRID_MED_FILE=${GRID_MED_FILE:-$OS_FILES_ROOT/DiscoDivas/med.g1000.4pop.tsv}
        export GRID_PCA_WEIGHT=${GRID_PCA_WEIGHT:-$OS_FILES_ROOT/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight}
        export GRID_PCA_FILE=${GRID_PCA_FILE:-$DIR_DATA/gen/ukb/pca_proj/ukb.discodivas.pca.tsv.gz}
        export GRID_ARG_HAP_DIR=${GRID_ARG_HAP_DIR:-$DIR_DATA/gen/ukb/37/hap}
        export GRID_PRSFORMER_UPSTREAM=${GRID_PRSFORMER_UPSTREAM:-$DIR_SOFTWARE/PRSformer}
        export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-$GRID_OUTPUT_ROOT/logs}
        if [[ $submit == hpc ]]; then
            export GRID_CONDA_ENV=${GRID_CONDA_ENV:-$DIR_WORK/.conda/envs/grid}
            export GRID_MEMORY_CAP_GB=${GRID_MEMORY_CAP_GB:-$MEMORY_DEFAULT_GIB}
        fi
    fi
    # Caller chooses its scientific runtime. Planning itself needs only python3.
    export WORKFLOW_PYTHON=${WORKFLOW_PYTHON:-${SHARED_PYTHON:-python3}}
}
