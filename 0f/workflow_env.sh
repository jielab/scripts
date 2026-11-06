#!/usr/bin/env bash
# Lightweight path/resource initialization. Never imports numerical libraries,
# reads cohort matrices, activates Conda, or executes an analysis on a login node.
# Call workflow_init PROJECT "$@"; consume the returned WORKFLOW_ARGS array.

workflow_init() {
    local project=${1:?project required}
    shift
    local shared_root
    shared_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
    source "$shared_root/0os.sh" || return
    local -a args=()
    # This wrapper historically accepts a bare --wait as TRUE.
    while (($#)); do
        if [[ $1 == --wait && ${2:-} != TRUE && ${2:-} != FALSE && ${2:-} != true && ${2:-} != false ]]; then
            args+=(--wait TRUE); shift
        else args+=("$1"); shift; fi
    done
    workflow_args "${args[@]}" || return
    export WORKFLOW_PROJECT=$project
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    if [[ $project == grid ]]; then
        export GRID_OUTPUT_ROOT=${GRID_OUTPUT_ROOT:-$ANALYSIS_ROOT/grid}
        export GRID_ANALYSIS_ROOT=${GRID_ANALYSIS_ROOT:-$GRID_OUTPUT_ROOT}
        export GRID_DATA_ROOT=${GRID_DATA_ROOT:-$DIR_WORK}
        export GRID_GWAS_DIR=${GRID_GWAS_DIR:-$DIR_GWAS/4grid/common}
        if [[ -z ${GRID_TARGET_DIR:-} ]]; then
            GRID_TARGET_DIR=$DIR_GEN/ukb/37/typ
            if [[ ! -s $GRID_TARGET_DIR/chr1.pgen && ! -s $GRID_TARGET_DIR/chr1.bed && -s $DIR_GEN/ukb/37/hap/chr1.pgen ]]; then
                GRID_TARGET_DIR=$DIR_GEN/ukb/37/hap
            fi
            export GRID_TARGET_DIR
        fi
        export GRID_IMP_DIR=${GRID_IMP_DIR:-$DIR_GEN/ukb/37/imp}
        export GRID_PHE_FILE=${GRID_PHE_FILE:-$DIR_UKB/Rdata/all.rds}
        export GRID_SCORE_DIR=${GRID_SCORE_DIR:-$DIR_GEN/ukb/pgs}
        export GRID_CSX_REF_DIR=${GRID_CSX_REF_DIR:-$DIR_LD/csx}
        # Sample exclusion is opt-in; no default withdrawal/ID file is required.
        export GRID_REMOVE=${GRID_REMOVE:-}
        export GRID_REL_FILE=${GRID_REL_FILE:-$DIR_GEN/ukb/37/typ/ukb_rel.dat}
        export GRID_MED_FILE=${GRID_MED_FILE:-$OS_FILES_ROOT/DiscoDivas/med.g1000.4pop.tsv}
        export GRID_PCA_WEIGHT=${GRID_PCA_WEIGHT:-$OS_FILES_ROOT/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight}
        export GRID_PCA_FILE=${GRID_PCA_FILE:-$DIR_GEN/ukb/pca_proj/ukb.discodivas.pca.tsv.gz}
        export GRID_ARG_HAP_DIR=${GRID_ARG_HAP_DIR:-$DIR_GEN/ukb/37/hap}
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
