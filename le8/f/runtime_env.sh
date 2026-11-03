#!/usr/bin/env bash
# LE8 paths derive from the shared OS table; explicit caller values win.
# This file is lightweight and safe on a login node.
export LE8_ANALYSIS_ROOT=${LE8_ANALYSIS_ROOT:-$ANALYSIS_ROOT/le8}
export UKB_PHE=${UKB_PHE:-$DIR_UKB}
export UKB_COMMON=${UKB_COMMON:-$UKB_PHE/common}
export UKB_BGEN_DIR=${UKB_BGEN_DIR:-$DIR_GEN/ukb/37/imp}
export LE8_UKB_REL_FILE=${LE8_UKB_REL_FILE:-$DIR_GEN/ukb/37/typ/ukb_rel.dat}
export LE8_GWAS_DIR=${LE8_GWAS_DIR:-$DIR_GWAS/main}
export LE8_PQTL_IV_DIR=${LE8_PQTL_IV_DIR:-$DIR_GWAS/prot}
export LE8_MQTL_IV_DIR=${LE8_MQTL_IV_DIR:-$DIR_GWAS/met}
export LE8_REFGEN_ROOT=${LE8_REFGEN_ROOT:-$DIR_GEN/1kg}
export C2_PROT_HERITABILITY_FILE=${C2_PROT_HERITABILITY_FILE:-$LE8_PQTL_IV_DIR/ppp.h2.csv}
export C5_CELL_ATLAS=${C5_CELL_ATLAS:-$DIR_ANNOT/cellage/atlas.csv}
export LE8_GO_DIR=${LE8_GO_DIR:-$DIR_ANNOT/le8/go}
export LE8_GGSEG_DIR=${LE8_GGSEG_DIR:-$DIR_ANNOT/le8/ggseg}
export LE8_SHARED_HELPERS=${LE8_SHARED_HELPERS:-$LE8_ROOT/../0f}
export LE8_RESULTS_R=${LE8_RESULTS_R:-$LE8_ROOT/../0f/results.R}
export PHE_F=${PHE_F:-$LE8_ROOT/../0f/phenotype.sh}
export CIS_INDEX_FUNCTIONS=${CIS_INDEX_FUNCTIONS:-$LE8_ROOT/../gwas/f/format.f.sh}
export CIS_INDEX_DATA_ROOT=${CIS_INDEX_DATA_ROOT:-$DIR_GWAS}
export LE8_MEMORY_LIMIT_GB=${LE8_MEMORY_LIMIT_GB:-$MEMORY_DEFAULT_GIB}
export RUN_GPU_COLOC=${RUN_GPU_COLOC:-FALSE}
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-$LE8_ANALYSIS_ROOT/logs/runs}
if [[ $submit == hpc ]]; then
    # Login profiles may export a TaiYi/SAIGE library; keep this runtime isolated.
    unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE
    export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
    export PATH="$(dirname -- "$LE8_RSCRIPT"):$PATH"
    export PYTHON_BIN=${PYTHON_BIN:-${LE8_PYTHON:-$DIR_WORK/.conda/envs/le8/bin/python3}}
    export R_BIN=${R_BIN:-$LE8_RSCRIPT}
    export C5_CIGMA_PYTHON=${C5_CIGMA_PYTHON:-$DIR_WORK/.conda/envs/le8-cigma/bin/python3}
    export CONDA_SH=${CONDA_SH:-$OS_GU_CONDA/etc/profile.d/conda.sh}
fi
if [[ ${WORKFLOW_WORKER:-0} == 1 ]]; then
    if [[ $submit == hpc && -z ${LSB_JOBID:-} ]]; then
        echo 'LE8 HPC execution requires an LSF allocation; use --plan-only to inspect jobs.' >&2
        return 2
    fi
    mkdir -p "${TMPDIR:-/tmp}"
    export TMP=${TMPDIR:-/tmp} TEMP=${TMPDIR:-/tmp}
    export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/python-cache
    export LE8_RESOURCE_POOL=${TMPDIR:-/tmp}/le8-resources.json
    export LE8_GLOBAL_RESOURCE_POOL=$LE8_RESOURCE_POOL
    export MRLINK2_EIGH_CACHE=${MRLINK2_EIGH_CACHE:-${LE8_CACHE_ROOT:-${TMPDIR:-/tmp}}/mrlink2-eigh}
    export MRLINK2_PREPARE_CACHE=${MRLINK2_PREPARE_CACHE:-${LE8_CACHE_ROOT:-${TMPDIR:-/tmp}}/mrlink2-prepared}
    export MRLINK2_VALIDATED_CACHE=${MRLINK2_VALIDATED_CACHE:-${LE8_CACHE_ROOT:-${TMPDIR:-/tmp}}/mrlink2-validated}
fi
