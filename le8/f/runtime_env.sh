#!/usr/bin/env bash
# LE8 paths derive from the shared OS table; explicit caller values win.
# This file is lightweight and safe on a login node.
export LE8_ANALYSIS_ROOT=${LE8_ANALYSIS_ROOT:-$ANALYSIS_ROOT/le8}
export UKB_PHE=${UKB_PHE:-$DIR_WORK/data/ukb}
export UKB_COMMON=${UKB_COMMON:-$UKB_PHE/common}
export UKB_BGEN_DIR=${UKB_BGEN_DIR:-$DIR_DATA/gen/ukb/37/imp}
export LE8_UKB_REL_FILE=${LE8_UKB_REL_FILE:-$DIR_DATA/gen/ukb/37/typ/ukb_rel.dat}
export LE8_GWAS_DIR=${LE8_GWAS_DIR:-$DIR_DATA/gwas/main}
export LE8_PQTL_IV_DIR=${LE8_PQTL_IV_DIR:-$DIR_DATA/gwas/prot}
export LE8_MQTL_IV_DIR=${LE8_MQTL_IV_DIR:-$DIR_DATA/gwas/met}
export LE8_REFGEN_ROOT=${LE8_REFGEN_ROOT:-$DIR_DATA/gen/1kg}
export C2_PROT_HERITABILITY_FILE=${C2_PROT_HERITABILITY_FILE:-$LE8_PQTL_IV_DIR/ppp.h2.csv}
export C5_CELL_ATLAS=${C5_CELL_ATLAS:-$DIR_DATA/annot/cellage/atlas.csv}
export LE8_GO_DIR=${LE8_GO_DIR:-$DIR_DATA/annot/le8/go}
export LE8_GGSEG_DIR=${LE8_GGSEG_DIR:-$DIR_DATA/annot/le8/ggseg}
export LE8_SHARED_HELPERS=${LE8_SHARED_HELPERS:-$LE8_ROOT/../0f}
export LE8_RESULTS_R=${LE8_RESULTS_R:-$LE8_ROOT/../0f/results.R}
export PHE_F=${PHE_F:-$LE8_ROOT/../0f/phenotype.sh}
export CIS_INDEX_FUNCTIONS=${CIS_INDEX_FUNCTIONS:-$LE8_ROOT/../gwas/f/format.f.sh}
export CIS_INDEX_DATA_ROOT=${CIS_INDEX_DATA_ROOT:-$DIR_DATA/gwas}
export LE8_MEMORY_LIMIT_GB=${LE8_MEMORY_LIMIT_GB:-$MEMORY_DEFAULT_GIB}
export RUN_GPU_COLOC=${RUN_GPU_COLOC:-FALSE}
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-$LE8_ANALYSIS_ROOT/logs/runs}
if [[ $submit == hpc ]]; then
    export PYTHON_BIN=${PYTHON_BIN:-${LE8_PYTHON:-$SHARED_PYTHON}}
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
