#!/usr/bin/env bash
# Resolve CPU-worker settings; no installs, numerical work, or path migration.
LD_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
export LD_CODE_DIR="$LD_ROOT/f"
if [[ ${LD_WORKER:-0} != 1 && ${LD_CONFIG_LOADED:-0} != 1 && -f ${LD_CONFIG:-$LD_ROOT/ld.config.sh} ]]; then
    source "${LD_CONFIG:-$LD_ROOT/ld.config.sh}"
fi
export SUBMIT=${SUBMIT:-auto}
if [[ $SUBMIT == auto ]]; then
    if command -v bsub >/dev/null 2>&1 || [[ -n ${LSB_JOBID:-} ]]; then SUBMIT=hpc; else SUBMIT=local; fi
fi
[[ $SUBMIT == hpc || $SUBMIT == local ]] || { echo 'SUBMIT must be hpc/local/auto' >&2; return 2; }
export GRCH=${GRCH:-38} GENOME_BUILD=${GRCH:-38}
[[ $GRCH == 37 || $GRCH == 38 ]] || { echo 'GRCH must be 37/38' >&2; return 2; }
export DIR_WORK=${DIR_WORK:-$HOME} DIR_GEN=${DIR_GEN:-${DIR_WORK:-$HOME}/data/gen}
export ANALYSIS_ROOT=${ANALYSIS_ROOT:-$DIR_WORK/analysis}
if [[ $GRCH == 38 ]]; then export OUT=${OUT:-$ANALYSIS_ROOT/ld}
else export OUT=${OUT:-$ANALYSIS_ROOT/ld/37}; fi
export PFILE_DIR=${PFILE_DIR:-$DIR_GEN/1kg/$GRCH/pfile}
export VCF_DIR=${VCF_DIR:-$DIR_GEN/1kg/$GRCH/vcf}
export PANEL=${PANEL:-$DIR_GEN/1kg/$GRCH/samples.txt}
export INPUT_FORMAT=${INPUT_FORMAT:-pgen}
export PYTHON=${LD_PYTHON:-python3} PLINK2=${PLINK2:-plink2} LD_CXX=${LD_CXX:-g++}
export CXX="$LD_ROOT/f/cxx.sh"
export CHRS=${CHRS:-1-22} SUBMIT_WAIT=${SUBMIT_WAIT:-FALSE} LD_REPLACE=${LD_REPLACE:-FALSE}
# Retained for local execution and CLI compatibility; HPC workers are independent.
export LD_JOBS=${LD_JOBS:-1}
export THREADS=${THREADS:-8} MEMORY_MB=${MEMORY_MB:-8192} LD_MEMORY_GIB=${LD_MEMORY_GIB:-64}
export LD_QUEUE=${LD_QUEUE:-33} LD_HOURS=${LD_HOURS:-120}
export LD_SLOTS=${LD_SLOTS:-$((THREADS+1))}
export LD_TMP_ROOT=${LD_TMP_ROOT:-$OUT} LD_TMP_DIR=${LD_TMP_DIR:-$OUT/tmp}
export LD_STORAGE=${LD_STORAGE:-blocked-zstd} MIN_FREE_GB=${MIN_FREE_GB:-100}
export D_GRID_KB=${D_GRID_KB:-250,500,1000,2000,5000,10000} MAX_D_KB=${MAX_D_KB:-80000}
export STABLE_CONFIRMATIONS=${STABLE_CONFIRMATIONS:-2}
export REFERENCE_COUNT_POLICY=${REFERENCE_COUNT_POLICY:-off} REFERENCE_MAX_REL_DEVIATION=${REFERENCE_MAX_REL_DEVIATION:-0.10}
export COINCIDENT_POLICY=${COINCIDENT_POLICY:-group} LD_TARGETS=${LD_TARGETS:-}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES= WORKFLOW_DEVICE=cpu
[[ $SUBMIT_WAIT == TRUE || $SUBMIT_WAIT == FALSE ]] || { echo '--wait must be TRUE/FALSE' >&2; return 2; }
[[ $LD_REPLACE == TRUE || $LD_REPLACE == FALSE ]] || { echo '--replace must be TRUE/FALSE' >&2; return 2; }
for key in LD_JOBS THREADS MEMORY_MB LD_MEMORY_GIB LD_SLOTS LD_HOURS; do
    [[ ${!key} =~ ^[1-9][0-9]*$ ]] || { echo "$key must be a positive integer" >&2; return 2; }
done
((THREADS <= 1024)) || { echo 'threads must be 1..1024' >&2; return 2; }
if [[ $SUBMIT == local ]]; then
    ((LD_JOBS <= 8)) || { echo 'Local jobs-N must be 1..8' >&2; return 2; }
else
    ((LD_JOBS <= 10000)) || { echo 'HPC jobs-N must be 1..10000' >&2; return 2; }
fi
((LD_SLOTS >= THREADS+1)) || { echo 'job-cores must cover PLINK threads plus one compressor' >&2; return 2; }
((MEMORY_MB <= LD_MEMORY_GIB*1024)) || { echo 'PLINK MEMORY_MB exceeds job memory' >&2; return 2; }
[[ $LD_QUEUE =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo 'Invalid CPU queue' >&2; return 2; }
for key in OUT PFILE_DIR VCF_DIR PANEL LD_TMP_ROOT LD_TMP_DIR; do
    printf -v "$key" '%s' "$(realpath -m -- "${!key}")"
done
