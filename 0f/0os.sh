#!/usr/bin/env bash
# One OS-specific path table. Shell names use underscores, not hyphens.
# Override submit, jobs_N or DIR_* before calling a project entry.
if [[ ${submit:-${SUBMIT:-auto}} == auto ]]; then
    if [[ -n ${LSB_JOBID:-} ]] || command -v bsub >/dev/null 2>&1; then submit=hpc; else submit=local; fi
else
    submit=${submit:-$SUBMIT}
fi
case "$submit" in local|hpc) ;; *) echo 'submit must be local or hpc' >&2; return 2 ;; esac
export submit SUBMIT="$submit"
if [[ $submit == local ]]; then jobs_N=${jobs_N:-1}; else jobs_N=${jobs_N:-8}; fi
export jobs_N
[[ $jobs_N =~ ^[1-9][0-9]*$ ]] && (( jobs_N <= 10000 )) || { echo 'jobs_N must be an integer in 1..10000' >&2; return 2; }
export DATE_SCRATCH=${DATE_SCRATCH:-2026-10-31}
export SCRIPTS_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
if [[ $submit == hpc ]]; then
    MEMORY_DEFAULT_GIB=160
    export DIR_WORK=${DIR_WORK:-/work/sph-huangj}
    export DIR_DATA=${DIR_DATA:-/data/sph-huangj}
    export DIR_SCRATCH=${DIR_SCRATCH:-/scratch/$DATE_SCRATCH/sph-huangj}
    export SHARED_PYTHON=${SHARED_PYTHON:-$DIR_WORK/.conda/envs/py311/bin/python3}
    OS_GU_CONDA=$DIR_WORK/software/gu/conda
    OS_AXT_ROOT=$DIR_WORK/data/annot/axt/37
    OS_FILES_ROOT=$DIR_DATA/files
    OS_SNP_LIST=$DIR_DATA/ukb/common/snp.lst
else
    MEMORY_DEFAULT_GIB=32
    export DIR_WORK=${DIR_WORK:-/mnt/d}
    export DIR_DATA=${DIR_DATA:-/mnt/f}
    export DIR_SCRATCH=${DIR_SCRATCH:-$DIR_WORK}
    export SHARED_PYTHON=${SHARED_PYTHON:-python3}
    OS_GU_CONDA=$DIR_WORK/software/anaconda3
    [[ -x $OS_GU_CONDA/bin/conda || ! -x $HOME/anaconda3/bin/conda ]] || OS_GU_CONDA=$HOME/anaconda3
    OS_AXT_ROOT=$DIR_DATA/annot/axt/37
    OS_FILES_ROOT=$DIR_WORK/files
    OS_SNP_LIST=$DIR_WORK/data/ukb/common/snp.lst
fi
# Default per-job budget, NOT detected RAM or a maximum; stages may request more.
export MEMORY_DEFAULT_GIB
export dir_work=$DIR_WORK dir_data=$DIR_DATA dir_scratch=$DIR_SCRATCH
export DIR_SOFTWARE=${DIR_SOFTWARE:-$DIR_WORK/software}
export ANALYSIS_ROOT=$DIR_SCRATCH/analysis METAG_ROOT=$DIR_SCRATCH/metag
export SCRATCH_ANALYSIS_ROOT=$ANALYSIS_ROOT SUBMIT_WAIT=${SUBMIT_WAIT:-FALSE}
export OS_GU_CONDA OS_AXT_ROOT OS_FILES_ROOT OS_SNP_LIST

workflow_args() {
    WORKFLOW_ARGS=()
    while (($#)); do
        case "$1" in
            --submit)
                if [[ ${2:-} == local || ${2:-} == hpc ]]; then
                    if [[ $2 != "$submit" ]]; then
                        unset DIR_WORK DIR_DATA DIR_SCRATCH DIR_SOFTWARE SHARED_PYTHON
                    fi
                    export submit=$2 SUBMIT=$2; shift 2
                else
                    # metag uses bare --submit to dispatch a generated plan.
                    WORKFLOW_ARGS+=("$1"); shift
                fi ;;
            --jobs-N|--jobs_N|--wait|--date-scratch)
                [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "$1 requires a value" >&2; return 2; }
                case "$1" in
                    --jobs-N|--jobs_N) export jobs_N=$2 ;;
                    --wait) export SUBMIT_WAIT=${2^^} ;;
                    --date-scratch)
                        if [[ $submit == hpc && ${DIR_SCRATCH:-} == /scratch/$DATE_SCRATCH/sph-huangj ]]; then unset DIR_SCRATCH; fi
                        export DATE_SCRATCH=$2 ;;
                esac
                shift 2 ;;
            *) WORKFLOW_ARGS+=("$1"); shift ;;
        esac
    done
    [[ $SUBMIT_WAIT == TRUE || $SUBMIT_WAIT == FALSE ]] || { echo '--wait requires TRUE or FALSE' >&2; return 2; }
    [[ $DATE_SCRATCH =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] && [[ $(date -d "$DATE_SCRATCH" +%F 2>/dev/null) == "$DATE_SCRATCH" ]] || { echo 'Invalid scratch date' >&2; return 2; }
    source "$SCRIPTS_ROOT/0f/0os.sh"
}
