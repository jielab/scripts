#!/usr/bin/env bash
# One OS-specific path table. Shell names use underscores, not hyphens.
# Main projects default to HPC; local runs must opt in. 0LOCAL forces local.
# Explicit auto remains available for callers that want LSF detection.
# Override submit, jobs_N or DIR_* before calling a project entry.
submit=${submit:-${SUBMIT:-hpc}}
if [[ $submit == auto ]]; then
    if [[ -n ${LSB_JOBID:-} ]] || command -v bsub >/dev/null 2>&1; then submit=hpc; else submit=local; fi
fi
case "$submit" in local|hpc) ;; *) echo 'submit must be local or hpc' >&2; return 2 ;; esac
export submit SUBMIT="$submit"
if [[ $submit == local ]]; then jobs_N=${jobs_N:-1}; else jobs_N=${jobs_N:-8}; fi
export jobs_N
[[ $jobs_N =~ ^[1-9][0-9]*$ ]] && (( jobs_N <= 10000 )) || { echo 'jobs_N must be an integer in 1..10000' >&2; return 2; }
# Path defaults are tracked so re-sourcing after a mode/root change recomputes
# derived paths, while explicit custom overrides survive. No directories created.
_os_default() {
    local name=$1 value=$2 marker=_OS_DEFAULT_$1
    if [[ ! -v $name || ${!name} == "${!marker-}" ]]; then
        printf -v "$name" '%s' "$value"
    fi
    export "$name"
    printf -v "$marker" '%s' "$value"
    export "$marker"
}
export DATE_SCRATCH=${DATE_SCRATCH:-2026-11-03}
[[ $DATE_SCRATCH =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] && [[ $(date -d "$DATE_SCRATCH" +%F 2>/dev/null) == "$DATE_SCRATCH" ]] || { echo 'Invalid scratch date' >&2; return 2; }
export SCRIPTS_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
# Two roots: DIR_DAT holds small files/code/software/analysis; DIR_BIG holds data.
# DIR_DATA and DIR_SCRATCH are compatibility aliases for DIR_BIG and DIR_TMP.
# TaiYi /data -> QiMing /work/sph-huangj/data (never /data on QiMing).
if [[ $submit == hpc ]]; then
    _os_default MEMORY_DEFAULT_GIB 160
    _os_default DIR_WORK /work/sph-huangj
    _os_default DIR_DAT "$DIR_WORK"
    _os_default DIR_BIG "$DIR_WORK/data"
    _os_default DIR_TMP "$DIR_DAT"
    _os_default OS_GU_CONDA "$DIR_WORK/software/gu/conda"
    _os_default SHARED_PYTHON "$DIR_WORK/.conda/envs/py311/bin/python3"
else
    _os_default MEMORY_DEFAULT_GIB 32
    _os_default DIR_WORK /mnt/d
    _os_default DIR_DAT /mnt/d
    _os_default DIR_BIG /mnt/f
    _os_default DIR_TMP "$DIR_DAT"
    local_conda=$DIR_WORK/software/anaconda3
    [[ -x $local_conda/bin/conda || ! -x $HOME/anaconda3/bin/conda ]] || local_conda=$HOME/anaconda3
    _os_default OS_GU_CONDA "$local_conda"
    unset local_conda
    _os_default SHARED_PYTHON python3
fi
_os_default DIR_DATA "$DIR_BIG"
_os_default DIR_SCRATCH "$DIR_TMP"
# DIR_TMP/DIR_SCRATCH remain optional legacy aliases; no dated scratch is needed.
# Most persistent inputs follow the large-data root; local UKB phenotypes use DIR_DAT/phe.
_os_default DIR_ANNOT "$DIR_DATA/annot"
_os_default DIR_GEN "$DIR_DATA/gen"
_os_default DIR_GWAS "$DIR_DATA/gwas"
_os_default DIR_LD "$DIR_DATA/LD"
_os_default DIR_METAG "$DIR_DATA/metag"
_os_default DIR_SOFTWARE "$DIR_DAT/software"
_os_default DIR_FILES "$DIR_DAT/files"
if [[ $submit == local ]]; then
    _os_default DIR_UKB "$DIR_DAT/phe/ukb"
else
    _os_default DIR_UKB "$DIR_WORK/phe/ukb"
fi
# Historical local phenotype root. Keep migration policy in the path table.
os_resolve_ukb_path() {
    local requested=${1%/}
    if [[ $submit == local && $requested == /mnt/f/ukb && ! -d $requested && -d $DIR_UKB ]]; then
        printf 'UKB path migrated: %s -> %s (DIR_UKB from 0f/0os.sh)\n' "$1" "$DIR_UKB" >&2
        printf '%s\n' "$DIR_UKB"
    else
        printf '%s\n' "$1"
    fi
}
# Additional resource roots shared by Shell, Python, and R.
_os_default DIR_HOME "$HOME"
_os_default DIR_RUNTIME_TMP /tmp
_os_default DIR_WINDOWS_DRIVE /mnt/c
_os_default DIR_WINDOWS "$DIR_WINDOWS_DRIVE/Windows"
_os_default DIR_WINDOWS_PROGRAMDATA "$DIR_WINDOWS_DRIVE/ProgramData"
_os_default DIR_WINDOWS_PROGRAMFILES "$DIR_WINDOWS_DRIVE/Program Files"
_os_default DIR_WINDOWS_PROGRAMFILES_X86 "$DIR_WINDOWS_DRIVE/Program Files (x86)"
_os_default DIR_WINDOWS_USERS "$DIR_WINDOWS_DRIVE/Users"
_os_default OS_RSCRIPT_WINDOWS "$DIR_SOFTWARE/win/R-4.5.1/bin/Rscript.exe"
_os_default OS_RTOOLS_WINDOWS "$DIR_SOFTWARE/win/rtools45"
# Adapters query this list; root names are defined only in this table.
OS_PATH_KEYS=(SCRIPTS_ROOT DIR_WORK DIR_DAT DIR_BIG DIR_DATA DIR_TMP DIR_SCRATCH
    DIR_SOFTWARE DIR_FILES DIR_UKB DIR_ANNOT DIR_GEN DIR_GWAS DIR_LD DIR_METAG
    DIR_RUNTIME_TMP DIR_HOME DIR_WINDOWS_DRIVE DIR_WINDOWS DIR_WINDOWS_PROGRAMDATA
    DIR_WINDOWS_PROGRAMFILES DIR_WINDOWS_PROGRAMFILES_X86 DIR_WINDOWS_USERS
    OS_RSCRIPT_WINDOWS OS_RTOOLS_WINDOWS ANALYSIS_ROOT)
if [[ $submit == hpc ]]; then
    for _os_desktop_name in DIR_WINDOWS_DRIVE DIR_WINDOWS DIR_WINDOWS_PROGRAMDATA DIR_WINDOWS_PROGRAMFILES DIR_WINDOWS_PROGRAMFILES_X86 DIR_WINDOWS_USERS OS_RSCRIPT_WINDOWS OS_RTOOLS_WINDOWS; do
        export -n "$_os_desktop_name" "_OS_DEFAULT_$_os_desktop_name"
    done
    unset _os_desktop_name
fi
_os_default OS_FILES_ROOT "$DIR_FILES"
_os_default OS_SNP_LIST "$DIR_UKB/common/snp.lst"
_os_default OS_AXT_ROOT "$DIR_ANNOT/axt/37"
# All new analysis outputs are persistent and independent of --dir-tmp.
_os_default ANALYSIS_ROOT "$DIR_DAT/analysis"
_os_default DIR_ANALYSIS "$ANALYSIS_ROOT"
_os_default SCRATCH_ANALYSIS_ROOT "$ANALYSIS_ROOT"
_os_default PERSISTENT_ANALYSIS_ROOT "$ANALYSIS_ROOT"
_os_default METAG_ROOT "$DIR_METAG"
for _os_project in ALLELE GRID GWAS LE8; do
    case $_os_project in ALLELE) _os_leaf=allele3 ;; *) _os_leaf=${_os_project,,} ;; esac
    _os_default "${_os_project}_ANALYSIS_ROOT" "$ANALYSIS_ROOT/$_os_leaf"
done
_os_default GU_ANALYSIS_BASE "$ANALYSIS_ROOT/gu"
_os_default GRID_OUTPUT_ROOT "$GRID_ANALYSIS_ROOT"
_os_default ALLELE_ARCHAIC_ROOT "$DIR_GEN"
_os_default LD_ANALYSIS_ROOT "$ANALYSIS_ROOT/ld"
# Raw reads stay under DIR_BIG/metag/raw; generated state/results use analysis.
_os_default METAG_STATE_ROOT "$ANALYSIS_ROOT/metag"
_os_default METAG_REFERENCE_RAW "$DIR_ANNOT/metag/raw"
_os_default METAG_READS_RAW "$DIR_METAG/raw"
_os_default METAG_TOOLS_DIR "$DIR_SOFTWARE/metag_parallel"
# QiMing queue defaults: inspect current availability with queueinfo.
if [[ $submit == hpc ]]; then
    _os_default HPC_CPU_QUEUE 33
    _os_default HPC_GPU_QUEUE 4a100-40
    _os_default HPC_MEDIUM_QUEUE 38
    _os_default HPC_HIGHMEM_QUEUE 1t81c
    _os_default HPC_HUGEMEM_QUEUE v3-6t
    # Host RAM, not VRAM. Request one GPU unless a backend really uses more.
    # 4a100-40: at least 384 GiB host RAM in the supplied QiMing inventory.
    _os_default HPC_GPU_MEMORY_GIB 160
    _os_default HPC_WALLTIME 120:00
    _os_default ALLELE_JOB_QUEUE "$HPC_CPU_QUEUE"
    _os_default GU_QUEUE "$HPC_CPU_QUEUE"
    # dispatch.py selects a CPU queue that fits the summed worker budget.
    _os_default GWAS_QUEUE auto
    _os_default LD_QUEUE "$HPC_CPU_QUEUE"
    _os_default LE8_QUEUE "$HPC_CPU_QUEUE"
    _os_default LE8_GPU_QUEUE "$HPC_GPU_QUEUE"
    _os_default LE8_PYTHON "$DIR_WORK/.conda/envs/le8/bin/python3"
    _os_default LE8_RSCRIPT "$DIR_WORK/.conda/envs/le8/bin/Rscript"
    _os_default SAIGE_ENV_DIR "$DIR_WORK/.conda/envs/saige"
    _os_default GWAS_RUNTIME_ENV "$SCRIPTS_ROOT/gwas/f/runtime.env"
    _os_default GRID_HPC_CONFIG "$SCRIPTS_ROOT/grid/config/hpc.qiming.json"
fi
export MEMORY_DEFAULT_GIB
export dir_work=$DIR_WORK dir_data=$DIR_DATA dir_scratch=$DIR_SCRATCH
export dir_dat=$DIR_DAT dir_big=$DIR_BIG dir_tmp=$DIR_TMP
export SUBMIT_WAIT=${SUBMIT_WAIT:-FALSE}
unset _os_project _os_leaf

workflow_args() {
    WORKFLOW_ARGS=()
    while (($#)); do
        case "$1" in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
        case "$1" in
            --dir-dat|--dir-big|--dir-tmp)
                [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "$1 requires a directory" >&2; return 2; }
                [[ $2 == /* ]] || { echo "$1 requires an absolute path: $2" >&2; return 2; }
                case "$1" in
                    --dir-dat) export DIR_DAT=${2%/} ;;
                    --dir-big) export DIR_BIG=${2%/} DIR_DATA=${2%/} ;;
                    --dir-tmp) export DIR_TMP=${2%/} DIR_SCRATCH=${2%/} ;;
                esac
                shift 2 ;;
            --submit)
                if [[ ${2:-} == local || ${2:-} == hpc || ${2:-} == auto ]]; then
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
