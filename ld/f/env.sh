#!/usr/bin/env bash
export DATE_SCRATCH=${DATE_SCRATCH:-2026-11-05}
# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
source "$_paths_directory/0f/0os.sh"
unset _paths_directory

# Source this file. No conda activation or network installation is required.
HPC_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export LD_CODE_DIR="$HPC_SCRIPT_DIR/f"
source "$HPC_SCRIPT_DIR/../0f/0os.sh"
PS1=${PS1-}
if [[ $submit == hpc ]]; then
    if [[ ${CONDA_PREFIX:-} != "$DIR_WORK/.conda/envs/py311" ]]; then
        set +u
        source "${LD_CONDA_SH:-$OS_GU_CONDA/etc/profile.d/conda.sh}"
        conda activate "$DIR_WORK/.conda/envs/py311"
        set -u
    fi
    export PYTHON=${LD_PYTHON:-$DIR_WORK/.conda/envs/py311/bin/python3}
    export LD_ZSTD_PREFIX=${LD_ZSTD_PREFIX:-$CONDA_PREFIX}
    export CXX="$HPC_SCRIPT_DIR/f/cxx.sh"
else
    export PYTHON=${LD_PYTHON:-python3}
    export CXX=${CXX:-g++}
fi
export PATH="$DIR_SOFTWARE/bin:$PATH"
export PLINK2=${PLINK2:-$DIR_SOFTWARE/bin/plink2}
export GRCH="${GRCH:-${GENOME_BUILD:-38}}"
if [[ "$GRCH" != 37 && "$GRCH" != 38 ]]; then
    echo 'GRCH must be 37 or 38' >&2
    return 2
fi
export GENOME_BUILD="$GRCH"
export INPUT_FORMAT="${INPUT_FORMAT:-pgen}"
# When switching builds in a shell that sourced this file before, retarget
# only our canonical defaults. Custom paths remain explicit overrides.
if [[ "$GRCH" == 38 ]]; then
    export REFERENCE_COUNT_POLICY=off
else
    export REFERENCE_COUNT_POLICY="${REFERENCE_COUNT_POLICY:-warn}"
fi
if [[ $submit == hpc ]]; then
    export SCRATCH_ANALYSIS_ROOT="/scratch/$DATE_SCRATCH/$(id -un)/analysis"
    if [[ $GRCH == 38 ]]; then
        _os_default OUT "$SCRATCH_ANALYSIS_ROOT/ld"
    else
        _os_default OUT "$SCRATCH_ANALYSIS_ROOT/ld/$GRCH"
    fi
else
    export SCRATCH_ANALYSIS_ROOT="$ANALYSIS_ROOT"
    _os_default OUT "$SCRATCH_ANALYSIS_ROOT/ld/$GRCH"
fi
_os_default PFILE_DIR "$DIR_GEN/1kg/$GRCH/pfile"
_os_default PANEL "$DIR_GEN/1kg/$GRCH/samples.txt"
_os_default LD_TMP_DIR "$OUT/tmp"
export LD_RESUME_DIR="${LD_RESUME_DIR-}"
export LD_STORAGE="${LD_STORAGE:-blocked-zstd}"
export MAX_D_KB="${MAX_D_KB:-${PILOT_MAX_D_KB:-40000}}"
export PARALLEL_CHROMS="${PARALLEL_CHROMS:-1}"
export THREADS="${THREADS:-4}"
export MEMORY_MB="${MEMORY_MB:-8192}"
export PYTHONDONTWRITEBYTECODE=1
# Reset obsolete canonical manifests inherited before the code directory moved.
case "${LD_INPUT_MANIFEST:-}" in
    ${SCRIPTS_ROOT}/ld/f/input_sizes.*.tsv|${SCRIPTS_ROOT}/ld/hpc/input_sizes.*.tsv)
        unset LD_INPUT_MANIFEST ;;
esac
export LD_QUEUE="${LD_QUEUE:-${HPC_CPU_QUEUE:-33}}"
export LD_MEMORY_GIB="${LD_MEMORY_GIB:-$MEMORY_DEFAULT_GIB}"
export LD_WALLTIME="${LD_WALLTIME:-${HPC_WALLTIME:-120:00}}"
export LD_INPUT_MANIFEST="${LD_INPUT_MANIFEST:-$HPC_SCRIPT_DIR/f/input_sizes.$GRCH.tsv}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

export AUTO_EXTEND_D=${AUTO_EXTEND_D:-TRUE}
