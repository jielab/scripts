#!/usr/bin/env bash
# Source this file. No conda activation or network installation is required.
HPC_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export LD_CODE_DIR="$HPC_SCRIPT_DIR/f"
# This HPC's older conda initialization reads PS1 even in batch shells.
PS1="${PS1-}"
if [[ "${CONDA_PREFIX:-}" != /work/sph-huangj/.conda/envs/py311 ]]; then
    source /share/apps/anaconda3/2020.7/etc/profile.d/conda.sh
    conda activate /work/sph-huangj/.conda/envs/py311
fi
if ! type module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh
fi
module load gcc/11.2.0
export PYTHON="/work/sph-huangj/.conda/envs/py311/bin/python3"
export PATH="/work/sph-huangj/.conda/envs/py311/bin:/work/sph-huangj/software/bin:$PATH"
# Use the existing Zstandard static library, avoiding a conflicting conda
# libstdc++ on LD_LIBRARY_PATH. cxx.sh substitutes only the -lzstd argument.
export LD_ZSTD_PREFIX="${LD_ZSTD_PREFIX:-/share/apps/anaconda3/2020.7}"
export CXX="$HPC_SCRIPT_DIR/f/cxx.sh"
export PLINK2="${PLINK2:-/work/sph-huangj/software/bin/plink2}"
export GRCH="${GRCH:-${GENOME_BUILD:-38}}"
if [[ "$GRCH" != 37 && "$GRCH" != 38 ]]; then
    echo 'GRCH must be 37 or 38' >&2
    return 2
fi
export GENOME_BUILD="$GRCH"
export INPUT_FORMAT="${INPUT_FORMAT:-pgen}"
# When switching builds in a shell that sourced this file before, retarget
# only our canonical defaults. Custom paths remain explicit overrides.
for old_build in 37 38; do
    [[ "${OUT:-}" != "/data/sph-huangj/analysis/ld/$old_build" ]] || unset OUT
    [[ "${PFILE_DIR:-}" != "/data/sph-huangj/gen/1kg/$old_build/pfile" ]] || unset PFILE_DIR
    [[ "${PANEL:-}" != "/data/sph-huangj/gen/1kg/$old_build/samples.txt" ]] || unset PANEL
    [[ "${LD_TMP_DIR:-}" != "/data/sph-huangj/analysis/ld/$old_build/tmp" ]] || unset LD_TMP_DIR
    [[ "${LD_INPUT_MANIFEST:-}" != "$HPC_SCRIPT_DIR/f/input_sizes.$old_build.tsv" ]] || unset LD_INPUT_MANIFEST
    [[ "${LD_INPUT_MANIFEST:-}" != "$HPC_SCRIPT_DIR/hpc/input_sizes.$old_build.tsv" ]] || unset LD_INPUT_MANIFEST
done
if [[ "$GRCH" == 38 ]]; then
    export REFERENCE_COUNT_POLICY=off
else
    export REFERENCE_COUNT_POLICY="${REFERENCE_COUNT_POLICY:-warn}"
fi
export OUT="${OUT:-/data/sph-huangj/analysis/ld/$GRCH}"
export PFILE_DIR="${PFILE_DIR:-/data/sph-huangj/gen/1kg/$GRCH/pfile}"
export PANEL="${PANEL:-/data/sph-huangj/gen/1kg/$GRCH/samples.txt}"
export LD_TMP_DIR="${LD_TMP_DIR:-$OUT/tmp}"
export LD_RESUME_DIR="${LD_RESUME_DIR-}"
export LD_STORAGE="${LD_STORAGE:-blocked-zstd}"
export MAX_D_KB="${MAX_D_KB:-${PILOT_MAX_D_KB:-40000}}"
export PARALLEL_CHROMS="${PARALLEL_CHROMS:-1}"
export THREADS="${THREADS:-4}"
export MEMORY_MB="${MEMORY_MB:-8192}"
export PYTHONDONTWRITEBYTECODE=1
export LD_INPUT_MANIFEST="${LD_INPUT_MANIFEST:-$HPC_SCRIPT_DIR/f/input_sizes.$GRCH.tsv}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
