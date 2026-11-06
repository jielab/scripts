#!/usr/bin/env bash
set -euo pipefail
source /work/sph-huangj/scripts/grid/f/0.common.sh
grid_activate_environment
export GRID_CSX_BACKEND=cpu GRID_CSX_GPU=0 GRID_CSX_MODULE=csx
export GRID_CSX_GPU_MIN_BLOCK=2048 GRID_CSX_BLAS_THREADS=1
export GRID_CSX_LD_CACHE=/work/sph-huangj/scripts/grid/.preflight/20261105/cache/csx-ld GRID_CSX_LD_CACHE_GIB=16
