#!/usr/bin/env bash
set -euo pipefail
source /work/sph-huangj/scripts/grid/f/0.common.sh
grid_activate_environment
export GRID_CSX_BACKEND=cpu GRID_CSX_GPU=0 GRID_CSX_MODULE=csx
