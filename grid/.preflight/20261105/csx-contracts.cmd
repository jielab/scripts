#!/usr/bin/env bash
set -uo pipefail
[[ -n ${LSB_JOBID:-} ]] || exit 2
cd /work/sph-huangj/scripts/grid
source ../0f/workflow_env.sh
workflow_init grid
AUDIT=/work/sph-huangj/scripts/grid/.preflight/20261105
export WORKFLOW_WORKER=1 WORKFLOW_RESOURCE_MANAGER=lsf WORKFLOW_MEMORY_GIB=16 GRID_MEMORY_CAP_GB=16 GRID_THREADS=2 GRID_JOBS=1 GRID_CACHE_ROOT="$AUDIT/cache"
export PATH="$GRID_CONDA_ENV/bin:$PATH" PYTHONPATH="$PWD/f:$PWD/../0f" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
bash 2.csx.sh --check --output-root "$AUDIT/native" > "$AUDIT/2.csx.sh.recheck.log" 2>&1
printf 'csx\t%s\n' "$?" > "$AUDIT/contracts-exits.tsv"
python3 "$AUDIT/data-contracts.py" > "$AUDIT/data-contracts.log" 2>&1
printf 'data-contracts\t%s\n' "$?" >> "$AUDIT/contracts-exits.tsv"
cat "$AUDIT/contracts-exits.tsv"
