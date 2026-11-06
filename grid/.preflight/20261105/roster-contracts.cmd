#!/usr/bin/env bash
set -uo pipefail
[[ -n ${LSB_JOBID:-} ]] || exit 2
cd /work/sph-huangj/scripts/grid
source ../0f/workflow_env.sh
workflow_init grid
AUDIT=/work/sph-huangj/scripts/grid/.preflight/20261105
export PATH="$GRID_CONDA_ENV/bin:$PATH" PYTHONPATH="$PWD/f:$PWD/2test:$PWD/../0f" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 GRID_CACHE_ROOT="$AUDIT/cache"
python3 "$AUDIT/data-contracts.py" > "$AUDIT/data-contracts-fixed.log" 2>&1
printf 'data-contracts\t%s\n' "$?" > "$AUDIT/roster-exits.tsv"
python3 "$AUDIT/cohort-alignment.py" > "$AUDIT/cohort-alignment.log" 2>&1
printf 'cohort-alignment\t%s\n' "$?" >> "$AUDIT/roster-exits.tsv"
cat "$AUDIT/roster-exits.tsv"
