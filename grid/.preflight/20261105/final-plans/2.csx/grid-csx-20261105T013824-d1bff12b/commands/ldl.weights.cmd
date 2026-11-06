#!/usr/bin/env bash
#BSUB -J grid-csx.ldl.weights
#BSUB -q 38
#BSUB -n 32
#BSUB -W 120:00
#BSUB -M 256GB
#BSUB -R "span[hosts=1] rusage[mem=256GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/2.csx/grid-csx-20261105T013824-d1bff12b/logs/ldl.weights.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/2.csx/grid-csx-20261105T013824-d1bff12b/logs/ldl.weights.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 32 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: (none)
# Scientific command: bash /work/sph-huangj/scripts/grid/2.csx.sh weights --traits ldl --output-root /work/sph-huangj/analysis/grid
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/2.csx/grid-csx-20261105T013824-d1bff12b/plan.json --key ldl.weights
