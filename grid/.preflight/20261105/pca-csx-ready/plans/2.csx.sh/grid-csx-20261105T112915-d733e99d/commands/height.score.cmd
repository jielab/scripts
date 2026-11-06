#!/usr/bin/env bash
#BSUB -J grid-csx.height.score
#BSUB -q v4-3t
#BSUB -n 32
#BSUB -W 120:00
#BSUB -M 256GB
#BSUB -R "span[hosts=1] rusage[mem=256GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/2.csx.sh/grid-csx-20261105T112915-d733e99d/logs/height.score.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/2.csx.sh/grid-csx-20261105T112915-d733e99d/logs/height.score.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 32 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: height.weights
# Scientific command: bash /work/sph-huangj/scripts/grid/2.csx.sh score --traits height --output-root /work/sph-huangj/analysis/grid --jobs 4
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/2.csx.sh/grid-csx-20261105T112915-d733e99d/plan.json --key height.score
