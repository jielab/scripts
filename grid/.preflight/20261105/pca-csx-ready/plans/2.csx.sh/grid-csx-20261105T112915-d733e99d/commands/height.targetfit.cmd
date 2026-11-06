#!/usr/bin/env bash
#BSUB -J grid-csx.height.targetfit
#BSUB -q v4-3t
#BSUB -n 8
#BSUB -W 24:00
#BSUB -M 32GB
#BSUB -R "span[hosts=1] rusage[mem=32GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/2.csx.sh/grid-csx-20261105T112915-d733e99d/logs/height.targetfit.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/2.csx.sh/grid-csx-20261105T112915-d733e99d/logs/height.targetfit.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 8 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: shared.prepare
# Scientific command: bash /work/sph-huangj/scripts/grid/2.csx.sh targetfit --traits height --analysis-root /work/sph-huangj/analysis/grid
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/2.csx.sh/grid-csx-20261105T112915-d733e99d/plan.json --key height.targetfit
