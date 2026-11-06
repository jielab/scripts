#!/usr/bin/env bash
#BSUB -J grid-pca.pca
#BSUB -q 33
#BSUB -n 10
#BSUB -W 24:00
#BSUB -M 32GB
#BSUB -R "span[hosts=1] rusage[mem=32GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/1.pca.sh/grid-pca-20261105T112913-090d4d60/logs/pca.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/1.pca.sh/grid-pca-20261105T112913-090d4d60/logs/pca.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 10 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: (none)
# Scientific command: bash /work/sph-huangj/scripts/grid/1.pca.sh --traits height,ldl,t2dm --output-root /work/sph-huangj/analysis/grid --jobs 4
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/1.pca.sh/grid-pca-20261105T112913-090d4d60/plan.json --key pca
