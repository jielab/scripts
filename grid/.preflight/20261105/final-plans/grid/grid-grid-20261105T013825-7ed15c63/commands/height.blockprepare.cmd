#!/usr/bin/env bash
#BSUB -J grid-grid.height.blockprepare
#BSUB -q 38
#BSUB -n 32
#BSUB -W 120:00
#BSUB -M 256GB
#BSUB -R "span[hosts=1] rusage[mem=256GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/grid/grid-grid-20261105T013825-7ed15c63/logs/height.blockprepare.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/grid/grid-grid-20261105T013825-7ed15c63/logs/height.blockprepare.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 32 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: shared.prepare
# Scientific command: bash /work/sph-huangj/scripts/grid/2.csx.sh blockprepare --traits height --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.json --device cpu --threads 8
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/grid/grid-grid-20261105T013825-7ed15c63/plan.json --key height.blockprepare
