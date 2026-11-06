#!/usr/bin/env bash
#BSUB -J grid-grid.t2dm.twins
#BSUB -q 8v100-32-sc
#BSUB -n 8
#BSUB -W 120:00
#BSUB -M 160GB
#BSUB -R "span[hosts=1] rusage[mem=160GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/grid/grid-grid-20261105T013825-7ed15c63/logs/t2dm.twins.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/grid/grid-grid-20261105T013825-7ed15c63/logs/t2dm.twins.%J.err
#BSUB -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 8 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: t2dm.prepare
# Scientific command: bash /work/sph-huangj/scripts/grid/grid.sh twins --traits t2dm --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.json --device cuda --threads 8
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/final-plans/grid/grid-grid-20261105T013825-7ed15c63/plan.json --key t2dm.twins
