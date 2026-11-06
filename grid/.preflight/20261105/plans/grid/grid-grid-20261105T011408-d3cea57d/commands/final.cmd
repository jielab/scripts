#!/usr/bin/env bash
#BSUB -J grid-grid.final
#BSUB -q 33
#BSUB -n 4
#BSUB -W 24:00
#BSUB -M 16GB
#BSUB -R "span[hosts=1] rusage[mem=16GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/plans/grid/grid-grid-20261105T011408-d3cea57d/logs/final.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/plans/grid/grid-grid-20261105T011408-d3cea57d/logs/final.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 4 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: architecture
# Scientific command: bash /work/sph-huangj/scripts/grid/grid.sh final --traits height,ldl,t2dm --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.json
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/plans/grid/grid-grid-20261105T011408-d3cea57d/plan.json --key final
