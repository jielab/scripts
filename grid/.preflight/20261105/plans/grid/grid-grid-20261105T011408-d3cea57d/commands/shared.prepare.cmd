#!/usr/bin/env bash
#BSUB -J grid-grid.shared.prepare
#BSUB -q 38
#BSUB -n 32
#BSUB -W 120:00
#BSUB -M 256GB
#BSUB -R "span[hosts=1] rusage[mem=256GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/plans/grid/grid-grid-20261105T011408-d3cea57d/logs/shared.prepare.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/plans/grid/grid-grid-20261105T011408-d3cea57d/logs/shared.prepare.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 32 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: (none)
# Scientific command: bash /work/sph-huangj/scripts/grid/2.csx.sh prepare --traits height,ldl,t2dm --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.json
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/plans/grid/grid-grid-20261105T011408-d3cea57d/plan.json --key shared.prepare
