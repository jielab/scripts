#!/usr/bin/env bash
#BSUB -J grid-prsformer.t2dm.all
#BSUB -q 8v100-32-sc
#BSUB -n 32
#BSUB -W 120:00
#BSUB -M 256GB
#BSUB -R "span[hosts=1] rusage[mem=256GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261106/queue-prsformer/plans/grid-prsformer-20261105T234343-b215b7e0/logs/t2dm.all.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261106/queue-prsformer/plans/grid-prsformer-20261105T234343-b215b7e0/logs/t2dm.all.%J.err
#BSUB -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 32 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: (none)
# Scientific command: bash /work/sph-huangj/scripts/grid/4.prsformer.sh all --traits t2dm --analysis-root /work/sph-huangj/analysis/grid --config config/abm.qiming-v100.json
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261106/queue-prsformer/plans/grid-prsformer-20261105T234343-b215b7e0/plan.json --key t2dm.all
