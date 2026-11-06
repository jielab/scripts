#!/usr/bin/env bash
#BSUB -J grid-prsformer.ldl.prsformer2
#BSUB -q 8v100-32-sc
#BSUB -n 8
#BSUB -W 120:00
#BSUB -M 160GB
#BSUB -R "span[hosts=1] rusage[mem=160GB/host]"
#BSUB -o /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/4.prsformer.sh/grid-prsformer-20261105T112917-d7bd76f2/logs/ldl.prsformer2.%J.out
#BSUB -e /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/4.prsformer.sh/grid-prsformer-20261105T112917-d7bd76f2/logs/ldl.prsformer2.%J.err
#BSUB -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 8 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: shared.prepare
# Scientific command: bash /work/sph-huangj/scripts/grid/4.prsformer.sh prsformer2 --traits ldl --analysis-root /work/sph-huangj/analysis/grid --device cuda --threads 8
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/scripts/grid/.preflight/20261105/pca-csx-ready/plans/4.prsformer.sh/grid-prsformer-20261105T112917-d7bd76f2/plan.json --key ldl.prsformer2
