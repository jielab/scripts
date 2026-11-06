#!/bin/bash -l
#BSUB -q 33
#BSUB -n 10
#BSUB -M 16GB
#BSUB -W 0:20
#BSUB -R "span[hosts=1] rusage[mem=16GB/host]"
#BSUB -J ukb.schema-check
#BSUB -o "/work/sph-huangj/scripts/ukb/.validation-preflight.HxkHlK/schema.out"
#BSUB -e "/work/sph-huangj/scripts/ukb/.validation-preflight.HxkHlK/schema.err"
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec /work/sph-huangj/.conda/envs/le8/bin/Rscript "/work/sph-huangj/scripts/ukb/.validation-preflight.HxkHlK/schema.R"
