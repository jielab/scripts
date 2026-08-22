#!/bin/sh

## submit file for LOFO

#SBATCH --partition=compute
#SBATCH --job-name=LOFO
#SBATCH --account=sc-users
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --output=<path>/%x-%A.out

## change directory
cd <path>

echo "[LOG] Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "[LOG] Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "[LOG] Job ID: $SLURM_JOB_ID"
echo "[LOG] Node ID: $SLURM_NODEID"
echo "[LOG] Node List: $SLURM_NODELIST"
echo "[LOG] Job ID: $SLURM_ARRAY_TASK_ID"

echo "[LOG] Running script with ${input} as input and ${output} as output"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "[LOG] Starting script at: $date"

R_CONTAINER='<path>'
R_SCRIPT='<path>/03_LOFO.R'
BIND_DIR="<path>,<path>"

singularity exec \
  --bind $BIND_DIR \
  $R_CONTAINER Rscript $R_SCRIPT ${input} ${output}

printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "[LOG] Finishing script at: $date"
echo "[LOG] Done!"
exit $(echo $?)
