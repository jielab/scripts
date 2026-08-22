#!/bin/sh

## Submit file for imputation

#SBATCH --partition=compute
#SBATCH --job-name=Imputation
#SBATCH --account=sc-users
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=115
#SBATCH --mem-per-cpu=3G
#SBATCH --output=<path>%x-%A_%a.out
#SBATCH --array=1-6%6

## Change directory
cd <path>

## Logging
echo "[LOG] Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "[LOG] Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "[LOG] Job ID: $SLURM_JOB_ID"
echo "[LOG] Node ID: $SLURM_NODEID"
echo "[LOG] Node List: $SLURM_NODELIST"
echo "[LOG] Job ID: $SLURM_ARRAY_TASK_ID"

## Container and Directories
R_CONTAINER='<path>'
R_SCRIPT='<path>/01_imputation.R'
BIND_DIR="<path>,<path>"

## Run the Array Job
singularity exec --bind $BIND_DIR $R_CONTAINER Rscript $R_SCRIPT