#!/bin/sh

## This script submits the job for 02.2_feature_selection_all_individual.R
## it creates 18 array jobs, running feature selection for each outcome

#SBATCH --partition=compute
#SBATCH --job-name=IndividualDomainTop20
#SBATCH --account=sc-users
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-18%18
#SBATCH --output=<path>/%x-%A-%2a.out

## define the outcome and followup
outc_values=(<"outcomes">)
cdat_values=(<"followup">)

## use SLURM_ARRAY_TASK_ID to select the appropriate values
index=$((SLURM_ARRAY_TASK_ID - 1))
outc="${outc_values[$index]}"
cdat="${cdat_values[$index]}"

## do some logging
echo "[LOG] Running script with outcome: ${outc} and follow-up: ${cdat}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "[LOG] Starting script at: $date"

## define container and script
R_CONTAINER='<path>'
R_SCRIPT='<path>/02.2_feature_selection_all_individual.R'
BIND_DIR="<path>,<path>"

## execute the script in the container
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT ${outc} ${cdat}

## some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
