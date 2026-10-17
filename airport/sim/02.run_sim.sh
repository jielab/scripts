#!/bin/bash
source /mnt/d/scripts/0f/console.sh
set -e
background=FALSE
case "${1:-}" in
  --background) background=TRUE ;;
  '') ;;
  *) echo 'Usage: 02.run_sim.sh [--background]' >&2; exit 2 ;;
esac


# 🚩 02.run_sim

label=run_sim
outdir=/work/sph-huangj/analysis/avia
mkdir -p $outdir

echo "#!/bin/bash
module load python/anaconda3/2020.7
source activate
conda activate R4.4.2
cp /work/sph-huangj/scripts/avia/02.run_sim.R ./
R CMD BATCH 02.run_sim.R
" >$outdir/02.run_sim.cmd
cd $outdir
if [[ $background == TRUE ]]; then
  bsub -q short -n 40 -J avia.$label -o $label.LOG -e $label.ERR <02.run_sim.cmd
else
  bash 02.run_sim.cmd
fi
