#!/bin/bash
source /mnt/d/scripts/0f/console.sh
set -e
background=FALSE
case "${1:-}" in
  --background) background=TRUE ;;
  '') ;;
  *) echo 'Usage: 01.cfr.sh [--background]' >&2; exit 2 ;;
esac


# 🚩 01.cfr

dir=/work/sph-huangj

for X in $(
	cd $dir/data/avia
	cat iso3c.txt
); do

	label=$X
	outdir=$dir/analysis/avia/cfr/$X
	mkdir -p $outdir
	if [[ -f $outdir/$label.Rout ]]; then
		echo $label already run
		continue
	fi

	echo "#!/bin/bash
	module load python/anaconda3/2020.7
	source activate
	conda activate R4.4.2
	cat /work/sph-huangj/scripts/avia/01.cfr.R | sed -e '7 s/?/$X/' > $label.R
	R CMD BATCH $label.R
	" >$outdir/$label.cmd
	cd $outdir
	if [[ $background == TRUE ]]; then
		bsub -q short -n 40 -J avia.$label -o $label.LOG -e $label.ERR <$label.cmd
	else
		bash "$label.cmd"
	fi

done
