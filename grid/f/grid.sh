#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
if [[ $GRID_DRY_RUN == TRUE ]]; then echo 'PLAN GRID: arg check -> bridge 1csx inputs -> LD -> transport -> blocked fit -> weights -> scores';exit 0;fi
# Internal GRID steps run together; evaluation uses the separate Yeval.sh entry.
exec {grid_lock}>"$gdir/run.lock"
flock -n "$grid_lock" || _grid_die "Another GRID run is active: $gdir"
rm -f "$gdir/GRID_RUN.txt"
bash "$ROOT/f/arg.sh" "$@"
grid_run python3 "$ROOT/f/grid_inputs.py" --trait "$trait" --gwas-dir "$GRID_GWAS_DIR" --snpinfo "$GRID_CSX_SNPINFO" --chrs "${CHRS[*]}" --out "$out"
steps=(ld transport fit weights score)
for step in "${steps[@]}"; do bash "$ROOT/f/$step.sh" "$@"; done
cat > "$gdir/GRID_RUN.txt" <<META
created=$(date -Is)
trait=$trait
chromosomes=${CHRS[*]}
method=genealogy-informed shrinkage of PRS-CSx population effects toward a shared effect
validation=blocked out-of-fold transportability model
local_genealogy=$GRID_ARG_TREES_DIR/chrCHR.variants.tsv.gz
ld_source=PRS-CSx HDF5 reference panels
participant_phenotype_used_for_weights=FALSE
META
echo "GRID completed: $out/scores/grid.tsv.gz"
