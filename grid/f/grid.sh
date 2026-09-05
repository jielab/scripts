#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
# Internal steps run together; only grid and eval are public GRID modules.
steps=(arg ld transport fit weights score)
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
