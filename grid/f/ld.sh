#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
ld_hdf5_for(){
  local p=${1,,} c=$2
  find "$GRID_CSX_REF_DIR" -type f -iname "*chr${c}*.hdf5" -ipath "*${p}*" -print | sort | head -n1
}
make_ld(){
  local p c f dest
  for p in "${POPS[@]}"; do p=${p^^}; mkdir -p "$GRID_LDSCORE_DIR/$p"
    for c in "${CHRS[@]}"; do
      dest="$GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz"; [[ -s $dest && $GRID_REPLACE == FALSE ]] && continue
      f=$(ld_hdf5_for "$p" "$c"); [[ -s $f ]] || _grid_die "No PRS-CSx HDF5 LD file for $p chr$c below $GRID_CSX_REF_DIR"
      grid_run_logged "$out/log/grid.ld.$p.chr$c.log" python3 "$ROOT/f/extract_ld_scores.py" --hdf5 "$f" --pop "$p" --chr "$c" --out "$dest"
    done
  done
}

make_ld
echo "GRID ld completed: $gdir"
