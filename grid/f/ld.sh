#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
ld_hdf5_for(){
  local p=${1,,} c=$2
  local ref_type=1kg d f
  [[ $(basename -- "$GRID_CSX_SNPINFO") != snpinfo_mult_ukbb_hm3 ]] || ref_type=ukbb
  d="$GRID_CSX_REF_DIR/ldblk_${ref_type}_$p"
  [[ -d $d ]] || d="$GRID_CSX_REF_DIR/ldblk_${ref_type}_${p^^}"
  f="$d/ldblk_${ref_type}_chr$c.hdf5";[[ -s $f ]] || return 1;printf '%s\n' "$f"
}
make_ld(){
  local p c f dest
  for p in "${POPS[@]}"; do p=${p^^}; mkdir -p "$GRID_LDSCORE_DIR/$p"
    for c in "${CHRS[@]}"; do
      dest="$GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz"
      f=$(ld_hdf5_for "$p" "$c"); [[ -s $f ]] || _grid_die "No PRS-CSx HDF5 LD file for $p chr$c below $GRID_CSX_REF_DIR"
      key=$(python3 "$ROOT/f/pipeline_io.py" signature "$p" "$c" --files "$f" "$ROOT/f/extract_ld_scores.py")
      [[ -s $dest && -f $dest.signature && $(cat "$dest.signature") == "$key" && $GRID_REPLACE == FALSE ]] && continue
      grid_run_logged "$out/log/grid.ld.$p.chr$c.log" python3 "$ROOT/f/extract_ld_scores.py" --hdf5 "$f" --pop "$p" --chr "$c" --out "$dest.tmp.gz"
      mv -f -- "$dest.tmp.gz" "$dest";printf '%s\n' "$key" > "$dest.signature"
    done
  done
}

make_ld
echo "GRID ld completed: $gdir"
