#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
make_weights(){
  [[ -s $gdir/model/variant_conservation.tsv.gz ]] || _grid_die 'GRID model output missing; rerun ./grid.sh grid'
  [[ -s $out/csx/manifest.tsv ]] || _grid_die 'Run ./csx.sh first'
  grid_run_logged "$out/log/grid.weights.log" python3 "$ROOT/f/make_grid_weights.py" --weights-dir "$out/csx/weights" --conservation "$gdir/model/variant_conservation.tsv.gz" --manifest "$out/csx/manifest.tsv" --out-dir "$gdir/weights" --pops "$GRID_POPS" --chrs "${CHRS[*]}"
}

make_weights
echo "GRID weights completed: $gdir"
