#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
make_fit(){
  local t="$gdir/transport.tsv.gz"
  [[ -s $t ]] || _grid_die 'GRID transport table missing; rerun ./grid.sh grid'
  grid_run_logged "$out/log/grid.transport.fit.log" python3 "$ROOT/f/fit_transport_model.py" --input "$t" --out-dir "$gdir/model" --ridge-alpha "$GRID_RIDGE_ALPHA" --conservation-min "$GRID_CONSERVATION_MIN" --conservation-max "$GRID_CONSERVATION_MAX" --seed "$GRID_SEED"
}

make_fit
echo "GRID fit completed: $gdir"
