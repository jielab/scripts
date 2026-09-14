#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
make_transport(){
  for c in "${CHRS[@]}"; do
    [[ -s $GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz ]] || _grid_die "Missing ARG feature file $GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
    if grid_is_true "$GRID_REQUIRE_LD"; then for p in "${POPS[@]}"; do [[ -s $GRID_LDSCORE_DIR/${p^^}/chr$c.ldscore.tsv.gz ]] || _grid_die "Missing LD score for ${p^^} chr$c"; done; fi
  done
  t="$gdir/transport.tsv.gz"
  cmd=(python3 "$ROOT/f/build_transport_table.py" --trait "$trait" --pops "$GRID_POPS" --chrs "${CHRS[*]}" --sumstats-dir "$out/sumstats/bychr" --arg-dir "$GRID_ARG_TREES_DIR" --ldscore-dir "$GRID_LDSCORE_DIR" --centers "$GRID_MED_FILE" --max-snps-per-chr "$GRID_MAX_SNPS_PER_CHR" --out "$t")
  [[ -z $GRID_EXTERNAL_AGE ]] || cmd+=(--external-age "$GRID_EXTERNAL_AGE")
  [[ -s $t && $GRID_REPLACE == FALSE ]] || grid_run_logged "$out/log/grid.transport.log" "${cmd[@]}"
}

make_transport
echo "GRID transport completed: $gdir"
