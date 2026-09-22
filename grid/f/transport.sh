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
  inputs=("$out/grid/inputs.signature" "$GRID_MED_FILE" "$ROOT/f/build_transport_table.py")
  [[ -z $GRID_EXTERNAL_AGE ]] || inputs+=("$GRID_EXTERNAL_AGE")
  for c in "${CHRS[@]}";do
    inputs+=("$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz")
    for p in "${POPS[@]}";do
      inputs+=("$out/sumstats/bychr/$trait.$p.chr$c.tsv.gz")
      [[ ! -f $GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz ]] || inputs+=("$GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz")
    done
  done
  key=$(python3 "$ROOT/f/pipeline_io.py" signature "$GRID_MAX_SNPS_PER_CHR" "$GRID_REQUIRE_LD" --files "${inputs[@]}")
  if [[ -s $t && -f $t.signature && $(cat "$t.signature") == "$key" && $GRID_REPLACE == FALSE ]];then echo 'SKIP matching GRID transport';return;fi
  rm -f -- "$t.signature"
  grid_run_logged "$out/log/grid.transport.log" "${cmd[@]}"
  printf '%s\n' "$key" > "$t.signature"
}

make_transport
echo "GRID transport completed: $gdir"
