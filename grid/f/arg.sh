#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
# shellcheck source=f/common.sh
source "$ROOT/f/common.sh" "$@"
grid_parse_args "$@"

read -r -a CHRS <<< "$(grid_expand_chrs "$GRID_CHRS")"
(( ${#CHRS[@]} )) || _grid_die 'No chromosomes selected'
case "$GRID_ARG_ACTION" in
  all|check) ;;
  prepare|infer|convert|features|affinity)
    echo "ERROR: ARG construction moved to refGen.sh data preparation." >&2
    echo "Run: bash /mnt/d/scripts/0data/refGen.sh make-arg --method needle --dir-gen $(dirname -- "$GRID_ARG_HAP_DIR") --dir-pfile $GRID_ARG_HAP_DIR --map-dir ${GRID_ARG_MAP_DIR:-MAP_DIR} --ancestry-file $GRID_ANCESTRY_FILE --chr $GRID_CHRS" >&2
    exit 2
    ;;
  *) _grid_die "bad --arg-action=$GRID_ARG_ACTION (GRID now supports check/all only)";;
esac

missing=0
for c in "${CHRS[@]}"; do
  files=(
    "$GRID_ARG_OUT/argn/chr$c.argn"
    "$GRID_ARG_TREES_DIR/chr$c.trees"
    "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"
    "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv"
    "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
  )
  for f in "${files[@]}"; do
    if [[ ! -s $f ]]; then echo "ERROR: missing GRID ARG-Needle data for chr$c: $f" >&2; missing=1; fi
  done
done
if ((missing)); then
  echo "Prepare the reusable ARG data with:" >&2
  echo "  bash /mnt/d/scripts/0data/refGen.sh make-arg --method needle --dir-gen $(dirname -- "$GRID_ARG_HAP_DIR") --dir-pfile $GRID_ARG_HAP_DIR --map-dir ${GRID_ARG_MAP_DIR:-MAP_DIR} --ancestry-file $GRID_ANCESTRY_FILE --chr $GRID_CHRS" >&2
  exit 2
fi

for c in "${CHRS[@]}"; do
  gzip -t "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
  python3 - "$GRID_ARG_TREES_DIR/chr$c.trees" "$c" <<'PY'
import sys,tskit
ts=tskit.load(sys.argv[1])
if min(ts.num_trees,ts.num_samples,ts.num_sites,ts.num_mutations) < 1:
    raise SystemExit(
        f"ERROR: invalid GRID ARG-Needle tree for chr{sys.argv[2]}: "
        f"trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}"
    )
print(
    f"GRID ARG CHECK PASS chr{sys.argv[2]} trees={ts.num_trees} "
    f"samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}"
)
PY
done
echo "GRID ARG-Needle data ready: $GRID_ARG_OUT"
