#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
args=("$@")
run(){ local m=$1; shift; echo "===== GRID $m $(date -Is) =====" >&2; bash "$ROOT/f/$m.sh" "$@"; }
# Fail before expensive work when required software/data are absent.
run preflight --scope core "${args[@]}"
grid_is_true "$GRID_RUN_PCA" && run pca "${args[@]}"
grid_is_true "$GRID_RUN_ANCESTRY" && run ancestry "${args[@]}"
if grid_is_true "$GRID_RUN_ARG"; then run preflight --scope arg "${args[@]}"; run arg "${args[@]}"; fi
for trait in $(grid_csv_words "$GRID_TRAITS"); do
  targs=("${args[@]}" --trait "$trait")
  grid_is_true "$GRID_RUN_CSX" && run csx "${targs[@]}"
  grid_is_true "$GRID_RUN_DISCO" && run disco "${targs[@]}"
  if grid_is_true "$GRID_RUN_GRID"; then run preflight --scope grid "${targs[@]}"; run grid "${targs[@]}"; fi
  grid_is_true "$GRID_RUN_EVAL" && run eval "${targs[@]}"
done
echo "GRID all completed: $GRID_OUTPUT_ROOT"
