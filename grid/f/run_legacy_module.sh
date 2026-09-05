#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
module=$1; shift
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
legacy="$ROOT/grid.legacy.sh"; [[ -s $legacy ]] || { echo "ERROR: missing $legacy" >&2; exit 2; }
help=$(bash "$legacy" help 2>&1 || true)
cmd=(bash "$legacy" "$module")
add(){ local k=$1 v=$2; grep -Fq -- "$k" <<<"$help" && cmd+=("$k" "$v"); }
add --data-root "$GRID_DATA_ROOT"
add --dir-gen "$GRID_TARGET_DIR"
add --dir-imp "$GRID_IMP_DIR"
add --dir-gwas "$GRID_GWAS_DIR"
add --phe-file "$GRID_PHE_FILE"
add --output-dir "$GRID_OUTPUT_ROOT"
add --jobs "$GRID_JOBS"
add --threads "$GRID_THREADS"
add --med-file "$GRID_MED_FILE"
add --pca-weight "$GRID_PCA_WEIGHT"
add --cov-pcs "$GRID_COV_PCS"
add --distance-pcs "$GRID_DISTANCE_PCS"
add --auto-ancestry "$GRID_AUTO_ANCESTRY"
add --auto-ancestry-file "$GRID_ANCESTRY_FILE"
add --ancestry-prob-threshold "$GRID_ANCESTRY_PROB_MIN"
add --replace "$GRID_REPLACE"
printf '[GRID legacy-%s]' "$module" >&2; printf ' %q' "${cmd[@]}" >&2; echo >&2
exec "${cmd[@]}"
