#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
trait=${GRID_TRAIT,,}; POPS=($(grid_csv_words "$GRID_POPS")); CHRS=($(grid_expand_chrs "$GRID_CHRS"))
out="$GRID_OUTPUT_ROOT/$trait"; gdir="$out/grid"; mkdir -p "$gdir" "$out/log" "$out/scores" "$gdir/weights" "$gdir/model"
