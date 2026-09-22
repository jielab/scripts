#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
trait=${GRID_TRAIT,,}; POPS=($(grid_csv_words "$GRID_POPS")); CHRS=($(grid_expand_chrs "$GRID_CHRS"))
suffix=''; [[ ${CHRS[*]} == "$(seq -s ' ' 1 22)" ]] || suffix="/chr$(IFS=,; echo "${CHRS[*]}")"
[[ $trait == height || $trait == ldl || $trait == t2dm ]] || _grid_die 'Unknown trait'
for c in "${CHRS[@]}";do [[ $c =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || _grid_die "Invalid autosome: $c";done
out="$GRID_OUTPUT_ROOT/$trait$suffix"; gdir="$out/grid"; mkdir -p "$gdir" "$out/log" "$out/scores" "$gdir/weights" "$gdir/model"
