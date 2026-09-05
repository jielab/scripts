#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
trait=${GRID_TRAIT,,}; out="$GRID_OUTPUT_ROOT/$trait"; csx="$out/scores/csx.tsv.gz"; [[ -s $csx ]] || _grid_die "Run grid.sh csx --trait $trait first"
# Locate the PCA projection produced by the preserved PCA module.
pca=''
for f in "$GRID_OUTPUT_ROOT/reference/pca/ukb.projected.tsv.gz" /mnt/d/data/ukb/phe/pca/ukb.pca.projected.tsv.gz /mnt/d/data/ukb/phe/pca/ukb_pca_projected.tsv.gz /mnt/d/data/ukb/phe/pca/ukb.pca.tsv.gz /mnt/d/data/ukb/phe/pca/ukb.pca.rds; do [[ -s $f ]] && { pca=$f; break; }; done
mkdir -p "$out/scores" "$out/log"
cmd=(Rscript "$ROOT/f/disco_zero.R" --csx "$csx" --out "$out/scores/disco_zero.tsv.gz" --centers "$GRID_MED_FILE" --ancestry "$GRID_ANCESTRY_FILE" --n-pcs "$GRID_DISTANCE_PCS")
[[ -z $pca ]] || cmd+=(--pca "$pca")
grid_run_logged "$out/log/disco.log" "${cmd[@]}"
echo "DiscoDivas zero-shot comparator completed: $out/scores/disco_zero.tsv.gz"
