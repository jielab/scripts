#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
trait=${GRID_TRAIT,,}; out="$GRID_OUTPUT_ROOT/$trait"; [[ -s $out/scores/grid.tsv.gz ]] || _grid_die "Run grid.sh grid --trait $trait first"
pca=''
for f in "$GRID_PCA_FILE" "$GRID_OUTPUT_ROOT/reference/pca/ukb.projected.tsv.gz" /mnt/d/data/ukb/phe/pca/ukb.pca.projected.tsv.gz /mnt/d/data/ukb/phe/pca/ukb_pca_projected.tsv.gz /mnt/d/data/ukb/phe/pca/ukb.pca.tsv.gz /mnt/d/data/ukb/phe/pca/ukb.pca.rds; do [[ -s $f ]] && { pca=$f; break; }; done
cmd=(Rscript "$ROOT/f/eval_grid.R" --trait "$trait" --phe "$GRID_PHE_FILE" --ancestry "$GRID_ANCESTRY_FILE" --score-dir "$out/scores" --out-dir "$out/eval" --phenotype-col "$GRID_PHENOTYPE_COL" --covariates "$GRID_EVAL_COVARIATES" --bootstrap "$GRID_EVAL_REPEATS" --folds "$GRID_EVAL_FOLDS" --min-group "$GRID_MIN_GROUP" --tuned "$GRID_EVAL_TUNED" --write-predictions "$GRID_WRITE_PREDICTIONS" --seed "$GRID_SEED")
[[ -z $pca ]] || cmd+=(--pca "$pca")
mkdir -p "$out/log"; grid_run_logged "$out/log/eval.log" "${cmd[@]}"
echo "Evaluation completed: $out/eval/performance.tsv"
