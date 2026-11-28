#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
case "${1:-}" in
  ""|help|-h|--help) cat <<'HELP'
GRID learns from one PRS-CSx META coefficient set per trait on the shared BED.
  ./grid.sh prepare --grch 38 --traits height,ldl,t2dm
  ./grid.sh run --traits height --method block.glmnet,block.glmnet.disco,blocknet.disco
  ./grid.sh run --traits t2dm --method block.glmnet.disco,blocknet.disco
  ./grid.sh run    # three traits, both META models, four-channel comparators, ABM/twins, paper
  ./grid.sh run --method universal.auto,universal.fixed    # only the single-weight models
Methods: universal.auto, universal.fixed, block.glmnet, block.lightgbm, block.glmnet.disco, blocknet.disco, blocknet.disco.global, block.stack,
         abm, abm.deep, twins, architecture; robustness is an optional sensitivity analysis.
         direct.linear, direct.global_mlp, direct.regional_ridge, direct.blocknet, blocknet.twins
Options: --method CSV, --blocknet-run cpu|gpu [cpu], --cores N, --wait TRUE|FALSE,
         --dir-ld-bed DIR --ld-bed uni.bed (or four ancestry BED filenames),
         --config FILE, --analysis-root DIR, --replace TRUE|FALSE, --plan-only.
./grid.sh predict --model MODEL.joblib --blocks DIR --covariates AGE_SEX.tsv --out scores.tsv.gz
One-weight blocks: main/TRAIT/universal.prepare.{auto,fixed}; four-channel scores: main/TRAIT/regions.
Use ./grid.paper.sh run for a report-only refresh.
HELP
  exit 0;;
esac
exec bash "$root/f/run_entry.sh" grid "$@"
