#!/usr/bin/env bash
# Compare PRS models and display PRS-CSx across ancestry groups and genetic distance.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
usage(){ cat <<'HELP'
Yeval — paired out-of-fold PRS comparisons by genetic ancestry.

  ./Yeval.sh --trait height --type ct --covar-name age,sex,PC1,PC2
  ./Yeval.sh --trait ldl --type ct --covar-name age,sex,PC1,PC2,drug.lipid
  ./Yeval.sh --trait t2dm --type t2e --covar-name age,sex,PC1,PC2,drug.dm,drug.htn
  ./Yeval.sh --trait t2dm --type dt --phenotype-col t2dm

Methods: COJO, PRS-CSx-auto-meta, PRS-CSx (four-score regression),
DiscoDivas-tuned. Saved DiscoDivas-untuned and fixed-meta are supplements.
Tuning is performed inside each outer training fold; no phi-grid selection.

Inputs:
  --score-dir DIR            /mnt/d/data/ukb/pgs; trait/csx.pgs.gz
  --pgs-file FILE            Explicit CSx table; --pt-file / --disco-file also supported.
  --pheno-file FILE          /mnt/d/data/ukb/phe/Rdata/all.rds (RDS or TSV/gzip)
  --ancestry-file FILE       /mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz
  --group-col NAME           genetic_ancestry
  --pca-file FILE            Same folder/ukb.discodivas.pca.tsv.gz
  --med-file FILE            /mnt/d/files/DiscoDivas/med.g1000.4pop.tsv
  --grid-file FILE           Optional GRID table; adds GRID-tuned and saved GRID scores.
  --covar-name LIST          age,sex,PC1,PC2; comma-separated, or none
  --phenotype-col NAME       Default trait; default t2dm dt derives baseline Yr2e/Yt2e.
  --event-col / --time-col   Default trait.Yt2e / trait.t2e
  --remove FILE             /mnt/d/files/ukb.exclude.id; negative IDs also excluded

Evaluation:
  --type ct|dt|t2e           ct/dt: OOF partial R2; dt also AUC/Brier; t2e: Harrell C.
                            Binary R2 is OBSERVED SCALE, not liability R2.
                            t2e also reports covariates-only C and paired delta C.
  --disco-tune TRUE|FALSE    TRUE; FALSE evaluates only the saved untuned Disco score.
  --disco-a LIST             1,1,1,1 in AFR,EAS,EUR,SAS order (same as 2disco.sh).
  --min-anchor N            100 training people per ancestry and fold
  --distance-pcs N          10; distance to the 1KG EUR reference center
  --distance-bins N         Up to 10 quantile bins per original ancestry group
  --min-bin-events N        20 events/cases AND non-events/controls per bin
                            Sparse groups use fewer bins; only PRS-CSx is plotted.
  --folds N / --seed N      5 / 20260904
  --bootstrap N             200; 0 disables intervals for a smoke run
  --min-n N                 100 per target/bin
  --write-predictions TRUE|FALSE  FALSE; TRUE writes one compressed OOF table
  --prevalence SPEC          cohort / K / EUR=...,AFR=...; descriptive context only
  --out-root DIR            /mnt/d/analysis/grid/Yeval
  --allow-missing-scores     Explicit partial report
  --check                   Validate inputs; no fitting or COJO scoring

COJO scoring if --pt-file is absent:
  --dir-gwas DIR            /mnt/e/gwas/4grid/common
  --dir-gen DIR             /mnt/e/ukbGen/37/imp
  --pt-effect bJ|b           bJ; --threads N: 4

Outputs stay in <out-root>/<trait>/; changing outcome type overwrites this report.
The four-panel distance figure pairs categorical ancestry with continuous distance.
distance_performance.tsv contains the plotted bin estimates and sample/event counts.
A completed run has SUCCESS. Failed runs do not replace the previous report.
Set GRID_RSCRIPT to choose an R executable; the activated grid environment is preferred.
HELP
}
case "${1:-}" in -h|--help|help|'') usage; exit 0;; esac
args=("$@")
trait='' type='' outroot=/mnt/d/analysis/grid/Yeval score_dir=/mnt/d/data/ukb/pgs
pt_file='' gwas_dir=/mnt/e/gwas/4grid/common gen_dir=/mnt/e/ukbGen/37/imp
pt_effect=bJ threads=4 remove=/mnt/d/files/ukb.exclude.id check=FALSE
while (($#)); do
  case "$1" in
    --check) check=TRUE;shift;;
    --allow-missing-scores) shift;;
    *) [[ $# -ge 2 && $2 != --* ]] || { echo "Missing value: $1" >&2; exit 2; }
       case "$1" in
         --trait) trait=$2;; --type) type=$2;; --out-root) outroot=$2;;
         --score-dir) score_dir=$2;; --pt-file) pt_file=$2;; --dir-gwas) gwas_dir=$2;;
         --dir-gen) gen_dir=$2;; --pt-effect) pt_effect=$2;; --threads) threads=$2;;
         --remove) remove=$2;;
         --method|--pgs-file|--disco-file|--pheno-file|--ancestry-file|--group-col|--covar-name|--phenotype-col|--event-col|--time-col|--prevalence|--pca-file|--med-file|--distance-pcs|--distance-bins|--min-bin-events|--folds|--seed|--bootstrap|--min-n|--disco-tune|--disco-a|--min-anchor|--write-predictions|--grid-file) :;;
         *) echo "Unknown option: $1" >&2; exit 2;;
       esac;shift 2;;
  esac
done
[[ $trait =~ ^[A-Za-z0-9_.-]+$ && $trait != . && $trait != .. ]] || { echo 'Invalid/missing --trait' >&2; exit 2; }
[[ $type == ct || $type == dt || $type == t2e ]] || { echo '--type must be ct, dt or t2e' >&2; exit 2; }
out="$outroot/$trait"; mkdir -p "$out"
exec {lock}>"$out/run.lock"; flock -n "$lock" || { echo "Evaluation already running: $out" >&2; exit 1; }
source "$ROOT/f/environment.sh"
source "$ROOT/f/r_runtime.sh"
grid_select_r data.table,ggplot2,survival,pROC,patchwork
r=("${GRID_R[@]}")
{ printf 'Command: ';printf '%q ' "$ROOT/Yeval.sh" "${args[@]}";printf '\n'; } > "$out/eval.log"
echo "Evaluation output: $out"
[[ $check == TRUE ]] || rm -f -- "$out/SUCCESS"
if [[ $check == FALSE && -z $pt_file ]]; then
  source "$ROOT/f/environment.sh"
  python3 "$ROOT/f/yeval_pt.py" --trait "$trait" --dir-gwas "$gwas_dir" --dir-gen "$gen_dir" \
    --output "$score_dir/$trait/pt.pgs.gz" --effect "$pt_effect" \
    --threads "$threads" --remove "$remove" 2>&1 | tee -a "$out/eval.log"
fi
# Rscript reads its source incrementally. Snapshot all R files so editing the
# workspace while a long evaluation runs cannot change the executing program.
runtime=$(mktemp -d "${TMPDIR:-/tmp}/yeval-runtime.XXXXXX")
trap 'rm -rf -- "$runtime"' EXIT
cp "$ROOT/f/Yeval.R" "$ROOT/f/yeval_plots.R" "$ROOT/f/yeval_disco.R" "$runtime/"
mkdir -p "$runtime/report"
"${r[@]}" "$runtime/Yeval.R" "${args[@]}" --run-dir "$runtime/report" 2>&1 | tee -a "$out/eval.log"
if [[ $check == FALSE ]]; then
  for file in "$runtime/report"/*; do
    cp -- "$file" "$out/$(basename -- "$file").tmp"
    mv -f -- "$out/$(basename -- "$file").tmp" "$out/$(basename -- "$file")"
  done
  [[ -f $runtime/report/predictions.tsv.gz ]] || rm -f -- "$out/predictions.tsv.gz"
  date -Is > "$out/SUCCESS"
fi
# Remove only known obsolete Yeval outputs, after a successful full report.
# Keep the lock inode stable: unlinking it could allow concurrent evaluations.
if [[ $check == FALSE ]]; then
  obsolete=(command.sh comparison.pdf combined_scores.pdf paired_improvement.pdf
    genetic_landscape.pdf genetic_landscape.png distance_performance.pdf
    fold_coefficients.tsv folds.tsv.gz genetic_distance.tsv.gz manifest.tsv
    methods.tsv paired_comparison.tsv prevalence.tsv skipped.tsv
    pt.commands.jsonl pt.log pt.matches.tsv pt.plink.log pt.variants.tsv)
  for name in "${obsolete[@]}"; do rm -f -- "$out/$name"; done
fi
