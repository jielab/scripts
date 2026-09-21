#!/usr/bin/env bash
# Compare ancestry-matched PT, csx.auto, phenotype-tuned PRS-CSX and saved DiscoDivas.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
usage(){ cat <<'HELP'
Yeval — paired, held-out UKB prediction comparisons by genetic ancestry.

  ./Yeval.sh --trait height --type ct
  ./Yeval.sh --trait ldl --type ct --covar-name age,sex,PC1,PC2,drug.lipid
  ./Yeval.sh --trait t2dm --type t2e --covar-name age,sex,PC1,PC2
  ./Yeval.sh --trait t2dm --type dt --covar-name age,sex,PC1,PC2 --prevalence cohort

All evaluation outputs: /mnt/d/analysis/grid/Yeval/<trait>/ (no method/type folders).
Open report.html for all figures, performance and methods. Its PNG images and
linked downloads (plots.pdf, performance.tsv, cohort.tsv, methods.md) are kept,
along with eval.log and run.lock. Intermediate evaluation tables are not exported.

Methods:
  PT            Same-ancestry .jma.cojo SNP/refA/bJ, PLINK2 --score; cached pt.pgs.gz.
  PRS-CS-multi  csx.auto (requested display name; see methods.md for provenance).
  PRS-CSX       Joint regression of csx.AFR/EAS/EUR/SAS, trained within each fold.
  DiscoDivas    Saved disco.pgs.gz, evaluated on exactly the same held-out people.
  Supplement    csx.meta and each of the four individual csx.* scores.

Options:
  --trait Y --type ct|dt|t2e  Required. dt: 0/1 disease; t2e: censored survival.
  --score-dir DIR            /mnt/d/data/ukb/pgs
  --pgs-file FILE            Override csx.pgs.gz; --disco-file / --pt-file also supported.
  --pheno-file FILE          /mnt/d/data/ukb/phe/Rdata/all.rds (full ancestry cohort).
  --ancestry-file FILE       /mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz
  --group-col NAME           genetic_ancestry; column in phenotype or ancestry file.
  --covar-name LIST          age,sex,PC1,PC2; comma-separated, or none.
  --phenotype-col NAME       Default Y. Default t2dm dt is baseline ICD10 T2D status
                            from dated prevalent/incident status (see methods.md).
  --event-col / --time-col   Defaults Y.Yt2e / Y.t2e; used only for t2e.
  --prevalence K|cohort|EUR=K,AFR=K,EAS=K,SAS=K
                            dt only. Default cohort: ancestry-specific prevalence
                            before score/covariate filtering; cohort-based assumption.
                            Not used for t2e (Cox/Harrell C) or ct.
  --pca-file FILE            Reference-projected PCs; same default as 2disco.sh.
  --med-file FILE            /mnt/d/files/DiscoDivas/med.g1000.4pop.tsv
  --distance-pcs N           10; must match saved DiscoDivas PCA provenance.
  --distance-bins N          5 quantile bins within each ancestry.
  --folds N / --seed N       5 / 20260904; folds shared by all methods.
  --bootstrap N             200 paired subject resamples of fixed OOF predictions.
  --min-n N                 100; insufficient groups/bins are explicitly recorded.
  --dir-gwas DIR            /mnt/e/gwas/4grid/common; AFR also accepts t2dm.AFA.
  --dir-gen DIR             /mnt/e/ukbGen/37/imp; chr1..22 pfiles or bfiles.
  --pt-effect NAME          bJ (joint COJO effect); may explicitly select b.
  --threads N               4 for PLINK2. Scoring uses selected COJO variants.
  --remove FILE             /mnt/d/files/ukb.exclude.id; negative IDs also excluded.
  --out-root DIR            /mnt/d/analysis/grid/Yeval
  --check                   Validate evaluation inputs, without fitting or PT scoring.
  --allow-missing-scores    Explicit partial report; missing methods never zero-filled.

Legacy --method csx|disco|all is accepted; the report always compares all methods.
--type dt reports Liability R2 plus AUC/Brier; --type t2e reports Harrell C, not
Liability R2. New runs overwrite same-named files in the trait root, including
when switching between dt and t2e; outputs have no type prefix.
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
         --method|--pgs-file|--disco-file|--pheno-file|--ancestry-file|--group-col|--covar-name|--phenotype-col|--event-col|--time-col|--prevalence|--pca-file|--med-file|--distance-pcs|--distance-bins|--folds|--seed|--bootstrap|--min-n) :;;
         *) echo "Unknown option: $1" >&2; exit 2;;
       esac;shift 2;;
  esac
done
[[ $trait =~ ^[A-Za-z0-9_.-]+$ && $trait != . && $trait != .. ]] || { echo 'Invalid/missing --trait' >&2; exit 2; }
[[ $type == ct || $type == dt || $type == t2e ]] || { echo '--type must be ct, dt or t2e' >&2; exit 2; }
out="$outroot/$trait"; mkdir -p "$out"
exec {lock}>"$out/run.lock"; flock -n "$lock" || { echo "Evaluation already running: $out" >&2; exit 1; }
r=(env -u R_ENVIRON_USER -u R_LIBS_USER /usr/bin/Rscript)
{ printf 'Command: ';printf '%q ' "$ROOT/Yeval.sh" "${args[@]}";printf '\n'; } > "$out/eval.log"
echo "Evaluation output: $out"
if [[ $check == FALSE && -z $pt_file ]]; then
  source "$ROOT/f/environment.sh"
  python3 "$ROOT/f/yeval_pt.py" --trait "$trait" --dir-gwas "$gwas_dir" --dir-gen "$gen_dir" \
    --output "$score_dir/$trait/pt.pgs.gz" --effect "$pt_effect" \
    --threads "$threads" --remove "$remove" 2>&1 | tee -a "$out/eval.log"
fi
# Rscript reads its source incrementally. Snapshot both files so editing the
# workspace while a long evaluation runs cannot change the executing program.
runtime=$(mktemp -d "${TMPDIR:-/tmp}/yeval-runtime.XXXXXX")
trap 'rm -rf -- "$runtime"' EXIT
cp "$ROOT/f/Yeval.R" "$ROOT/f/yeval_plots.R" "$runtime/"
"${r[@]}" "$runtime/Yeval.R" "${args[@]}" 2>&1 | tee -a "$out/eval.log"
# Remove only known obsolete Yeval outputs, after a successful full report.
# Keep the lock inode stable: unlinking it could allow concurrent evaluations.
if [[ $check == FALSE ]]; then
  obsolete=(command.sh comparison.pdf combined_scores.pdf paired_improvement.pdf
    genetic_landscape.pdf distance_performance.pdf distance_performance.tsv
    fold_coefficients.tsv folds.tsv.gz genetic_distance.tsv.gz manifest.tsv
    methods.tsv paired_comparison.tsv predictions.tsv.gz prevalence.tsv skipped.tsv
    pt.commands.jsonl pt.log pt.matches.tsv pt.plink.log pt.variants.tsv)
  for name in "${obsolete[@]}"; do rm -f -- "$out/$name"; done
fi
