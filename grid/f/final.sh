#!/usr/bin/env bash
# Final performance comparison; retains the original Yeval fitting and plotting engine.
set -euo pipefail
umask 077
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
usage() {
	cat <<'HELP'
GRID final — compare saved CSx, Disco, PRSformer and ABM results.

Example usage:
  cd /mnt/d/scripts/grid
  bash grid.sh final
  # Optional: one trait or validation only.
  bash grid.sh final --trait height
  bash grid.sh final --check

Defaults: height,ldl,t2dm; ct for height/ldl, baseline dt for t2dm.
Saved learned predictions are evaluated directly, without recalibration or refitting.
All methods use the same common test participants; baseline regressions are fitted
only in common development participants. Different saved splits use their intersection.
Without learned results, the original Yeval out-of-fold baseline comparison is used.

Inputs:
  --traits LIST / --trait NAME  Default height,ldl,t2dm.
  --score-dir DIR              /mnt/d/data/ukb/pgs (CSx, Disco, existing COJO).
  --grid-root DIR              /mnt/d/analysis/grid/grid (saved ABM outputs).
  --prsformer-root DIR         Prefer grid/benchmark; otherwise /mnt/d/analysis/grid.
  --grid-file FILE|none        Override/disable ABM; requires --grid-split FILE.
  --prsformer-file FILE|none   Override/disable PRSformer; requires --prsformer-split FILE.
  --pheno-file FILE            /mnt/d/data/ukb/phe/Rdata/all.rds.
  --ancestry-file FILE         /mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz.
  --pca-file / --med-file      Original Yeval projected PCs / reference centres.
  --covar-name LIST            age,sex,PC1,PC2; or none.
  --remove FILE               /mnt/d/files/ukb.exclude.id.
  --pgs-file / --disco-file / --pt-file FILE  Explicit baseline score tables.

Evaluation:
  --type ct|dt|t2e             Optional override for one trait. Existing t2dm learned
                              models predict baseline disease, so they require dt.
  --disco-tune TRUE|FALSE      FALSE: evaluate saved Disco; TRUE: development-only tuning.
  --posterior-mode auto|required|off  auto: use posterior if available, otherwise no SD panel.
  --bootstrap N               200 paired bootstrap resamples; 0 disables intervals.
  --min-n N                   100 per ancestry/bin.
  --write-predictions TRUE|FALSE  FALSE; private individual predictions use RDS.
  --out-root DIR              /mnt/d/analysis/grid/final; one subdirectory per trait.
  --check                     Validate inputs/splits only; no fitting or formal outputs.
  --dry-run                   Print evaluation commands only.
  --score-cojo                Explicitly generate missing COJO scores using final.py.
                              Default final uses existing results only.

Original Yeval options for distance bins, posterior uncertainty, covariates, folds,
COJO inputs and prevalence remain supported. Each figure has a matching XLSX table.
Missing optional methods are reported; malformed or incompatible learned results fail.
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0 ;; esac
traits=height,ldl,t2dm
type=
outroot=/mnt/d/analysis/grid/final
score_dir=/mnt/d/data/ukb/pgs
check=FALSE
dry_run=FALSE
score_cojo=FALSE
gwas_dir=/mnt/f/gwas/4grid/common
gen_dir=/mnt/f/gen/ukb/37/imp
pt_effect=bJ
threads=4
remove=/mnt/d/files/ukb.exclude.id
pt_file=
forward=()
while (($#)); do
	case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
	case $1 in
		--check) check=TRUE; forward+=("$1"); shift ;;
		--dry-run) dry_run=TRUE; shift ;;
		--score-cojo) score_cojo=TRUE; shift ;;
		--allow-missing-scores) forward+=("$1"); shift ;;
		*)
			[[ $# -ge 2 && $2 != --* ]] || { echo "Missing value: $1" >&2; exit 2; }
			case $1 in
				--trait|--traits) traits=$2; shift 2; continue ;;
				--type) type=$2; shift 2; continue ;;
				--out-root) outroot=$2; shift 2; continue ;;
				--score-dir) score_dir=$2 ;;
				--pt-file) pt_file=$2 ;;
				--dir-gwas) gwas_dir=$2 ;;
				--dir-gen) gen_dir=$2 ;;
				--pt-effect) pt_effect=$2 ;;
				--threads) threads=$2 ;;
				--remove) remove=$2 ;;
				--method|--pgs-file|--disco-file|--pheno-file|--ancestry-file|--group-col|--covar-name|--phenotype-col|--event-col|--time-col|--prevalence|--pca-file|--med-file|--distance-pcs|--distance-bins|--min-bin-events|--folds|--seed|--bootstrap|--min-n|--disco-tune|--disco-a|--min-anchor|--write-predictions|--grid-root|--grid-file|--grid-split|--prsformer-root|--prsformer-file|--prsformer-split|--posterior-file|--posterior-mode|--genetic-variance-file|--training-centers|--pca-space|--allow-chromosome-subset|--individual-metric|--individual-max-points|--distance-source) ;;
				*) echo "Unknown option: $1" >&2; exit 2 ;;
			esac
			forward+=("$1" "$2"); shift 2 ;;
	esac
done
[[ $traits =~ ^(height|ldl|t2dm)(,(height|ldl|t2dm))*$ ]] || { echo 'Use --traits height,ldl,t2dm' >&2; exit 2; }
IFS=, read -r -a trait_list <<<"$traits"
[[ -z $type || ( ${#trait_list[@]} == 1 && $type =~ ^(ct|dt|t2e)$ ) ]] || { echo '--type requires one trait and ct|dt|t2e' >&2; exit 2; }
runtime=
trap '[[ -z $runtime ]] || rm -rf -- "$runtime"' EXIT
if [[ $dry_run == FALSE ]]; then
	source "$ROOT/f/0.common.sh"
	grid_activate_environment
	grid_select_r data.table,ggplot2,survival,pROC,patchwork,openxlsx,jsonlite,digest,zip
	r=("${GRID_R[@]}")
else
	r=(Rscript)
fi
run_trait() {
	local trait=$1 kind=$type out="$outroot/$1" cache log file
	[[ -n $kind ]] || { if [[ $trait == t2dm ]]; then kind=dt; else kind=ct; fi; }
	local command=("${r[@]}" "$ROOT/f/final.R" --trait "$trait" --type "$kind" --out-root "$outroot" "${forward[@]}")
	printf 'RUN '; printf '%q ' "${command[@]}"; printf '\n'
	[[ $dry_run == FALSE ]] || return 0
	cache="/tmp/grid/final/$(printf '%s' "$out" | sha256sum | cut -c1-16)"
	mkdir -p "$cache"
	exec {lock}>"$cache/run.lock"
	flock -n "$lock" || { echo "Final evaluation already running: $out" >&2; exit 1; }
	log="$cache/final.log"
	if [[ $check == FALSE && $score_cojo == TRUE && -z $pt_file ]]; then
		python3 "$ROOT/f/final.py" --trait "$trait" --dir-gwas "$gwas_dir" --dir-gen "$gen_dir" \
			--output "$score_dir/$trait/Yeval.cojo.rds" --effect "$pt_effect" --threads "$threads" --remove "$remove" 2>&1 | tee "$log"
	fi
	runtime=$(mktemp -d /tmp/grid-final.XXXXXX)
	cp "$ROOT/f/final.R" "$runtime/final.R"
	cp "$ROOT/../0f/results.R" "$runtime/results.R"
	GRID_RESULTS_R="$runtime/results.R" "${r[@]}" "$runtime/final.R" --trait "$trait" --type "$kind" \
		--out-root "$outroot" "${forward[@]}" --run-dir "$runtime/report" 2>&1 | tee -a "$log"
	if [[ $check == FALSE ]]; then
		mkdir -p "$out"
		for file in "$runtime/report"/*; do
			cp -p -- "$file" "$out/$(basename -- "$file")"
			cmp -s -- "$file" "$out/$(basename -- "$file")"
		done
		for file in final.paired_improvement.png final.paired_improvement.xlsx final.distance_bins.png final.distance_bins.xlsx final.individuals.rds final.predictions.rds; do
			[[ -f $runtime/report/$file ]] || rm -f -- "$out/$file"
		done
	fi
	rm -rf -- "$runtime"
	runtime=
	flock -u "$lock"
	exec {lock}>&-
}
for trait in "${trait_list[@]}"; do run_trait "$trait"; done
