#!/usr/bin/env bash
# Shared environment, configuration, resource caps, checks and pipeline runner.


# 🚩 grid_activate_environment
grid_activate_environment() {
	export PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/tmp/python-cache
	export TMPDIR=/tmp TMP=/tmp TEMP=/tmp
	# Make the supported environment usable for non-interactive GRID calls without
	# requiring `conda activate` first.
	GRID_CONDA_ENV=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
	if [[ -x "$GRID_CONDA_ENV/bin/python3" && ${CONDA_PREFIX:-} != "$GRID_CONDA_ENV" ]]; then
		grid_conda_init="$(dirname -- "$(dirname -- "$GRID_CONDA_ENV")")/etc/profile.d/conda.sh"
		if [[ -f "$grid_conda_init" ]]; then
			# Third-party activation/deactivation hooks may read unset variables (e.g.
			# Anaconda's GeoTIFF hook). Relax nounset only while Conda runs its hooks,
			# then restore the caller's setting on both success and failure.
			grid_conda_nounset=FALSE
			[[ $- != *u* ]] || grid_conda_nounset=TRUE
			set +u
			if source "$grid_conda_init" && conda activate "$GRID_CONDA_ENV"; then
				grid_conda_status=0
			else
				grid_conda_status=$?
			fi
			if [[ $grid_conda_nounset == TRUE ]]; then set -u; fi
			unset grid_conda_nounset
			if ((grid_conda_status != 0)); then
				printf 'ERROR: Could not activate Conda environment: %s\n' "$GRID_CONDA_ENV" >&2
				unset grid_conda_init
				return "$grid_conda_status"
			fi
			unset grid_conda_status
		else
			export PATH="$GRID_CONDA_ENV/bin:$PATH"
		fi
		unset grid_conda_init
	fi
	if [[ -d "$GRID_CONDA_ENV/lib/R/library" ]]; then
		export R_ENVIRON_USER=/dev/null
		export R_LIBS_USER="$GRID_CONDA_ENV/lib/R/library"
	fi
	export PYTHONNOUSERSITE=1
	if [[ -z ${ARG_NEEDLE_HOME:-} && -x "$GRID_CONDA_ENV/bin/python3" ]]; then
		ARG_NEEDLE_HOME=$(
			"$GRID_CONDA_ENV/bin/python3" - <<'PY' 2>/dev/null || true
import importlib.util
from pathlib import Path
spec = importlib.util.find_spec("arg_needle.scripts.infer_args_advanced")
if spec is not None and spec.origin is not None:
    print(Path(spec.origin).resolve().parents[2])
PY
		)
		export ARG_NEEDLE_HOME
	fi
}
# Prefer the activated environment; explicitly allow GRID_RSCRIPT overrides.
grid_select_r() {
	local packages=$1 candidate
	local -a runtime
	if [[ -n ${GRID_RSCRIPT:-} ]]; then
		GRID_R=("$GRID_RSCRIPT")
		# The system R must use its own startup configuration and package libraries.
		if [[ $GRID_RSCRIPT -ef /usr/bin/Rscript ]]; then
			GRID_R=(env -u R_ENVIRON_USER -u R_LIBS_USER "$GRID_RSCRIPT")
		fi
		"${GRID_R[@]}" -e "p<-strsplit('$packages',',',fixed=TRUE)[[1]];stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))" || return
		return
	fi
	for candidate in "$(command -v Rscript || true)" /usr/bin/Rscript; do
		[[ -n $candidate && -x $candidate ]] || continue
		runtime=("$candidate")
		if [[ $candidate -ef /usr/bin/Rscript ]]; then
			runtime=(env -u R_ENVIRON_USER -u R_LIBS_USER "$candidate")
		fi
		if "${runtime[@]}" -e "p<-strsplit('$packages',',',fixed=TRUE)[[1]];stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))" >/dev/null 2>&1; then
			GRID_R=("${runtime[@]}")
			return
		fi
	done
	echo "ERROR: No usable R runtime; install packages: $packages, or set GRID_RSCRIPT" >&2
	return 1
}

source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../0f" && pwd -P)/memory_cap.sh"


# 🚩 grid_configure
grid_configure() {
	# Shared configuration and option parsing for GRID modules.
	set -euo pipefail
	export LC_ALL=C
	export PYTHONDONTWRITEBYTECODE=1
	export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
	export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
	export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
	export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

	GRID_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
	GRID_F="$GRID_ROOT/f"
	GRID_ORIGINAL_ARGS=("$@")

	# Core defaults.
	GRID_DATA_ROOT=${GRID_DATA_ROOT:-/mnt/d}
	GRID_TRAIT=${GRID_TRAIT:-height}
	GRID_TRAITS=${GRID_TRAITS:-height,ldl,t2dm}
	GRID_POPS=${GRID_POPS:-AFR,EAS,EUR,SAS}
	GRID_CHRS=${GRID_CHRS:-1-22}
	GRID_JOBS=${GRID_JOBS:-4}
	GRID_THREADS=${GRID_THREADS:-8}
	GRID_REPLACE=${GRID_REPLACE:-FALSE}
	GRID_DRY_RUN=${GRID_DRY_RUN:-FALSE}
	GRID_GWAS_DIR=${GRID_GWAS_DIR:-/mnt/f/gwas/4grid/common}
	if [[ -z ${GRID_TARGET_DIR+x} ]]; then
		GRID_TARGET_DIR=/mnt/d/data/ukb/gen/typ
		if [[ ! -s $GRID_TARGET_DIR/chr1.pgen && ! -s $GRID_TARGET_DIR/chr1.bed && -s /mnt/f/gen/ukb/37/hap/chr1.pgen ]]; then
			GRID_TARGET_DIR=/mnt/f/gen/ukb/37/hap
		fi
	fi
	GRID_IMP_DIR=${GRID_IMP_DIR:-/mnt/f/gen/ukb/37/imp}
	GRID_PHE_FILE=${GRID_PHE_FILE:-/mnt/d/data/ukb/phe/Rdata/phe.rds}
	GRID_OUTPUT_ROOT=${GRID_OUTPUT_ROOT:-/mnt/d/analysis/grid}
	GRID_KEEP=${GRID_KEEP:-}
	GRID_REMOVE=${GRID_REMOVE:-/mnt/d/files/ukb.exclude.id}
	GRID_N_GWAS=${GRID_N_GWAS:-}

	# PRS-CSx.
	GRID_BIM_EXPLICIT=${GRID_CSX_BIM_PREFIX:+TRUE}
	GRID_BIM_EXPLICIT=${GRID_BIM_EXPLICIT:-FALSE}
	GRID_SNPINFO_EXPLICIT=${GRID_CSX_SNPINFO:+TRUE}
	GRID_SNPINFO_EXPLICIT=${GRID_SNPINFO_EXPLICIT:-FALSE}
	GRID_CSX_REF_DIR=${GRID_CSX_REF_DIR:-/mnt/f/refLD/csx}
	GRID_CSX_SNPINFO=${GRID_CSX_SNPINFO:-$GRID_CSX_REF_DIR/snpinfo_mult_1kg_hm3}
	GRID_CSX_BIM_PREFIX=${GRID_CSX_BIM_PREFIX:-$GRID_TARGET_DIR/ukb_array}
	GRID_PHI=${GRID_PHI:-1e-2}
	GRID_MCMC_ITER=${GRID_MCMC_ITER:-4000}
	GRID_MCMC_BURNIN=${GRID_MCMC_BURNIN:-2000}
	GRID_MCMC_THIN=${GRID_MCMC_THIN:-5}
	GRID_SEED=${GRID_SEED:-20260904}
	GRID_SUMSTATS_CHUNK=${GRID_SUMSTATS_CHUNK:-500000}

	# PCA/ancestry/DiscoDivas.
	GRID_MED_FILE=${GRID_MED_FILE:-/mnt/d/files/DiscoDivas/med.g1000.4pop.tsv}
	GRID_PCA_WEIGHT=${GRID_PCA_WEIGHT:-/mnt/d/files/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight}
	GRID_COV_PCS=${GRID_COV_PCS:-20}
	GRID_DISTANCE_PCS=${GRID_DISTANCE_PCS:-10}
	GRID_PCA_FILE=${GRID_PCA_FILE:-}
	GRID_ANCESTRY_FILE=${GRID_ANCESTRY_FILE:-}
	GRID_ANCESTRY_PROB_MIN=${GRID_ANCESTRY_PROB_MIN:-0.90}
	GRID_ANCHOR_PROB_MIN=${GRID_ANCHOR_PROB_MIN:-0.999}
	GRID_ANCHOR_MAX_PER_GROUP=${GRID_ANCHOR_MAX_PER_GROUP:-10000}

	# ARG-Needle. The phased BGENs are UKB Field 22438, GRCh37.
	GRID_ARG_ACTION=${GRID_ARG_ACTION:-check}
	GRID_ARG_HAP_DIR=${GRID_ARG_HAP_DIR:-/mnt/f/gen/ukb/37/hap}
	GRID_ARG_MAP_DIR=${GRID_ARG_MAP_DIR:-}
	GRID_ARG_MAP_PATTERN=${GRID_ARG_MAP_PATTERN:-}
	if [[ -z ${GRID_ARG_OUT+x} ]]; then
		GRID_ARG_OUT="$(dirname -- "$GRID_ARG_HAP_DIR")/arg.needle"
		GRID_ARG_OUT_AUTO=TRUE
	else
		GRID_ARG_OUT_AUTO=FALSE
	fi
	GRID_ARG_SCRATCH=${GRID_ARG_SCRATCH:-}
	GRID_ARG_FULL=${GRID_ARG_FULL:-FALSE}
	GRID_ARG_MAX_INDIVIDUALS=${GRID_ARG_MAX_INDIVIDUALS:-20000}
	GRID_ARG_SEED_HAPLOTYPES=${GRID_ARG_SEED_HAPLOTYPES:-4000}
	GRID_ARG_THREADS=${GRID_ARG_THREADS:-8}
	GRID_ARG_JOBS=${GRID_ARG_JOBS:-1}
	GRID_ARG_ANCHORS_PER_POP=${GRID_ARG_ANCHORS_PER_POP:-1000}
	GRID_ARG_WINDOW_BP=${GRID_ARG_WINDOW_BP:-1000000}
	GRID_ARG_NEEDLE_HOME=${GRID_ARG_NEEDLE_HOME:-${ARG_NEEDLE_HOME:-}}
	GRID_ARG_NORMALIZE=${GRID_ARG_NORMALIZE:-TRUE}
	GRID_ARG_MAF_MIN=${GRID_ARG_MAF_MIN:-0.001}
	GRID_ARG_GENO_MAX=${GRID_ARG_GENO_MAX:-0.05}
	GRID_ARG_KEEP=${GRID_ARG_KEEP:-}
	GRID_ARG_TREES_DIR=${GRID_ARG_TREES_DIR:-$GRID_ARG_OUT/trees}
	GRID_ARG_AFFINITY=${GRID_ARG_AFFINITY:-FALSE}

	# Evolutionary transport model and GRID scoring.
	GRID_LDSCORE_DIR=${GRID_LDSCORE_DIR:-}
	GRID_REQUIRE_LD=${GRID_REQUIRE_LD:-TRUE}
	GRID_EXTERNAL_AGE=${GRID_EXTERNAL_AGE:-}
	GRID_MAX_SNPS_PER_CHR=${GRID_MAX_SNPS_PER_CHR:-0}
	GRID_RIDGE_ALPHA=${GRID_RIDGE_ALPHA:-10}
	GRID_TRANSPORT_MODEL=${GRID_TRANSPORT_MODEL:-evolutionary_full}
	GRID_CONSERVATION_MIN=${GRID_CONSERVATION_MIN:-0.05}
	GRID_CONSERVATION_MAX=${GRID_CONSERVATION_MAX:-0.995}
	GRID_LOCAL=${GRID_LOCAL:-FALSE}

	# Evaluation.
	GRID_EVAL_REPEATS=${GRID_EVAL_REPEATS:-100}
	GRID_EVAL_FOLDS=${GRID_EVAL_FOLDS:-5}
	GRID_MIN_GROUP=${GRID_MIN_GROUP:-100}
	GRID_EVAL_TUNED=${GRID_EVAL_TUNED:-TRUE}
	GRID_EVAL_ZERO_SHOT=${GRID_EVAL_ZERO_SHOT:-TRUE}
	GRID_WRITE_PREDICTIONS=${GRID_WRITE_PREDICTIONS:-FALSE}
	GRID_ETHNICITY_COL=${GRID_ETHNICITY_COL:-auto}
	GRID_PHENOTYPE_COL=${GRID_PHENOTYPE_COL:-auto}
	GRID_EVAL_COVARIATES=${GRID_EVAL_COVARIATES:-auto}

	_grid_die() {
		echo "ERROR: $*" >&2
		exit 2
	}
	_grid_need_value() { [[ $# -ge 2 && -n ${2:-} && ${2:-} != --* ]] || _grid_die "$1 requires a value"; }
	grid_bool() { case "${1^^}" in TRUE | T | 1 | YES | Y | ON) echo TRUE ;; FALSE | F | 0 | NO | N | OFF) echo FALSE ;; *) return 2 ;; esac }
	grid_is_true() { [[ $(grid_bool "$1") == TRUE ]]; }
	grid_require_bool() {
		local x
		x=$(grid_bool "$2") || _grid_die "$1 must be TRUE/FALSE"
		printf '%s\n' "$x"
	}

	grid_parse_args() {
		while (($#)); do
			case "$1" in
				--trait)
					_grid_need_value "$@"
					GRID_TRAIT=${2,,}
					shift 2
					;;
				--traits)
					_grid_need_value "$@"
					GRID_TRAITS=${2,,}
					shift 2
					;;
				--pops)
					_grid_need_value "$@"
					GRID_POPS=${2^^}
					shift 2
					;;
				--chrs | --chr | --arg-chrs)
					_grid_need_value "$@"
					GRID_CHRS=$2
					shift 2
					;;
				--jobs)
					_grid_need_value "$@"
					GRID_JOBS=$2
					shift 2
					;;
				--threads)
					_grid_need_value "$@"
					GRID_THREADS=$2
					shift 2
					;;
				--data-root)
					_grid_need_value "$@"
					GRID_DATA_ROOT=${2%/}
					shift 2
					;;
				--dir-gwas)
					_grid_need_value "$@"
					GRID_GWAS_DIR=${2%/}
					shift 2
					;;
				--dir-gen)
					_grid_need_value "$@"
					GRID_TARGET_DIR=${2%/}
					shift 2
					;;
				--dir-imp)
					_grid_need_value "$@"
					GRID_IMP_DIR=${2%/}
					shift 2
					;;
				--phe-file | --phe)
					_grid_need_value "$@"
					GRID_PHE_FILE=$2
					shift 2
					;;
				--output-root | --output-dir)
					_grid_need_value "$@"
					GRID_OUTPUT_ROOT=${2%/}
					shift 2
					;;
				--keep | --target-ids)
					_grid_need_value "$@"
					GRID_KEEP=$2
					shift 2
					;;
				--remove | --exclude-ids)
					_grid_need_value "$@"
					GRID_REMOVE=$2
					shift 2
					;;
				--n-gwas)
					_grid_need_value "$@"
					GRID_N_GWAS=$2
					shift 2
					;;
				--replace)
					_grid_need_value "$@"
					GRID_REPLACE=$(grid_require_bool "$1" "$2")
					shift 2
					;;
				--dry-run)
					_grid_need_value "$@"
					GRID_DRY_RUN=$(grid_require_bool "$1" "$2")
					shift 2
					;;

				--csx-ref-dir | --ref-dir)
					_grid_need_value "$@"
					GRID_CSX_REF_DIR=${2%/}
					shift 2
					;;
				--csx-snpinfo)
					_grid_need_value "$@"
					GRID_CSX_SNPINFO=$2
					GRID_SNPINFO_EXPLICIT=TRUE
					shift 2
					;;
				--csx-bim-prefix | --bim-prefix)
					_grid_need_value "$@"
					GRID_CSX_BIM_PREFIX=$2
					GRID_BIM_EXPLICIT=TRUE
					shift 2
					;;
				--phi)
					_grid_need_value "$@"
					GRID_PHI=$2
					shift 2
					;;
				--mcmc-iter)
					_grid_need_value "$@"
					GRID_MCMC_ITER=$2
					shift 2
					;;
				--mcmc-burnin)
					_grid_need_value "$@"
					GRID_MCMC_BURNIN=$2
					shift 2
					;;
				--mcmc-thin)
					_grid_need_value "$@"
					GRID_MCMC_THIN=$2
					shift 2
					;;
				--seed)
					_grid_need_value "$@"
					GRID_SEED=$2
					shift 2
					;;
				--sumstats-chunk)
					_grid_need_value "$@"
					GRID_SUMSTATS_CHUNK=$2
					shift 2
					;;

				--med-file)
					_grid_need_value "$@"
					GRID_MED_FILE=$2
					shift 2
					;;
				--pca-file)
					_grid_need_value "$@"
					GRID_PCA_FILE=$2
					shift 2
					;;
				--pca-weight)
					_grid_need_value "$@"
					GRID_PCA_WEIGHT=$2
					shift 2
					;;
				--cov-pcs)
					_grid_need_value "$@"
					GRID_COV_PCS=$2
					shift 2
					;;
				--distance-pcs)
					_grid_need_value "$@"
					GRID_DISTANCE_PCS=$2
					shift 2
					;;
				--ancestry-file | --auto-ancestry-file)
					_grid_need_value "$@"
					GRID_ANCESTRY_FILE=$2
					shift 2
					;;
				--ancestry-prob-min | --ancestry-prob-threshold)
					_grid_need_value "$@"
					GRID_ANCESTRY_PROB_MIN=$2
					shift 2
					;;
				--anchor-prob-min)
					_grid_need_value "$@"
					GRID_ANCHOR_PROB_MIN=$2
					shift 2
					;;
				--anchor-max-per-group)
					_grid_need_value "$@"
					GRID_ANCHOR_MAX_PER_GROUP=$2
					shift 2
					;;
				--ethnicity-col)
					_grid_need_value "$@"
					GRID_ETHNICITY_COL=$2
					shift 2
					;;
				--phenotype-col)
					_grid_need_value "$@"
					GRID_PHENOTYPE_COL=$2
					shift 2
					;;
				--eval-covariates)
					_grid_need_value "$@"
					GRID_EVAL_COVARIATES=$2
					shift 2
					;;

				--arg-action)
					_grid_need_value "$@"
					GRID_ARG_ACTION=${2,,}
					shift 2
					;;
				--arg-hap-dir | --arg-hap-root | --ukb-hap-root)
					_grid_need_value "$@"
					GRID_ARG_HAP_DIR=${2%/}
					if [[ $GRID_ARG_OUT_AUTO == TRUE ]]; then
						GRID_ARG_OUT="$(dirname -- "$GRID_ARG_HAP_DIR")/arg.needle"
						GRID_ARG_TREES_DIR="$GRID_ARG_OUT/trees"
					fi
					shift 2
					;;
				--arg-map-dir)
					_grid_need_value "$@"
					GRID_ARG_MAP_DIR=${2%/}
					shift 2
					;;
				--arg-map-pattern)
					_grid_need_value "$@"
					GRID_ARG_MAP_PATTERN=$2
					shift 2
					;;
				--arg-out | --arg-dir)
					_grid_need_value "$@"
					GRID_ARG_OUT=${2%/}
					GRID_ARG_OUT_AUTO=FALSE
					GRID_ARG_TREES_DIR="$GRID_ARG_OUT/trees"
					shift 2
					;;
				--arg-trees-dir)
					_grid_need_value "$@"
					GRID_ARG_TREES_DIR=${2%/}
					shift 2
					;;
				--arg-scratch)
					_grid_need_value "$@"
					GRID_ARG_SCRATCH=${2%/}
					shift 2
					;;
				--arg-full)
					_grid_need_value "$@"
					GRID_ARG_FULL=$(grid_require_bool "$1" "$2")
					shift 2
					;;
				--arg-max-individuals | --arg-max-samples)
					_grid_need_value "$@"
					GRID_ARG_MAX_INDIVIDUALS=$2
					shift 2
					;;
				--arg-seed-haplotypes)
					_grid_need_value "$@"
					GRID_ARG_SEED_HAPLOTYPES=$2
					shift 2
					;;
				--arg-threads)
					_grid_need_value "$@"
					GRID_ARG_THREADS=$2
					shift 2
					;;
				--arg-jobs)
					_grid_need_value "$@"
					GRID_ARG_JOBS=$2
					shift 2
					;;
				--arg-anchors-per-pop)
					_grid_need_value "$@"
					GRID_ARG_ANCHORS_PER_POP=$2
					shift 2
					;;
				--arg-window-bp)
					_grid_need_value "$@"
					GRID_ARG_WINDOW_BP=$2
					shift 2
					;;
				--arg-needle-home)
					_grid_need_value "$@"
					GRID_ARG_NEEDLE_HOME=${2%/}
					shift 2
					;;
				--arg-normalize)
					_grid_need_value "$@"
					GRID_ARG_NORMALIZE=$(grid_require_bool "$1" "$2")
					shift 2
					;;
				--arg-maf-min)
					_grid_need_value "$@"
					GRID_ARG_MAF_MIN=$2
					shift 2
					;;
				--arg-geno-max)
					_grid_need_value "$@"
					GRID_ARG_GENO_MAX=$2
					shift 2
					;;
				--arg-keep)
					_grid_need_value "$@"
					GRID_ARG_KEEP=$2
					shift 2
					;;
				--arg-affinity)
					_grid_need_value "$@"
					GRID_ARG_AFFINITY=$(grid_require_bool "$1" "$2")
					shift 2
					;;

				--ldscore-dir)
					_grid_need_value "$@"
					GRID_LDSCORE_DIR=${2%/}
					shift 2
					;;
				--require-ld)
					_grid_need_value "$@"
					GRID_REQUIRE_LD=$(grid_require_bool "$1" "$2")
					shift 2
					;;
				--external-age)
					_grid_need_value "$@"
					GRID_EXTERNAL_AGE=$2
					shift 2
					;;
				--grid-max-snps-per-chr)
					_grid_need_value "$@"
					GRID_MAX_SNPS_PER_CHR=$2
					shift 2
					;;
				--grid-ridge-alpha)
					_grid_need_value "$@"
					GRID_RIDGE_ALPHA=$2
					shift 2
					;;
				--grid-transport-model)
					_grid_need_value "$@"
					GRID_TRANSPORT_MODEL=$2
					shift 2
					;;
				--conservation-min)
					_grid_need_value "$@"
					GRID_CONSERVATION_MIN=$2
					shift 2
					;;
				--conservation-max)
					_grid_need_value "$@"
					GRID_CONSERVATION_MAX=$2
					shift 2
					;;
				--grid-local)
					_grid_need_value "$@"
					GRID_LOCAL=$(grid_require_bool "$1" "$2")
					shift 2
					;;

				--eval-repeats)
					_grid_need_value "$@"
					GRID_EVAL_REPEATS=$2
					shift 2
					;;
				--eval-folds)
					_grid_need_value "$@"
					GRID_EVAL_FOLDS=$2
					shift 2
					;;
				--min-group)
					_grid_need_value "$@"
					GRID_MIN_GROUP=$2
					shift 2
					;;
				--eval-tuned)
					_grid_need_value "$@"
					GRID_EVAL_TUNED=$(grid_require_bool "$1" "$2")
					shift 2
					;;
				--eval-zero-shot)
					_grid_need_value "$@"
					GRID_EVAL_ZERO_SHOT=$(grid_require_bool "$1" "$2")
					shift 2
					;;
				--write-predictions)
					_grid_need_value "$@"
					GRID_WRITE_PREDICTIONS=$(grid_require_bool "$1" "$2")
					shift 2
					;;

				-h | --help)
					bash "$GRID_ROOT/grid.sh" --help
					exit 0
					;;
				*) _grid_die "unknown option '$1'" ;;
			esac
		done
		GRID_OUTPUT_ROOT=${GRID_OUTPUT_ROOT%/}
		[[ $GRID_BIM_EXPLICIT == TRUE ]] || GRID_CSX_BIM_PREFIX="$GRID_TARGET_DIR/ukb_array"
		[[ $GRID_SNPINFO_EXPLICIT == TRUE ]] || GRID_CSX_SNPINFO="$GRID_CSX_REF_DIR/snpinfo_mult_1kg_hm3"
		[[ -n $GRID_LDSCORE_DIR ]] || GRID_LDSCORE_DIR="$GRID_OUTPUT_ROOT/reference/ldscore"
		[[ -n $GRID_PCA_FILE ]] || GRID_PCA_FILE="$GRID_DATA_ROOT/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz"
		[[ -n $GRID_ANCESTRY_FILE ]] || GRID_ANCESTRY_FILE="$(dirname -- "$GRID_PCA_FILE")/ukb.ancestry.auto.tsv.gz"
		[[ $GRID_JOBS =~ ^[1-9][0-9]*$ ]] || _grid_die "--jobs must be a positive integer"
		[[ $GRID_THREADS =~ ^[1-9][0-9]*$ ]] || _grid_die "--threads must be a positive integer"
		[[ $GRID_ARG_THREADS =~ ^[1-9][0-9]*$ ]] || _grid_die "--arg-threads must be a positive integer"
		[[ $GRID_ARG_JOBS =~ ^[1-9][0-9]*$ ]] || _grid_die "--arg-jobs must be a positive integer"
	}

	grid_csv_words() { printf '%s' "$1" | tr ',;' '  '; }
	grid_expand_chrs() {
		python3 - "$1" <<'PY'
import re,sys
s=sys.argv[1].replace(',', ' ').replace(';',' ')
out=[]
for tok in s.split():
    m=re.fullmatch(r'(\d+)-(\d+)',tok)
    if m:
        a,b=map(int,m.groups()); out.extend(map(str,range(a,b+1)))
    else: out.append(tok.replace('chr','').replace('CHR',''))
seen=[]
for x in out:
    if x not in seen: seen.append(x)
print(' '.join(seen))
PY
	}

	grid_find_gwas() {
		local trait=${1,,} pop=${2^^} f tag
		tag="$trait.$pop"
		f="$GRID_GWAS_DIR/$tag/gwas/$tag.gz"
		if [[ -s $f ]]; then
			printf '%s\n' "$f"
			return
		fi
		# AFA is the source study's African-American label; use the AFR LD panel.
		if [[ $trait == t2dm && $pop == AFR ]]; then
			f="$GRID_GWAS_DIR/t2dm.AFA/gwas/t2dm.AFA.gz"
			if [[ -s $f ]]; then
				printf '%s\n' "$f"
				return
			fi
		fi
		f=$(find -L "$GRID_GWAS_DIR" -maxdepth 1 -type f -iname "${trait}.${pop}.gz" -print -quit 2>/dev/null || true)
		[[ -n $f ]] || return 1
		printf '%s\n' "$f"
	}

	grid_bgen_for() { find "$GRID_ARG_HAP_DIR" -maxdepth 1 -type f -name "ukb22438_c${1}_b0_v2.bgen" -print -quit; }
	grid_pfile_for() {
		local p="$GRID_ARG_HAP_DIR/chr${1}"
		[[ -s $p.pgen && (-s $p.pvar || -s $p.pvar.zst) && -s $p.psam ]] && printf '%s\n' "$p"
	}
	grid_sample_for() {
		local p
		p=$(find "$GRID_ARG_HAP_DIR" -maxdepth 1 -type f -name "ukb22438_c${1}_b0_v2_s*.sample" -print -quit)
		if [[ -n $p ]]; then
			printf '%s\n' "$p"
		elif [[ -s $GRID_ARG_HAP_DIR/chr${1}.psam ]]; then
			printf '%s\n' "$GRID_ARG_HAP_DIR/chr${1}.psam"
		fi
	}

	grid_map_for() {
		local c=$1 p
		if [[ -n $GRID_ARG_MAP_PATTERN ]]; then
			p=${GRID_ARG_MAP_PATTERN//\{chr\}/$c}
			p=${p//%CHR%/$c}
			[[ -s $p ]] && {
				printf '%s\n' "$p"
				return
			}
		fi
		[[ -n $GRID_ARG_MAP_DIR ]] || return 1
		for p in \
			"$GRID_ARG_MAP_DIR/chr${c}.map" "$GRID_ARG_MAP_DIR/chr${c}.txt" "$GRID_ARG_MAP_DIR/chr${c}.map.gz" \
			"$GRID_ARG_MAP_DIR/chr${c}.b37.gmap.gz" \
			"$GRID_ARG_MAP_DIR/genetic_map_GRCh37_chr${c}.txt" "$GRID_ARG_MAP_DIR/genetic_map_GRCh37_chr${c}.txt.gz" \
			"$GRID_ARG_MAP_DIR/genetic_map_chr${c}_combined_b37.txt" "$GRID_ARG_MAP_DIR/genetic_map_chr${c}_combined_b37.txt.gz" \
			"$GRID_ARG_MAP_DIR/genetic_map_chr${c}_combined_b37.txt.map"; do
			[[ -s $p ]] && {
				printf '%s\n' "$p"
				return
			}
		done
		return 1
	}

	grid_target_mode() {
		local c=$1 d=$GRID_TARGET_DIR
		if [[ -s $d/chr$c.pgen && (-s $d/chr$c.pvar || -s $d/chr$c.pvar.zst) && -s $d/chr$c.psam ]]; then
			echo pfile
		elif [[ -s $d/chr$c.bed && -s $d/chr$c.bim && -s $d/chr$c.fam ]]; then
			echo bfile
		else return 1; fi
	}

	grid_run() {
		printf '[GRID]'
		printf ' %q' "$@"
		printf '\n' >&2
		grid_is_true "$GRID_DRY_RUN" || "$@"
	}

	grid_run_logged() {
		local log=$1
		shift
		mkdir -p "$(dirname "$log")"
		printf '[GRID]'
		printf ' %q' "$@"
		printf '\n' | tee -a "$log" >&2
		grid_is_true "$GRID_DRY_RUN" || "$@" >>"$log" 2>&1
	}

	grid_checksum_record() {
		local out=$1
		shift
		mkdir -p "$(dirname "$out")"
		{
			date -Is
			for f in "$@"; do [[ -e $f ]] && stat -c '%n\t%s\t%Y' "$f"; done
		} >"$out"
	}
}

grid_enable_logging() {
	# Record exact tool invocations; verbose tool output stays in per-step logs.
	grid_run_logged() {
		local logfile=$1
		shift
		mkdir -p "$(dirname -- "$logfile")"
		(
			flock 9
			printf '%q ' "$@" >&9
			printf ' > %q 2>&1\n' "$logfile" >&9
		) 9>>"$GRID_COMMAND_FILE"
		[[ $GRID_DRY_RUN == FALSE ]] || return 0
		if "$@" >"$logfile" 2>&1; then return 0; else
			local rc=$?
			echo "ERROR: command failed (exit $rc): $logfile" >&2
			tail -n 12 "$logfile" >&2
			return "$rc"
		fi
	}
	grid_run() { grid_run_logged "$logdir/step.$(date +%s%N).$BASHPID.log" "$@"; }
	need() { [[ -s $1 ]] || _grid_die "Missing/empty file: $1"; }
	join_comma() {
		local IFS=,
		echo "$*"
	}
	# Publish verified results; staging and rollback copies stay in /tmp.
	publish() {
		local src=$1 dst=$2 previous
		need "$src"
		mkdir -p "$(dirname -- "$dst")"
		previous=$(mktemp /tmp/grid-publish.XXXXXX)
		if [[ -f $dst ]]; then cp -p -- "$dst" "$previous"; fi
		if ! cp -- "$src" "$dst" || ! cmp -s -- "$src" "$dst"; then
			if [[ -s $previous ]]; then cp -p -- "$previous" "$dst"; else rm -f -- "$dst"; fi
			rm -f -- "$previous"
			return 1
		fi
		rm -f -- "$previous"
	}
	# A subset run must never overwrite a genome-wide result.
	suffix=''
	if [[ ${CHRS[*]} != "$(seq -s ' ' 1 22)" ]]; then suffix=".chr$(join_comma "${CHRS[@]}")"; fi
}


# 🚩 grid_pipeline
grid_pipeline() (
	# Shared launcher for PCA, PRS-CSx and official DiscoDivas.
	set -euo pipefail
	ROOT=$GRID_COMMON_ROOT
	method=$1
	shift
	grid_activate_environment
	grid_configure "$@"
	GRID_CHECK=FALSE
	GRID_STAGE=all
	GRID_SCORE_DIR=${GRID_SCORE_DIR:-/mnt/d/data/ukb/pgs}
	GRID_DISCO_A=1,1,1,1
	GRID_REGRESS_PCA=TRUE
	GRID_CSX_MODELS=populations,auto,meta
	GRID_POSTERIOR=TRUE
	GRID_POSTERIOR_FREQ_DIR=
	GRID_POSTERIOR_MEMORY=${GRID_POSTERIOR_MEMORY:-8192}
	GRID_SCORE_MEMORY=${GRID_SCORE_MEMORY:-2048}
	args=()
	while (($#)); do
		case "$1" in
			--check)
				GRID_CHECK=TRUE
				shift
				;;
			--posterior)
				_grid_need_value "$@"
				GRID_POSTERIOR=$(grid_require_bool "$1" "$2")
				shift 2
				;;
			--posterior-frequency-dir)
				_grid_need_value "$@"
				GRID_POSTERIOR_FREQ_DIR=$2
				shift 2
				;;
			--posterior-memory)
				_grid_need_value "$@"
				GRID_POSTERIOR_MEMORY=$2
				shift 2
				;;
			--score-memory)
				_grid_need_value "$@"
				GRID_SCORE_MEMORY=$2
				shift 2
				;;
			--models)
				_grid_need_value "$@"
				GRID_CSX_MODELS=$2
				shift 2
				;;
			--stage)
				_grid_need_value "$@"
				GRID_STAGE=$2
				shift 2
				;;
			--score-dir)
				_grid_need_value "$@"
				GRID_SCORE_DIR=${2%/}
				shift 2
				;;
			--a-list)
				_grid_need_value "$@"
				GRID_DISCO_A=$2
				shift 2
				;;
			--regress-pca)
				_grid_need_value "$@"
				GRID_REGRESS_PCA=$(grid_require_bool "$1" "$2")
				shift 2
				;;
			--trait)
				_grid_need_value "$@"
				GRID_TRAITS=${2,,}
				shift 2
				;;
			--traits)
				_grid_need_value "$@"
				GRID_TRAITS=${2,,}
				shift 2
				;;
			*)
				args+=("$1")
				shift
				;;
		esac
	done
	grid_parse_args "${args[@]}"
	if [[ $method == csx ]]; then
		for mem in "$GRID_SCORE_MEMORY" "$GRID_POSTERIOR_MEMORY"; do
			[[ $mem =~ ^[1-9][0-9]{0,6}$ ]] && ((mem >= 640)) || _grid_die 'PLINK memory must be an integer >= 640 MiB'
		done
		# Leave 25% of the task cap for Python, PLINK overhead and file cache.
		# The cgroup is the hard aggregate limit; --memory only sizes PLINK's workspace.
		GRID_SCORE_JOBS=$GRID_JOBS
		if [[ -n ${GRID_MEMORY_CAP_GB:-} && $GRID_STAGE != weights ]]; then
			score_budget=$((GRID_MEMORY_CAP_GB * 1024 * 3 / 4))
			((GRID_SCORE_MEMORY <= score_budget)) || _grid_die '--score-memory exceeds 75% of the total memory cap'
			if [[ $GRID_POSTERIOR == TRUE && ,$GRID_CSX_MODELS, == *,populations,* ]]; then
				((GRID_POSTERIOR_MEMORY <= score_budget)) || _grid_die '--posterior-memory exceeds 75% of the total memory cap'
			fi
			max_score_jobs=$((score_budget / GRID_SCORE_MEMORY))
			if ((GRID_SCORE_JOBS > max_score_jobs)); then GRID_SCORE_JOBS=$max_score_jobs; fi
		fi
		echo "[csx] Scoring: $GRID_SCORE_JOBS parallel chromosomes, $GRID_SCORE_MEMORY MiB PLINK workspace each; inference jobs=$GRID_JOBS"
		python3 -c 'import numpy,pandas,scipy,h5py' || _grid_die 'Missing PRS-CSx Python dependencies; update the grid environment'
		[[ -n $GRID_CSX_MODELS ]] || _grid_die '--models cannot be empty'
		for model in $(grid_csv_words "$GRID_CSX_MODELS"); do
			[[ $model == populations || $model == auto || $model == meta ]] || _grid_die '--models accepts populations,auto,meta'
			[[ $model != meta || $GRID_PHI != auto ]] || _grid_die 'Use --models populations,auto with --phi auto; meta requires a fixed phi'
		done
	fi
	[[ $GRID_STAGE == all || $GRID_STAGE == weights || $GRID_STAGE == score ]] || _grid_die '--stage must be all, weights or score'
	[[ $method == csx || $GRID_STAGE == all ]] || _grid_die '--stage applies only to CSx'
	# Validate before creating any output or launching computation.
	[[ $GRID_POPS == AFR,EAS,EUR,SAS ]] || _grid_die 'This workflow requires --pops AFR,EAS,EUR,SAS (order matches the supplied reference centers)'
	read -r -a CHRS <<<"$(grid_expand_chrs "$GRID_CHRS")"
	((${#CHRS[@]})) || _grid_die 'Empty chromosome list'
	for c in "${CHRS[@]}"; do [[ $c =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || _grid_die "Invalid autosome: $c"; done
	[[ $method != pca || ${CHRS[*]} == "$(seq -s ' ' 1 22)" ]] || _grid_die 'PCA requires all 22 autosomes'
	POPS=(AFR EAS EUR SAS)
	for t in $(grid_csv_words "$GRID_TRAITS"); do [[ $t == height || $t == ldl || $t == t2dm ]] || _grid_die "Unknown trait: $t"; done
	[[ -n $GRID_TRAITS ]] || _grid_die 'Empty trait list'
	work="$GRID_OUTPUT_ROOT/$method"
	if [[ $method == pca || $method == disco ]]; then
		work=$(python3 "$ROOT/f/0.common.py" cache-path "$work")
	fi
	logdir=$(python3 "$ROOT/f/0.common.py" cache-path "$GRID_OUTPUT_ROOT/$method/log")
	cmddir=$(python3 "$ROOT/f/0.common.py" cache-path "$GRID_OUTPUT_ROOT/$method/cmd")
	mkdir -p "$work" "$cmddir" "$logdir"
	run_id=$(date +%Y%m%dT%H%M%S).$$
	GRID_COMMAND_FILE="$cmddir/$method.$run_id.sh"
	GRID_RUN_LOG="$logdir/$method.$run_id.log"
	printf '#!/usr/bin/env bash\nset -euo pipefail\nsource %q\ngrid_activate_environment\n' "$ROOT/f/0.common.sh" >"$GRID_COMMAND_FILE"
	exec > >(tee -a "$GRID_RUN_LOG") 2>&1
	trap 'rc=$?; if ((rc)); then echo "ERROR: exit=$rc; details: $GRID_RUN_LOG"; fi' EXIT
	grid_enable_logging
	echo "[$method] work directory: $work"
	echo "[$method] command file: $GRID_COMMAND_FILE"
	echo "[$method] log file: $GRID_RUN_LOG"
	if [[ $method == pca ]]; then
		pca_run
	else
		for trait in $(grid_csv_words "$GRID_TRAITS"); do
			if [[ $method == csx ]]; then
				for csx_model in $(grid_csv_words "$GRID_CSX_MODELS"); do
					case "$csx_model" in
						populations) (
							GRID_TRAIT=$trait
							csx_population_run
						) ;;
						auto | meta) (
							GRID_TRAIT=$trait
							csx_combined_run
						) ;;
						*) _grid_die '--models accepts populations,auto,meta' ;;
					esac
				done
			else
				(
					GRID_TRAIT=$trait
					disco_run
				)
			fi
		done
	fi
)


# 🚩 grid_preflight
grid_preflight() (
	ROOT=$GRID_COMMON_ROOT
	set -euo pipefail
	# --scope belongs only to preflight; remove it before the common parser.
	scope=all
	args=()
	while (($#)); do
		case "$1" in
			--scope)
				[[ $# -ge 2 ]] || {
					echo 'ERROR: --scope requires core|arg|grid|all' >&2
					exit 2
				}
				scope=${2,,}
				shift 2
				;;
			*)
				args+=("$1")
				shift
				;;
		esac
	done
	# shellcheck source=f/0.common.sh
	grid_activate_environment
	grid_configure "${args[@]}"
	grid_parse_args "${args[@]}"
	case "$scope" in core | arg | grid | all) ;; *)
		echo "ERROR: bad --scope=$scope" >&2
		exit 2
		;;
	esac

	report="$(python3 "$ROOT/f/0.common.py" cache-path "$GRID_OUTPUT_ROOT/preflight")/${scope}.tsv"
	mkdir -p "$(dirname "$report")"
	mkdir -p "$(dirname "$report")"
	printf 'scope\tstatus\titem\tdetail\n' >"$report"
	fail=0
	record() {
		local st=$1 item=$2 detail=${3:-}
		printf '%s\t%s\t%s\t%s\n' "$scope" "$st" "$item" "$detail" | tee -a "$report"
		[[ $st != FAIL ]] || fail=1
	}
	check_cmd() { if command -v "$1" >/dev/null 2>&1; then record PASS "command:$1" "$(command -v "$1")"; else record FAIL "command:$1" missing; fi; }
	check_file() { if [[ -s $1 ]]; then record PASS "$2" "$1"; else record FAIL "$2" "$1"; fi; }
	check_dir() { if [[ -d $1 ]]; then record PASS "$2" "$1"; else record FAIL "$2" "$1"; fi; }
	check_py() {
		local mod=$1
		if python3 - "$mod" <<'PY' >/dev/null 2>&1
import importlib,sys
importlib.import_module(sys.argv[1])
PY
		then record PASS "python:$mod" importable; else record FAIL "python:$mod" not_importable; fi
	}
	check_r() {
		local pkg=$1
		if Rscript -e "quit(status=ifelse(requireNamespace('$pkg',quietly=TRUE),0,1))" >/dev/null 2>&1; then record PASS "R:$pkg" installed; else record FAIL "R:$pkg" missing; fi
	}

	if [[ $scope == core || $scope == all ]]; then
		for x in bash awk sed grep sort gzip python3 Rscript plink2; do check_cmd "$x"; done
		for pvar in "$GRID_TARGET_DIR"/*.pvar.zst; do
			[[ -s "$pvar" && ! -s "${pvar%.zst}" ]] || continue
			check_cmd zstdcat
			break
		done
		for m in numpy pandas scipy h5py; do check_py "$m"; done
		for p in data.table ggplot2 patchwork openxlsx; do check_r "$p"; done
		check_file "$GRID_PHE_FILE" phenotype_rds
		check_file "$GRID_PCA_WEIGHT" DiscoDivas_pca_weight
		check_file "$GRID_MED_FILE" DiscoDivas_reference_centers
		check_dir "$GRID_CSX_REF_DIR" PRSCSx_reference_root
		check_file "$GRID_CSX_SNPINFO" PRSCSx_snpinfo
		check_file "$GRID_CSX_BIM_PREFIX.bim" target_validation_bim
		prscx=$(find "$GRID_ROOT/f/csx" -maxdepth 3 -type f \( -iname 'PRScsx.py' -o -iname 'prscsx.py' \) -print -quit 2>/dev/null || true)
		[[ -n $prscx ]] && record PASS PRSCSx_program "$prscx" || record FAIL PRSCSx_program "$GRID_ROOT/f/csx"
		for c in $(grid_expand_chrs "$GRID_CHRS" | awk '{print $1}'); do
			if grid_target_mode "$c" >/dev/null 2>&1; then record PASS "target_genotype_chr$c" "$GRID_TARGET_DIR"; else record FAIL "target_genotype_chr$c" "$GRID_TARGET_DIR"; fi
		done
		for t in $(grid_csv_words "$GRID_TRAITS"); do
			for p in $(grid_csv_words "$GRID_POPS"); do
				if f=$(grid_find_gwas "$t" "$p"); then record PASS "gwas:${t}.${p}" "$f"; else record FAIL "gwas:${t}.${p}" "$GRID_GWAS_DIR/${t}.${p}.gz"; fi
			done
		done
	fi

	if [[ $scope == arg || $scope == all ]]; then
		for x in gzip python3; do check_cmd "$x"; done
		check_py tskit
		for c in $(grid_expand_chrs "$GRID_CHRS"); do
			check_file "$GRID_ARG_OUT/argn/chr$c.argn" "ARG_Needle_argn_chr$c"
			check_file "$GRID_ARG_TREES_DIR/chr$c.trees" "ARG_Needle_trees_chr$c"
			check_file "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" "ARG_Needle_sample_map_chr$c"
			check_file "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" "ARG_Needle_anchors_chr$c"
			check_file "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz" "ARG_Needle_features_chr$c"
		done
		if ((fail)); then
			record WARN ARG_build_command "bash /mnt/d/scripts/gu/arg.sh build --method needle --dir-gen $(dirname -- "$GRID_ARG_HAP_DIR") --chr $GRID_CHRS"
		fi
	fi

	if [[ $scope == grid || $scope == all ]]; then
		for m in numpy pandas scipy h5py tskit; do check_py "$m"; done
		check_file "$GRID_MED_FILE" global_PCA_centers
		# The LD score files are generated from the already downloaded PRS-CSx HDF5 panels.
		for p in $(grid_csv_words "$GRID_POPS"); do
			low=${p,,}
			d=$(find "$GRID_CSX_REF_DIR" -maxdepth 1 -type d -iname "ldblk*${low}*" -print -quit 2>/dev/null || true)
			[[ -n $d ]] && record PASS "PRSCSx_HDF5_LD:$p" "$d" || record FAIL "PRSCSx_HDF5_LD:$p" "$GRID_CSX_REF_DIR"
		done
		for c in $(grid_expand_chrs "$GRID_CHRS"); do
			f="$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
			[[ -s $f ]] && record PASS "ARG_features_chr$c" "$f" || record WARN "ARG_features_chr$c" 'created by arg.sh build --method needle'
		done
	fi

	if ((fail)); then
		echo "PREFLIGHT FAILED: $report" >&2
		exit 1
	fi
	echo "PREFLIGHT PASS (warnings may remain): $report"
)

GRID_COMMON_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
	case "${1:-help}" in
		-h | --help | help) echo 'Usage: bash f/0.common.sh preflight [--scope core|arg|grid|all] [shared options]' ;;
		preflight)
			shift
			grid_preflight "$@"
			;;
		*)
			echo "Unknown common command: $1" >&2
			exit 2
			;;
	esac
fi
