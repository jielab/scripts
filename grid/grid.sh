#!/usr/bin/env bash
# GRID: GPU individual-reference matching and fixed held-out evaluation.
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1


# 🚩 Command-line help
usage() {
	cat <<'HELP'
GRID — Genetic Risk based on Individual Distance

Usage:
  cd /mnt/d/scripts/grid
  bash grid.sh --traits height,ldl,t2dm --check
  bash grid.sh --traits height,ldl,t2dm
  bash grid.sh --check-device
  bash grid.sh --stage fit --trait height
  bash grid.sh --stage predict --trait height --model-file MODEL --data-file NEW_PEOPLE
  bash grid.sh all --traits height,ldl,t2dm --run-prsformer \
    --prsformer-python ~/.venvs/grid-prsformer/bin/python --prsformer-device cuda:0

Stages (--stage):
  prepare    Align phenotype, CSx scores and PCs; freeze the common 50/50 split.
  evolution  Validate score lineage; build features only for explicit annotations.
  fit        Tune within training; freeze predictors, matching and screening rules.
  predict    Apply a frozen --model-file to new individuals in --data-file.
  report     Evaluate the frozen test half; publish figures, workbooks and explanations.
  all        Run prepare -> evolution -> fit -> report (default).
  check      Validate inputs/configuration without fitting or publishing.

Optional modules:
  pca        Forward to the existing 0.pca.sh; use 'grid.sh pca --help'.
  csx        Forward to the existing 1.csx.sh.
  disco      Forward to the existing 2.disco.sh.
  prsformer  Use the common cohort and test set prepared by GRID.
  all        Run all GRID stages, optionally also the PRSformer benchmark.
  check      Same as --stage check.

Common data options:
  --trait NAME / --traits LIST  height,ldl,t2dm.
  --pheno-file FILE            /mnt/d/data/ukb/phe/Rdata/all.rds.
  --score-dir DIR              /mnt/d/data/ukb/pgs.
  --pca-file FILE              /mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz.
  --ancestry-file FILE         /mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz.
  --dir-gen DIR               /mnt/f/gen/ukb/37/hap.
  --gwas-dir DIR              /mnt/f/gwas/4grid/common.
  --snpinfo FILE              /mnt/f/refLD/csx/snpinfo_mult_1kg_hm3.
  --weights-file FILE         Optional canonical weights; {trait} is expanded by Python.
  --data-file FILE            Optional prepared eid,split,y,csx.* input table.
  --split-file FILE           Optional fixed eid,split=train|test roster.
  --split-group-file FILE     Optional eid,family_id (or group) to separate relatives.
  --annotation-file FILE      Optional canonical annotations; default uses CSx + PCs.
  --age-permutation-mode MODE auto (labeled), chr-maf (strict), or chr-only.
  --min-age-variants N         Minimum nonzero-weight, confidently assigned ages.
  --build NAME                GRCh37.
  --allow-proxy-only          For explicit annotations: permit inadequate age coverage.
  --out-root DIR              /mnt/d/analysis/grid/GRID.
  --cache-dir DIR             /tmp/grid-cache/grid; sensitive temporary data.

Controls:
  --abm-backend METHOD         selective_attention (default) or reference (legacy).
  --device DEVICE              GRID attention: cuda (default); cuda:0 selects GPU 0.
  --retrieval METHOD           cuda (default), kd_tree or hnsw.
  --check-device               Test real GRID GPU forward/backward without cohort IO.
  --attention-epochs N         Maximum neural training epochs (20); tune_model selects.
  --python PATH               GRID_PYTHON, then python3 from the automatically
                              activated grid Conda environment. An explicit
                              interpreter uses its own environment without activation.
  --feature-manifest FILE     predict: saved feature-contract and input-file hashes.
  --check                     Same as --stage check.
  --dry-run [TRUE|FALSE]       Print commands only; do not read data or import models.
  --python-help               Show every method/data option implemented by f/grid.py.
  --run-prsformer [TRUE|FALSE] With stage all: prepare -> evolution -> PRSformer -> fit -> report.
  --prsformer-OPTION VALUE    Pass --OPTION VALUE to 3.prsformer.sh (e.g. --prsformer-epochs 10).
  --prsformer-arg ARG         Pass one exact PRSformer token; repeat for a flag and value.
  --prsformer-root DIR        Default OUTROOT/benchmark; separate from earlier runs.

Separate PRSformer benchmark:
  bash grid.sh --stage prepare --traits height,ldl,t2dm
  bash grid.sh prsformer --traits height,ldl,t2dm --prsformer-device cuda:0
  bash grid.sh --stage fit --traits height,ldl,t2dm
  bash grid.sh --stage report --traits height,ldl,t2dm

prepare writes common.keep, families.tsv.gz and prsformer.split.tsv.gz in CACHE.
The original test half remains test; the training half is divided internally into
80% fit and 20% validation (approximately 40/10/50 of the original cohort).
PRSformer results go to PRSFORMER_ROOT/prsformer and PRSFORMER_ROOT/scores/<trait>.
Use identical phenotype, covariates, traits and cohort options across stage commands.

Other options are forwarded unchanged to f/grid.py. CSx weights/scores are reused;
PCA, CSx and Disco are not rebuilt automatically. Without a family roster, random
splitting does not establish unrelatedness. Matched people are empirical reference
individuals, not biological twins or guarantees of individual prediction accuracy.
HELP
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
need_value() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || die "Missing value for $1"; }
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
MODULE=grid
if (($#)); then
	case $1 in
		-h | --help | help) usage; exit 0 ;;
		--*) ;;
		*) MODULE=$1; shift ;;
	esac
fi


# 🚩 Existing standalone methods
case $MODULE in
	pca) exec bash "$ROOT/0.pca.sh" "$@" ;;
	csx) exec bash "$ROOT/1.csx.sh" "$@" ;;
	disco) exec bash "$ROOT/2.disco.sh" "$@" ;;
	grid | all | check | prsformer) ;;
	*) die "Unknown module: $MODULE (see --help)" ;;
esac


# 🚩 Runtime, inputs and argument arrays
PYTHON=${GRID_PYTHON:-python3}
PYTHON_EXPLICIT=FALSE
[[ -z ${GRID_PYTHON:-} ]] || PYTHON_EXPLICIT=TRUE
CACHE_DIR=/tmp/grid-cache/grid
OUT_ROOT=/mnt/d/analysis/grid/GRID
PRSFORMER_ROOT=
STAGE=all
[[ $MODULE != check ]] || STAGE=check
TRAITS=height,ldl,t2dm
RUN_PRSFORMER=FALSE
DRY_RUN=FALSE
PYTHON_HELP=FALSE
RUN_LOG=
GRID_ARGS=()
PRS_SHARED=()
PRS_ARGS=()
while (($#)); do
	case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
	case $1 in
		-h | --help | help) usage; exit 0 ;;
		--python-help) PYTHON_HELP=TRUE; shift ;;
		--python) need_value "$@"; PYTHON=$2; PYTHON_EXPLICIT=TRUE; shift 2 ;;
		--check) STAGE=check; shift ;;
		--check-device) STAGE=check; GRID_ARGS+=(--check-device); shift ;;
		--stage) need_value "$@"; STAGE=$2; shift 2 ;;
		--cache-dir) need_value "$@"; CACHE_DIR=$2; shift 2 ;;
		--out-root) need_value "$@"; OUT_ROOT=$2; shift 2 ;;
		--trait | --traits) need_value "$@"; TRAITS=$2; shift 2 ;;
		--prsformer-root) need_value "$@"; PRSFORMER_ROOT=$2; shift 2 ;;
		--run-prsformer | --dry-run)
			flag=$1
			flag_value=TRUE
			shift
			if (($#)); then
				case ${1^^} in TRUE | FALSE) flag_value=${1^^}; shift ;; esac
			fi
			if [[ $flag == --run-prsformer ]]; then RUN_PRSFORMER=$flag_value; else DRY_RUN=$flag_value; fi
			;;
		--prsformer-arg)
			[[ $# -ge 2 && -n $2 ]] || die '--prsformer-arg requires one exact token'
			PRS_ARGS+=("$2")
			shift 2
			;;
		--prsformer-*)
			need_value "$@"
			PRS_ARGS+=("--${1#--prsformer-}" "$2")
			shift 2
			;;
		--pheno-file | --ancestry-file | --group-col | --dir-gen | --covariates | --height-col | --ldl-col | --t2dm-col | --chrs | --remove | --seed)
			need_value "$@"; GRID_ARGS+=("$1" "$2"); PRS_SHARED+=("$1" "$2"); shift 2 ;;
		--snpinfo)
			need_value "$@"; GRID_ARGS+=("$1" "$2"); PRS_SHARED+=(--snp-list "$2"); shift 2 ;;
		*) GRID_ARGS+=("$1"); shift ;;
	esac
done
case $STAGE in prepare | evolution | fit | predict | report | all | check) ;; *) die "Unknown --stage: $STAGE" ;; esac
[[ $MODULE != check || $STAGE == check ]] || die "The check module requires stage check"
[[ $MODULE != prsformer || $STAGE == all || $STAGE == check ]] || die "For PRSformer stages, use --prsformer-mode prepare/train/predict/report"
[[ -n $CACHE_DIR && $CACHE_DIR != / ]] || die 'Use a dedicated non-root --cache-dir'
[[ -n $OUT_ROOT && $OUT_ROOT != / ]] || die 'Use a dedicated non-root --out-root'
[[ $RUN_PRSFORMER == FALSE || $STAGE == all ]] || die '--run-prsformer requires --stage all'
[[ -n $PRSFORMER_ROOT ]] || PRSFORMER_ROOT="$OUT_ROOT/benchmark"

# Keep PLINK and other executables in the same environment as the selected Python.
# Preserve the supplied environment's bin directory even when Python is a symlink.
if [[ $DRY_RUN == FALSE && ( $MODULE != prsformer || $PYTHON_HELP == TRUE ) ]]; then
	if [[ $PYTHON_EXPLICIT == FALSE ]]; then
	CONDA_INIT="$HOME/miniforge3/etc/profile.d/conda.sh"
	[[ -f $CONDA_INIT ]] || CONDA_INIT="$HOME/anaconda3/etc/profile.d/conda.sh"
	[[ -f $CONDA_INIT ]] || die 'Conda initialization script unavailable; install Miniforge/Conda first'
	# Activation/deactivation hooks may read unset variables in another environment.
	set +u
	# shellcheck disable=SC1090
	if source "$CONDA_INIT" && conda activate grid; then
		CONDA_STATUS=0
	else
		CONDA_STATUS=$?
	fi
	set -u
	((CONDA_STATUS == 0)) || die 'Could not activate Conda environment grid; run install_grid.sh'
	fi

	PYTHON_COMMAND=$(command -v "$PYTHON") || die "Python unavailable: $PYTHON; set --python or run install_grid.sh"
	PYTHON_DIRECTORY=$(cd -- "$(dirname -- "$PYTHON_COMMAND")" && pwd -P)
	PYTHON="$PYTHON_DIRECTORY/$(basename -- "$PYTHON_COMMAND")"
	export PATH="$PYTHON_DIRECTORY:$PATH"
fi

run_command() {
	printf 'RUN '
	printf '%q ' "$@"
	printf '\n'
	[[ $DRY_RUN == FALSE ]] || return 0
	if [[ -n $RUN_LOG ]]; then
		# pipefail preserves the worker's failure status while retaining the
		# diagnostics of unattended/delayed runs after the terminal is closed.
		"$@" 2>&1 | tee -a "$RUN_LOG"
	else
		"$@"
	fi
}
if [[ $PYTHON_HELP == TRUE ]]; then
	run_command "$PYTHON" "$ROOT/f/grid.py" --help
	exit 0
fi


# 🚩 One shared cohort and one fixed outer test set
run_grid_stage() {
	run_command "$PYTHON" "$ROOT/f/grid.py" --stage "$1" \
		--cache-dir "$CACHE_DIR" --out-root "$OUT_ROOT" --traits "$TRAITS" \
		--prsformer-root "$PRSFORMER_ROOT" "${GRID_ARGS[@]}"
}
run_prsformer() {
	local path
	local check_args=()
	[[ $STAGE != check ]] || check_args=(--check)
	if [[ $DRY_RUN == FALSE ]]; then
		for path in common.keep families.tsv.gz prsformer.split.tsv.gz; do
			[[ -s $CACHE_DIR/$path ]] || die "Missing shared roster $CACHE_DIR/$path; run grid.sh --stage prepare with these cohort settings first"
		done
	fi
	# Common cohort/output arguments follow model options so an optional benchmark
	# cannot silently move test people to training or overwrite an earlier split.
	run_command bash "$ROOT/3.prsformer.sh" \
		--train-fraction 0.4 --validation-fraction 0.1 --cache-dir "$CACHE_DIR/prsformer" \
		"${PRS_ARGS[@]}" "${PRS_SHARED[@]}" "${check_args[@]}" --traits "$TRAITS" \
		--keep "$CACHE_DIR/common.keep" \
		--split-file "$CACHE_DIR/prsformer.split.tsv.gz" \
		--split-group-file "$CACHE_DIR/families.tsv.gz" \
		--output-root "$PRSFORMER_ROOT" --score-dir "$PRSFORMER_ROOT/scores"
}


# 🚩 Prevent concurrent stages from changing the shared split
if [[ $DRY_RUN == FALSE ]]; then
	if [[ $MODULE != prsformer ]]; then
		command -v "$PYTHON" >/dev/null 2>&1 || die "Python unavailable: $PYTHON; set --python or run install_grid.sh"
		[[ -f $ROOT/f/grid.py ]] || die "Missing GRID implementation: $ROOT/f/grid.py"
	fi
	if [[ $STAGE != check ]]; then
		command -v flock >/dev/null 2>&1 || die 'flock is required for concurrent-run protection'
		mkdir -p -- "$CACHE_DIR"
		exec {grid_lock}>"$CACHE_DIR/grid.lock"
		flock -n "$grid_lock" || die "Another GRID run is using $CACHE_DIR"
		mkdir -p -- "$CACHE_DIR/logs"
		RUN_LOG=$(mktemp "$CACHE_DIR/logs/grid.$(date -u +%Y%m%dT%H%M%SZ).XXXXXX.log")
		printf 'LOG %s\n' "$RUN_LOG"
	fi
fi
if [[ $MODULE == prsformer ]]; then
	run_prsformer
elif [[ $RUN_PRSFORMER == FALSE ]]; then
	run_grid_stage "$STAGE"
else
	run_grid_stage prepare
	run_grid_stage evolution
	run_prsformer
	run_grid_stage fit
	run_grid_stage report
fi
