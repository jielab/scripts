#!/usr/bin/env bash
# GRID: unified entry for benchmark methods, evaluation, and individual-matching ABM.
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
need_value() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || die "Missing value for $1"; }

usage() {
	cat <<'HELP'
GRID — Genetic Risk based on Individual Distance

Modules (only these three):
  prsformer  Reproduce PRSformer; automatically prepare its shared GRID roster.
  abm        Run the innovative ABM method: prepare, features, fit and publish.
  final      Evaluate available CSx, Disco, PRSformer and ABM results and plot performance.

Example usage (WSL; default traits: height,ldl,t2dm):
  cd /mnt/d/scripts/grid

  # Upstream: skip completed methods when their inputs/settings have not changed.
  bash 0.pca.sh
  bash 1.csx.sh
  bash 2.disco.sh

  # Optional: train PRSformer only when no completed result is available.
  bash grid.sh prsformer

  # Run ABM, then compare all available saved results.
  bash grid.sh abm,final

Controls:
  --traits LIST      Default height,ldl,t2dm.
  --check            Validate inputs without fitting.
  --dry-run          Print commands without running analysis.
  --help             Show module help.
  --python-help      Show advanced ABM options.

Notes:
  * ABM outputs: /mnt/d/analysis/grid/grid; final reports: /mnt/d/analysis/grid/final.
  * Comma-separated modules run in the listed order and stop on the first error.
    Shared sequence options: --traits, --check, --dry-run. For method-specific
    options, run the modules separately.
  * ABM runs its complete workflow with no --stage option needed.
  * final uses saved PRSformer/ABM predictions; it never retrains those models.
    Learned methods are compared on their common held-out test participants.
    Baseline-only evaluation retains the original Yeval cross-validation workflow.
  * Skip prsformer when the desired completed benchmark already exists.
  * CSx cache reuse still republishes scores, which may trigger Disco recalculation.
  * ABM refits and overwrites its outputs by default; --no-replace changes this.

HELP
}

MODULE=abm
if (($#)); then
	case $1 in
		-h|--help|help) usage; exit 0 ;;
		--*) ;;
		*) MODULE=$1; shift ;;
	esac
fi

# Validate the complete sequence before starting any expensive model run.
if [[ $MODULE == *,* ]]; then
	[[ $MODULE =~ ^(prsformer|abm|final)(,(prsformer|abm|final))+$ ]] || die "Invalid module sequence: $MODULE (use prsformer, abm, final)"
	IFS=, read -r -a MODULES <<<"$MODULE"
	declare -A SEEN_MODULES=()
	for module in "${MODULES[@]}"; do
		[[ -z ${SEEN_MODULES[$module]+present} ]] || die "Repeated module: $module"
		SEEN_MODULES[$module]=1
	done
	SEQUENCE_ARGS=()
	SEQUENCE_DRY_RUN=FALSE
	while (($#)); do
		case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
		case $1 in
			-h|--help|help) usage; exit 0 ;;
			--trait|--traits) need_value "$@"; SEQUENCE_ARGS+=(--traits "$2"); shift 2 ;;
			--check) SEQUENCE_ARGS+=(--check); shift ;;
			--dry-run)
				shift
				SEQUENCE_DRY_RUN=TRUE
				if (($#)); then
					case ${1^^} in TRUE|FALSE) SEQUENCE_DRY_RUN=${1^^}; shift ;; esac
				fi
				;;
			*) die "Option $1 is not shared by a module sequence; run modules separately for method-specific options" ;;
		esac
	done
	[[ $SEQUENCE_DRY_RUN == FALSE ]] || SEQUENCE_ARGS+=(--dry-run)
	for module in "${MODULES[@]}"; do
		printf '\nMODULE %s\n' "$module"
		bash "$ROOT/grid.sh" "$module" "${SEQUENCE_ARGS[@]}" || exit "$?"
	done
	exit 0
fi

case "$MODULE" in
	-h|--help|help) usage; exit 0 ;;
	final)
		[[ -f "$ROOT/f/final.sh" ]] || die "Missing $ROOT/f/final.sh; apply the update package first"
		exec bash "$ROOT/f/final.sh" "$@"
		;;
	abm|prsformer) ;;
	*) die "Unknown module: $MODULE (see: bash grid.sh help)" ;;
esac

# One public PRSformer entry. Independent-cohort mode preserves existing
# training caches/checkpoints and runtime checks without creating a second script.
if [[ $MODULE == prsformer ]]; then
	PRS_STANDALONE=FALSE
	PRS_MODULE_ARGS=()
	for arg in "$@"; do
		if [[ $arg == --standalone ]]; then
			PRS_STANDALONE=TRUE
		else
			PRS_MODULE_ARGS+=("$arg")
		fi
	done
	if [[ $PRS_STANDALONE == TRUE ]]; then
		exec bash "$ROOT/f/prsformer.sh" "${PRS_MODULE_ARGS[@]}"
	fi
	for arg in "$@"; do
		case $arg in -h|--help|help)
			cat <<'HELP'
GRID PRSformer — one public entry: bash grid.sh prsformer

Default: automatically prepare the common GRID cohort and frozen outer test half.
Default traits: height,ldl,t2dm; default device: cuda.
  cd /mnt/d/scripts/grid
  bash grid.sh prsformer

Optional: republish reports from an existing shared-split run without training.
  bash grid.sh prsformer --prsformer-mode report --prsformer-replace TRUE

Common options:
  --cache-dir DIR       GRID cache containing common.keep, families.tsv.gz and
                        prsformer.split.tsv.gz (default /tmp/grid-cache/grid).
  --out-root DIR        GRID report root (default /mnt/d/analysis/grid/grid).
  --prsformer-root DIR  Benchmark output root (default OUTROOT/benchmark).
  --prsformer-OPTION VALUE  Model options, e.g. --prsformer-python PATH,
                            --prsformer-epochs 30, --prsformer-mode train.
  --prsformer-arg ARG   One exact model token, e.g. --prsformer-arg --resume.
  --check / --dry-run   Validate or print the shared-cohort invocation.

Existing independent 60/20/20 runs, custom cohorts and runtime-only checks:
  bash grid.sh prsformer --standalone --help
  bash grid.sh prsformer --standalone --check-runtime
  bash grid.sh prsformer --standalone --mode report --replace
In --standalone mode, use ordinary model flags (--python, --device, --mode, etc.).
The explicit mode preserves the original cache/output defaults for independent runs.
HELP
			exit 0 ;;
		esac
	done
fi

PYTHON=${GRID_PYTHON:-python3}
PYTHON_EXPLICIT=FALSE
[[ -z ${GRID_PYTHON:-} ]] || PYTHON_EXPLICIT=TRUE
CACHE_DIR=/tmp/grid-cache/grid
OUT_ROOT=/mnt/d/analysis/grid/grid
PRSFORMER_ROOT=
STAGE=all
TRAITS=height,ldl,t2dm
REPLACE=TRUE
DRY_RUN=FALSE
PYTHON_HELP=FALSE
RUN_LOG=
GRID_ARGS=()
PRS_SHARED=()
PRS_ARGS=()

while (($#)); do
	case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
	case $1 in
		-h|--help|help) usage; exit 0 ;;
		--python-help) PYTHON_HELP=TRUE; shift ;;
		--python) need_value "$@"; PYTHON=$2; PYTHON_EXPLICIT=TRUE; shift 2 ;;
		--replace) REPLACE=TRUE; shift ;;
		--no-replace) REPLACE=FALSE; shift ;;
		--check) STAGE=check; shift ;;
		--check-device) STAGE=check; GRID_ARGS+=(--check-device); shift ;;
		--stage)
			need_value "$@"
			STAGE=$2
			[[ $STAGE != features ]] || STAGE=evolution
			shift 2
			;;
		--cache-dir) need_value "$@"; CACHE_DIR=$2; shift 2 ;;
		--out-root) need_value "$@"; OUT_ROOT=$2; shift 2 ;;
		--trait|--traits) need_value "$@"; TRAITS=$2; shift 2 ;;
		--prsformer-root) need_value "$@"; PRSFORMER_ROOT=$2; shift 2 ;;
		--dry-run)
			flag=$1
			flag_value=TRUE
			shift
			if (($#)); then
				case ${1^^} in TRUE|FALSE) flag_value=${1^^}; shift ;; esac
			fi
			DRY_RUN=$flag_value
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
		--pheno-file|--ancestry-file|--group-col|--dir-gen|--covariates|--height-col|--ldl-col|--t2dm-col|--chrs|--remove|--seed)
			need_value "$@"; GRID_ARGS+=("$1" "$2"); PRS_SHARED+=("$1" "$2"); shift 2 ;;
		--snpinfo)
			need_value "$@"; GRID_ARGS+=("$1" "$2"); PRS_SHARED+=(--snp-list "$2"); shift 2 ;;
		*) GRID_ARGS+=("$1"); shift ;;
	esac
done

case $STAGE in prepare|evolution|fit|predict|report|all|check) ;;
	*) die "Unknown ABM --stage: $STAGE" ;;
esac
[[ $MODULE != prsformer || $STAGE == all || $STAGE == check ]] || die 'For PRSformer stages, use --prsformer-mode prepare/train/predict/report'
[[ -n $CACHE_DIR && $CACHE_DIR != / ]] || die 'Use a dedicated non-root --cache-dir'
[[ -n $OUT_ROOT && $OUT_ROOT != / ]] || die 'Use a dedicated non-root --out-root'
[[ -n $PRSFORMER_ROOT ]] || PRSFORMER_ROOT="$OUT_ROOT/benchmark"

if [[ $DRY_RUN == FALSE ]]; then
	if [[ $PYTHON_EXPLICIT == FALSE ]]; then
		CONDA_INIT="$HOME/miniforge3/etc/profile.d/conda.sh"
		[[ -f $CONDA_INIT ]] || CONDA_INIT="$HOME/anaconda3/etc/profile.d/conda.sh"
		[[ -f $CONDA_INIT ]] || die 'Conda initialization script unavailable; install Miniforge/Conda first'
		set +u
		if source "$CONDA_INIT" && conda activate "${GRID_ENV_NAME:-grid}"; then CONDA_STATUS=0; else CONDA_STATUS=$?; fi
		set -u
		((CONDA_STATUS == 0)) || die 'Could not activate Conda environment grid; run install_grid.sh'
	fi
	PYTHON_COMMAND=$(command -v "$PYTHON") || die "Python unavailable: $PYTHON"
	PYTHON_DIRECTORY=$(cd -- "$(dirname -- "$PYTHON_COMMAND")" && pwd -P)
	PYTHON="$PYTHON_DIRECTORY/$(basename -- "$PYTHON_COMMAND")"
	export PATH="$PYTHON_DIRECTORY:$PATH"
fi

run_command() {
	printf 'RUN '; printf '%q ' "$@"; printf '\n'
	[[ $DRY_RUN == FALSE ]] || return 0
	if [[ -n $RUN_LOG ]]; then "$@" 2>&1 | tee -a "$RUN_LOG"; else "$@"; fi
}

if [[ $PYTHON_HELP == TRUE ]]; then
	run_command "$PYTHON" "$ROOT/f/main.py" --help
	exit 0
fi

run_abm_stage() {
	local replace_args=()
	[[ $REPLACE == FALSE ]] || replace_args=(--replace)
	run_command "$PYTHON" "$ROOT/f/main.py" --stage "$1" \
		--cache-dir "$CACHE_DIR" --out-root "$OUT_ROOT" --traits "$TRAITS" \
		--prsformer-root "$PRSFORMER_ROOT" --publish-only "${replace_args[@]}" "${GRID_ARGS[@]}"
}

run_prsformer() {
	local path
	local check_args=()
	[[ $STAGE != check ]] || check_args=(--check)
	if [[ $DRY_RUN == FALSE ]]; then
		for path in common.keep families.tsv.gz prsformer.split.tsv.gz; do
			[[ -s $CACHE_DIR/$path ]] || die "Missing shared roster $CACHE_DIR/$path; run 'grid.sh prsformer' to prepare it"
		done
	fi
	run_command bash "$ROOT/f/prsformer.sh" \
		--train-fraction 0.4 --validation-fraction 0.1 --cache-dir "$CACHE_DIR/prsformer" \
		"${PRS_ARGS[@]}" "${PRS_SHARED[@]}" "${check_args[@]}" --traits "$TRAITS" \
		--keep "$CACHE_DIR/common.keep" \
		--split-file "$CACHE_DIR/prsformer.split.tsv.gz" \
		--split-group-file "$CACHE_DIR/families.tsv.gz" \
		--output-root "$PRSFORMER_ROOT" --score-dir "$PRSFORMER_ROOT/scores"
}

if [[ $DRY_RUN == FALSE ]]; then
	if [[ $MODULE != prsformer ]]; then
		command -v "$PYTHON" >/dev/null 2>&1 || die "Python unavailable: $PYTHON"
		[[ -f $ROOT/f/main.py ]] || die "Missing GRID implementation: $ROOT/f/main.py"
	fi
	if [[ $STAGE != check ]]; then
		command -v flock >/dev/null 2>&1 || die 'flock is required'
		mkdir -p -- "$CACHE_DIR"
		exec {grid_lock}>"$CACHE_DIR/grid.lock"
		flock -n "$grid_lock" || die "Another GRID run is using $CACHE_DIR"
		mkdir -p -- "$CACHE_DIR/logs"
		RUN_LOG=$(mktemp "$CACHE_DIR/logs/grid.$(date -u +%Y%m%dT%H%M%SZ).XXXXXX.log")
		printf 'LOG %s\n' "$RUN_LOG"
	fi
fi

if [[ $MODULE == prsformer ]]; then
	# Preparation is shared, but model training remains entirely within PRSformer.
	# Read-only checks and report/predict-only actions do not rewrite the roster.
	PRS_MODE=all
	for ((i=0; i<${#PRS_ARGS[@]}; i++)); do
		[[ ${PRS_ARGS[$i]} != --mode ]] || PRS_MODE=${PRS_ARGS[$((i+1))]}
	done
	if [[ $STAGE != check && ( $PRS_MODE == all || $PRS_MODE == prepare ) ]]; then
		run_abm_stage prepare
	fi
	run_prsformer
else
	run_abm_stage "$STAGE"
fi
