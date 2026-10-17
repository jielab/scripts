#!/usr/bin/env bash
# PRSformer shared-cohort entry; public command: ./grid.sh prsformer.
# Parse the complete script before long jobs so edits cannot replay a dispatch.
{
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)

# The public Python entry normally owns logging; also support direct invocation.
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-/tmp/grid-logs/runs}
source "$ROOT/../0f/console.sh"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
need_value() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || die "Missing value for $1"; }
validate_traits() {
	[[ $1 =~ ^(height|ldl|t2dm)(,(height|ldl|t2dm))*$ ]] || die 'Use unique traits from height,ldl,t2dm'
	local trait
	local -a traits
	local -A seen=()
	IFS=, read -r -a traits <<<"$1"
	for trait in "${traits[@]}"; do
		[[ -z ${seen[$trait]+present} ]] || die "Repeated trait: $trait"
		seen[$trait]=1
	done
}

usage() {
    cat <<'HELP'
GRID PRSformer — one public entry: bash grid.sh prsformer

Default: reuse the frozen GRID cohort; prepare it only when no roster exists.
Default traits: height,ldl,t2dm; default device: cuda.
  cd /mnt/d/scripts/grid
  bash grid.sh prsformer

Optional: republish reports from an existing shared-split run without training.
  bash grid.sh prsformer --prsformer-mode report --prsformer-replace TRUE

Common options:
  --cache-dir DIR       GRID cache containing split.rds (default /tmp/grid-cache/grid).
                        A published OUTROOT/grid.split.rds can also supply the roster.
  --prsformer-cache-dir DIR  Prepared inputs (default PRSFORMER_ROOT/prsformer/cache).
  --prsformer-run-dir DIR    Checkpoints/history (default PRSFORMER_ROOT/prsformer/run).
  --out-root DIR        ABM report root (default /mnt/d/analysis/grid/abm).
  --prsformer-root DIR  Parent of the prsformer output folder (default parent of OUTROOT).
  --wait-lock          Wait for another run of the same module/cache.
  --prsformer-OPTION VALUE  Model options, e.g. --prsformer-python PATH,
                            --prsformer-epochs 30, --prsformer-mode train.
  --prsformer-arg ARG   One exact model token, e.g. --prsformer-arg --resume.
  --check / --dry-run   Validate or print the shared-cohort invocation.

Screen output shows stages, GPU checks and concise progress. Full logs are saved
under /tmp/grid-logs/runs (override SCRIPT_LOG_DIR); SCRIPT_VERBOSE=1 shows all output.

PRSformer stores an unchanged roster snapshot under its own CACHE/roster directory.
It requires the entire GRID cohort and will fail if phenotype/genotype alignment
drops participants. ABM features and caches are preserved. After training:
  bash grid.sh final
Or train and compare in order:
  bash grid.sh prsformer,final

Existing independent 60/20/20 runs, custom cohorts and runtime-only checks:
  bash grid.sh prsformer --standalone --help
  bash grid.sh prsformer --standalone --check-runtime
  bash grid.sh prsformer --standalone --mode report --replace
In --standalone mode, use ordinary model flags (--python, --device, --mode, etc.).
Both modes save prepared inputs and checkpoints under the prsformer output folder.
HELP
}

# Standalone runtime checks and independent-cohort runs use ordinary model flags.
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
    case $arg in -h|--help|help) usage; exit 0 ;; esac
done

PYTHON=${GRID_PYTHON:-python3}
PYTHON_EXPLICIT=FALSE
[[ -z ${GRID_PYTHON:-} ]] || PYTHON_EXPLICIT=TRUE
CACHE_DIR=/tmp/grid-cache/grid
OUT_ROOT=/mnt/d/analysis/grid/abm
PRSFORMER_ROOT=
PRS_CACHE_DIR=
STAGE=all
TRAITS=height,ldl,t2dm
REPLACE=TRUE
WAIT_LOCK=FALSE
SHARED_SEED=20260904
DRY_RUN=FALSE
RUN_LOG=
GRID_ARGS=()
PRS_SHARED=()
PRS_ARGS=()

while (($#)); do
	case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
	case $1 in
		-h|--help|help) usage; exit 0 ;;
		--python) need_value "$@"; PYTHON=$2; PYTHON_EXPLICIT=TRUE; shift 2 ;;
		--replace) REPLACE=TRUE; shift ;;
		--no-replace) REPLACE=FALSE; shift ;;
		--wait-lock) WAIT_LOCK=TRUE; shift ;;
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
		--prsformer-cache-dir) need_value "$@"; PRS_CACHE_DIR=$2; shift 2 ;;
		--dry-run)
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
		--seed) need_value "$@"; SHARED_SEED=$2; GRID_ARGS+=("$1" "$2"); PRS_SHARED+=("$1" "$2"); shift 2 ;;
		--pheno-file|--ancestry-file|--group-col|--dir-gen|--covariates|--height-col|--ldl-col|--t2dm-col|--chrs|--remove)
			need_value "$@"; GRID_ARGS+=("$1" "$2"); PRS_SHARED+=("$1" "$2"); shift 2 ;;
		--snpinfo)
			need_value "$@"; GRID_ARGS+=("$1" "$2"); PRS_SHARED+=(--snp-list "$2"); shift 2 ;;
		*) GRID_ARGS+=("$1"); shift ;;
	esac
done

validate_traits "$TRAITS"
[[ $STAGE == all || $STAGE == check ]] || die 'For PRSformer stages, use --prsformer-mode prepare/train/predict/report'
[[ -n $CACHE_DIR && $CACHE_DIR != / ]] || die 'Use a dedicated non-root --cache-dir'
[[ -n $OUT_ROOT && $OUT_ROOT != / ]] || die 'Use a dedicated non-root --out-root'
[[ -n $PRSFORMER_ROOT ]] || PRSFORMER_ROOT="$(dirname -- "$OUT_ROOT")"
[[ -n $PRS_CACHE_DIR ]] || PRS_CACHE_DIR="$PRSFORMER_ROOT/prsformer/cache"
PRS_ROSTER_DIR="$PRS_CACHE_DIR/roster"

if [[ $DRY_RUN == FALSE ]]; then
	if [[ $PYTHON_EXPLICIT == FALSE ]]; then
		CONDA_INIT="$HOME/miniforge3/etc/profile.d/conda.sh"
		[[ -f $CONDA_INIT ]] || CONDA_INIT="$HOME/anaconda3/etc/profile.d/conda.sh"
		[[ -f $CONDA_INIT ]] || die 'Conda initialization script unavailable; install Miniforge/Conda first'
		set +u
		# shellcheck disable=SC1090 # The installed Conda root is resolved above.
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

prepare_shared_cohort() {
	local replace_args=()
	[[ $REPLACE == FALSE ]] || replace_args=(--replace)
	run_command "$PYTHON" "$ROOT/f/main.py" --stage prepare \
		--cache-dir "$CACHE_DIR" --out-root "$OUT_ROOT" --traits "$TRAITS" \
		--prsformer-root "$PRSFORMER_ROOT" --publish-only "${replace_args[@]}" "${GRID_ARGS[@]}"
}

run_prsformer() {
	local path
	local check_args=()
	[[ $STAGE != check ]] || check_args=(--check)
	if [[ $DRY_RUN == FALSE ]]; then
		for path in common.keep families.tsv.gz prsformer.split.tsv.gz; do
			[[ -s $PRS_ROSTER_DIR/$path ]] || die "Missing shared roster $PRS_ROSTER_DIR/$path; run 'grid.sh prsformer' to prepare it"
		done
	fi
	run_command bash "$ROOT/f/prsformer.sh" \
		--train-fraction 0.4 --validation-fraction 0.1 \
		"${PRS_ARGS[@]}" "${PRS_SHARED[@]}" "${check_args[@]}" --traits "$TRAITS" \
		--cache-dir "$PRS_CACHE_DIR" --keep "$PRS_ROSTER_DIR/common.keep" \
		--split-file "$PRS_ROSTER_DIR/prsformer.split.tsv.gz" \
		--split-group-file "$PRS_ROSTER_DIR/families.tsv.gz" \
		--require-shared-cohort \
		--output-root "$PRSFORMER_ROOT" --score-dir "$PRSFORMER_ROOT/prsformer/scores"
}

acquire_lock() {
	local path=$1
	exec {module_lock}>"$path"
	if [[ $WAIT_LOCK == TRUE ]]; then
		printf 'WAIT for exclusive module cache: %s\n' "$path"
		flock "$module_lock"
	else
		flock -n "$module_lock" || die "Another run is using $path; use --wait-lock to queue this run"
	fi
}

if [[ $DRY_RUN == FALSE && $STAGE != check ]]; then
    command -v flock >/dev/null 2>&1 || die 'flock is required'
    mkdir -p -- "$CACHE_DIR" "$PRS_CACHE_DIR"
    acquire_lock "$PRS_CACHE_DIR/grid.wrapper.lock"
    mkdir -p -- "$CACHE_DIR/logs"
    RUN_LOG=$(mktemp "$CACHE_DIR/logs/grid.$(date -u +%Y%m%dT%H%M%SZ).XXXXXX.log")
    printf 'LOG %s\n' "$RUN_LOG"
fi

# Reuse frozen IDs and preserve all ABM features and input timestamps. Rewriting
# the roster on every invocation invalidates PRSformer cache/resume signatures.
PRS_MODE=all
for ((i=0; i<${#PRS_ARGS[@]}; i++)); do
	[[ ${PRS_ARGS[$i]} != --mode ]] || PRS_MODE=${PRS_ARGS[$((i+1))]}
done
if [[ $STAGE != check && ( $PRS_MODE == all || $PRS_MODE == prepare ) && ! -s $CACHE_DIR/split.rds && ! -s $OUT_ROOT/grid.split.rds ]]; then
	if [[ $DRY_RUN == FALSE ]]; then
		exec {prepare_lock}>"$CACHE_DIR/grid.lock"
		flock -n "$prepare_lock" || die "ABM is preparing the initial roster; retry after its preparation finishes"
	fi
	if [[ ! -s $CACHE_DIR/split.rds && ! -s $OUT_ROOT/grid.split.rds ]]; then prepare_shared_cohort; fi
	if [[ $DRY_RUN == FALSE ]]; then flock -u "$prepare_lock"; exec {prepare_lock}>&-; fi
fi
ROSTER_ARGS=(--cache-dir "$CACHE_DIR" --grid-split "$OUT_ROOT/grid.split.rds" --seed "$SHARED_SEED")
if [[ $STAGE == check ]]; then
	if [[ -d $PRS_ROSTER_DIR ]]; then
		ROSTER_ARGS+=(--check)
	elif [[ $DRY_RUN == FALSE ]]; then
		PRS_ROSTER_DIR=$(mktemp -d /tmp/grid-prsformer-roster-check.XXXXXX)
		trap 'rm -rf -- "$PRS_ROSTER_DIR"' EXIT
	fi
fi
ROSTER_ARGS+=(--roster-dir "$PRS_ROSTER_DIR")
run_command "$PYTHON" "$ROOT/f/prsformer.py" shared-roster "${ROSTER_ARGS[@]}"
run_prsformer

# Never read an edited script after a job returns.
exit 0
}

