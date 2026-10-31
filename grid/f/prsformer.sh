#!/usr/bin/env bash


# PRSformer: UKB genotype/phenotype alignment -> supervised training -> held-out scores.
# Official architecture: https://github.com/23andMe/PRSformer (developed by 23andMe, Inc.).
# Internal PRSformer module implementation; public entry: 4.prsformer.sh prsformer; helper: f/prsformer.py.
# Parse every stage before starting training so in-place script edits cannot
# replay training/publication when a long-running command returns.
{
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1

usage() {
	cat <<'HELP'
Usage:
    bash f/prsformer.sh --mode prepare|train|predict|report

Internal PRSformer runtime helper.
Public options: ./4.prsformer.sh prsformer --help
HELP
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
value_required() {
	[[ $# -ge 2 && -n $2 && $2 != --* ]] || die "Missing value for $1"
}
boolean_value() {
	FLAG_VALUE=TRUE
	FLAG_SHIFT=1
	if [[ $# -ge 2 ]]; then
		case ${2^^} in
			TRUE | FALSE) FLAG_VALUE=${2^^}; FLAG_SHIFT=2 ;;
		esac
	fi
}

# 🚩 Shared environment and PRSformer defaults
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
PYTHON=${GRID_PRSFORMER_PYTHON:-${GRID_PYTHON:-python3}}
PYTHON_EXPLICIT=FALSE
[[ -z ${GRID_PRSFORMER_PYTHON:-${GRID_PYTHON:-}} ]] || PYTHON_EXPLICIT=TRUE
UPSTREAM_DIR=${PRSFORMER_UPSTREAM_DIR:-/mnt/f/software/PRSformer}
MODE=all
TRAITS=height,ldl,t2dm
DIR_GEN=/mnt/f/gen/ukb/37/hap
PHENO_FILE=/mnt/d/data/ukb/Rdata/all.rds
ANCESTRY_FILE=/mnt/f/gen/ukb/pca_proj/ukb.ancestry.auto.tsv.gz
GROUP_COL=genetic_ancestry
COVARIATES=age,sex,PC1,PC2
HEIGHT_COL=height
LDL_COL=ldl
T2DM_COL=auto
SNP_LIST=/mnt/f/LD/csx/snpinfo_mult_1kg_hm3
CHRS=1-22
KEEP=none
REMOVE=/mnt/d/files/ukb.exclude.id
SPLIT_FILE=none
SPLIT_GROUP_FILE=none
REQUIRE_SHARED_COHORT=FALSE
TRAIN_FRACTION=0.6
VALIDATION_FRACTION=0.2
MAF=0.01
CALL_RATE=0.98
CHUNK_VARIANTS=32
CACHE_DIR=
RUN_DIR=
OUTPUT_ROOT=/mnt/d/analysis/grid
PUBLICATION_DIR=
SCORE_DIR=
MAX_SAMPLES=
MAX_VARIANTS=
DEVICE=cuda
ATTENTION=neighborhood
EMBED_DIM=64
HEADS=4
LAYERS=2
FF_DIM=128
KERNEL_SIZE=385
DILATION=1
# Keep 64 subjects/update. RTX 5090 Laptop benchmark: 16x4 keeps substantially
# more VRAM headroom than 32x2 at nearly the same throughput.
BATCH_SIZE=${PRSFORMER_BATCH_SIZE:-16}
ACCUMULATION_STEPS=${PRSFORMER_ACCUMULATION_STEPS:-4}
EPOCHS=30
PATIENCE=3
MIN_DELTA=1e-4
LR=5e-4
WEIGHT_DECAY=0.05
CLIP_GRAD=1.0
AMP=fp16
GRADIENT_CHECKPOINTING=TRUE
WORKERS=${PRSFORMER_WORKERS:-2}
THREADS=2
SEED=20260904
CHECKPOINT=
CHECK=FALSE
CHECK_RUNTIME=FALSE
DRY_RUN=FALSE
RESUME=FALSE
REPLACE=FALSE

# 🚩 Command-line options; preserve arguments as arrays, never eval shell text
while (($#)); do
	case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
	case $1 in
		-h | --help | help) usage; exit 0 ;;
		--check) CHECK=TRUE; shift ;;
		--check-runtime) CHECK_RUNTIME=TRUE; shift ;;
		--require-shared-cohort) REQUIRE_SHARED_COHORT=TRUE; shift ;;
		--dry-run | --resume | --replace)
			option=$1
			boolean_value "$@"
			case $option in
				--dry-run) DRY_RUN=$FLAG_VALUE ;;
				--resume) RESUME=$FLAG_VALUE ;;
				--replace) REPLACE=$FLAG_VALUE ;;
			esac
			shift "$FLAG_SHIFT"
			;;
		--gradient-checkpointing) GRADIENT_CHECKPOINTING=TRUE; shift ;;
		--no-gradient-checkpointing) GRADIENT_CHECKPOINTING=FALSE; shift ;;
		*)
			value_required "$@"
			case $1 in
				--mode | --stage) MODE=$2 ;;
				--trait | --traits) TRAITS=$2 ;;
				--python) PYTHON=$2; PYTHON_EXPLICIT=TRUE ;;
				--upstream-dir) UPSTREAM_DIR=$2 ;;
				--dir-gen) DIR_GEN=$2 ;;
				--pheno-file | --phe-file) PHENO_FILE=$2 ;;
				--ancestry-file) ANCESTRY_FILE=$2 ;;
				--group-col) GROUP_COL=$2 ;;
				--covariates) COVARIATES=$2 ;;
				--height-col) HEIGHT_COL=$2 ;;
				--ldl-col) LDL_COL=$2 ;;
				--t2dm-col) T2DM_COL=$2 ;;
				--snp-list) SNP_LIST=$2 ;;
				--chrs) CHRS=$2 ;;
				--keep) KEEP=$2 ;;
				--remove) REMOVE=$2 ;;
				--split-file) SPLIT_FILE=$2 ;;
				--split-group-file) SPLIT_GROUP_FILE=$2 ;;
				--train-fraction) TRAIN_FRACTION=$2 ;;
				--validation-fraction) VALIDATION_FRACTION=$2 ;;
				--maf) MAF=$2 ;;
				--call-rate) CALL_RATE=$2 ;;
				--chunk-variants) CHUNK_VARIANTS=$2 ;;
				--cache-dir) CACHE_DIR=$2 ;;
				--run-dir) RUN_DIR=$2 ;;
				--output-root) OUTPUT_ROOT=$2 ;;
				--publication-dir) PUBLICATION_DIR=$2 ;;
				--score-dir) SCORE_DIR=$2 ;;
				--max-samples) MAX_SAMPLES=$2 ;;
				--max-variants) MAX_VARIANTS=$2 ;;
				--device) DEVICE=$2 ;;
				--attention) ATTENTION=$2 ;;
				--embed-dim) EMBED_DIM=$2 ;;
				--heads) HEADS=$2 ;;
				--layers) LAYERS=$2 ;;
				--ff-dim) FF_DIM=$2 ;;
				--kernel-size) KERNEL_SIZE=$2 ;;
				--dilation) DILATION=$2 ;;
				--batch-size) BATCH_SIZE=$2 ;;
				--accumulation-steps) ACCUMULATION_STEPS=$2 ;;
				--epochs) EPOCHS=$2 ;;
				--patience) PATIENCE=$2 ;;
				--min-delta) MIN_DELTA=$2 ;;
				--lr) LR=$2 ;;
				--weight-decay) WEIGHT_DECAY=$2 ;;
				--clip-grad) CLIP_GRAD=$2 ;;
				--amp) AMP=$2 ;;
				--workers) WORKERS=$2 ;;
				--threads) THREADS=$2 ;;
				--seed) SEED=$2 ;;
				--checkpoint) CHECKPOINT=$2 ;;
				*) die "Unknown option: $1 (see --help)" ;;
			esac
			shift 2
			;;
	esac
done
case $MODE in all | prepare | train | predict | report) ;; *) die "Unknown mode: $MODE" ;; esac
case $ATTENTION in neighborhood | global) ;; *) die '--attention must be neighborhood or global' ;; esac
case $AMP in fp16 | bf16 | off) ;; *) die '--amp must be fp16, bf16, or off' ;; esac
[[ $RESUME == FALSE || $REPLACE == FALSE ]] || die '--resume and --replace are incompatible; resuming must reuse identical input caches'
[[ $RESUME == FALSE || $MODE == train || $MODE == all ]] || die '--resume applies to all/train only'
[[ -n $CACHE_DIR ]] || CACHE_DIR="$OUTPUT_ROOT/prsformer/cache"
[[ -n $RUN_DIR ]] || RUN_DIR="$OUTPUT_ROOT/prsformer/run"
[[ -n $SCORE_DIR ]] || SCORE_DIR="$OUTPUT_ROOT/prsformer/scores"
[[ -n $CHECKPOINT ]] || CHECKPOINT="$RUN_DIR/model.pt"
MODEL_COVARIATES=$COVARIATES
[[ ${COVARIATES,,} != none ]] || MODEL_COVARIATES=
# Shared preparation fixes the full roster independently of the training tasks.
# A trait subset must reuse/rebuild the same cache and frozen input signature.
PREPARE_TRAITS=$TRAITS
[[ $REQUIRE_SHARED_COHORT == FALSE ]] || PREPARE_TRAITS=height,ldl,t2dm

# 🚩 The helper interface is shared by standalone and combined modes
PREPARE=("$PYTHON" "$ROOT/f/prsformer.py" prepare
	--dir-gen "$DIR_GEN" --pheno-file "$PHENO_FILE" --ancestry-file "$ANCESTRY_FILE"
	--group-col "$GROUP_COL" --traits "$PREPARE_TRAITS" --covariates "$COVARIATES"
	--height-col "$HEIGHT_COL" --ldl-col "$LDL_COL" --t2dm-col "$T2DM_COL"
	--snp-list "$SNP_LIST" --chrs "$CHRS" --keep "$KEEP" --remove "$REMOVE"
	--split-file "$SPLIT_FILE" --split-group-file "$SPLIT_GROUP_FILE"
	--train-fraction "$TRAIN_FRACTION" --validation-fraction "$VALIDATION_FRACTION"
	--seed "$SEED" --maf "$MAF" --call-rate "$CALL_RATE"
	--chunk-variants "$CHUNK_VARIANTS" --cache-dir "$CACHE_DIR")
[[ -z $MAX_SAMPLES ]] || PREPARE+=(--max-samples "$MAX_SAMPLES")
[[ -z $MAX_VARIANTS ]] || PREPARE+=(--max-variants "$MAX_VARIANTS")
[[ $REQUIRE_SHARED_COHORT == FALSE ]] || PREPARE+=(--require-shared-cohort)
[[ $REPLACE == FALSE ]] || PREPARE+=(--replace)
MODEL_ARGS=(--genotypes "$CACHE_DIR/genotypes.npy" --data "$CACHE_DIR/data.tsv.gz"
	--variants "$CACHE_DIR/variants.tsv.gz" --upstream-dir "$UPSTREAM_DIR"
	--out-dir "$RUN_DIR" --traits "$TRAITS" --covariates "$MODEL_COVARIATES"
	--device "$DEVICE" --attention "$ATTENTION" --embed-dim "$EMBED_DIM"
	--heads "$HEADS" --layers "$LAYERS" --ff-dim "$FF_DIM"
	--kernel-size "$KERNEL_SIZE" --dilation "$DILATION" --batch-size "$BATCH_SIZE"
	--accumulation-steps "$ACCUMULATION_STEPS" --epochs "$EPOCHS" --patience "$PATIENCE"
	--min-delta "$MIN_DELTA" --lr "$LR" --weight-decay "$WEIGHT_DECAY"
	--clip-grad "$CLIP_GRAD" --amp "$AMP" --workers "$WORKERS" --threads "$THREADS" --seed "$SEED")
if [[ $GRADIENT_CHECKPOINTING == TRUE ]]; then MODEL_ARGS+=(--gradient-checkpointing); else MODEL_ARGS+=(--no-gradient-checkpointing); fi
RUNTIME=("$PYTHON" "$ROOT/f/prsformer.py" runtime "${MODEL_ARGS[@]}")
TRAIN=("$PYTHON" "$ROOT/f/prsformer.py" train "${MODEL_ARGS[@]}")
[[ $RESUME == FALSE ]] || TRAIN+=(--resume)
PREDICT=("$PYTHON" "$ROOT/f/prsformer.py" predict "${MODEL_ARGS[@]}" --checkpoint "$CHECKPOINT")
PUBLISH=("$PYTHON" "$ROOT/f/prsformer.py" publish --run-dir "$RUN_DIR"
	--cache-dir "$CACHE_DIR" --output-root "$OUTPUT_ROOT" --score-dir "$SCORE_DIR" --traits "$TRAITS")
[[ -z $PUBLICATION_DIR ]] || PUBLISH+=(--publication-dir "$PUBLICATION_DIR")
[[ $REPLACE == FALSE ]] || PUBLISH+=(--replace)

run_command() {
	printf '\nRUN '
	printf '%q ' "$@"
	printf '\n'
	[[ $DRY_RUN == FALSE ]] || return 0
	if [[ $CHECK == TRUE || $CHECK_RUNTIME == TRUE ]]; then
		"$@"
	else
		# Keep the execution record and checkpoints visible before epoch one ends.
		mkdir -p -- "$RUN_DIR"
		{ printf '\nRUN '; printf '%q ' "$@"; printf '\n'; } >> "$RUN_DIR/run.log"
		"$@" 2>&1 | tee -a "$RUN_DIR/run.log"
	fi
}
cache_ready() {
	[[ -s $CACHE_DIR/genotypes.npy && -s $CACHE_DIR/data.tsv.gz && -s $CACHE_DIR/variants.tsv.gz && -s $CACHE_DIR/prepare.json ]]
}
need_cache() {
	cache_ready || die "Prepared cache incomplete: $CACHE_DIR; run --mode prepare first"
	if [[ -n ${GRID_WORKFLOW_ANALYSIS_ROOT:-} && $CHECK == FALSE ]]; then
		"$PYTHON" "$ROOT/f/workflow_status.py" --analysis-root "$GRID_WORKFLOW_ANALYSIS_ROOT" --traits "$TRAITS" --prs-cache "$CACHE_DIR"
	fi
}

prepare_for_training() {
	if [[ $RESUME == TRUE ]]; then
		printf 'RESUME: reuse frozen prepared inputs in %s; training validates them against the saved checkpoint.\n' "$CACHE_DIR"
	else
		run_command "${PREPARE[@]}"
	fi
}

# 🚩 Dependencies: fail before decoding a production dataset when GPU setup is invalid
dependencies() {
	"$PYTHON" - "$MODE" "$ATTENTION" "$DEVICE" "$AMP" "$UPSTREAM_DIR" "$SPLIT_GROUP_FILE" "$ROOT" <<'PY'
import importlib
from pathlib import Path
import sys

mode, attention, device_name, amp, upstream, group_file, root = sys.argv[1:]
modules = {"numpy", "pandas"}
if mode in ("all", "prepare"):
	modules.update(("pgenlib", "pyreadr"))
	if group_file.lower() != "none":
		modules.add("sklearn")
if mode in ("all", "report"):
	modules.update(("matplotlib", "openpyxl", "pyreadr"))
if mode in ("all", "train", "predict"):
	modules.add("scipy")
for name in sorted(modules):
	try:
		importlib.import_module(name)
	except (ImportError, OSError) as exc:
		raise SystemExit(f"ERROR: Cannot import {name} in {sys.executable}: {exc}. Follow install.sh --help installation instructions.")
if mode in ("all", "train", "predict"):
	for name in ("model.py", "modules.py", "utils.py"):
		if not (Path(upstream) / "src" / name).is_file():
			raise SystemExit(f"ERROR: Missing official source {Path(upstream) / 'src' / name}; set --upstream-dir to the documented PRSformer checkout.")
	try:
		from types import SimpleNamespace
		spec = importlib.util.spec_from_file_location("prsformer_runtime", Path(root) / "f/prsformer.py")
		adapter = importlib.util.module_from_spec(spec)
		spec.loader.exec_module(adapter)
		torch, device, _ = adapter.runtime(SimpleNamespace(device=device_name, amp=amp, threads=1), attention)
		adapter.load_official_model(upstream, attention)
		print(f"DEPENDENCIES OK: python={sys.executable}; torch={torch.__version__}; device={device}; attention={attention}")
	except (ImportError, OSError, RuntimeError, ValueError, AssertionError) as exc:
		raise SystemExit(f"ERROR: PRSformer runtime preflight failed: {exc} See install.sh --help.")
else:
	print(f"DEPENDENCIES OK: python={sys.executable}; mode={mode}; GPU not required for this mode")
PY
	# Exercise compiled kernels before a full cohort is decoded or written.
	if [[ $MODE == all ]]; then run_command "${RUNTIME[@]}"; fi
}

# 🚩 Print-only execution plan
if [[ $DRY_RUN == TRUE ]]; then
	if [[ $CHECK_RUNTIME == TRUE ]]; then
		run_command "${RUNTIME[@]}"
	elif [[ $CHECK == TRUE ]]; then
		case $MODE in
			all)
				if [[ $RESUME == FALSE ]]; then run_command "${PREPARE[@]}" --check; fi
				run_command "${TRAIN[@]}" --check
				;;
			prepare) run_command "${PREPARE[@]}" --check ;;
			train) run_command "${TRAIN[@]}" --check ;;
			predict) run_command "${PREDICT[@]}" --check ;;
			report) printf 'PLAN: check prepared cache and staged reporting files.\n' ;;
		esac
	else
		case $MODE in
			all) prepare_for_training; run_command "${TRAIN[@]}"; run_command "${PUBLISH[@]}" ;;
			prepare) run_command "${PREPARE[@]}" ;;
			train) run_command "${TRAIN[@]}" ;;
			predict) run_command "${PREDICT[@]}" ;;
			report) run_command "${PUBLISH[@]}" ;;
		esac
	fi
	printf '\nPLAN ONLY: no dataset was read and no runtime compatibility was validated.\n'
	exit 0
fi
if [[ $PYTHON_EXPLICIT == FALSE ]]; then
	# Reuse installed dependencies. Consolidating installer files does not require
	# reinstalling a previously working PRSformer environment.
	PYTHON=
	for candidate in "$HOME/miniforge3/envs/${GRID_ENV_NAME:-grid}/bin/python" \
		"$HOME/anaconda3/envs/${GRID_ENV_NAME:-grid}/bin/python" \
		"$HOME/.venvs/grid-prsformer/bin/python" "$(command -v python3 || true)"; do
		[[ -x $candidate ]] || continue
		if "$candidate" - "$ATTENTION" <<'PY' >/dev/null 2>&1
import importlib.metadata as metadata
import importlib.util
import sys
packages = ('torch', 'numpy', 'pandas', 'scipy', 'pgenlib', 'pyreadr', 'matplotlib', 'openpyxl')
if not all(importlib.util.find_spec(name) is not None for name in packages):
	raise SystemExit(1)
if sys.argv[1] == 'neighborhood':
	stack = tuple(metadata.version(name).split('+')[0] for name in ('torch', 'natten'))
	if stack not in {('2.12.0', '0.21.7')}:
		raise SystemExit(1)
PY
		then PYTHON=$candidate; break; fi
	done
	[[ -n $PYTHON ]] || die 'No existing compatible PRSformer environment found; select --python, or use install.sh to install missing dependencies'
	printf 'PRSformer Python: %s (existing environment)\n' "$PYTHON"
fi
PYTHON_COMMAND=$(command -v "$PYTHON") || die "Python unavailable: $PYTHON"
PYTHON_DIRECTORY=$(cd -- "$(dirname -- "$PYTHON_COMMAND")" && pwd -P)
PYTHON="$PYTHON_DIRECTORY/$(basename -- "$PYTHON_COMMAND")"
export PATH="$PYTHON_DIRECTORY:$PATH"
# Bind every command to the selected interpreter, including an existing venv
# whose directory may provide python without a python3 alias.
PREPARE[0]=$PYTHON
RUNTIME[0]=$PYTHON
TRAIN[0]=$PYTHON
PREDICT[0]=$PYTHON
PUBLISH[0]=$PYTHON
command -v -- "$PYTHON" >/dev/null 2>&1 || die "Python not found: $PYTHON; run bash install.sh, or select an existing environment with --python"
command -v flock >/dev/null 2>&1 || die 'flock is missing; install the Linux/WSL util-linux package before running PRSformer'
[[ -f $ROOT/f/prsformer.py ]] || die "Missing PRSformer implementation helpers under $ROOT/f"
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS

# 🚩 Read-only checks; no inference, model fitting, or publication
if [[ $CHECK_RUNTIME == TRUE ]]; then
	run_command "${RUNTIME[@]}"
	exit 0
fi
if [[ $CHECK == TRUE ]]; then
	dependencies
	if [[ $MODE == prepare || ( $MODE == all && $RESUME == FALSE ) ]]; then run_command "${PREPARE[@]}" --check; fi
	case $MODE in
		all)
			if cache_ready; then run_command "${TRAIN[@]}" --check; else printf 'CHECK: raw inputs passed; prepare the cache to enable the official-model data/forward/backward check.\n'; fi
			;;
		train) need_cache; run_command "${TRAIN[@]}" --check ;;
		predict) need_cache; run_command "${PREDICT[@]}" --check ;;
		report)
			[[ -s $RUN_DIR/metrics.tsv && -s $RUN_DIR/test_predictions.tsv.gz ]] || die "Missing staged metrics/predictions under $RUN_DIR"
			printf 'CHECK: report dependencies and staged input files are present; no files published.\n'
			;;
	esac
	exit 0
fi
dependencies

# 🚩 Hold one canonical-cache lock across preparation, training, prediction and publication
# The model helper has its own independent output-directory lock. Do not acquire
# this cache key again inside a helper: this descriptor stays open until shell exit.
mkdir -p -- "$CACHE_DIR"
CACHE_LOCK_FILE="$CACHE_DIR/.pipeline.lock"
exec {CACHE_LOCK_FD}>"$CACHE_LOCK_FILE"
flock -n "$CACHE_LOCK_FD" || die "Another PRSformer invocation is using this cache: $CACHE_DIR. Wait for it to finish or choose a different --cache-dir."

if [[ $MODE == all || $MODE == train ]]; then
	if [[ $RESUME == TRUE && -f $RUN_DIR/evaluation.json ]]; then
		die "Run already complete: $RUN_DIR; use --mode report to publish it, or a new --run-dir to retrain."
	fi
	if [[ $RESUME == FALSE && ( -e $RUN_DIR/model.pt || -e $RUN_DIR/training.pt ) ]]; then
		die "A saved model already exists in $RUN_DIR. Use --resume for an interrupted run, --mode report for a completed run, or a new --run-dir to retrain; --replace does not delete models."
	fi
fi
case $MODE in
	all) prepare_for_training; need_cache; run_command "${TRAIN[@]}"; run_command "${PUBLISH[@]}" ;;
	prepare) run_command "${PREPARE[@]}" ;;
	train) need_cache; run_command "${TRAIN[@]}" ;;
	predict) need_cache; run_command "${PREDICT[@]}" ;;
	report) run_command "${PUBLISH[@]}" ;;
esac

# Keep unfinished runs resumable. Successful evaluations no longer need the
# large dosage matrix for reports; preserve its metadata and a verified receipt.
if [[ $MODE != prepare && ${PRSFORMER_KEEP_GENOTYPES:-1} != 1 ]]; then
	run_command "$PYTHON" "$ROOT/f/prsformer_cleanup.py" --cache-dir "$CACHE_DIR" --run-dir "$RUN_DIR" --cache-lock-held
fi

# Exit inside the parsed block without rereading the script from disk.
exit 0
}
