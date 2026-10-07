#!/usr/bin/env bash
# PRSformer: UKB genotype/phenotype alignment -> supervised training -> held-out scores.
# Official architecture: https://github.com/23andMe/PRSformer (developed by 23andMe, Inc.).
# Independent grid adapter; helpers: f/3.prsformer_data.py and f/3.prsformer.py.
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1

usage() {
	cat <<'HELP'
PRSformer — joint height / ldl / baseline t2dm training

Usage (WSL / Linux):
  cd /mnt/d/scripts/grid
  bash 3.prsformer.sh --traits height,ldl,t2dm --check
  bash 3.prsformer.sh --traits height,ldl,t2dm
  bash 3.prsformer.sh --mode train --resume
  bash 3.prsformer.sh --mode report

Modes:
  all       Prepare -> train -> publish (default; train also scores held-out test data).
  prepare   Align inputs, split subjects, training-only SNP QC, dosage memmap.
  train     Train from an existing cache and evaluate the validation-selected model.
  predict   Recalculate test predictions from --checkpoint (default RUN/model.pt).
  report    Publish staged metrics as XLSX/PNG and individual scores as RDS.

Data and paths:
  --trait NAME / --traits LIST  height,ldl,t2dm (joint masked multi-task training).
  --python PATH                ~/.venvs/grid-prsformer/bin/python; see README setup.
  --upstream-dir DIR           /mnt/f/software/PRSformer (official checkout).
  --dir-gen DIR                /mnt/f/gen/ukb/37/hap (chr*.pgen/pvar/psam).
  --pheno-file FILE            /mnt/d/data/ukb/phe/Rdata/all.rds.
  --ancestry-file FILE         /mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz.
  --group-col NAME             genetic_ancestry.
  --covariates LIST            age,sex,PC1,PC2; numeric columns, training-only fitting.
  --height-col / --ldl-col     height / ldl.
  --t2dm-col NAME              auto; use an explicitly defined 0/1 baseline endpoint.
  --snp-list FILE              /mnt/f/refLD/csx/snpinfo_mult_1kg_hm3.
  --chrs LIST                  1-22; genomic subsets require their own cache/run dirs.
  --keep FILE / --remove FILE  none / /mnt/d/files/ukb.exclude.id; 'none' disables.
  --split-file FILE            Explicit eid,split table; default none.
  --split-group-file FILE      Relatedness/family grouping; default none.
  --train-fraction N           0.6; --validation-fraction 0.2; remainder is test.
  --maf N / --call-rate N      0.01 / 0.98, evaluated using training genotypes only.
  --chunk-variants N           32; controls genotype preparation working memory.
  --cache-dir DIR              /tmp/grid-prsformer (sensitive temporary data).
  --run-dir DIR                CACHE/run (models and temporary test predictions).
  --output-root DIR            /mnt/d/analysis/grid (aggregate reports).
  --score-dir DIR              /mnt/d/data/ukb/pgs (individual RDS).
  --max-samples / --max-variants N  Explicit small-data development limits.

Model and training:
  --device DEVICE              cuda; e.g. cuda:0.
  --attention TYPE             neighborhood; global is limited to <=4096 SNPs.
  --embed-dim N / --heads N     64 / 4.
  --layers N / --ff-dim N       2 / 128.
  --kernel-size N              385; --dilation 1 (or one value per layer).
  --batch-size N               1; --accumulation-steps 64.
  --epochs N / --patience N     30 / 5; --min-delta 1e-4.
  --lr N / --weight-decay N     5e-4 / 0.05; --clip-grad 1.0.
  --amp fp16|bf16|off           fp16; CPU smoke tests require off.
  --no-gradient-checkpointing  Disable activation checkpointing (enabled by default).
  --workers N / --threads N    0 / 4.
  --seed N                     20260904; controls subject split and training.
  --checkpoint FILE            RUN/model.pt for --mode predict.

Controls:
  --check                      Validate without training or publishing. all/prepare
                               first checks raw data; an existing cache also enables
                               an actual small official-model forward/backward check.
  --dry-run [TRUE|FALSE]        Print commands only; do not read datasets or import GPU code.
  --resume [TRUE|FALSE]         Resume an interrupted all/train run with identical settings.
  --replace [TRUE|FALSE]        Replace preparation cache and/or published reports/scores.
                               Saved models are retained: choose a new --run-dir to retrain.
  --help                       Show this message.

No packages are installed automatically. Production neighborhood attention requires
the documented PyTorch 2.6.0 + CUDA 12.4 + NATTEN 0.17.5 environment and NVIDIA GPU.
This trains on individual genotypes/phenotypes; no CSx GWAS weights or pretraining
are implied. RDS scores contain held-out test individuals only, not all-subject OOF.
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

# 🚩 Independent environment and grid defaults
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PYTHON=${GRID_PRSFORMER_PYTHON:-${HOME}/.venvs/grid-prsformer/bin/python}
UPSTREAM_DIR=${PRSFORMER_UPSTREAM_DIR:-/mnt/f/software/PRSformer}
MODE=all
TRAITS=height,ldl,t2dm
DIR_GEN=/mnt/f/gen/ukb/37/hap
PHENO_FILE=/mnt/d/data/ukb/phe/Rdata/all.rds
ANCESTRY_FILE=/mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz
GROUP_COL=genetic_ancestry
COVARIATES=age,sex,PC1,PC2
HEIGHT_COL=height
LDL_COL=ldl
T2DM_COL=auto
SNP_LIST=/mnt/f/refLD/csx/snpinfo_mult_1kg_hm3
CHRS=1-22
KEEP=none
REMOVE=/mnt/d/files/ukb.exclude.id
SPLIT_FILE=none
SPLIT_GROUP_FILE=none
TRAIN_FRACTION=0.6
VALIDATION_FRACTION=0.2
MAF=0.01
CALL_RATE=0.98
CHUNK_VARIANTS=32
CACHE_DIR=/tmp/grid-prsformer
RUN_DIR=
OUTPUT_ROOT=/mnt/d/analysis/grid
SCORE_DIR=/mnt/d/data/ukb/pgs
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
BATCH_SIZE=1
ACCUMULATION_STEPS=64
EPOCHS=30
PATIENCE=5
MIN_DELTA=1e-4
LR=5e-4
WEIGHT_DECAY=0.05
CLIP_GRAD=1.0
AMP=fp16
GRADIENT_CHECKPOINTING=TRUE
WORKERS=0
THREADS=4
SEED=20260904
CHECKPOINT=
CHECK=FALSE
DRY_RUN=FALSE
RESUME=FALSE
REPLACE=FALSE

# 🚩 Command-line options; preserve arguments as arrays, never eval shell text
while (($#)); do
	case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
	case $1 in
		-h | --help | help) usage; exit 0 ;;
		--check) CHECK=TRUE; shift ;;
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
				--python) PYTHON=$2 ;;
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
[[ -n $RUN_DIR ]] || RUN_DIR="$CACHE_DIR/run"
[[ -n $CHECKPOINT ]] || CHECKPOINT="$RUN_DIR/model.pt"
MODEL_COVARIATES=$COVARIATES
[[ ${COVARIATES,,} != none ]] || MODEL_COVARIATES=

# 🚩 The helper interface is shared by standalone and combined modes
PREPARE=("$PYTHON" "$ROOT/f/3.prsformer_data.py" prepare
	--dir-gen "$DIR_GEN" --pheno-file "$PHENO_FILE" --ancestry-file "$ANCESTRY_FILE"
	--group-col "$GROUP_COL" --traits "$TRAITS" --covariates "$COVARIATES"
	--height-col "$HEIGHT_COL" --ldl-col "$LDL_COL" --t2dm-col "$T2DM_COL"
	--snp-list "$SNP_LIST" --chrs "$CHRS" --keep "$KEEP" --remove "$REMOVE"
	--split-file "$SPLIT_FILE" --split-group-file "$SPLIT_GROUP_FILE"
	--train-fraction "$TRAIN_FRACTION" --validation-fraction "$VALIDATION_FRACTION"
	--seed "$SEED" --maf "$MAF" --call-rate "$CALL_RATE"
	--chunk-variants "$CHUNK_VARIANTS" --cache-dir "$CACHE_DIR")
[[ -z $MAX_SAMPLES ]] || PREPARE+=(--max-samples "$MAX_SAMPLES")
[[ -z $MAX_VARIANTS ]] || PREPARE+=(--max-variants "$MAX_VARIANTS")
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
TRAIN=("$PYTHON" "$ROOT/f/3.prsformer.py" train "${MODEL_ARGS[@]}")
[[ $RESUME == FALSE ]] || TRAIN+=(--resume)
PREDICT=("$PYTHON" "$ROOT/f/3.prsformer.py" predict "${MODEL_ARGS[@]}" --checkpoint "$CHECKPOINT")
PUBLISH=("$PYTHON" "$ROOT/f/3.prsformer_data.py" publish --run-dir "$RUN_DIR"
	--cache-dir "$CACHE_DIR" --output-root "$OUTPUT_ROOT" --score-dir "$SCORE_DIR" --traits "$TRAITS")
[[ $REPLACE == FALSE ]] || PUBLISH+=(--replace)

run_command() {
	printf '\nRUN '
	printf '%q ' "$@"
	printf '\n'
	[[ $DRY_RUN == FALSE ]] || return 0
	"$@"
}
cache_ready() {
	[[ -s $CACHE_DIR/genotypes.npy && -s $CACHE_DIR/data.tsv.gz && -s $CACHE_DIR/variants.tsv.gz && -s $CACHE_DIR/prepare.json ]]
}
need_cache() { cache_ready || die "Prepared cache incomplete: $CACHE_DIR; run --mode prepare first"; }

# 🚩 Dependencies: fail before decoding a production dataset when GPU setup is invalid
dependencies() {
	"$PYTHON" - "$MODE" "$ATTENTION" "$DEVICE" "$AMP" "$UPSTREAM_DIR" "$SPLIT_GROUP_FILE" <<'PY'
import importlib
from pathlib import Path
import sys

mode, attention, device_name, amp, upstream, group_file = sys.argv[1:]
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
		raise SystemExit(f"ERROR: Cannot import {name} in {sys.executable}: {exc}. Follow README.prsformer.md installation instructions.")
if mode in ("all", "train", "predict"):
	for name in ("model.py", "modules.py", "utils.py"):
		if not (Path(upstream) / "src" / name).is_file():
			raise SystemExit(f"ERROR: Missing official source {Path(upstream) / 'src' / name}; set --upstream-dir to the documented PRSformer checkout.")
	try:
		torch = importlib.import_module("torch")
		device = torch.device(device_name)
		if device.type not in ("cuda", "cpu"):
			raise RuntimeError("Only cuda or cpu devices are supported.")
		if device.type == "cuda":
			if not torch.cuda.is_available():
				raise RuntimeError("CUDA is unavailable to this Python environment.")
			torch.cuda.get_device_properties(device)
		if device.type == "cpu" and (attention != "global" or amp != "off"):
			raise RuntimeError("CPU is only for an explicit small --attention global --amp off smoke test.")
		if attention == "neighborhood":
			if str(torch.__version__).split("+")[0] != "2.6.0" or torch.version.cuda != "12.4":
				raise RuntimeError("Use the pinned torch 2.6.0 CUDA 12.4 build for production neighborhood attention.")
			natten = importlib.import_module("natten")
			if str(getattr(natten, "__version__", "")).split("+")[0] != "0.17.5":
				raise RuntimeError("Use natten==0.17.5+torch260cu124; newer APIs are not interchangeable.")
			if not all(hasattr(natten, name) for name in ("NeighborhoodAttention1D", "use_fused_na", "is_fna_enabled")):
				raise RuntimeError("The NATTEN build lacks the official PRSformer attention API.")
		if amp == "bf16" and device.type == "cuda" and not torch.cuda.is_bf16_supported():
			raise RuntimeError("BF16 is not supported on this GPU; use fp16 or off.")
		print(f"DEPENDENCIES OK: python={sys.executable}; torch={torch.__version__}; device={device}; attention={attention}")
	except (ImportError, OSError, RuntimeError, ValueError, AssertionError) as exc:
		raise SystemExit(f"ERROR: PRSformer runtime preflight failed: {exc} See README.prsformer.md.")
else:
	print(f"DEPENDENCIES OK: python={sys.executable}; mode={mode}; GPU not required for this mode")
PY
}

# 🚩 Print-only execution plan
if [[ $DRY_RUN == TRUE ]]; then
	if [[ $CHECK == TRUE ]]; then
		case $MODE in
			all) run_command "${PREPARE[@]}" --check; run_command "${TRAIN[@]}" --check ;;
			prepare) run_command "${PREPARE[@]}" --check ;;
			train) run_command "${TRAIN[@]}" --check ;;
			predict) run_command "${PREDICT[@]}" --check ;;
			report) printf 'PLAN: check prepared cache and staged reporting files.\n' ;;
		esac
	else
		case $MODE in
			all) run_command "${PREPARE[@]}"; run_command "${TRAIN[@]}"; run_command "${PUBLISH[@]}" ;;
			prepare) run_command "${PREPARE[@]}" ;;
			train) run_command "${TRAIN[@]}" ;;
			predict) run_command "${PREDICT[@]}" ;;
			report) run_command "${PUBLISH[@]}" ;;
		esac
	fi
	printf '\nPLAN ONLY: no dataset was read and no runtime compatibility was validated.\n'
	exit 0
fi
command -v -- "$PYTHON" >/dev/null 2>&1 || die "Python not found: $PYTHON; use --python or create the environment described in README.prsformer.md"
command -v flock >/dev/null 2>&1 || die 'flock is missing; install the Linux/WSL util-linux package before running PRSformer'
[[ -f $ROOT/f/3.prsformer_data.py && -f $ROOT/f/3.prsformer.py ]] || die "Keep 3.prsformer.sh beside its f/ helper directory"
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS

# 🚩 Read-only checks; no inference, model fitting, or publication
if [[ $CHECK == TRUE ]]; then
	if [[ $MODE == all || $MODE == prepare ]]; then run_command "${PREPARE[@]}" --check; fi
	dependencies
	case $MODE in
		all)
			if cache_ready; then run_command "${TRAIN[@]}" --check; else printf 'CHECK: raw inputs passed; prepare the cache to enable the official-model data/forward/backward check.\n'; fi
			;;
		train) need_cache; run_command "${TRAIN[@]}" --check ;;
		predict) need_cache; run_command "${PREDICT[@]}" --check ;;
		report)
			need_cache
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
CACHE_LOCK_KEY=$("$PYTHON" - "$CACHE_DIR" <<'PY'
import hashlib
from pathlib import Path
import sys
canonical_cache = str(Path(sys.argv[1]).expanduser().resolve())
print(hashlib.sha256(canonical_cache.encode()).hexdigest())
PY
)
CACHE_LOCK_FILE="/tmp/grid-prsformer-cache-${CACHE_LOCK_KEY}.lock"
exec {CACHE_LOCK_FD}>"$CACHE_LOCK_FILE"
flock -n "$CACHE_LOCK_FD" || die "Another PRSformer invocation is using this cache: $CACHE_DIR. Wait for it to finish or choose a different --cache-dir."

if [[ $MODE == all || $MODE == train ]]; then
	if [[ $RESUME == FALSE && ( -e $RUN_DIR/model.pt || -e $RUN_DIR/training.pt ) ]]; then
		die "A saved model already exists in $RUN_DIR. Use --resume for an interrupted run, --mode report for a completed run, or a new --run-dir to retrain; --replace does not delete models."
	fi
fi
case $MODE in
	all) run_command "${PREPARE[@]}"; need_cache; run_command "${TRAIN[@]}"; run_command "${PUBLISH[@]}" ;;
	prepare) run_command "${PREPARE[@]}" ;;
	train) need_cache; run_command "${TRAIN[@]}" ;;
	predict) need_cache; run_command "${PREDICT[@]}" ;;
	report) need_cache; run_command "${PUBLISH[@]}" ;;
esac
