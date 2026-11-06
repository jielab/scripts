#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# This entry initializes LOCAL paths and is not a supported --submit hpc pipeline.
# The bsub example above applies only after its worker paths/runtime are ported.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
submit=${submit:-${SUBMIT:-local}}
source "$_paths_directory/0f/0os.sh"
unset _paths_directory


# Resolve --dir-dat/--dir-big/--dir-tmp through the shared path table.
source "$SCRIPTS_ROOT/0f/local_paths.sh"
local_workflow_init "$@" || exit $?
set -- "${WORKFLOW_ARGS[@]}"

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "${SCRIPTS_ROOT}/0f/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Default  Predict article scores

Main options:
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./09c_predict_pubmedbert_pfi.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 09c_predict_pubmedbert_pfi
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"
: "${PFI_PRED_BATCH_SIZE:=64}"
: "${PFI_PRED_LIMIT:=0}"
: "${PFI_PRED_FRESH_SAMPLE:=0}"
: "${PFI_PRED_FULL_TEXT_ONLY:=1}"
: "${PFI_PRED_PROGRESS_EVERY:=10000}"
: "${PFI_CPU_FULL_WARNING_ROWS:=1000000}"
sample_scope="alltext"
[[ "$PFI_PRED_FULL_TEXT_ONLY" == "1" ]] && sample_scope="fulltext"
: "${PFI_PRED_SAMPLE_PATH:=$DAT_DIR/pubmedbert_prediction_${sample_scope}_sample_${PFI_PRED_SAMPLE_ROWS}_seed${PFI_PRED_SAMPLE_SEED}.parquet}"
args=()
[[ "$PFI_PRED_FULL_TEXT_ONLY" == "1" ]] && args+=(--full_text_only)
[[ "$PFI_PRED_FRESH_SAMPLE" == "1" ]] && args+=(--fresh_sample)
[[ "$PFI_BERT_FP16" == "1" ]] && args+=(--fp16)
[[ "$PFI_ALLOW_CPU_FULL_PRED" == "1" ]] && args+=(--allow_cpu_full)
if [[ "$PFI_PRED_SAMPLE_ROWS" != "0" ]]; then
	args+=(--sample_rows "$PFI_PRED_SAMPLE_ROWS" --sample_seed "$PFI_PRED_SAMPLE_SEED" --sample_path "$PFI_PRED_SAMPLE_PATH")
fi
"$PFI_PYTHON_BIN" -u f/predict_pubmedbert_pfi.py \
	--blinded "$PFI_BLINDED_TEXT" \
	--model_dir "$MODEL_DIR/pubmedbert" \
	--out "$PFI_PUBMEDBERT_PREDICTIONS" \
	--batch_size "$PFI_PRED_BATCH_SIZE" \
	--max_length "$PFI_BERT_MAX_LENGTH" \
	--threads "${PFI_DUCKDB_THREADS:-4}" \
	--limit "$PFI_PRED_LIMIT" \
	--progress_every "$PFI_PRED_PROGRESS_EVERY" \
	--cpu_full_warning_rows "$PFI_CPU_FULL_WARNING_ROWS" \
	"${args[@]}" \
	2>&1 | tee "$LOG_DIR/09c_predict_pubmedbert_pfi.log"
