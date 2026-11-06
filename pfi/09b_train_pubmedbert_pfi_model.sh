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
    Default  Train PubMedBERT

Main options:
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./09b_train_pubmedbert_pfi_model.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 09b_train_pubmedbert_pfi_model
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$MODEL_DIR" "$DAT_DIR" "$LOG_DIR"
fp16_arg=()
[[ "$PFI_BERT_FP16" == "1" ]] && fp16_arg=(--fp16)
"$PFI_PYTHON_BIN" -u f/train_pubmedbert_pfi.py \
	--training "$PFI_TRAINING_DATASET" \
	--model_name "$PFI_PUBMEDBERT_MODEL" \
	--fallback_model_name "$PFI_PUBMEDBERT_FALLBACK" \
	--model_dir "$MODEL_DIR/pubmedbert" \
	--metrics_out "$DAT_DIR/pubmedbert_training_metrics.json" \
	--valid_pred_out "$DAT_DIR/pubmedbert_validation_predictions.csv" \
	--fit_summary_out "$DAT_DIR/pubmedbert_training_fit_summary.csv" \
	--epochs "$PFI_BERT_EPOCHS" \
	--batch_size "$PFI_BERT_BATCH_SIZE" \
	--grad_accum "$PFI_BERT_GRAD_ACCUM" \
	--max_length "$PFI_BERT_MAX_LENGTH" \
	--lr "$PFI_BERT_LR" \
	--seed "$PFI_REVIEW_SEED" \
	"${fp16_arg[@]}" \
	2>&1 | tee "$LOG_DIR/09b_train_pubmedbert_pfi_model.log"
