#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$MODEL_DIR" "$DAT_DIR" "$LOG_DIR"
fp16_arg=()
[[ "$PFI_BERT_FP16" == "1" ]] && fp16_arg=(--fp16)
"$PFI_PYTHON_BIN" -u f/train_pubmedbert_pfi_v2.py \
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
