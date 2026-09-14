#!/usr/bin/env bash
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
