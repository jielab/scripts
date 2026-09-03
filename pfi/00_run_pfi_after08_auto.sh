#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR" "$OUTPUT_DIR/analysis" "$PFI_ANALYSIS_DATA_DIR"

LOG_FILE="$LOG_DIR/00_run_pfi_after08_auto.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "===== 00_run_pfi_after08_auto.sh started at $(date -Is) ====="
echo "This starts after 06/07 have produced articles_master and blinded_text."
echo "If no expert scored CSV exists, weak-supervision pseudo labels will be generated automatically."
echo "If no real journal-rank file exists, a temporary proxy rank will be generated from articles_master."

if [[ ! -e "$PFI_ARTICLES_MASTER" ]]; then
  echo "ERROR: articles_master is missing: $PFI_ARTICLES_MASTER" >&2
  echo "Run 06_build_master_dataset.sh first." >&2
  exit 1
fi
if [[ ! -e "$PFI_BLINDED_TEXT" ]]; then
  echo "ERROR: blinded_text is missing: $PFI_BLINDED_TEXT" >&2
  echo "Run 07_make_blinded_text.sh first." >&2
  exit 1
fi

# Features are needed for review sampling and pseudo labels.
if [[ ! -e "$PFI_ARTICLE_FEATURES" ]]; then
  echo "INFO: article features are missing; running 08_extract_basic_features.sh"
  bash "$PFI_ROOT/08_extract_basic_features.sh"
fi

bash "$PFI_ROOT/08b_normalize_journal_rank.sh" "${1:-}"
bash "$PFI_ROOT/09a_make_pfi_review_inputs.sh"
bash "$PFI_ROOT/09a_auto_score_pfi_review_dataset.sh"
bash "$PFI_ROOT/09a_prepare_pfi_training_dataset.sh" "$PFI_JOURNAL_RANK_NORMALIZED" "$PFI_PAPER_SCORE_REVIEW_PSEUDO"
bash "$PFI_ROOT/09b_train_pubmedbert_pfi_model.sh"
bash "$PFI_ROOT/09c_predict_pubmedbert_pfi.sh"
bash "$PFI_ROOT/09d_compute_pfi_index.sh" "$PFI_JOURNAL_RANK_NORMALIZED"
bash "$PFI_ROOT/10_run_pfi_analysis_R.sh"
bash "$PFI_ROOT/12_audit_strong_outliers.sh"

echo "===== 00_run_pfi_after08_auto.sh finished at $(date -Is) ====="
echo "Main outputs:"
echo "  $PFI_TRAINING_DATASET"
echo "  $PFI_PUBMEDBERT_PREDICTIONS"
echo "  $PFI_FAIRNESS_RESULTS_CSV"
echo "  $OUTPUT_DIR/analysis"
