#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"

rank_file="${1:-$PFI_JOURNAL_RANK_NORMALIZED}"
paper_score_file="${2:-$PFI_PAPER_SCORE_REVIEW_SCORED}"

# Ensure journal rank exists. If no real rank was provided, 08b can create a temporary proxy rank.
if [[ ! -s "$rank_file" ]]; then
  echo "INFO: normalized journal-rank file is missing; running 08b_normalize_journal_rank.sh" >&2
  bash "$PFI_ROOT/08b_normalize_journal_rank.sh"
  rank_file="$PFI_JOURNAL_RANK_NORMALIZED"
fi

# Ensure review input exists.
if [[ ! -s "$PFI_PAPER_SCORE_REVIEW_CSV" || ! -s "$PFI_PAPER_SCORE_REVIEW_KEY" ]]; then
  echo "INFO: blinded review inputs are missing; running 09a_make_pfi_review_inputs.sh" >&2
  bash "$PFI_ROOT/09a_make_pfi_review_inputs.sh"
fi

# Prefer human/expert labels if present. If missing, automatically create weak-supervision labels.
if [[ ! -s "$paper_score_file" ]]; then
  pseudo_file="$PFI_PAPER_SCORE_REVIEW_PSEUDO"
  if [[ ! -s "$pseudo_file" ]]; then
    echo "INFO: scored review CSV is missing; creating weak-supervision pseudo labels automatically." >&2
    bash "$PFI_ROOT/09a_auto_score_pfi_review_dataset.sh" "$PFI_PAPER_SCORE_REVIEW_CSV" "$pseudo_file"
  fi
  paper_score_file="$pseudo_file"
  echo "INFO: using weak-supervision pseudo labels for training: $paper_score_file" >&2
fi

"$PFI_PYTHON_BIN" -u f/make_pfi_training_dataset_v2.py \
  --review_scores "$paper_score_file" \
  --review_key "$PFI_PAPER_SCORE_REVIEW_KEY" \
  --blinded "$PFI_BLINDED_TEXT" \
  --features "$PFI_ARTICLE_FEATURES" \
  --journal_rank "$rank_file" \
  --master "$PFI_ARTICLES_MASTER" \
  --out "$PFI_TRAINING_DATASET" \
  --summary_out "$DAT_DIR/training_dataset_summary.csv" \
  --fit_summary_out "$DAT_DIR/training_label_fit_summary.csv" \
  --audit_csv_out "$DAT_DIR/training_records_preview.csv" \
  --max_rows "$PFI_MAX_TRAIN_ROWS" \
  --min_text_chars "$PFI_MIN_TEXT_CHARS" \
  --seed "$PFI_REVIEW_SEED" \
  2>&1 | tee "$LOG_DIR/09a_prepare_pfi_training_dataset.log"
