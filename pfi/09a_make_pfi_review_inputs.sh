#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"
: "${PFI_HIGH_ANCHOR_CSV:=}"
: "${PFI_LOW_ANCHOR_CSV:=}"
extra_args=()
[[ -n "$PFI_HIGH_ANCHOR_CSV" && -s "$PFI_HIGH_ANCHOR_CSV" ]] && extra_args+=(--high_anchor_csv "$PFI_HIGH_ANCHOR_CSV")
[[ -n "$PFI_LOW_ANCHOR_CSV" && -s "$PFI_LOW_ANCHOR_CSV" ]] && extra_args+=(--low_anchor_csv "$PFI_LOW_ANCHOR_CSV")
if [[ -z "${PFI_REVIEW_MAX_SCAN_ROWS:-}" ]]; then
  PFI_REVIEW_MAX_SCAN_ROWS=$(( PFI_REVIEW_ROWS * 2000 ))
  if (( PFI_REVIEW_MAX_SCAN_ROWS < 200000 )); then
    PFI_REVIEW_MAX_SCAN_ROWS=200000
  fi
fi
"$PFI_PYTHON_BIN" -u f/make_pfi_review_inputs_v2.py \
  --master "$PFI_ARTICLES_MASTER" \
  --blinded "$PFI_BLINDED_TEXT" \
  --features "$PFI_ARTICLE_FEATURES" \
  --journal_rank "$PFI_JOURNAL_RANK_NORMALIZED" \
  --out_csv "$PFI_PAPER_SCORE_REVIEW_CSV" \
  --out_xlsx "$PFI_PAPER_SCORE_REVIEW_XLSX" \
  --key_out "$PFI_PAPER_SCORE_REVIEW_KEY" \
  --summary_out "$PFI_PAPER_SCORE_REVIEW_SUMMARY" \
  --n "$PFI_REVIEW_ROWS" \
  --min_text_chars "$PFI_MIN_TEXT_CHARS" \
  --max_scan_rows "$PFI_REVIEW_MAX_SCAN_ROWS" \
  --threads "${PFI_DUCKDB_THREADS:-4}" \
  --seed "$PFI_REVIEW_SEED" \
  "${extra_args[@]}" \
  2>&1 | tee "$LOG_DIR/09a_make_pfi_review_inputs.log"
