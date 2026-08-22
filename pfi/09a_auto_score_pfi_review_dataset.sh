#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"

: "${PFI_ARTICLE_METRICS_CSV:=}"
: "${PFI_HIGH_ANCHOR_CSV:=}"
: "${PFI_LOW_ANCHOR_CSV:=}"
: "${PFI_AUTO_SCORE_MODE:=weak_supervision}"

in_csv="${1:-$PFI_PAPER_SCORE_REVIEW_CSV}"
out_csv="${2:-$PFI_PAPER_SCORE_REVIEW_PSEUDO}"

if [[ ! -s "$in_csv" ]]; then
  echo "ERROR: blinded review CSV is missing: $in_csv" >&2
  echo "Run: ./09a_make_pfi_review_inputs.sh" >&2
  exit 1
fi

extra_args=()
[[ -n "$PFI_ARTICLE_METRICS_CSV" && -s "$PFI_ARTICLE_METRICS_CSV" ]] && extra_args+=(--article_metrics "$PFI_ARTICLE_METRICS_CSV")
[[ -n "$PFI_HIGH_ANCHOR_CSV" && -s "$PFI_HIGH_ANCHOR_CSV" ]] && extra_args+=(--high_anchor_csv "$PFI_HIGH_ANCHOR_CSV")
[[ -n "$PFI_LOW_ANCHOR_CSV" && -s "$PFI_LOW_ANCHOR_CSV" ]] && extra_args+=(--low_anchor_csv "$PFI_LOW_ANCHOR_CSV")

"$PFI_PYTHON_BIN" -u f/auto_score_pfi_review_dataset.py \
  --review_csv "$in_csv" \
  --review_key "$PFI_PAPER_SCORE_REVIEW_KEY" \
  --features "$PFI_ARTICLE_FEATURES" \
  --out_csv "$out_csv" \
  --audit_csv "$PFI_PAPER_SCORE_REVIEW_PSEUDO_AUDIT" \
  --mode "$PFI_AUTO_SCORE_MODE" \
  "${extra_args[@]}" \
  2>&1 | tee "$LOG_DIR/09a_auto_score_pfi_review_dataset.log"

cat >&2 <<EOF
OK: pseudo-scored review file written to:
  $out_csv
Use it for training if manual expert labels are unavailable:
  bash 09a_prepare_pfi_training_dataset.sh $PFI_JOURNAL_RANK_NORMALIZED $out_csv
EOF
