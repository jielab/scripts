#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$PFI_ANALYSIS_DATA_DIR" "$LOG_DIR"
rank_file="${1:-$PFI_JOURNAL_RANK_NORMALIZED}"
institution_rank_file="${2:-}"
args=()
[[ "$PFI_RESULTS_FULL_TEXT_ONLY" == "1" ]] && args+=(--full_text_only)
[[ -n "$institution_rank_file" && -s "$institution_rank_file" ]] && args+=(--institution_rank "$institution_rank_file")
"$PFI_PYTHON_BIN" -u f/compute_pfi_index.py \
  --master "$PFI_ARTICLES_MASTER" \
  --predictions "$PFI_PUBMEDBERT_PREDICTIONS" \
  --journal_rank "$rank_file" \
  --out_csv "$PFI_FAIRNESS_RESULTS_CSV" \
  --out_parquet "$PFI_FAIRNESS_RESULTS_PARQUET" \
  --summary_dir "$PFI_ANALYSIS_DATA_DIR" \
  "${args[@]}" \
  2>&1 | tee "$LOG_DIR/09d_compute_pfi_index.log"
