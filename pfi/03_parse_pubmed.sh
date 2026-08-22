#!/usr/bin/env bash
set -euo pipefail




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Shell options, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$PARQUET_DIR/pubmed/baseline" "$PARQUET_DIR/pubmed/updatefiles" "$LOG_DIR"

max_files=0
if [[ "$PFI_TEST_MODE" == "1" ]]; then
  max_files=3
  echo "PFI_TEST_MODE=1: parsing only first 3 baseline/update files."
fi

"$PFI_PYTHON_BIN" f/parse_pubmed.py \
  --in_dir "$PUBMED_BASELINE_DIR" \
  --out_dir "$PARQUET_DIR/pubmed/baseline" \
  --min_year "$MIN_YEAR" \
  --workers "$N_JOBS" \
  --max_files "$max_files" 2>&1 | tee "$LOG_DIR/03_parse_pubmed_baseline.log"

"$PFI_PYTHON_BIN" f/parse_pubmed.py \
  --in_dir "$PUBMED_UPDATE_DIR" \
  --out_dir "$PARQUET_DIR/pubmed/updatefiles" \
  --min_year "$MIN_YEAR" \
  --workers "$N_JOBS" \
  --max_files "$max_files" 2>&1 | tee "$LOG_DIR/03_parse_pubmed_updatefiles.log"
