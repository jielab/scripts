#!/usr/bin/env bash
set -euo pipefail


# 🚩 Shell options, paths, and shared inputs
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$PFI_WORK_PARQUET_DIR" "$LOG_DIR"

{
  echo "[$(date -Is)] 08_extract_basic_features.sh started"
  "$PFI_PYTHON_BIN" -u f/extract_basic_features.py \
    --blinded "$PFI_BLINDED_TEXT" \
    --out "$PFI_ARTICLE_FEATURES" \
    --workers "${PFI_FEATURE_WORKERS:-$N_JOBS}" \
    --progress_every_rows "${PFI_FEATURE_PROGRESS_EVERY_ROWS:-250000}"
  echo "[$(date -Is)] 08_extract_basic_features.sh finished"
} 2>&1 | tee "$LOG_DIR/08_extract_basic_features.log"
