#!/usr/bin/env bash
set -euo pipefail


# 🚩 Shell options, paths, and shared inputs
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$PFI_WORK_PARQUET_DIR" "$PFI_WORK_DUCKDB_TEMP_DIR" "$LOG_DIR"

full_text_arg=()
: "${PFI_BLINDED_FULL_TEXT_ONLY:=1}"
if [[ "$PFI_BLINDED_FULL_TEXT_ONLY" == "1" ]]; then
  full_text_arg=(--full_text_only)
fi
memory_limit_arg=()
if [[ -n "${PFI_DUCKDB_MEMORY_LIMIT:-}" ]]; then
  memory_limit_arg=(--memory_limit "$PFI_DUCKDB_MEMORY_LIMIT")
fi
resume_arg=()
: "${PFI_BLINDED_RESUME:=1}"
if [[ "$PFI_BLINDED_RESUME" == "1" ]]; then
  resume_arg=(--resume)
fi

{
  echo "[$(date -Is)] 07_make_blinded_text.sh started"
  "$PFI_PYTHON_BIN" -u f/make_blinded_text.py \
    --master "$PFI_ARTICLES_MASTER" \
    --out "$PFI_BLINDED_TEXT" \
    --temp_dir "$PFI_WORK_DUCKDB_TEMP_DIR" \
    --threads "${PFI_DUCKDB_THREADS:-4}" \
    "${memory_limit_arg[@]}" \
    --engine "${PFI_BLINDED_ENGINE:-pyarrow}" \
    --batch_size "${PFI_BLINDED_BATCH_SIZE:-1000}" \
    --rows_per_file "${PFI_BLINDED_ROWS_PER_FILE:-100000}" \
    --progress_every_rows "${PFI_BLINDED_PROGRESS_EVERY_ROWS:-250000}" \
    "${resume_arg[@]}" \
    "${full_text_arg[@]}" \
    --partition_by_year
  echo "[$(date -Is)] 07_make_blinded_text.sh finished"
} 2>&1 | tee "$LOG_DIR/07_make_blinded_text.log"
