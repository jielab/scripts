#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"
n="${1:-1000}"
"$PFI_PYTHON_BIN" -u f/sample_llm_review_batch.py \
  --blinded "$PFI_BLINDED_TEXT" \
  --features "$PFI_ARTICLE_FEATURES" \
  --out "$DAT_DIR/llm_review_sample.jsonl" \
  --n "$n" 2>&1 | tee "$LOG_DIR/11_sample_llm_review_batch.log"
