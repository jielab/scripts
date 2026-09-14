#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$PFI_ANALYSIS_DATA_DIR" "$LOG_DIR"
"$PFI_PYTHON_BIN" -u f/audit_outliers_author_institution.py \
  --fairness "$PFI_FAIRNESS_RESULTS_PARQUET" \
  --out_dir "$PFI_ANALYSIS_DATA_DIR/outlier_audit" \
  2>&1 | tee "$LOG_DIR/12_audit_strong_outliers.log"
