#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$OUTPUT_DIR/analysis" "$PFI_ANALYSIS_DATA_DIR" "$LOG_DIR"
"$PFI_RSCRIPT_BIN" 01_run_pfi_analysis.R \
  "$PFI_FAIRNESS_RESULTS_CSV" \
  "$OUTPUT_DIR/analysis" \
  "$PFI_ANALYSIS_DATA_DIR" \
  2>&1 | tee "$LOG_DIR/10_run_pfi_analysis_R.log"
