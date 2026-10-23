#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Default  Parse downloaded PubMed XML

Main options:
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./03_parse_pubmed.sh
HELP
			exit 0
			;;
	esac
fi
set -euo pipefail


# 🚩 Shell options, paths, and shared inputs
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$PARQUET_DIR/pubmed/baseline" "$PARQUET_DIR/pubmed/updatefiles" "$LOG_DIR"

max_files=0
if [[ "$PFI_TEST_MODE" == "1" ]]; then
	max_files=3
	echo "Test mode: parsing only the first 3 baseline/update files."
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
