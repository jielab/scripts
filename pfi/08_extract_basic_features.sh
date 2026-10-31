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
    Default  Extract article features

Main options:
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./08_extract_basic_features.sh
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
