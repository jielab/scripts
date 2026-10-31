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
    Default  Run the complete PFI workflow

Main options:
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./00_run_pfi_all.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 00_run_pfi_all
set -euo pipefail

bash 02_download_pubmed.sh
bash 03_parse_pubmed.sh
bash 04_download_pmc_oa_bulk_xml.sh
bash 05_parse_pmc_oa_bulk_xml.sh
bash 06_build_master_dataset.sh
bash 07_make_blinded_text.sh
bash 08_extract_basic_features.sh
rank_args=()
if [[ $# -gt 0 && -n "${1:-}" ]]; then
	rank_args=("$1")
fi
bash 08b_normalize_journal_rank.sh "${rank_args[@]}"
bash 09a_make_pfi_review_inputs.sh
cat <<'MSG'
STOP: Fill the blinded review XLSX, save it as paper_score_review_1000_scored.csv, then run:
  ./pfi.sh prepare-training
  ./pfi.sh train
  ./pfi.sh predict
  ./pfi.sh compute
  ./pfi.sh analysis
MSG
