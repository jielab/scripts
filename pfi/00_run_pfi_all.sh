#!/usr/bin/env bash
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
