#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR" "$PFI_ROOT/conf"

: "${PFI_ALLOW_PROXY_JOURNAL_RANK:=0}"
: "${PFI_PROXY_JOURNAL_RANK_FULL_TEXT_ONLY:=1}"

rank_file="${1:-${PFI_JOURNAL_RANK_REAL:-}}"
if [[ -z "${rank_file:-}" ]]; then
  rank_file="$PFI_ROOT/journal_rank_real.csv"
fi

schema_file="$PFI_ROOT/conf/journal_rank_schema.csv"
if [[ -s "$PFI_ROOT/journal_rank_schema.csv" && ! -s "$schema_file" ]]; then
  cp -p "$PFI_ROOT/journal_rank_schema.csv" "$schema_file"
fi

looks_like_schema=0
if [[ -s "$rank_file" && "$(basename "$rank_file")" == "journal_rank_schema.csv" ]]; then
  looks_like_schema=1
fi

need_proxy=0
if [[ ! -s "$rank_file" ]]; then
  need_proxy=1
elif [[ "$looks_like_schema" == "1" ]]; then
  need_proxy=1
fi

if [[ "$need_proxy" == "1" ]]; then
  if [[ "$PFI_ALLOW_PROXY_JOURNAL_RANK" != "1" ]]; then
    echo "ERROR: journal rank input is missing or is only journal_rank_schema.csv: $rank_file" >&2
    echo "Use ./pfi.sh normalize-rank --allow-proxy-journal-rank TRUE only for a temporary test rank." >&2
    exit 1
  fi
  if [[ ! -e "$PFI_ARTICLES_MASTER" ]]; then
    echo "ERROR: journal rank input is missing and articles_master is also missing: $PFI_ARTICLES_MASTER" >&2
    echo "Run 06_build_master_dataset.sh first, or provide a real JCR/SJR rank file." >&2
    exit 1
  fi
  proxy_csv="$DAT_DIR/journal_rank_proxy_from_master.csv"
  proxy_audit="$DAT_DIR/journal_rank_proxy_from_master_audit.csv"
  echo "WARNING: real journal rank input is missing or invalid: $rank_file" >&2
  echo "WARNING: building TEMPORARY proxy journal rank from articles_master so the pipeline can run through." >&2
  echo "WARNING: use this proxy only for debugging/feasibility, not for final scientific interpretation." >&2
  proxy_args=()
  [[ "$PFI_PROXY_JOURNAL_RANK_FULL_TEXT_ONLY" == "1" ]] && proxy_args+=(--full_text_only)
  "$PFI_PYTHON_BIN" -u f/build_proxy_journal_rank_from_master.py \
    --master "$PFI_ARTICLES_MASTER" \
    --out "$proxy_csv" \
    --audit_csv "$proxy_audit" \
    --min_year "${MIN_YEAR:-2000}" \
    "${proxy_args[@]}" \
    2>&1 | tee "$LOG_DIR/08b_build_proxy_journal_rank.log"
  rank_file="$proxy_csv"
fi

# If the user passes a real file, save a durable copy under the configured location.
if [[ -n "${PFI_JOURNAL_RANK_REAL:-}" && -s "$rank_file" && "$rank_file" != "$PFI_JOURNAL_RANK_REAL" && "$rank_file" != "$DAT_DIR/journal_rank_proxy_from_master.csv" ]]; then
  mkdir -p "$(dirname "$PFI_JOURNAL_RANK_REAL")"
  cp -p "$rank_file" "$PFI_JOURNAL_RANK_REAL"
  echo "INFO: copied journal rank source to stable location: $PFI_JOURNAL_RANK_REAL" >&2
  rank_file="$PFI_JOURNAL_RANK_REAL"
fi

"$PFI_PYTHON_BIN" -u f/normalize_journal_rank.py \
  --rank "$rank_file" \
  --schema "$schema_file" \
  --out "$PFI_JOURNAL_RANK_NORMALIZED" \
  --audit_csv "$DAT_DIR/journal_rank_normalized_audit.csv" \
  2>&1 | tee "$LOG_DIR/08b_normalize_journal_rank.log"
