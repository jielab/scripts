#!/usr/bin/env bash
set -euo pipefail


# 🚩 Shell options, paths, and shared inputs
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
check_raw_mounts
mkdir -p "$PMC_OA_RAW_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/04_download_pmc_oa_bulk_xml.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "===== 04_download_pmc_oa_bulk_xml.sh started at $(date -Is) ====="
echo "Logging to: $LOG_FILE"

ARIA2_MAX_CONCURRENT_DOWNLOADS="${ARIA2_MAX_CONCURRENT_DOWNLOADS:-6}"
ARIA2_MAX_CONNECTION_PER_SERVER="${ARIA2_MAX_CONNECTION_PER_SERVER:-4}"
ARIA2_SPLIT="${ARIA2_SPLIT:-4}"
ARIA2_MIN_SPLIT_SIZE="${ARIA2_MIN_SPLIT_SIZE:-20M}"

# PMC OA bulk XML baseline files. Groups: commercial, non-commercial, and other licenses.
# NCBI moved archived PMC FTP files under deprecated/ on 2026-04-13 and plans to
# remove these bulk packages in 2026-08. After that, switch this step to AWS.
BASE="${PMC_OA_BULK_BASE:-https://ftp.ncbi.nlm.nih.gov/pub/pmc/deprecated/oa_bulk}"
: > "$LOG_DIR/pmc_oa_bulk_xml_urls.txt"
MISSING_URLS="$LOG_DIR/pmc_oa_bulk_xml_urls_missing.txt"

limit=0
if [[ "$PFI_TEST_MODE" == "1" ]]; then
  limit=10
  echo "Test mode: downloading only the first 10 PMC OA bulk XML files."
fi

for group in oa_comm oa_noncomm oa_other; do
  url="$BASE/$group/xml/"
  echo "Creating PMC OA URL list from: $url"
  "$PFI_PYTHON_BIN" f/make_url_list_from_ncbi_dir.py \
    --url "$url" \
    --pattern '.*\.tar\.gz$' \
    --out "$LOG_DIR/pmc_${group}_xml_urls.txt" \
    --limit "$limit"
  cat "$LOG_DIR/pmc_${group}_xml_urls.txt" >> "$LOG_DIR/pmc_oa_bulk_xml_urls.txt"
done

n=$(wc -l < "$LOG_DIR/pmc_oa_bulk_xml_urls.txt")
if [[ "$n" -eq 0 ]]; then
  echo "ERROR: no PMC OA XML tar.gz URLs found." >&2
  exit 1
fi

"$PFI_PYTHON_BIN" - "$LOG_DIR/pmc_oa_bulk_xml_urls.txt" "$PMC_OA_RAW_DIR" "$MISSING_URLS" <<'PY'
import os
import sys
from urllib.parse import urlparse

listfile, outdir, missing_listfile = sys.argv[1:]
missing = []

with open(listfile, encoding="utf-8") as f:
    for line in f:
        url = line.strip()
        if not url:
            continue
        name = os.path.basename(urlparse(url).path)
        target = os.path.join(outdir, name)
        aria2 = target + ".aria2"
        if not os.path.exists(target) or os.path.exists(aria2):
            missing.append(url)

with open(missing_listfile, "w", encoding="utf-8") as f:
    for url in missing:
        f.write(url + "\n")
PY

n_missing=$(wc -l < "$MISSING_URLS")
downloaded_before=$(find "$PMC_OA_RAW_DIR" -maxdepth 1 -name '*.tar.gz' | wc -l)
incomplete_before=$(find "$PMC_OA_RAW_DIR" -maxdepth 1 -name '*.aria2' | wc -l)

echo "PMC OA pre-download status:"
echo "  URL list: $n"
echo "  Directory has $downloaded_before .tar.gz files"
echo "  Directory has $incomplete_before .aria2 resume files"
echo "  Need download/resume now: $n_missing"

if [[ "$n_missing" -eq 0 ]]; then
  rm -f "$MISSING_URLS"
  echo "OK: PMC OA download complete for $PMC_OA_RAW_DIR."
  exit 0
fi

echo "Downloading/resuming $n_missing PMC OA XML tar.gz files to: $PMC_OA_RAW_DIR"
echo "aria2 parallel settings:"
echo "  max concurrent downloads: $ARIA2_MAX_CONCURRENT_DOWNLOADS"
echo "  max connections per server: $ARIA2_MAX_CONNECTION_PER_SERVER"
echo "  split per file: $ARIA2_SPLIT"
echo "  min split size: $ARIA2_MIN_SPLIT_SIZE"
download_status=0
if command -v aria2c >/dev/null 2>&1; then
  aria2_args=(
    --continue=true
    --max-concurrent-downloads="$ARIA2_MAX_CONCURRENT_DOWNLOADS"
    --max-connection-per-server="$ARIA2_MAX_CONNECTION_PER_SERVER"
    --split="$ARIA2_SPLIT"
    --min-split-size="$ARIA2_MIN_SPLIT_SIZE"
    --dir="$PMC_OA_RAW_DIR"
    --input-file="$MISSING_URLS"
    --allow-overwrite=false
    --auto-file-renaming=false
    --summary-interval=60
    --console-log-level=warn
    --max-tries=0
    --retry-wait=60
    --timeout=60
    --connect-timeout=30
    --log="$LOG_DIR/04_download_pmc_oa_bulk_xml.aria2.log"
    --log-level=notice
  )
  set +e
  aria2c "${aria2_args[@]}"
  download_status=$?
  set -e
else
  wget_args=(
    --continue
    --no-use-server-timestamps
    --directory-prefix="$PMC_OA_RAW_DIR"
    --input-file="$MISSING_URLS"
  )
  set +e
  wget "${wget_args[@]}"
  download_status=$?
  set -e
fi

echo "Downloaded PMC OA tarballs:"
downloaded=$(find "$PMC_OA_RAW_DIR" -maxdepth 1 -name '*.tar.gz' | wc -l)
incomplete=$(find "$PMC_OA_RAW_DIR" -maxdepth 1 -name '*.aria2' | wc -l)
"$PFI_PYTHON_BIN" - "$LOG_DIR/pmc_oa_bulk_xml_urls.txt" "$PMC_OA_RAW_DIR" "$MISSING_URLS" <<'PY'
import os
import sys
from urllib.parse import urlparse

listfile, outdir, missing_listfile = sys.argv[1:]
missing = []

with open(listfile, encoding="utf-8") as f:
    for line in f:
        url = line.strip()
        if not url:
            continue
        name = os.path.basename(urlparse(url).path)
        target = os.path.join(outdir, name)
        aria2 = target + ".aria2"
        if not os.path.exists(target) or os.path.exists(aria2):
            missing.append(url)

with open(missing_listfile, "w", encoding="utf-8") as f:
    for url in missing:
        f.write(url + "\n")
PY
remaining=$(wc -l < "$MISSING_URLS")
echo "$downloaded"
echo "Download status for $PMC_OA_RAW_DIR:"
echo "  URL list: $n"
echo "  Directory has $downloaded .tar.gz files"
echo "  Directory has $incomplete .aria2 resume files"
echo "  Still need download/resume: $remaining"

if [[ "$download_status" -ne 0 || "$remaining" -ne 0 || "$incomplete" -ne 0 ]]; then
  echo "ERROR: PMC OA download is incomplete." >&2
  echo "URL list has $n files, directory has $downloaded .tar.gz files, but still has $incomplete .aria2 resume files and $remaining files needing download/resume." >&2
  echo "Re-run this script to resume: ./04_download_pmc_oa_bulk_xml.sh" >&2
  if [[ "$incomplete" -gt 0 ]]; then
    echo "First incomplete files:" >&2
    find "$PMC_OA_RAW_DIR" -maxdepth 1 -name '*.aria2' -printf '%f\n' | sort | head -20 >&2
  fi
  echo "First missing/incomplete URLs are listed in: $MISSING_URLS" >&2
  head -20 "$MISSING_URLS" >&2
  if [[ "$download_status" -eq 0 ]]; then
    exit 1
  fi
  exit "$download_status"
fi

rm -f "$MISSING_URLS"
echo "OK: PMC OA download complete for $PMC_OA_RAW_DIR."
