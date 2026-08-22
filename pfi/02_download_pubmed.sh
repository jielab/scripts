#!/usr/bin/env bash
set -euo pipefail




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Shell options, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
check_raw_mounts

mkdir -p "$PUBMED_BASELINE_DIR" "$PUBMED_UPDATE_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/02_download_pubmed.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "===== 02_download_pubmed.sh started at $(date -Is) ====="
echo "Logging to: $LOG_FILE"

# NCBI FTP/HTTPS can intermittently close TLS handshakes when many files are
# opened in parallel. Use conservative defaults; override from the shell if needed.
ARIA2_MAX_CONCURRENT_DOWNLOADS="${ARIA2_MAX_CONCURRENT_DOWNLOADS:-2}"
ARIA2_MAX_CONNECTION_PER_SERVER="${ARIA2_MAX_CONNECTION_PER_SERVER:-1}"
ARIA2_SPLIT="${ARIA2_SPLIT:-1}"
ARIA2_MIN_SPLIT_SIZE="${ARIA2_MIN_SPLIT_SIZE:-20M}"
PFI_DOWNLOAD_TOOL="${PFI_DOWNLOAD_TOOL:-aria2}"
PFI_DOWNLOAD_MAX_TRIES="${PFI_DOWNLOAD_MAX_TRIES:-20}"
PFI_DOWNLOAD_RETRY_WAIT="${PFI_DOWNLOAD_RETRY_WAIT:-30}"
PFI_DOWNLOAD_TIMEOUT="${PFI_DOWNLOAD_TIMEOUT:-60}"
PFI_DOWNLOAD_CONNECT_TIMEOUT="${PFI_DOWNLOAD_CONNECT_TIMEOUT:-30}"

download_dir () {
  local url="$1"
  local outdir="$2"
  local listfile="$3"
  local limit="$4"
  local missing_listfile="${listfile%.txt}_missing.txt"
  local pattern='pubmed*.xml.gz'

  mkdir -p "$outdir"
  echo "Creating URL list from: $url"

  "$PFI_PYTHON_BIN" f/make_url_list_from_ncbi_dir.py \
    --url "$url" \
    --pattern 'pubmed.*\.xml\.gz$' \
    --out "$listfile" \
    --limit "$limit"

  local n
  n=$(wc -l < "$listfile")
  if [[ "$n" -eq 0 ]]; then
    echo "ERROR: no pubmed*.xml.gz URLs found from $url" >&2
    exit 1
  fi

  "$PFI_PYTHON_BIN" - "$listfile" "$outdir" "$missing_listfile" <<'PY'
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

  local n_missing
  n_missing=$(wc -l < "$missing_listfile")
  echo "Found $n URLs; $n_missing need download/resume in $outdir"
  if [[ "$n_missing" -eq 0 ]]; then
    local downloaded
    local incomplete
    downloaded=$(find "$outdir" -maxdepth 1 -name "$pattern" | wc -l)
    incomplete=$(find "$outdir" -maxdepth 1 -name '*.aria2' | wc -l)
    rm -f "$missing_listfile"
    echo "Download status for $outdir:"
    echo "  URL list: $n"
    echo "  Directory has $downloaded .xml.gz files"
    echo "  Directory has $incomplete .aria2 resume files"
    echo "  Still need download/resume: 0"
    echo "OK: PubMed download complete for $outdir."
    return 0
  fi

  echo "Downloading/resuming $n_missing files to $outdir"
  echo "aria2 parallel settings:"
  echo "  max concurrent downloads: $ARIA2_MAX_CONCURRENT_DOWNLOADS"
  echo "  max connections per server: $ARIA2_MAX_CONNECTION_PER_SERVER"
  echo "  split per file: $ARIA2_SPLIT"
  echo "  min split size: $ARIA2_MIN_SPLIT_SIZE"
  local download_status=0
  if [[ "$PFI_DOWNLOAD_TOOL" == "aria2" && -n "$(command -v aria2c || true)" ]]; then
    local aria2_args=(
      --continue=true
      --max-concurrent-downloads="$ARIA2_MAX_CONCURRENT_DOWNLOADS"
      --max-connection-per-server="$ARIA2_MAX_CONNECTION_PER_SERVER"
      --split="$ARIA2_SPLIT"
      --min-split-size="$ARIA2_MIN_SPLIT_SIZE"
      --dir="$outdir"
      --input-file="$missing_listfile"
      --allow-overwrite=false
      --auto-file-renaming=false
      --file-allocation=none
      --max-tries="$PFI_DOWNLOAD_MAX_TRIES"
      --retry-wait="$PFI_DOWNLOAD_RETRY_WAIT"
      --timeout="$PFI_DOWNLOAD_TIMEOUT"
      --connect-timeout="$PFI_DOWNLOAD_CONNECT_TIMEOUT"
      --lowest-speed-limit=1K
      --summary-interval=60
      --console-log-level=warn
      --log="${missing_listfile%.txt}.aria2.log"
      --log-level=notice
    )
    set +e
    aria2c "${aria2_args[@]}"
    download_status=$?
    set -e
  else
    echo "Using wget downloader. To force this mode, run: PFI_DOWNLOAD_TOOL=wget ./02_download_pubmed.sh"
    local wget_args=(
      --continue
      --tries="$PFI_DOWNLOAD_MAX_TRIES"
      --waitretry="$PFI_DOWNLOAD_RETRY_WAIT"
      --timeout="$PFI_DOWNLOAD_TIMEOUT"
      --no-use-server-timestamps
      --directory-prefix="$outdir"
      --input-file="$missing_listfile"
    )
    set +e
    wget "${wget_args[@]}"
    download_status=$?
    set -e
  fi

  "$PFI_PYTHON_BIN" - "$listfile" "$outdir" "$missing_listfile" <<'PY'
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

  local downloaded
  local incomplete
  local remaining
  downloaded=$(find "$outdir" -maxdepth 1 -name "$pattern" | wc -l)
  incomplete=$(find "$outdir" -maxdepth 1 -name '*.aria2' | wc -l)
  remaining=$(wc -l < "$missing_listfile")

  echo "Download status for $outdir:"
  echo "  URL list: $n"
  echo "  Directory has $downloaded .xml.gz files"
  echo "  Directory has $incomplete .aria2 resume files"
  echo "  Still need download/resume: $remaining"

  if [[ "$download_status" -ne 0 || "$remaining" -ne 0 || "$incomplete" -ne 0 ]]; then
    echo "ERROR: PubMed download is incomplete for $outdir." >&2
    echo "Re-run this script to resume: ./02_download_pubmed.sh" >&2
    if [[ "$remaining" -gt 0 ]]; then
      echo "First missing/incomplete URLs are listed in: $missing_listfile" >&2
      head -20 "$missing_listfile" >&2
    fi
    if [[ "$download_status" -eq 0 ]]; then
      return 1
    fi
    return "$download_status"
  fi

  rm -f "$missing_listfile"
  echo "OK: PubMed download complete for $outdir."
}

limit=0
if [[ "$PFI_TEST_MODE" == "1" ]]; then
  limit=3
  echo "PFI_TEST_MODE=1: downloading only first 3 baseline/update files."
fi

echo "Downloading PubMed baseline XML files to: $PUBMED_BASELINE_DIR"
download_dir "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/" \
  "$PUBMED_BASELINE_DIR" "$LOG_DIR/pubmed_baseline_urls.txt" "$limit"

echo "Downloading PubMed daily update XML files to: $PUBMED_UPDATE_DIR"
download_dir "https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/" \
  "$PUBMED_UPDATE_DIR" "$LOG_DIR/pubmed_update_urls.txt" "$limit"

echo "Downloaded baseline files:"
find "$PUBMED_BASELINE_DIR" -name 'pubmed*.xml.gz' | wc -l

echo "Downloaded update files:"
find "$PUBMED_UPDATE_DIR" -name 'pubmed*.xml.gz' | wc -l
