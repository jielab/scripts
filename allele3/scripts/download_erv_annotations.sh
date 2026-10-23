#!/usr/bin/env bash
# Shared all-repeat joined annotation. Atomic publication; safe concurrent reuse.
set -euo pipefail
OUT=${1:-/mnt/f/annot/repeats/hg38}
mkdir -p "$OUT"
exec 9>"$OUT/.joined-download.lock"
flock -x 9
FILE="$OUT/rmskJoinedCurrent.txt.gz"
URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmskJoinedCurrent.txt.gz
if [[ -s $FILE ]]; then gzip -t "$FILE"; printf '%s\n' "$FILE"; exit 0; fi
# Each attempt re-detects the current length, retaining progress after truncation.
for attempt in {1..9}; do
    if curl --fail --location --connect-timeout 30 --speed-time 60 --speed-limit 1024 --continue-at - -o "$FILE.partial" "$URL"; then break; fi
    if (( attempt == 9 )); then exit 1; fi
    sleep 2
done
gzip -t "$FILE.partial"
mv "$FILE.partial" "$FILE"
sha256sum "$FILE" > "$FILE.sha256"
printf '%s\n' "$FILE"
