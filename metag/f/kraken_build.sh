#!/usr/bin/env bash
set -euo pipefail
MAX_GB="${1:-72}"
ROOT="$(python3 -c 'from common import DIR_REF; print(DIR_REF)')"
OUT="$(python3 -c 'from common import DIR_LIB; print(DIR_LIB / "kraken2")')"
DB="$OUT/db"
mkdir -p "$DB"

echo "Kraken2 comparator only. Reference FASTA will NOT be deleted."
echo "max-db-size=${MAX_GB} GiB"

command -v kraken2-build >/dev/null || { echo "kraken2-build missing" >&2; exit 1; }

for metadata in nodes.dmp names.dmp nucl_gb.accession2taxid nucl_wgs.accession2taxid; do
  if [[ ! -f "$DB/taxonomy/$metadata" ]]; then
    echo "Missing local Kraken2 metadata: $DB/taxonomy/$metadata" >&2
    echo "Download Kraken2 taxonomy manually using the command in 0.README.md, then rerun." >&2
    exit 1
  fi
done

# Add downloaded reference FASTA directly. This can take a long time.
# kraken2-build accepts FASTA; for gz files we stream-decompress to temporary named pipes/files is unsafe.
# Build a stable symlink/decompression library once under derived analysis storage, never touch source refs.
LIB="$OUT/library"
mkdir -p "$LIB"
find "$ROOT" -type f -name '*_genomic.fna.gz' -print0 | while IFS= read -r -d '' f; do
  bn="$(basename "$f" .gz)"
  dst="$LIB/$bn"
  if [[ ! -s "$dst" ]]; then
    pigz -dc "$f" > "$dst.partial"
    mv "$dst.partial" "$dst"
  fi
done

find "$LIB" -type f -name '*.fna' -print0 | while IFS= read -r -d '' f; do
  kraken2-build --add-to-library "$f" --db "$DB"
done

ARGS=(--build --db "$DB" --threads 16)
if [[ "$MAX_GB" != 0 ]]; then
  ARGS+=(--max-db-size "$((MAX_GB*1024*1024*1024))")
fi
kraken2-build "${ARGS[@]}"

echo "PASS Kraken2 comparator DB: $DB"
echo "SOURCE REFERENCE RETAINED: $ROOT"
