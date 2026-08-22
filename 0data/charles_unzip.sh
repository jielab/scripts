#!/usr/bin/env bash
set -euo pipefail



usage() {
  cat <<'USAGE'
Usage:
  ./charles_unzip.sh
  RAW=[raw_dir] CLEAN=[clean_dir] ./charles_unzip.sh
  ./charles_unzip.sh -h

What it does:
  1) Find all .rar and .zip files under RAW.
  2) Keep the first-level folder name, such as 2008, 2011, 2020, Harmonized.
  3) Extract each archive to a temporary folder.
  4) Copy only the extracted files into CLEAN/<wave>/, removing long internal folders.

Example:
  RAW/2008/xxx.rar containing CHARLS 20101118/Main/Healthcare_20101118.dta
  -> CLEAN/2008/Healthcare_20101118.dta

Duplicate rule:
  - Never overwrite.
  - Same filename + identical content: skip and record in _extract_qc/skipped_identical.tsv
  - Same filename + different content: keep both by adding archive-name prefix and record in
    _extract_qc/renamed_conflicts.tsv

Defaults:
  RAW=/mnt/d/data/charls/raw
  CLEAN=/mnt/d/data/charls/clean

Install tools if needed:
  sudo apt-get update && sudo apt-get install -y unrar unzip
USAGE
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
fi

RAW="${RAW:-/mnt/d/data/charls/raw}"
CLEAN="${CLEAN:-/mnt/d/data/charls/clean}"
KEEP_TMP="${KEEP_TMP:-0}"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    echo "Install with: sudo apt-get update && sudo apt-get install -y unrar unzip" >&2
    exit 1
  }
}

for x in unrar unzip find sort cmp cp basename dirname sed awk wc date mkdir rm tr; do
  need_cmd "$x"
done

[[ -d "$RAW" ]] || { echo "ERROR: RAW directory does not exist: $RAW" >&2; exit 1; }

# Avoid accidental deletion of a broad folder.
case "$CLEAN" in
  /|/mnt|/mnt/d|/mnt/d/data|"")
    echo "ERROR: unsafe CLEAN path: $CLEAN" >&2
    exit 1
    ;;
esac

sanitize() {
  printf '%s' "$1" | sed 's#[/[:space:]]#_#g; s#[^A-Za-z0-9._-]#_#g; s#_\{2,\}#_#g; s#^_##; s#_$##'
}

archive_stub() {
  local x
  x="$(basename -- "$1")"
  x="${x%.*}"
  sanitize "$x"
}

unique_name() {
  local dir="$1" stub="$2" base="$3"
  local stem ext out n
  if [[ "$base" == *.* ]]; then
    stem="${base%.*}"
    ext=".${base##*.}"
  else
    stem="$base"
    ext=""
  fi

  out="$dir/${stub}__${stem}${ext}"
  n=2
  while [[ -e "$out" ]]; do
    out="$dir/${stub}__${stem}.${n}${ext}"
    n=$((n + 1))
  done
  printf '%s' "$out"
}

extract_one_archive() {
  local archive="$1" tmp="$2"
  case "${archive,,}" in
    *.rar) unrar x -o+ "$archive" "$tmp/" >> "$LOG" 2>&1 ;;
    *.zip) unzip -qq -o "$archive" -d "$tmp" >> "$LOG" 2>&1 ;;
    *) echo "Unsupported archive: $archive" >> "$LOG"; return 2 ;;
  esac
}

printf 'RAW   = %s\n' "$RAW"
printf 'CLEAN = %s\n' "$CLEAN"
printf 'Removing old CLEAN directory...\n'
rm -rf -- "$CLEAN"
mkdir -p -- "$CLEAN"

QC="$CLEAN/_extract_qc"
TMP="$CLEAN/.tmp_extract"
mkdir -p -- "$QC" "$TMP"

LOG="$QC/extract.log"
ARCHIVES="$QC/archive_list.tsv"
MANIFEST="$QC/extracted_manifest.tsv"
SKIPPED="$QC/skipped_identical.tsv"
RENAMED="$QC/renamed_conflicts.tsv"
FAILED="$QC/failed_archives.tsv"
FILELIST="$QC/file_list.txt"

: > "$LOG"
printf 'wave\tarchive_rel\tarchive_abs\n' > "$ARCHIVES"
printf 'wave\tarchive_rel\tinternal_file\toutput_file\taction\tsize_bytes\n' > "$MANIFEST"
printf 'wave\tarchive_rel\tinternal_file\texisting_output\n' > "$SKIPPED"
printf 'wave\tarchive_rel\tinternal_file\toriginal_output\trenamed_output\n' > "$RENAMED"
printf 'wave\tarchive_rel\tarchive_abs\terror\n' > "$FAILED"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG" >&2; }

cleanup() {
  if [[ "$KEEP_TMP" == "1" ]]; then
    log "KEEP_TMP=1, temporary folder kept: $TMP"
  else
    rm -rf -- "$TMP"
  fi
}
trap cleanup EXIT

# Create the same first-level folders as RAW.
while IFS= read -r -d '' d; do
  mkdir -p -- "$CLEAN/$(basename -- "$d")"
done < <(find "$RAW" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

mapfile -d '' archives < <(
  find "$RAW" -type f \( -iname '*.rar' -o -iname '*.zip' \) -print0 | sort -z
)

[[ ${#archives[@]} -gt 0 ]] || { echo "ERROR: no .rar or .zip files found under $RAW" >&2; exit 1; }
log "Found ${#archives[@]} archives."

i=0
for archive in "${archives[@]}"; do
  i=$((i + 1))
  rel="${archive#$RAW/}"
  wave="${rel%%/*}"
  [[ "$wave" != "$rel" ]] || wave="unknown_wave"
  outdir="$CLEAN/$wave"
  mkdir -p -- "$outdir"

  printf '%s\t%s\t%s\n' "$wave" "$rel" "$archive" >> "$ARCHIVES"

  tmpdir="$TMP/${i}_$(sanitize "$rel")"
  rm -rf -- "$tmpdir"
  mkdir -p -- "$tmpdir"

  log "EXTRACT [$i/${#archives[@]}]: $rel -> $CLEAN/$wave/"

  if ! extract_one_archive "$archive" "$tmpdir"; then
    printf '%s\t%s\t%s\t%s\n' "$wave" "$rel" "$archive" "extract_failed" >> "$FAILED"
    echo "ERROR: extraction failed: $archive" >&2
    echo "See log: $LOG" >&2
    exit 1
  fi

  n_files=0
  while IFS= read -r -d '' f; do
    internal="${f#$tmpdir/}"
    case "$internal" in
      __MACOSX/*|*/__MACOSX/*|.DS_Store|*/.DS_Store|Thumbs.db|*/Thumbs.db) continue ;;
    esac

    base="$(basename -- "$f")"
    target="$outdir/$base"
    target_rel="${target#$CLEAN/}"
    size="$(wc -c < "$f" | tr -d ' ')"

    if [[ ! -e "$target" ]]; then
      cp -p -- "$f" "$target"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$wave" "$rel" "$internal" "$target_rel" "copied" "$size" >> "$MANIFEST"
    elif cmp -s -- "$f" "$target"; then
      printf '%s\t%s\t%s\t%s\n' "$wave" "$rel" "$internal" "$target_rel" >> "$SKIPPED"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$wave" "$rel" "$internal" "$target_rel" "skipped_identical" "$size" >> "$MANIFEST"
    else
      stub="$(archive_stub "$rel")"
      new_target="$(unique_name "$outdir" "$stub" "$base")"
      new_rel="${new_target#$CLEAN/}"
      cp -p -- "$f" "$new_target"
      printf '%s\t%s\t%s\t%s\t%s\n' "$wave" "$rel" "$internal" "$target_rel" "$new_rel" >> "$RENAMED"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$wave" "$rel" "$internal" "$new_rel" "renamed_conflict" "$size" >> "$MANIFEST"
    fi

    n_files=$((n_files + 1))
  done < <(find "$tmpdir" -type f -print0 | sort -z)

  if [[ $n_files -eq 0 ]]; then
    log "WARNING: no regular files found inside archive: $rel"
  fi

  [[ "$KEEP_TMP" == "1" ]] || rm -rf -- "$tmpdir"
done

find "$CLEAN" -type f ! -path "$QC/*" ! -path "$TMP/*" -printf '%P\n' | sort > "$FILELIST"

n_out="$(wc -l < "$FILELIST" | tr -d ' ')"
n_skip="$(awk 'NR > 1 {n++} END {print n+0}' "$SKIPPED")"
n_rename="$(awk 'NR > 1 {n++} END {print n+0}' "$RENAMED")"

log "DONE. Final output files: $n_out"
log "Skipped identical duplicates: $n_skip"
log "Renamed different-content conflicts: $n_rename"

echo
echo "Done. Main output folders:"
find "$CLEAN" -mindepth 1 -maxdepth 1 -type d ! -name '_extract_qc' ! -name '.tmp_extract' -printf '  %p\n' | sort

echo
echo "QC files:"
echo "  $FILELIST"
echo "  $MANIFEST"
echo "  $SKIPPED"
echo "  $RENAMED"
echo "  $LOG"
