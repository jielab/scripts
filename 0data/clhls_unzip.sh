#!/usr/bin/env bash
set -eo pipefail



usage() {
  cat <<'USAGE'
Usage:
  ./clhls_unzip.sh
  RAW=[raw_dir] CLEAN=[clean_dir] ./clhls_unzip.sh
  RAW=[raw_dir] CLEAN=[clean_dir] QC=[qc_dir] ./clhls_unzip.sh
  ./clhls_unzip.sh -h

What it does:
  1) Find all .zip and .rar files under RAW.
  2) Extract each dataverse_files.zip from the long Chinese-named folder.
  3) Put extracted files into CLEAN/<short_folder>/.
  4) Put QC files in _extract_qc at the same level as raw/ and clean/.
  5) Flatten internal folders to avoid very long paths.
  6) Never overwrite files:
       - same filename + identical content: skip and record in skipped_identical.tsv
       - same filename + different content: rename and record in renamed_conflicts.tsv

Default paths:
  RAW=/mnt/d/data/clhls杜克/raw
  CLEAN=/mnt/d/data/clhls杜克/clean
  QC=/mnt/d/data/clhls杜克/_extract_qc

Install tools if needed:
  sudo apt-get update && sudo apt-get install -y unzip unrar
USAGE
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
fi

RAW="${RAW:-/mnt/d/data/clhls杜克/raw}"
CLEAN="${CLEAN:-/mnt/d/data/clhls杜克/clean}"
QC="${QC:-$(dirname "$RAW")/_extract_qc}"
KEEP_TMP="${KEEP_TMP:-0}"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    echo "Install with: sudo apt-get update && sudo apt-get install -y unzip unrar" >&2
    exit 1
  }
}

for x in unzip find sort cmp cp basename dirname sed awk wc date mkdir rm tr grep paste chmod; do
  need_cmd "$x"
done

[[ -d "$RAW" ]] || { echo "ERROR: RAW directory does not exist: $RAW" >&2; exit 1; }

# Avoid accidental deletion of broad or important folders.
unsafe_path() {
  case "$1" in
    /|/mnt|/mnt/d|/mnt/d/data|""|.) return 0 ;;
    *) return 1 ;;
  esac
}

if unsafe_path "$CLEAN"; then
  echo "ERROR: unsafe CLEAN path: $CLEAN" >&2
  exit 1
fi

if unsafe_path "$QC"; then
  echo "ERROR: unsafe QC path: $QC" >&2
  exit 1
fi

if [[ "$QC" == "$RAW" || "$QC" == "$CLEAN" ]]; then
  echo "ERROR: QC cannot be the same as RAW or CLEAN." >&2
  echo "RAW=$RAW" >&2
  echo "CLEAN=$CLEAN" >&2
  echo "QC=$QC" >&2
  exit 1
fi

sanitize() {
  printf '%s' "$1" |
    sed 's#[/[:space:]]#_#g; s#[^A-Za-z0-9._-]#_#g; s#_\{2,\}#_#g; s#^_##; s#_$##'
}

years_from_name() {
  printf '%s\n' "$1" |
    grep -oE '(19|20)[0-9]{2}' |
    awk '!seen[$0]++' |
    paste -sd '_' -
}

clhls_folder_label() {
  local name="$1"
  local years label

  years="$(years_from_name "$name" || true)"
  [[ -n "$years" ]] || years="unknown"

  if [[ "$name" == *"截面数据"* ]]; then
    label="cross_${years}"
  elif [[ "$name" == *"社区数据"* ]]; then
    label="community_${years}"
  elif [[ "$name" == *"生物医学指标"* ]]; then
    label="biomarker_${years}"
  elif [[ "$name" == *"成年子女"* || "$name" == *"配对样本"* ]]; then
    label="paired_child_${years}"
  elif [[ "$name" == *"追踪数据"* || "$name" == *"跟踪调查"* ]]; then
    label="tracking_${years}"
  else
    label="$(sanitize "$name")"
  fi

  printf '%s' "$label"
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
  local archive="$1"
  local tmp="$2"

  case "${archive,,}" in
    *.zip)
      unzip -qq -o "$archive" -d "$tmp" >> "$LOG" 2>&1
      ;;
    *.rar)
      need_cmd unrar
      unrar x -o+ "$archive" "$tmp/" >> "$LOG" 2>&1
      ;;
    *)
      echo "Unsupported archive: $archive" >> "$LOG"
      return 2
      ;;
  esac
}

printf 'RAW   = %s\n' "$RAW"
printf 'CLEAN = %s\n' "$CLEAN"
printf 'QC    = %s\n' "$QC"
printf 'Removing old CLEAN and QC directories...\n'

rm -rf -- "$CLEAN" "$QC"
mkdir -p -- "$CLEAN" "$QC"

TMP="$QC/.tmp_extract"
mkdir -p -- "$TMP"

LOG="$QC/extract.log"
ARCHIVES="$QC/archive_list.tsv"
MAP="$QC/folder_map.tsv"
MANIFEST="$QC/extracted_manifest.tsv"
SKIPPED="$QC/skipped_identical.tsv"
RENAMED="$QC/renamed_conflicts.tsv"
FAILED="$QC/failed_archives.tsv"
FILELIST="$QC/file_list.txt"
LEFT_ARCHIVES="$QC/remaining_archives_in_clean.txt"
SUMMARY="$QC/summary.txt"

: > "$LOG"

printf 'raw_folder\toutput_folder\n' > "$MAP"
printf 'output_folder\traw_folder\tarchive_rel\tarchive_abs\n' > "$ARCHIVES"
printf 'output_folder\traw_folder\tarchive_rel\tinternal_file\toutput_file\taction\tsize_bytes\n' > "$MANIFEST"
printf 'output_folder\traw_folder\tarchive_rel\tinternal_file\texisting_output\n' > "$SKIPPED"
printf 'output_folder\traw_folder\tarchive_rel\tinternal_file\toriginal_output\trenamed_output\n' > "$RENAMED"
printf 'output_folder\traw_folder\tarchive_rel\tarchive_abs\terror\n' > "$FAILED"

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "$LOG" >&2
}

cleanup() {
  if [[ "$KEEP_TMP" == "1" ]]; then
    log "KEEP_TMP=1, temporary folder kept: $TMP"
  else
    rm -rf -- "$TMP"
  fi
}
trap cleanup EXIT

declare -A WAVE_MAP
declare -A LABEL_OWNER

while IFS= read -r -d '' d; do
  raw_folder="$(basename -- "$d")"
  base_label="$(clhls_folder_label "$raw_folder")"
  label="$base_label"
  n=2

  while [[ -n "${LABEL_OWNER[$label]:-}" && "${LABEL_OWNER[$label]}" != "$raw_folder" ]]; do
    label="${base_label}__${n}"
    n=$((n + 1))
  done

  LABEL_OWNER["$label"]="$raw_folder"
  WAVE_MAP["$raw_folder"]="$label"

  mkdir -p -- "$CLEAN/$label"
  printf '%s\t%s\n' "$raw_folder" "$label" >> "$MAP"
done < <(find "$RAW" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

mapfile -d '' archives < <(
  find "$RAW" -type f \( -iname '*.zip' -o -iname '*.rar' \) -print0 | sort -z
)

[[ ${#archives[@]} -gt 0 ]] || {
  echo "ERROR: no .zip or .rar files found under $RAW" >&2
  exit 1
}

log "Found ${#archives[@]} archives."

i=0
for archive in "${archives[@]}"; do
  i=$((i + 1))

  rel="${archive#$RAW/}"
  raw_folder="${rel%%/*}"

  if [[ "$raw_folder" == "$rel" ]]; then
    raw_folder="unknown_folder"
    output_folder="unknown_folder"
  else
    output_folder="${WAVE_MAP[$raw_folder]:-$(sanitize "$raw_folder")}"
  fi

  outdir="$CLEAN/$output_folder"
  mkdir -p -- "$outdir"

  printf '%s\t%s\t%s\t%s\n' "$output_folder" "$raw_folder" "$rel" "$archive" >> "$ARCHIVES"

  tmpdir="$TMP/${i}_$(sanitize "$output_folder")"
  rm -rf -- "$tmpdir"
  mkdir -p -- "$tmpdir"

  log "EXTRACT [$i/${#archives[@]}]: $rel -> $CLEAN/$output_folder/"

  if ! extract_one_archive "$archive" "$tmpdir"; then
    printf '%s\t%s\t%s\t%s\t%s\n' "$output_folder" "$raw_folder" "$rel" "$archive" "extract_failed" >> "$FAILED"
    echo "ERROR: extraction failed: $archive" >&2
    echo "See log: $LOG" >&2
    exit 1
  fi

  n_files=0

  while IFS= read -r -d '' f; do
    internal="${f#$tmpdir/}"

    case "$internal" in
      __MACOSX/*|*/__MACOSX/*|.DS_Store|*/.DS_Store|Thumbs.db|*/Thumbs.db)
        continue
        ;;
    esac

    base="$(basename -- "$f")"
    target="$outdir/$base"
    target_rel="${target#$CLEAN/}"
    size="$(wc -c < "$f" | tr -d ' ')"

    if [[ ! -e "$target" ]]; then
      cp -p -- "$f" "$target"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$output_folder" "$raw_folder" "$rel" "$internal" "$target_rel" "copied" "$size" >> "$MANIFEST"

    elif cmp -s -- "$f" "$target"; then
      printf '%s\t%s\t%s\t%s\t%s\n' \
        "$output_folder" "$raw_folder" "$rel" "$internal" "$target_rel" >> "$SKIPPED"

      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$output_folder" "$raw_folder" "$rel" "$internal" "$target_rel" "skipped_identical" "$size" >> "$MANIFEST"

    else
      stub="$(sanitize "$output_folder")"
      new_target="$(unique_name "$outdir" "$stub" "$base")"
      new_rel="${new_target#$CLEAN/}"

      cp -p -- "$f" "$new_target"

      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$output_folder" "$raw_folder" "$rel" "$internal" "$target_rel" "$new_rel" >> "$RENAMED"

      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$output_folder" "$raw_folder" "$rel" "$internal" "$new_rel" "renamed_conflict" "$size" >> "$MANIFEST"
    fi

    n_files=$((n_files + 1))
  done < <(find "$tmpdir" -type f -print0 | sort -z)

  if [[ $n_files -eq 0 ]]; then
    log "WARNING: no regular files found inside archive: $rel"
  fi

  [[ "$KEEP_TMP" == "1" ]] || rm -rf -- "$tmpdir"
done

find "$CLEAN" -type f -printf '%P\n' | sort > "$FILELIST"
find "$CLEAN" -type f \( -iname '*.zip' -o -iname '*.rar' -o -iname '*.7z' \) -printf '%P\n' | sort > "$LEFT_ARCHIVES"

n_out="$(wc -l < "$FILELIST" | tr -d ' ')"
n_skip="$(awk 'NR > 1 {n++} END {print n+0}' "$SKIPPED")"
n_rename="$(awk 'NR > 1 {n++} END {print n+0}' "$RENAMED")"
n_failed="$(awk 'NR > 1 {n++} END {print n+0}' "$FAILED")"
n_left_archive="$(wc -l < "$LEFT_ARCHIVES" | tr -d ' ')"

{
  echo "RAW=$RAW"
  echo "CLEAN=$CLEAN"
  echo "QC=$QC"
  echo "archives_found=${#archives[@]}"
  echo "final_output_files=$n_out"
  echo "skipped_identical_duplicates=$n_skip"
  echo "renamed_different_content_conflicts=$n_rename"
  echo "failed_archives=$n_failed"
  echo "remaining_archives_in_clean=$n_left_archive"
} > "$SUMMARY"

log "DONE. Final output files: $n_out"
log "Skipped identical duplicates: $n_skip"
log "Renamed different-content conflicts: $n_rename"
log "Failed archives: $n_failed"
log "Remaining compressed files inside CLEAN: $n_left_archive"

if [[ "$n_left_archive" -gt 0 ]]; then
  log "NOTE: some compressed files remain inside CLEAN. See: $LEFT_ARCHIVES"
fi

echo
echo "Done. Main output folders:"
find "$CLEAN" -mindepth 1 -maxdepth 1 -type d -printf '  %p\n' | sort

echo
echo "QC files:"
echo "  $SUMMARY"
echo "  $MAP"
echo "  $FILELIST"
echo "  $MANIFEST"
echo "  $SKIPPED"
echo "  $RENAMED"
echo "  $FAILED"
echo "  $LEFT_ARCHIVES"
echo "  $LOG"

echo
echo "Useful checks:"
echo "  column -t -s \$'\t' $MAP"
echo "  column -t -s \$'\t' $RENAMED | less -S"
echo "  column -t -s \$'\t' $SKIPPED | less -S"
echo "  cat $SUMMARY"
echo "  find $CLEAN -maxdepth 2 -type f | head -50"
