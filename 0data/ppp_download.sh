#!/usr/bin/env bash
set -euo pipefail


# 🚩 Paths and defaults
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PY="$SCRIPT_DIR/f/ppp_download.py"
PYTHON_BIN=""

TOKEN_FILE_DEFAULT="/mnt/d/data/ukb/prot/authToken.txt"
OUT_DEFAULT="/mnt/d/data.BIG/gwas/prot/.project/download"
LOG_DEFAULT="/mnt/d/data.BIG/gwas/prot/.project/log"
PROT_DIR_DEFAULT="/mnt/d/data.BIG/gwas/prot"
PROT_TSV_DEFAULT="/mnt/d/files/ppp_3k.38.tsv"
PROT_BED_DEFAULT="/mnt/d/files/ppp_3k.38.bed"
FOLDER_DEFAULT="syn51365303"                    # UKB proteomics European discovery folder
NAME_REGEX_DEFAULT='\.tar$'
MANIFEST_DRYRUN="European_discovery.dryrun.tsv"
MANIFEST_DOWNLOAD="European_discovery.download_manifest.tsv"
LOG_FILE="ppp_download.European_discovery.log"
RENAME_LOG_FILE="ppp_rename.log"
RENAME_MANIFEST_FILE="ppp_rename_manifest.tsv"
RENAME_BADNAME_FILE="ppp_rename_badname.tsv"
RENAME_DUPLICATE_FILE="ppp_rename_duplicate_shortname.tsv"
RENAME_MISSING_FILE="ppp_rename_missing_protein.tsv"
UNZIP_LOG_FILE="ppp_unzip.log"
UNZIP_MANIFEST_FILE="ppp_unzip_manifest.tsv"
UNZIP_FAILED_FILE="ppp_unzip_failed.tsv"
RSID_MAP_ALL_FILE="olink_rsid_map_mac5_info03_b0_7_ALL_patched_v2.tsv"
RSID_MAP_LOOKUP_FILE="olink_rsid_map_mac5_info03_b0_7_ALL_patched_v2.sorted.tsv"
RSID_MAP_LOG_FILE="ppp_rsid_map.log"
RSID_MAP_MANIFEST_FILE="ppp_rsid_map_manifest.tsv"
RSID_MAP_FAILED_FILE="ppp_rsid_map_failed.tsv"


# 🚩 Command-line help
usage() {
  cat <<EOF
Usage:
  ./ppp_download.sh [mode] [options]

Modes:
  test dryrun foreground download rename unzip rsid_map

Mode options:
  dryrun      [extra ppp_download.py args]
  foreground  [--workers N] [--limit N] [extra args]
  download    [--workers N] [--limit N] [extra args]
  rename      [--run-rename TRUE|FALSE] [--replace-bed TRUE|FALSE]
  unzip       [--jobs N] [--overwrite TRUE|FALSE]
  rsid_map    [--jobs N|auto] [--threads N|auto] [--overwrite TRUE|FALSE] [--delete-raw TRUE|FALSE] [--tmpdir DIR] [--sort-memory SIZE]
              defaults: --jobs 4 --threads 4

Common options:
  --token-file FILE       Synapse token file.
  --output-dir DIR        Download/raw directory.
  --log-dir DIR           Log directory.
  --prot-dir DIR           PROT project directory.
  --prot-tsv FILE          Protein metadata TSV.
  --prot-bed FILE          Protein BED output.
  --folder ID             Synapse folder ID.
  --name-regex REGEX      Download filename filter.
  --python FILE           Python interpreter.
  --clean-dir DIR         Rename output directory.
  --bad-name-dir DIR      Bad-name quarantine directory.
  --rsid-map-dir DIR      rsID map directory.
  --rsid-map-all FILE     Combined rsID map.
  --rsid-map-lookup FILE  Sorted rsID lookup.
  --cache-lookup TRUE|FALSE Cache lookup data on the temporary filesystem.

Examples:
  ./ppp_download.sh test
  ./ppp_download.sh dryrun
  ./ppp_download.sh foreground --workers 2 --limit 5
  ./ppp_download.sh download --workers 4 #  不要设置太高
  ./ppp_download.sh rename
  ./ppp_download.sh rename --run-rename TRUE
  ./ppp_download.sh unzip
  ./ppp_download.sh unzip --jobs 8 # background; terminal can be closed
  tail -f $LOG_DEFAULT/$UNZIP_LOG_FILE
  ./ppp_download.sh rsid_map # auto CPU allocation; temp files under WSL ext4 /tmp/ppp_rsid_map
  tail -f $LOG_DEFAULT/$RSID_MAP_LOG_FILE
  pgrep -af 'ppp_download.py|ppp_download.sh' # 

Defaults:
  token file : $TOKEN_FILE_DEFAULT
  output dir : $OUT_DEFAULT
  log dir    : $LOG_DEFAULT
  unzip out  : next to each source .tar file
  PROT TSV    : $PROT_TSV_DEFAULT
  protein BED    : $PROT_BED_DEFAULT
  folder     : $FOLDER_DEFAULT
  file regex : $NAME_REGEX_DEFAULT

Path examples:
  ./ppp_download.sh download --token-file /path/authToken.txt --output-dir /path/raw --log-dir /path/log --folder syn... --workers 8
  ./ppp_download.sh rename --prot-dir /path/prot --prot-tsv /path/ppp_3k.38.tsv --prot-bed /path/ppp_3k.38.bed
  ./ppp_download.sh unzip --output-dir /path/raw
  unzip writes each [prot].raw.gz next to its source [prot].tar
  rsid_map writes common/[prot]/raw/[prot].gz from [prot].raw.gz, then deletes [prot].raw.gz by default
  rsid_map prefers pigz when installed; otherwise it falls back to gzip
  rsid_map caches the sorted lookup under its temp directory; use --cache-lookup FALSE to disable
EOF
}


# 🚩 Shared helpers
has_arg() {
  local target="$1"; shift || true
  for x in "$@"; do
    [[ "$x" == "$target" ]] && return 0
  done
  return 1
}

add_default_workers_if_missing() {
  local default_workers="$1"; shift
  if has_arg "--workers" "$@"; then
    printf '%s\0' "$@"
  else
    printf '%s\0' --workers "$default_workers" "$@"
  fi
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 2
  }
}

cpu_count() {
  local n
  n="$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')"
  [[ "$n" =~ ^[1-9][0-9]*$ ]] || n=1
  printf '%s\n' "$n"
}

gzip_decompress() {
  local file="$1" threads="${2:-1}"
  if command -v pigz >/dev/null 2>&1; then
    pigz -dc -p "$threads" -- "$file"
  else
    gzip -dc -- "$file"
  fi
}

gzip_compress() {
  local threads="${1:-1}"
  if command -v pigz >/dev/null 2>&1; then
    pigz -c -p "$threads"
  else
    gzip -c
  fi
}

bool_is_true() {
  case "${1^^}" in
    TRUE|T|YES|Y|1) return 0 ;;
    FALSE|F|NO|N|0|"") return 1 ;;
    *)
      echo "ERROR: expected TRUE or FALSE, got: $1" >&2
      exit 2
      ;;
  esac
}

parse_key_bool_arg() {
  local key="$1" default_value="$2"; shift 2
  local value="$default_value"
  while [[ $# -gt 0 ]]; do
    case "$1" in
      "$key")
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after $key" >&2; exit 2; }
        value="$2"
        shift 2
        ;;
      "$key"=*)
        value="${1#*=}"
        shift
        ;;
      *)
        echo "ERROR: unknown argument for this mode: $1" >&2
        exit 2
        ;;
    esac
  done
  case "${value^^}" in
    TRUE|T|YES|Y|1|FALSE|F|NO|N|0) printf '%s\n' "$value" ;;
    *) echo "ERROR: expected TRUE or FALSE for $key, got: $value" >&2; exit 2 ;;
  esac
}

task_log() {
  local line
  line="[$(date '+%F %T')] $*"
  printf '%s\n' "$line" | tee -a "$TASK_LOG" >&2
}

running_prot_tasks() {
  local self="${1:-$$}" bash_self="${2:-${BASHPID:-$$}}" parent="${3:-${PPID:-}}"
  ps -eo pid=,ppid=,args= 2>/dev/null | PROT_TASK_PATTERN='ppp_download[.]sh|ppp_download[.]py' awk -v self="$self" -v bash_self="$bash_self" -v parent="$parent" '
    BEGIN { pattern = ENVIRON["PROT_TASK_PATTERN"] }
    {
      pid = $1
      ppid = $2
      $1 = ""
      $2 = ""
      sub(/^[[:space:]]+/, "")
      cmd = $0
    }
    pid == self { next }
    pid == bash_self { next }
    pid == parent { next }
    ppid == self { next }
    ppid == bash_self { next }
    cmd ~ pattern { print pid "\t" cmd }
  '
}

check_existing_prot_tasks() {
  local old pids self bash_self parent
  self="$$"
  bash_self="${BASHPID:-$$}"
  parent="${PPID:-}"
  old="$(running_prot_tasks "$self" "$bash_self" "$parent" || true)"
  [[ -z "$old" ]] && return 0

  pids="$(printf '%s\n' "$old" | awk '{print $1}' | paste -sd ' ' -)"
  {
    echo "ERROR: existing ppp_download task(s) are still running:"
    printf '%s\n' "$old" | sed 's/^/  /'
    echo
    echo "Stop them first, then rerun this command:"
    echo "  kill $pids"
    echo
    echo "If they do not exit after a few seconds, use:"
    echo "  kill -9 $pids"
  } >&2
  return 1
}

unique_path() {
  local dir="$1" base="$2"
  local stem ext out n
  if [[ "$base" == *.* ]]; then
    stem="${base%.*}"
    ext=".${base##*.}"
  else
    stem="$base"
    ext=""
  fi

  out="$dir/$base"
  n=2
  while [[ -e "$out" ]]; do
    out="$dir/${stem}.${n}${ext}"
    n=$((n + 1))
  done
  printf '%s\n' "$out"
}


# 🚩 Protein annotation
generate_prot_bed() {
  local tsv="$1" bed="$2"
  local bed_dir tmp unsorted issues dup rc

  [[ -s "$tsv" ]] || { task_log "ERROR: PROT TSV not found or empty: $tsv"; return 2; }
  bed_dir="$(dirname -- "$bed")"
  mkdir -p -- "$bed_dir"
  tmp="$bed.tmp.$$"
  unsorted="$bed.tmp.$$.unsorted"
  issues="$bed.tmp.$$.issues"
  dup="$bed.tmp.$$.duplicates"
  rm -f -- "$tmp" "$unsorted" "$issues" "$dup"

  set +e
  awk -F '\t' -v OFS='\t' -v out="$unsorted" -v issues="$issues" -v dupfile="$dup" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        sub(/\r$/, "", $i)
      }
      if (NF < 13 || $5 != "Assay" || $6 != "Panel" || $11 != "chr" || $12 != "gene_start" || $13 != "gene_end") {
        print "expected columns 5=Assay, 6=Panel, 11=chr, 12=gene_start, 13=gene_end" > issues
        fatal = 1
        exit 3
      }
      next
    }
    {
      for (i = 1; i <= NF; i++) sub(/\r$/, "", $i)
      chr = $11
      start = $12
      end = $13
      prot = $5
      panel = $6

      if (prot == "" || panel == "" || chr == "" || start == "" || end == "") {
        print "row=" NR "\tprotein=" prot "\tpanel=" panel "\tchr=" chr "\tstart=" start "\tend=" end "\terror=missing_required_value" > issues
        bad = 1
      } else if (start !~ /^[0-9]+$/ || end !~ /^[0-9]+$/) {
        print "row=" NR "\tprotein=" prot "\tchr=" chr "\tstart=" start "\tend=" end "\terror=non_integer_coordinate" > issues
        bad = 1
      } else if ((end + 0) <= (start + 0)) {
        print "row=" NR "\tprotein=" prot "\tchr=" chr "\tstart=" start "\tend=" end "\terror=gene_end_not_greater_than_gene_start" > issues
        bad = 1
      }

      key = prot SUBSEP panel
      if (!(key in key_seen)) {
        key_seen[key] = 1
        key_order[++n_key] = key
        key_prot[key] = prot
        key_panel[key] = panel
        key_chr[key] = chr
        key_start[key] = start
        key_end[key] = end
        panel_count[prot]++
      } else {
        key_dup_count[key]++
        if (key_dup_count[key] == 1) {
          dup_order[++n_dup] = prot "\t" panel
        }
      }
    }
    END {
      if (fatal) exit 3
      for (i = 1; i <= n_dup; i++) print dup_order[i] > dupfile
      if (bad) exit 1
      for (i = 1; i <= n_key; i++) {
        key = key_order[i]
        prot = key_prot[key]
        panel = key_panel[key]
        out_prot = prot
        if (panel_count[prot] > 1) out_prot = prot "." panel
        if (out_prot in final_seen) {
          print "protein=" out_prot "\terror=duplicate_final_bed_name" > issues
          final_bad = 1
        } else {
          final_seen[out_prot] = 1
          print key_chr[key], key_start[key], key_end[key], out_prot, panel >> out
        }
      }
      if (final_bad) exit 1
    }
  ' "$tsv"
  rc=$?
  set -e

  if [[ -s "$dup" ]]; then
    task_log "WARNING duplicate Assay+Panel rows in TSV; keeping first record only:"
    sed 's/^/  /' "$dup" | tee -a "$TASK_LOG" >&2
  fi

  if (( rc != 0 )); then
    task_log "WARNING ppp_3k.38.bed was not generated because TSV validation failed."
    [[ -s "$issues" ]] && sed 's/^/  /' "$issues" | tee -a "$TASK_LOG" >&2
    rm -f -- "$tmp" "$unsorted" "$issues" "$dup"
    return "$rc"
  fi

  if ! awk -F '\t' -v OFS='\t' '
    function chr_key(chr, c) {
      c = toupper(chr)
      sub(/^CHR/, "", c)
      if (c ~ /^[0-9]+$/) return c + 0
      if (c == "X") return 23
      if (c == "Y") return 24
      if (c == "M" || c == "MT") return 25
      return 1000
    }
    {
      print chr_key($1), $2 + 0, $3 + 0, NR, $0
    }
  ' "$unsorted" | sort -t $'\t' -k1,1n -k2,2n -k3,3n -k4,4n | cut -f5- > "$tmp"; then
    task_log "ERROR: failed to sort BED by chr and gene_start."
    rm -f -- "$tmp" "$unsorted" "$issues" "$dup"
    return 2
  fi

  mv -- "$tmp" "$bed"
  rm -f -- "$unsorted" "$issues" "$dup"
  task_log "Generated sorted BED: $bed"
}

validate_prot_bed() {
  local bed="$1"
  local issues dup rc
  [[ -s "$bed" ]] || { task_log "ERROR: BED not found or empty: $bed"; return 2; }
  issues="$bed.validate.$$.issues"
  dup="$bed.validate.$$.duplicates"
  rm -f -- "$issues" "$dup"

  set +e
  awk -F '\t' '
    NF < 5 {
      print "row=" NR "\terror=expected_at_least_5_tab_separated_columns" > issues
      bad = 1
      next
    }
    {
      sub(/\r$/, "", $5)
      if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) {
        print "row=" NR "\tprotein=" $4 "\tstart=" $2 "\tend=" $3 "\terror=non_integer_coordinate" > issues
        bad = 1
      } else if (($3 + 0) <= ($2 + 0)) {
        print "row=" NR "\tprotein=" $4 "\tstart=" $2 "\tend=" $3 "\terror=gene_end_not_greater_than_gene_start" > issues
        bad = 1
      }
      seen[$4]++
      if (seen[$4] == 2) dup_order[++n_dup] = $4
    }
    END {
      for (i = 1; i <= n_dup; i++) print dup_order[i] > dupfile
      if (n_dup > 0 || bad) exit 1
    }
  ' issues="$issues" dupfile="$dup" "$bed"
  rc=$?
  set -e

  if [[ -s "$dup" ]]; then
    task_log "WARNING duplicate protein names in BED; BED must have unique column 4 values:"
    sed 's/^/  /' "$dup" | tee -a "$TASK_LOG" >&2
  fi
  if [[ -s "$issues" ]]; then
    task_log "WARNING invalid BED coordinates or shape:"
    sed 's/^/  /' "$issues" | tee -a "$TASK_LOG" >&2
  fi
  rm -f -- "$issues" "$dup"
  return "$rc"
}

ensure_prot_bed() {
  local replace_bed="${1:-FALSE}"
  if [[ ! -s "$PROT_BED" ]] || bool_is_true "$replace_bed"; then
    task_log "Preparing BED from TSV: $PROT_TSV"
    generate_prot_bed "$PROT_TSV" "$PROT_BED"
  else
    task_log "BED already exists; validating: $PROT_BED"
  fi
  validate_prot_bed "$PROT_BED"
}

derive_tar_protein() {
  local stem="$1"
  if [[ "$stem" =~ ^(.+)_([A-Za-z0-9]{6,10}(-[0-9]+)?(_[A-Za-z0-9]{6,10}(-[0-9]+)?)*)_[Oo][Ii][Dd][0-9]+_ ]]; then
    printf '%s\n' "${BASH_REMATCH[1]}"
    return 0
  fi
  return 1
}

derive_tar_panel() {
  local stem="$1"
  if [[ "$stem" =~ _[Oo][Ii][Dd][0-9]+_[Vv][^_]+_(.+)$ ]]; then
    printf '%s\n' "${BASH_REMATCH[1]}"
    return 0
  fi
  return 1
}

bed_protein_to_assay() {
  local protein="$1" panel="$2"
  if [[ -n "$panel" && "$protein" == *".$panel" ]]; then
    printf '%s\n' "${protein%.$panel}"
  else
    printf '%s\n' "$protein"
  fi
}


# 🚩 Archive renaming
run_rename_task() {
  local run_rename="FALSE" replace_bed="FALSE"
  local arg
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --run-rename)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --run-rename" >&2; exit 2; }
        run_rename="$2"
        shift 2
        ;;
      --run-rename=*)
        run_rename="${1#*=}"
        shift
        ;;
      --replace-bed)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --replace-bed" >&2; exit 2; }
        replace_bed="$2"
        shift 2
        ;;
      --replace-bed=*)
        replace_bed="${1#*=}"
        shift
        ;;
      *)
        echo "ERROR: unknown rename argument: $1" >&2
        exit 2
        ;;
    esac
  done
  bool_is_true "$run_rename" || true
  bool_is_true "$replace_bed" || true

  mkdir -p -- "$LOGDIR"
  TASK_LOG="$LOGDIR/$RENAME_LOG_FILE"
  local manifest="$LOGDIR/$RENAME_MANIFEST_FILE"
  local bad_manifest="$LOGDIR/$RENAME_BADNAME_FILE"
  local duplicate_manifest="$LOGDIR/$RENAME_DUPLICATE_FILE"
  local missing_manifest="$LOGDIR/$RENAME_MISSING_FILE"
  : > "$TASK_LOG"
  printf 'old_name\tnew_name\taction\tmatch_count\tmatches\n' > "$manifest"
  printf 'old_name\tbadname_path\taction\tmatch_count\tmatches\n' > "$bad_manifest"
  printf 'protein\tbed_panel\tselected_tar\ttar_file\ttar_panel\taction\n' > "$duplicate_manifest"
  printf 'protein\tbed_panel\taction\n' > "$missing_manifest"

  task_log "rename mode started; run_rename=$run_rename; raw=$OUT; badName=$BADNAME"
  ensure_prot_bed "$replace_bed"

  if [[ ! -d "$OUT" ]]; then
    task_log "WARNING raw directory does not exist, so only BED validation was completed: $OUT"
    return 0
  fi

  local -a proteins tars matches files
  local -A panel_by_protein assay_by_protein bed_by_name_key bed_by_assay_panel_key
  local bed_chr bed_start bed_end bed_panel assay key value
  while IFS=$'\t' read -r bed_chr bed_start bed_end prot bed_panel; do
    proteins+=("$prot")
    panel_by_protein["$prot"]="$bed_panel"
    assay="$(bed_protein_to_assay "$prot" "$bed_panel")"
    assay_by_protein["$prot"]="$assay"

    key="${prot,,}"
    bed_by_name_key["$key"]+="${bed_by_name_key["$key"]:+$'\t'}$prot"

    key="${assay,,}"$'\t'"${bed_panel,,}"
    bed_by_assay_panel_key["$key"]+="${bed_by_assay_panel_key["$key"]:+$'\t'}$prot"
  done < "$PROT_BED"
  mapfile -d '' tars < <(find "$OUT" -maxdepth 1 -type f -name '*.tar' -print0 | sort -z)
  if [[ ${#tars[@]} -eq 0 ]]; then
    task_log "WARNING no .tar files found under raw directory: $OUT"
    return 0
  fi

  local -A tar_base tar_protein tar_panel tar_assay tar_match_text tar_match_count target_count target_files selected_tar
  local tar base stem stem_key tar_assay_value tar_panel_value match_text target bad_target action n_renamed=0 n_bad=0 n_skip=0 n_would=0
  for tar in "${tars[@]}"; do
    base="$(basename -- "$tar")"
    stem="${base%.tar}"
    stem_key="${stem,,}"
    matches=()

    if [[ -n "${bed_by_name_key[$stem_key]+x}" ]]; then
      IFS=$'\t' read -r -a matches <<< "${bed_by_name_key[$stem_key]}"
      tar_assay_value="${assay_by_protein["${matches[0]}"]}"
      tar_panel_value="${panel_by_protein["${matches[0]}"]}"
    else
      tar_assay_value="$(derive_tar_protein "$stem" || true)"
      tar_panel_value="$(derive_tar_panel "$stem" || true)"
      if [[ -n "$tar_assay_value" && -n "$tar_panel_value" ]]; then
        key="${tar_assay_value,,}"$'\t'"${tar_panel_value,,}"
        if [[ -n "${bed_by_assay_panel_key[$key]+x}" ]]; then
          IFS=$'\t' read -r -a matches <<< "${bed_by_assay_panel_key[$key]}"
        fi
      fi
      if [[ ${#matches[@]} -eq 0 && -n "$tar_panel_value" ]]; then
        for prot in "${proteins[@]}"; do
          assay="${assay_by_protein["$prot"]}"
          bed_panel="${panel_by_protein["$prot"]}"
          if [[ "${tar_panel_value,,}" == "${bed_panel,,}" && "$stem_key" == "${assay,,}"_* ]]; then
            matches+=("$prot")
          fi
        done
      fi
    fi
    match_text="$(IFS=,; printf '%s' "${matches[*]-}")"
    tar_base["$tar"]="$base"
    tar_assay["$tar"]="$tar_assay_value"
    tar_panel["$tar"]="$tar_panel_value"
    tar_match_text["$tar"]="$match_text"
    tar_match_count["$tar"]="${#matches[@]}"

    if [[ ${#matches[@]} -eq 1 ]]; then
      prot="${matches[0]}"
      tar_protein["$tar"]="$prot"
      target_count["$prot"]=$(( ${target_count["$prot"]:-0} + 1 ))
      target_files["$prot"]+="${target_files["$prot"]:+$'\t'}$tar"
    fi
  done

  local duplicate_groups=0 duplicate_files=0 duplicate_selected=0 duplicate_to_bad=0
  local selected selected_reason item_status
  for prot in "${proteins[@]}"; do
    [[ ${target_count["$prot"]:-0} -gt 1 ]] || continue
    duplicate_groups=$((duplicate_groups + 1))
    duplicate_files=$((duplicate_files + target_count["$prot"]))
    IFS=$'\t' read -r -a files <<< "${target_files["$prot"]}"

    selected=""
    selected_reason=""
    for tar in "${files[@]}"; do
      if [[ "${tar_base["$tar"]}" == "$prot.tar" ]]; then
        selected="$tar"
        selected_reason="already_renamed"
        break
      fi
    done
    if [[ -z "$selected" ]]; then
      selected="${files[0]}"
      selected_reason="first_sorted_duplicate"
    fi

    if [[ -n "$selected" ]]; then
      selected_tar["$prot"]="$selected"
      duplicate_selected=$((duplicate_selected + 1))
    fi

    task_log "DUPLICATE BED name target: $prot; assay=${assay_by_protein["$prot"]}; bed_panel=${panel_by_protein["$prot"]}; files=${target_count["$prot"]}; selected=$(basename -- "${selected:-NONE}"); reason=$selected_reason"
    for tar in "${files[@]}"; do
      if [[ -n "$selected" && "$tar" == "$selected" ]]; then
        item_status="selected_${selected_reason}"
      else
        item_status="duplicate_to_badName"
        duplicate_to_bad=$((duplicate_to_bad + 1))
      fi
      task_log "  ${tar_base["$tar"]} panel=${tar_panel["$tar"]:-NA} action=$item_status"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$prot" \
        "${panel_by_protein["$prot"]}" \
        "$(basename -- "${selected:-NONE}")" \
        "${tar_base["$tar"]}" \
        "${tar_panel["$tar"]:-NA}" \
        "$item_status" >> "$duplicate_manifest"
    done
  done
  if (( duplicate_groups > 0 )); then
    task_log "duplicate BED-name targets: groups=$duplicate_groups; files=$duplicate_files; selected=$duplicate_selected; move_to_badName=$duplicate_to_bad"
    task_log "duplicate manifest: $duplicate_manifest"
  fi

  local missing_proteins=0
  for prot in "${proteins[@]}"; do
    [[ ${target_count["$prot"]:-0} -eq 0 ]] || continue
    missing_proteins=$((missing_proteins + 1))
    task_log "MISSING tar for BED protein: $prot; bed_panel=${panel_by_protein["$prot"]}"
    printf '%s\t%s\t%s\n' "$prot" "${panel_by_protein["$prot"]}" "no_matching_tar" >> "$missing_manifest"
  done
  if (( missing_proteins > 0 )); then
    task_log "missing BED proteins without tar: $missing_proteins"
    task_log "missing manifest: $missing_manifest"
  fi

  local duplicate_conflict
  for tar in "${tars[@]}"; do
    base="${tar_base["$tar"]}"
    match_text="${tar_match_text["$tar"]}"
    if [[ "${tar_match_count["$tar"]}" -eq 1 ]]; then
      prot="${tar_protein["$tar"]}"
      duplicate_conflict=0
      if [[ ${target_count["$prot"]:-0} -gt 1 && "${selected_tar["$prot"]:-}" != "$tar" ]]; then
        duplicate_conflict=1
      fi
      target="$OUT/$prot.tar"
      if (( duplicate_conflict )); then
        bad_target="$(unique_path "$BADNAME" "$base")"
        if bool_is_true "$run_rename"; then
          mkdir -p -- "$BADNAME"
          mv -- "$tar" "$bad_target"
          action="duplicate_short_name_moved_badName"
        else
          action="duplicate_short_name_would_move_badName"
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$base" "$bad_target" "$action" 1 "$match_text" >> "$bad_manifest"
        n_bad=$((n_bad + 1))
      elif [[ "$base" == "$prot.tar" ]]; then
        action="already_renamed"
        n_skip=$((n_skip + 1))
      elif [[ -e "$target" ]]; then
        bad_target="$(unique_path "$BADNAME" "$base")"
        if bool_is_true "$run_rename"; then
          mkdir -p -- "$BADNAME"
          mv -- "$tar" "$bad_target"
          action="target_exists_moved_badName"
        else
          action="target_exists_would_move_badName"
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$base" "$bad_target" "$action" "${tar_match_count["$tar"]}" "$match_text" >> "$bad_manifest"
        n_bad=$((n_bad + 1))
      else
        if bool_is_true "$run_rename"; then
          mv -- "$tar" "$target"
          action="renamed"
          n_renamed=$((n_renamed + 1))
        else
          action="would_rename"
          n_would=$((n_would + 1))
        fi
      fi
      printf '%s\t%s\t%s\t%s\t%s\n' "$base" "$(basename -- "$target")" "$action" "${tar_match_count["$tar"]}" "$match_text" >> "$manifest"
    else
      bad_target="$(unique_path "$BADNAME" "$base")"
      if [[ "${tar_match_count["$tar"]}" -eq 0 ]]; then
        action="no_match"
      else
        action="multiple_match"
      fi
      if bool_is_true "$run_rename"; then
        mkdir -p -- "$BADNAME"
        mv -- "$tar" "$bad_target"
        action="${action}_moved_badName"
      else
        action="${action}_would_move_badName"
      fi
      printf '%s\t%s\t%s\t%s\t%s\n' "$base" "$bad_target" "$action" "${tar_match_count["$tar"]}" "$match_text" >> "$bad_manifest"
      printf '%s\t%s\t%s\t%s\t%s\n' "$base" "$(basename -- "$bad_target")" "$action" "${tar_match_count["$tar"]}" "$match_text" >> "$manifest"
      n_bad=$((n_bad + 1))
    fi
  done

  task_log "rename scan done. renamed=$n_renamed; would_rename=$n_would; already_renamed=$n_skip; badName=$n_bad"
  task_log "manifest: $manifest"
  task_log "badName manifest: $bad_manifest"
}


# 🚩 Archive extraction
stream_tar_member_text() {
  local tar_path="$1" member="$2"
  if [[ "$member" == *.gz ]]; then
    tar -xOf "$tar_path" "$member" | gzip -dc
  else
    tar -xOf "$tar_path" "$member"
  fi
}

parse_unzip_args() {
  UNZIP_OVERWRITE="FALSE"
  UNZIP_JOBS="1"
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --overwrite)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --overwrite" >&2; exit 2; }
        UNZIP_OVERWRITE="$2"
        shift 2
        ;;
      --overwrite=*)
        UNZIP_OVERWRITE="${1#*=}"
        shift
        ;;
      --jobs)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --jobs" >&2; exit 2; }
        UNZIP_JOBS="$2"
        shift 2
        ;;
      --jobs=*)
        UNZIP_JOBS="${1#*=}"
        shift
        ;;
      *)
        echo "ERROR: unknown unzip argument: $1" >&2
        exit 2
        ;;
    esac
  done
  bool_is_true "$UNZIP_OVERWRITE" || true
  if ! [[ "$UNZIP_JOBS" =~ ^[0-9]+$ ]] || (( UNZIP_JOBS < 1 )); then
    echo "ERROR: --jobs must be a positive integer, got: $UNZIP_JOBS" >&2
    exit 2
  fi
}

unzip_one_tar() {
  local tar_path="$1" overwrite="$2" manifest="$3" failed="$4"
  local tar_dir base prot assay assay_l out_gz tmp_out chr chr_l member member_base member_base_l n_members size_bytes detail
  local -a members found output_chroms
  local -A chr_member

  tar_dir="$(dirname -- "$tar_path")"
  base="$(basename -- "$tar_path")"
  prot="${base%.tar}"
  assay="${assay_by_protein["$prot"]:-$prot}"
  assay_l="${assay,,}"
  out_gz="$tar_dir/$prot.raw.gz"

  if [[ -z "${known_protein[$prot]+x}" ]]; then
    detail="tar basename is not present in BED column 4; run rename first or inspect raw.badName"
    task_log "WARNING skip $base: $detail"
    printf '%s\t%s\t%s\t%s\n' "$prot" "$base" "unknown_protein_name" "$detail" >> "$failed"
    return 1
  fi

  if [[ -s "$out_gz" ]] && ! bool_is_true "$overwrite"; then
    size_bytes="$(stat -c '%s' -- "$out_gz")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$prot" "$base" "$out_gz" "skip_existing" "NA" "$size_bytes" >> "$manifest"
    return 0
  fi

  mapfile -t members < <(tar -tf "$tar_path" | grep -E '/?discovery_chr([0-9]+|X|Y)_' || true)
  detail=""
  output_chroms=()
  for chr in "${chroms[@]}"; do
    found=()
    chr_l="${chr,,}"
    for member in "${members[@]}"; do
      member_base="$(basename -- "$member")"
      member_base_l="${member_base,,}"
      if [[ "$member_base_l" == "discovery_chr${chr_l}_${assay_l}:"* || "$member_base_l" == "discovery_chr${chr_l}_${assay_l}_"* ]]; then
        found+=("$member")
      fi
    done
    if [[ ${#found[@]} -ne 1 ]]; then
      if [[ "$chr" == "Y" && ${#found[@]} -eq 0 ]]; then
        continue
      fi
      detail="${detail} chr${chr}=${#found[@]}"
    else
      chr_member["$chr"]="${found[0]}"
      output_chroms+=("$chr")
    fi
  done

  if [[ -n "$detail" ]]; then
    task_log "WARNING skip $base: expected exactly one discovery file per chr;$detail"
    printf '%s\t%s\t%s\t%s\n' "$prot" "$base" "bad_member_count" "$detail" >> "$failed"
    return 1
  fi

  tmp_out="$tar_dir/.${prot}.raw.gz.tmp.$$"
  rm -f -- "$tmp_out"
  n_members="${#output_chroms[@]}"
  task_log "UNZIP $base -> $out_gz"
  if ! {
    local first=1
    for chr in "${output_chroms[@]}"; do
      member="${chr_member[$chr]}"
      if [[ "$first" -eq 1 ]]; then
        stream_tar_member_text "$tar_path" "$member"
        first=0
      else
        stream_tar_member_text "$tar_path" "$member" | awk 'NR > 1'
      fi
    done
  } | gzip -c > "$tmp_out"; then
    rm -f -- "$tmp_out"
    task_log "WARNING failed to merge tar members: $base"
    printf '%s\t%s\t%s\t%s\n' "$prot" "$base" "merge_failed" "tar/gzip stream failed" >> "$failed"
    return 1
  fi

  mv -- "$tmp_out" "$out_gz"
  size_bytes="$(stat -c '%s' -- "$out_gz")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$prot" "$base" "$out_gz" "merged" "$n_members" "$size_bytes" >> "$manifest"
}

run_unzip_task() {
  parse_unzip_args "$@"

  for x in tar gzip awk find sort cut grep wc tr; do
    need_cmd "$x"
  done

  mkdir -p -- "$LOGDIR"
  TASK_LOG="$LOGDIR/$UNZIP_LOG_FILE"
  local manifest="$LOGDIR/$UNZIP_MANIFEST_FILE"
  local failed="$LOGDIR/$UNZIP_FAILED_FILE"
  : > "$TASK_LOG"
  printf 'protein\ttar_file\toutput_file\taction\tmember_count\tsize_bytes\n' > "$manifest"
  printf 'protein\ttar_file\terror\tdetail\n' > "$failed"

  task_log "unzip mode started; raw=$OUT; output=next_to_tar; overwrite=$UNZIP_OVERWRITE; jobs=$UNZIP_JOBS"
  ensure_prot_bed FALSE
  [[ -d "$OUT" ]] || { task_log "ERROR: raw directory does not exist: $OUT"; return 2; }

  local -A known_protein assay_by_protein
  local p bed_panel
  while IFS=$'\t' read -r _ _ _ p bed_panel; do
    known_protein["$p"]=1
    assay_by_protein["$p"]="$(bed_protein_to_assay "$p" "$bed_panel")"
  done < "$PROT_BED"

  local -a tars chroms
  chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
  mapfile -d '' tars < <(find "$OUT" -maxdepth 1 -type f -name '*.tar' -print0 | sort -z)
  if [[ ${#tars[@]} -eq 0 ]]; then
    task_log "WARNING no .tar files found under raw directory: $OUT"
    return 0
  fi

  local tar_path active=0 worker_fail=0 n_ok n_skip n_fail
  for tar_path in "${tars[@]}"; do
    unzip_one_tar "$tar_path" "$UNZIP_OVERWRITE" "$manifest" "$failed" &
    active=$((active + 1))
    if (( active >= UNZIP_JOBS )); then
      wait -n || worker_fail=$((worker_fail + 1))
      active=$((active - 1))
    fi
  done
  while (( active > 0 )); do
    wait -n || worker_fail=$((worker_fail + 1))
    active=$((active - 1))
  done

  n_ok="$(awk -F '\t' 'NR > 1 && $4 == "merged" {n++} END {print n + 0}' "$manifest")"
  n_skip="$(awk -F '\t' 'NR > 1 && $4 == "skip_existing" {n++} END {print n + 0}' "$manifest")"
  n_fail="$(awk 'NR > 1 {n++} END {print n + 0}' "$failed")"
  task_log "unzip done. merged=$n_ok; skipped_existing=$n_skip; failed_or_skipped_bad=$n_fail; worker_failures=$worker_fail"
  task_log "manifest: $manifest"
  task_log "failed: $failed"
  (( n_fail == 0 && worker_fail == 0 ))
}


# 🚩 rsID mapping
parse_rsid_map_args() {
  RSID_MAP_OVERWRITE="FALSE"
  RSID_MAP_DELETE_RAW="TRUE"
  RSID_MAP_JOBS=4
  RSID_MAP_THREADS=4
  # Keep external-sort scratch data on WSL's native ext4 filesystem. DrvFS
  # mounts such as /mnt/d and /mnt/h are much slower for this workload.
  RSID_MAP_TMPDIR=/tmp/ppp_rsid_map
  RSID_MAP_SORT_MEM=2G
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --overwrite)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --overwrite" >&2; exit 2; }
        RSID_MAP_OVERWRITE="$2"
        shift 2
        ;;
      --overwrite=*)
        RSID_MAP_OVERWRITE="${1#*=}"
        shift
        ;;
      --delete-raw)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --delete-raw" >&2; exit 2; }
        RSID_MAP_DELETE_RAW="$2"
        shift 2
        ;;
      --delete-raw=*)
        RSID_MAP_DELETE_RAW="${1#*=}"
        shift
        ;;
      --jobs)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --jobs" >&2; exit 2; }
        RSID_MAP_JOBS="$2"
        shift 2
        ;;
      --jobs=*)
        RSID_MAP_JOBS="${1#*=}"
        shift
        ;;
      --threads)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --threads" >&2; exit 2; }
        RSID_MAP_THREADS="$2"
        shift 2
        ;;
      --threads=*)
        RSID_MAP_THREADS="${1#*=}"
        shift
        ;;
      --tmpdir)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --tmpdir" >&2; exit 2; }
        RSID_MAP_TMPDIR="$2"
        shift 2
        ;;
      --tmpdir=*)
        RSID_MAP_TMPDIR="${1#*=}"
        shift
        ;;
      --sort-memory)
        [[ $# -ge 2 ]] || { echo "ERROR: missing value after --sort-memory" >&2; exit 2; }
        RSID_MAP_SORT_MEM="$2"
        shift 2
        ;;
      --sort-memory=*)
        RSID_MAP_SORT_MEM="${1#*=}"
        shift
        ;;
      *)
        echo "ERROR: unknown rsid_map argument: $1" >&2
        exit 2
        ;;
    esac
  done
  bool_is_true "$RSID_MAP_OVERWRITE" || true
  bool_is_true "$RSID_MAP_DELETE_RAW" || true
  bool_is_true "$RSID_MAP_CACHE_LOOKUP" || true
  if [[ "$RSID_MAP_JOBS" != "auto" ]] && { ! [[ "$RSID_MAP_JOBS" =~ ^[0-9]+$ ]] || (( RSID_MAP_JOBS < 1 )); }; then
    echo "ERROR: --jobs must be 'auto' or a positive integer, got: $RSID_MAP_JOBS" >&2
    exit 2
  fi
  if [[ "$RSID_MAP_THREADS" != "auto" ]] && { ! [[ "$RSID_MAP_THREADS" =~ ^[0-9]+$ ]] || (( RSID_MAP_THREADS < 1 )); }; then
    echo "ERROR: --threads must be 'auto' or a positive integer, got: $RSID_MAP_THREADS" >&2
    exit 2
  fi

  local cores
  cores="$(cpu_count)"
  if [[ "$RSID_MAP_JOBS" == "auto" ]]; then
    RSID_MAP_JOBS="$cores"
    if (( RSID_MAP_JOBS > 8 )); then
      RSID_MAP_JOBS=8
    fi
  fi
  if [[ "$RSID_MAP_THREADS" == "auto" ]]; then
    RSID_MAP_THREADS=$((cores / RSID_MAP_JOBS))
    if (( RSID_MAP_THREADS < 1 )); then
      RSID_MAP_THREADS=1
    fi
  fi
}

stage_rsid_map_lookup_cache() {
  bool_is_true "$RSID_MAP_CACHE_LOOKUP" || return 0

  local cache_dir cache_file cache_tmp
  cache_dir="$RSID_MAP_TMPDIR/lookup_cache"
  cache_file="$cache_dir/$(basename -- "$RSID_MAP_LOOKUP")"
  mkdir -p -- "$cache_dir"

  if [[ ! -s "$cache_file" ]]; then
    cache_tmp="$cache_file.tmp.$$"
    rm -f -- "$cache_tmp"
    task_log "caching rsid lookup on temp filesystem: $cache_file"
    if ! cp -- "$RSID_MAP_LOOKUP" "$cache_tmp"; then
      rm -f -- "$cache_tmp"
      task_log "ERROR: failed to cache rsid lookup: $cache_file"
      return 2
    fi
    mv -- "$cache_tmp" "$cache_file"
  else
    task_log "rsid lookup cache exists: $cache_file"
  fi

  RSID_MAP_LOOKUP="$cache_file"
}

rsid_map_chr_source_file() {
  local chr="$1"
  local old_nullglob
  local -a matches filtered

  old_nullglob="$(shopt -p nullglob || true)"
  shopt -s nullglob
  matches=( "$RSID_MAP_DIR"/*_chr${chr}_*.tsv )
  eval "$old_nullglob"

  filtered=()
  local f base
  for f in "${matches[@]}"; do
    base="$(basename -- "$f")"
    [[ "$base" == *"_ALL_"* ]] && continue
    filtered+=( "$f" )
  done

  if [[ ${#filtered[@]} -ne 1 ]]; then
    task_log "ERROR: expected exactly one rsid map file for chr${chr}, found ${#filtered[@]} under $RSID_MAP_DIR"
    return 2
  fi
  printf '%s\n' "${filtered[0]}"
}

ensure_rsid_map_all() {
  if [[ -s "$RSID_MAP_ALL" ]]; then
    task_log "rsid ALL map exists: $RSID_MAP_ALL"
    return 0
  fi

  local map_dir tmp chr src first
  local -a map_chroms
  map_dir="$(dirname -- "$RSID_MAP_ALL")"
  mkdir -p -- "$map_dir"
  tmp="$RSID_MAP_ALL.tmp.$$"
  rm -f -- "$tmp"

  map_chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X XY)
  task_log "building rsid ALL map: $RSID_MAP_ALL"
  first=1
  for chr in "${map_chroms[@]}"; do
    src="$(rsid_map_chr_source_file "$chr")" || { rm -f -- "$tmp"; return 2; }
    task_log "append chr${chr}: $src"
    if [[ "$first" -eq 1 ]]; then
      cat -- "$src" >> "$tmp"
      first=0
    else
      tail -n +2 -- "$src" >> "$tmp"
    fi
  done

  mv -- "$tmp" "$RSID_MAP_ALL"
  task_log "rsid ALL map built: $RSID_MAP_ALL"
}

ensure_rsid_map_lookup() {
  if [[ -s "$RSID_MAP_LOOKUP" ]]; then
    task_log "rsid lookup exists: $RSID_MAP_LOOKUP"
    return 0
  fi

  local tmp size_bytes
  mkdir -p -- "$RSID_MAP_TMPDIR"
  tmp="$RSID_MAP_LOOKUP.tmp.$$"
  rm -f -- "$tmp"

  task_log "building rsid lookup sorted by ID: $RSID_MAP_LOOKUP"
  if ! awk -F '\t' -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        sub(/\r$/, "", $i)
        col[$i] = i
      }
      if (!("ID" in col) || !("rsid" in col) || !("POS38" in col)) {
        print "ERROR: expected columns ID, rsid, POS38 in rsid map" > "/dev/stderr"
        exit 2
      }
      next
    }
    {
      id = $col["ID"]
      snp = $col["rsid"]
      pos = $col["POS38"]
      sub(/\r$/, "", id)
      sub(/\r$/, "", snp)
      sub(/\r$/, "", pos)
      if (id == "") next
      if (snp == "") snp = "NA"
      if (pos == "") pos = "NA"
      print id, snp, pos
    }
  ' "$RSID_MAP_ALL" | LC_ALL=C sort -T "$RSID_MAP_TMPDIR" -S "$RSID_MAP_SORT_MEM" -t $'\t' -k1,1 -u > "$tmp"; then
    rm -f -- "$tmp"
    task_log "ERROR: failed to build rsid lookup"
    return 2
  fi

  mv -- "$tmp" "$RSID_MAP_LOOKUP"
  size_bytes="$(stat -c '%s' -- "$RSID_MAP_LOOKUP")"
  task_log "rsid lookup built: $RSID_MAP_LOOKUP ($size_bytes bytes)"
}

rsid_map_one_raw() {
  local raw_gz="$1" overwrite="$2" delete_raw="$3" manifest="$4" failed="$5"
  local raw_dir raw_base prot out_gz tmp_out work_dir raw_sorted raw_count stats_file
  local header raw_rows mapped_rows na_rows size_bytes action

  raw_dir="$(dirname -- "$raw_gz")"
  raw_base="$(basename -- "$raw_gz")"
  prot="${raw_base%.raw.gz}"
  out_gz="$PROT_DIR/common/$prot/raw/$prot.gz"
  mkdir -p -- "$(dirname -- "$out_gz")"

  if [[ ! -e "$raw_gz" ]]; then
    action="skip_missing_raw"
    if [[ -s "$out_gz" ]]; then
      size_bytes="$(stat -c '%s' -- "$out_gz")"
    else
      size_bytes="NA"
    fi
    task_log "SKIP $raw_base: raw file missing; assuming rsid_map already completed"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$prot" "$raw_gz" "$out_gz" "$action" "NA" "NA" "NA" "$size_bytes" >> "$manifest"
    return 0
  fi

  if [[ -s "$out_gz" ]] && ! bool_is_true "$overwrite"; then
    action="skip_existing"
    if bool_is_true "$delete_raw"; then
      rm -f -- "$raw_gz"
      action="skip_existing_deleted_raw"
    fi
    size_bytes="$(stat -c '%s' -- "$out_gz")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$prot" "$raw_gz" "$out_gz" "$action" "NA" "NA" "NA" "$size_bytes" >> "$manifest"
    return 0
  fi

  header="$(gzip_decompress "$raw_gz" "$RSID_MAP_THREADS" 2>/dev/null | awk -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        sub(/\r$/, "", $i)
        printf "%s%s", (i == 1 ? "" : OFS), $i
      }
      print OFS "SNP" OFS "POS"
      exit
    }
  ' || true)"
  if [[ -z "$header" ]]; then
    task_log "WARNING: empty header in $raw_base"
    printf '%s\t%s\t%s\t%s\n' "$prot" "$raw_gz" "bad_raw_header" "empty header" >> "$failed"
    return 1
  fi

  work_dir="$(mktemp -d "$RSID_MAP_TMPDIR/rsid_map.XXXXXX")"
  raw_sorted="$work_dir/raw.idsort.tsv"
  raw_count="$work_dir/raw_rows.txt"
  stats_file="$work_dir/stats.tsv"
  tmp_out="$raw_dir/.${prot}.gz.tmp.$$"
  rm -f -- "$tmp_out"

  task_log "RSID_MAP $raw_base -> $out_gz"
  if ! gzip_decompress "$raw_gz" "$RSID_MAP_THREADS" | awk -v OFS='\t' -v count_file="$raw_count" '
    NR == 1 { next }
    {
      id = (NF >= 3 ? $3 : "")
      sub(/\r$/, "", id)
      printf "%s\t%d", id, NR - 1
      for (i = 1; i <= NF; i++) {
        sub(/\r$/, "", $i)
        printf "\t%s", $i
      }
      printf "\n"
      n++
    }
    END {
      print n + 0 > count_file
    }
  ' | LC_ALL=C sort -T "$work_dir" -S "$RSID_MAP_SORT_MEM" -t $'\t' -k1,1 > "$raw_sorted"; then
    rm -rf -- "$work_dir"
    rm -f -- "$tmp_out"
    task_log "WARNING: failed to sort raw file by ID: $raw_base"
    printf '%s\t%s\t%s\t%s\n' "$prot" "$raw_gz" "raw_sort_failed" "gzip/awk/sort failed" >> "$failed"
    return 1
  fi

  raw_rows="$(cat "$raw_count")"
  if ! {
    printf '%s\n' "$header"
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "$raw_sorted" "$RSID_MAP_LOOKUP" |
      awk -F '\t' -v OFS='\t' -v stats="$stats_file" '
        {
          if (NF >= 2 && $(NF - 1) == "NA" && $NF == "NA") {
            na++
          } else {
            mapped++
          }
          print
        }
        END {
          print mapped + 0, na + 0 > stats
        }
      ' |
      LC_ALL=C sort -T "$work_dir" -S "$RSID_MAP_SORT_MEM" -t $'\t' -k2,2n |
      cut -f3-
  } | gzip_compress "$RSID_MAP_THREADS" > "$tmp_out"; then
    rm -rf -- "$work_dir"
    rm -f -- "$tmp_out"
    task_log "WARNING: failed to merge rsid map: $raw_base"
    printf '%s\t%s\t%s\t%s\n' "$prot" "$raw_gz" "merge_failed" "join/sort/gzip failed" >> "$failed"
    return 1
  fi

  mapped_rows="NA"
  na_rows="NA"
  if [[ -s "$stats_file" ]]; then
    read -r mapped_rows na_rows < "$stats_file"
  fi

  mv -- "$tmp_out" "$out_gz"
  size_bytes="$(stat -c '%s' -- "$out_gz")"
  action="mapped"
  if bool_is_true "$delete_raw"; then
    rm -f -- "$raw_gz"
    action="mapped_deleted_raw"
  fi

  rm -rf -- "$work_dir"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$prot" "$raw_gz" "$out_gz" "$action" "$raw_rows" "$mapped_rows" "$na_rows" "$size_bytes" >> "$manifest"
}

run_rsid_map_task() {
  parse_rsid_map_args "$@"

  for x in gzip awk sort join cut wc tr find tail cat mktemp cp stat; do
    need_cmd "$x"
  done

  mkdir -p -- "$LOGDIR" "$RSID_MAP_TMPDIR"
  TASK_LOG="$LOGDIR/$RSID_MAP_LOG_FILE"
  local manifest="$LOGDIR/$RSID_MAP_MANIFEST_FILE"
  local failed="$LOGDIR/$RSID_MAP_FAILED_FILE"
  : > "$TASK_LOG"
  printf 'protein\traw_file\toutput_file\taction\traw_rows\tmapped_rows\tna_rows\tsize_bytes\n' > "$manifest"
  printf 'protein\traw_file\terror\tdetail\n' > "$failed"

  local compressor="gzip"
  if command -v pigz >/dev/null 2>&1; then
    compressor="pigz"
  fi
  task_log "rsid_map mode started; raw=$OUT; rsid_map_dir=$RSID_MAP_DIR; overwrite=$RSID_MAP_OVERWRITE; delete_raw=$RSID_MAP_DELETE_RAW; cpu=$(cpu_count); jobs=$RSID_MAP_JOBS; threads_per_job=$RSID_MAP_THREADS; compressor=$compressor; tmpdir=$RSID_MAP_TMPDIR; sort_memory=$RSID_MAP_SORT_MEM"
  [[ -d "$OUT" ]] || { task_log "ERROR: raw directory does not exist: $OUT"; return 2; }
  [[ -d "$RSID_MAP_DIR" ]] || { task_log "ERROR: rsid map directory does not exist: $RSID_MAP_DIR"; return 2; }

  ensure_rsid_map_all
  ensure_rsid_map_lookup
  stage_rsid_map_lookup_cache

  local -a raw_files
  mapfile -d '' raw_files < <(find "$OUT" -maxdepth 1 -type f -name '*.raw.gz' -print0 | sort -z)

  local raw_gz
  local active=0 worker_fail=0 n_fail n_done n_skip n_missing

  if [[ ${#raw_files[@]} -eq 0 ]]; then
    task_log "no .raw.gz files found under raw directory; nothing to process"
    task_log "rsid_map done. mapped=0; skipped_existing=0; skipped_missing_raw=0; failed=0; worker_failures=0"
    task_log "manifest: $manifest"
    task_log "failed: $failed"
    return 0
  fi

  for raw_gz in "${raw_files[@]}"; do
    rsid_map_one_raw "$raw_gz" "$RSID_MAP_OVERWRITE" "$RSID_MAP_DELETE_RAW" "$manifest" "$failed" &
    active=$((active + 1))
    if (( active >= RSID_MAP_JOBS )); then
      wait -n || worker_fail=$((worker_fail + 1))
      active=$((active - 1))
    fi
  done
  while (( active > 0 )); do
    wait -n || worker_fail=$((worker_fail + 1))
    active=$((active - 1))
  done

  n_done="$(awk -F '\t' 'NR > 1 && $4 ~ /^mapped/ {n++} END {print n + 0}' "$manifest")"
  n_skip="$(awk -F '\t' 'NR > 1 && $4 ~ /^skip_existing/ {n++} END {print n + 0}' "$manifest")"
  n_missing="$(awk -F '\t' 'NR > 1 && $4 == "skip_missing_raw" {n++} END {print n + 0}' "$manifest")"
  n_fail="$(awk 'NR > 1 {n++} END {print n + 0}' "$failed")"
  task_log "rsid_map done. mapped=$n_done; skipped_existing=$n_skip; skipped_missing_raw=$n_missing; failed=$n_fail; worker_failures=$worker_fail"
  task_log "manifest: $manifest"
  task_log "failed: $failed"
  (( n_fail == 0 && worker_fail == 0 ))
}


# 🚩 Dispatch
MODE="${1:-}"
if [[ -z "$MODE" || "$MODE" == "-h" || "$MODE" == "--help" ]]; then
  usage
  exit 0
fi
shift || true

TOKEN_FILE=$TOKEN_FILE_DEFAULT
OUT=$OUT_DEFAULT
LOGDIR=$LOG_DEFAULT
PROT_DIR=$PROT_DIR_DEFAULT
PROT_TSV=""
PROT_BED=""
CLEAN=""
BADNAME=""
PROT_FOLDER=$FOLDER_DEFAULT
NAME_REGEX=$NAME_REGEX_DEFAULT
RSID_MAP_DIR=""
RSID_MAP_ALL=""
RSID_MAP_LOOKUP=""
RSID_MAP_CACHE_LOOKUP=TRUE
mode_args=()

while (( $# )); do
  case "$1" in
    --token-file) TOKEN_FILE=${2:?ERROR: --token-file requires FILE}; shift 2 ;;
    --output-dir) OUT=${2:?ERROR: --output-dir requires DIR}; shift 2 ;;
    --log-dir) LOGDIR=${2:?ERROR: --log-dir requires DIR}; shift 2 ;;
    --prot-dir) PROT_DIR=${2:?ERROR: --prot-dir requires DIR}; shift 2 ;;
    --prot-tsv) PROT_TSV=${2:?ERROR: --prot-tsv requires FILE}; shift 2 ;;
    --prot-bed) PROT_BED=${2:?ERROR: --prot-bed requires FILE}; shift 2 ;;
    --folder) PROT_FOLDER=${2:?ERROR: --folder requires ID}; shift 2 ;;
    --name-regex) NAME_REGEX=${2:?ERROR: --name-regex requires REGEX}; shift 2 ;;
    --python) PYTHON_BIN=${2:?ERROR: --python requires FILE}; shift 2 ;;
    --clean-dir) CLEAN=${2:?ERROR: --clean-dir requires DIR}; shift 2 ;;
    --bad-name-dir) BADNAME=${2:?ERROR: --bad-name-dir requires DIR}; shift 2 ;;
    --rsid-map-dir) RSID_MAP_DIR=${2:?ERROR: --rsid-map-dir requires DIR}; shift 2 ;;
    --rsid-map-all) RSID_MAP_ALL=${2:?ERROR: --rsid-map-all requires FILE}; shift 2 ;;
    --rsid-map-lookup) RSID_MAP_LOOKUP=${2:?ERROR: --rsid-map-lookup requires FILE}; shift 2 ;;
    --cache-lookup)
      RSID_MAP_CACHE_LOOKUP=${2:?ERROR: --cache-lookup requires TRUE or FALSE}
      bool_is_true "$RSID_MAP_CACHE_LOOKUP" || true
      shift 2
      ;;
    *) mode_args+=("$1"); shift ;;
  esac
done
set -- "${mode_args[@]}"

PROT_TSV=${PROT_TSV:-$PROT_TSV_DEFAULT}
PROT_BED=${PROT_BED:-$(dirname -- "$PROT_TSV")/ppp_3k.38.bed}
CLEAN=${CLEAN:-$PROT_DIR}
BADNAME=${BADNAME:-${OUT%/}.badName}
RSID_MAP_DIR=${RSID_MAP_DIR:-$PROT_DIR/rsid_map}
RSID_MAP_ALL=${RSID_MAP_ALL:-$RSID_MAP_DIR/$RSID_MAP_ALL_FILE}
RSID_MAP_LOOKUP=${RSID_MAP_LOOKUP:-$RSID_MAP_DIR/$RSID_MAP_LOOKUP_FILE}
COMMON_CLI=(
  --token-file "$TOKEN_FILE" --output-dir "$OUT" --log-dir "$LOGDIR"
  --prot-dir "$PROT_DIR" --prot-tsv "$PROT_TSV" --prot-bed "$PROT_BED"
  --folder "$PROT_FOLDER" --name-regex "$NAME_REGEX" --clean-dir "$CLEAN"
  --bad-name-dir "$BADNAME" --rsid-map-dir "$RSID_MAP_DIR"
  --rsid-map-all "$RSID_MAP_ALL" --rsid-map-lookup "$RSID_MAP_LOOKUP"
  --cache-lookup "$RSID_MAP_CACHE_LOOKUP"
)
[[ -z $PYTHON_BIN ]] || COMMON_CLI+=(--python "$PYTHON_BIN")

mkdir -p "$LOGDIR"

COMMON_ARGS=(
  --out "$OUT"
  --folder "$PROT_FOLDER"
  --name-regex "$NAME_REGEX"
)

case "$MODE" in
  rename|unzip|rsid_map|dryrun|foreground|download)
    check_existing_prot_tasks
    ;;
esac

case "$MODE" in
  rename)
    run_rename_task "$@"
    ;;
  unzip)
    parse_unzip_args "$@"
    echo "Starting background unzip. Log: $LOGDIR/$UNZIP_LOG_FILE"
    nohup bash "$SCRIPT_DIR/ppp_download.sh" unzip-foreground "${COMMON_CLI[@]}" "$@" >/dev/null 2>&1 &
    echo "PID: $!"
    echo "Progress: tail -f $LOGDIR/$UNZIP_LOG_FILE"
    echo "Job: pgrep -af 'ppp_download.py|ppp_download.sh'"
    ;;
  unzip-foreground)
    run_unzip_task "$@"
    ;;
  rsid_map)
    parse_rsid_map_args "$@"
    echo "Starting background rsid_map. Log: $LOGDIR/$RSID_MAP_LOG_FILE"
    nohup bash "$SCRIPT_DIR/ppp_download.sh" rsid_map-foreground "${COMMON_CLI[@]}" "$@" >/dev/null 2>&1 &
    echo "PID: $!"
    echo "Progress: tail -f $LOGDIR/$RSID_MAP_LOG_FILE"
    echo "Job: pgrep -af 'ppp_download.py|ppp_download.sh'"
    ;;
  rsid_map-foreground)
    run_rsid_map_task "$@"
    ;;
  test)
    if [[ ! -s "$PY" ]]; then
      echo "ERROR: cannot find Python script: $PY" >&2
      exit 2
    fi
    if [[ -z "$PYTHON_BIN" ]]; then
      if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
      elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
      else
        echo "ERROR: cannot find python or python3 in PATH. Use --python /path/to/python." >&2
        exit 2
      fi
    fi
    if [[ ! -s "$TOKEN_FILE" ]]; then
      echo "ERROR: token file not found or empty: $TOKEN_FILE" >&2
      exit 2
    fi
    export SYNAPSE_AUTH_TOKEN="$(tr -d '\r\n\t '\''\"' < "$TOKEN_FILE")"
    "$PYTHON_BIN" "$PY" --check-login
    ;;
  dryrun)
    if [[ ! -s "$PY" ]]; then
      echo "ERROR: cannot find Python script: $PY" >&2
      exit 2
    fi
    if [[ -z "$PYTHON_BIN" ]]; then
      if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
      elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
      else
        echo "ERROR: cannot find python or python3 in PATH. Use --python /path/to/python." >&2
        exit 2
      fi
    fi
    if [[ ! -s "$TOKEN_FILE" ]]; then
      echo "ERROR: token file not found or empty: $TOKEN_FILE" >&2
      exit 2
    fi
    export SYNAPSE_AUTH_TOKEN="$(tr -d '\r\n\t '\''\"' < "$TOKEN_FILE")"
    "$PYTHON_BIN" "$PY" \
      "${COMMON_ARGS[@]}" \
      --manifest "$MANIFEST_DRYRUN" \
      --dry-run \
      "$@"
    ;;
  foreground)
    if [[ ! -s "$PY" ]]; then
      echo "ERROR: cannot find Python script: $PY" >&2
      exit 2
    fi
    if [[ -z "$PYTHON_BIN" ]]; then
      if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
      elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
      else
        echo "ERROR: cannot find python or python3 in PATH. Use --python /path/to/python." >&2
        exit 2
      fi
    fi
    if [[ ! -s "$TOKEN_FILE" ]]; then
      echo "ERROR: token file not found or empty: $TOKEN_FILE" >&2
      exit 2
    fi
    export SYNAPSE_AUTH_TOKEN="$(tr -d '\r\n\t '\''\"' < "$TOKEN_FILE")"
    # Default workers=4 for foreground testing unless user supplies --workers.
    mapfile -d '' EXTRA_ARGS < <(add_default_workers_if_missing 4 "$@")
    "$PYTHON_BIN" -u "$PY" \
      "${COMMON_ARGS[@]}" \
      --manifest "$MANIFEST_DOWNLOAD" \
      --retries 5 \
      "${EXTRA_ARGS[@]}"
    ;;
  download)
    if [[ ! -s "$PY" ]]; then
      echo "ERROR: cannot find Python script: $PY" >&2
      exit 2
    fi
    if [[ -z "$PYTHON_BIN" ]]; then
      if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
      elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
      else
        echo "ERROR: cannot find python or python3 in PATH. Use --python /path/to/python." >&2
        exit 2
      fi
    fi
    if [[ ! -s "$TOKEN_FILE" ]]; then
      echo "ERROR: token file not found or empty: $TOKEN_FILE" >&2
      exit 2
    fi
    export SYNAPSE_AUTH_TOKEN="$(tr -d '\r\n\t '\''\"' < "$TOKEN_FILE")"
    # Default workers=8 for full background download unless user supplies --workers.
    mapfile -d '' EXTRA_ARGS < <(add_default_workers_if_missing 8 "$@")
    echo "Starting background download. Log: $LOGDIR/$LOG_FILE"
    nohup "$PYTHON_BIN" -u "$PY" \
      "${COMMON_ARGS[@]}" \
      --manifest "$MANIFEST_DOWNLOAD" \
      --retries 5 \
      "${EXTRA_ARGS[@]}" \
      > "$LOGDIR/$LOG_FILE" 2>&1 &
    echo "PID: $!"
    echo "Progress: tail -f $LOGDIR/$LOG_FILE"
    echo "Job: pgrep -af 'ppp_download.py|ppp_download.sh'"
    ;;
  *)
    echo "ERROR: unknown mode: $MODE" >&2
    usage
    exit 2
    ;;
esac
