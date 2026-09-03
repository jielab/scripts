#!/usr/bin/env bash


# 🚩 Publish selected analysis results
#
# Usage:
#   ./github_copy.sh scan project,...    # rebuild DEST_ROOT/github_files.lst
#   ./github_copy.sh sync [project,...]  # copy approved files and remove stale ones

set -Eeuo pipefail
IFS=$'\n\t'

# These can be overridden for testing, for example:
#   GITHUB_COPY_SOURCE_ROOT=/tmp/source GITHUB_COPY_DEST_ROOT=/tmp/dest ./github_copy.sh scan le8
SOURCE_ROOT=${GITHUB_COPY_SOURCE_ROOT:-/mnt/d/analysis}
DEST_ROOT=${GITHUB_COPY_DEST_ROOT:-/mnt/d/github/analysis}
MANIFEST_NAME=github_files.lst

# GitHub rejects files above 100 MiB. Use a deliberately conservative limit.
MAX_FILE_BYTES=$((50 * 1024 * 1024))
# Compressed files are harder to inspect. Only small .gz files in gu/normalize
# are accepted, preserving the small normalized artifacts already in that tree.
MAX_GZIP_BYTES=$((10 * 1024 * 1024))

SRC_REAL=
DST_REAL=
REJECT_REASON=
FILE_SIZE=0
TEMP_FILE=
UKB_INPUT_FILE=
UKB_RESULT_FILE=
PROJECT_FILTER_ACTIVE=0
PROJECT_SPEC=
MANIFEST_PROJECT_SPEC=
declare -a SELECTED_PROJECT_LIST=()
declare -A SELECTED_PROJECT_SET=()

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

cleanup() {
    if [[ -n ${TEMP_FILE:-} && -e $TEMP_FILE ]]; then
        rm -f -- "$TEMP_FILE"
    fi
    if [[ -n ${UKB_INPUT_FILE:-} && -e $UKB_INPUT_FILE ]]; then
        rm -f -- "$UKB_INPUT_FILE"
    fi
    if [[ -n ${UKB_RESULT_FILE:-} && -e $UKB_RESULT_FILE ]]; then
        rm -f -- "$UKB_RESULT_FILE"
    fi
}
trap cleanup EXIT

usage() {
    cat <<'EOF'
Usage:
  github_copy.sh scan project,...
  github_copy.sh sync [project,...]

  scan  Examine only the explicitly named top-level datasets and atomically
        rebuild /mnt/d/github/analysis/github_files.lst. A dataset list is
        mandatory; scan never defaults to every folder.

  sync  Revalidate every manifest entry, copy changed files, and remove
        stale destination files within the selected dataset scope. Without a
        list, sync uses the scope recorded by the preceding scan. An explicit
        sync list must be a subset of that recorded scope.
        The destination .git directory and root-level repository files are kept.

Project selection:
  github_copy.sh scan le8,maha
  github_copy.sh sync le8,maha

  Only the named top-level dataset folders are scanned or synchronized. The ukb
  dataset cannot be selected and is never cleaned by sync.

Publish rules:
  le8/<group>/<analysis>/*.{csv,tsv,xlsx,png}
  le8/<group>/<analysis>/*/*.{csv,tsv,xlsx,png}
      Direct outputs and one child result directory (for example c1_correlate)
      are included below analysis roots such as cvd_cad/met and cvd_cad/prot.
  other datasets: the dataset root plus first- and second-level directories,
      limited to csv, tsv, xlsx and png result files.
  gu/normalize/** keeps its existing safe-file rule; small .gz files are allowed.
  nexstrain/* and nextstrain/* keep their existing root-file rule.

Hard safety rules are applied again during sync: every ukb directory is blocked,
symlinks and hidden paths are blocked, R data files/secrets/archives are blocked,
oversized files are blocked, and cache/temporary paths are blocked. Filenames
containing "ukb" are allowed when they are aggregate result files; only a path
component named exactly ukb is treated as the protected data directory.

After copying, sync checks non-exempt destination tables for high-confidence UKB
participant-data patterns. A detected file is logged, deleted from destination,
and removed from github_files.lst. The gu, grid and ems120 datasets skip this UKB
content check, but still receive every general safety check above.

Directories named mrlink2 or dandelion_network, plus cache and temporary
directories, are pruned during scanning and are not eligible for publication.
EOF
}

require_commands() {
    local command_name
    for command_name in find sort stat realpath mktemp mv cp cmp mkdir dirname rm chmod tr; do
        command -v "$command_name" >/dev/null 2>&1 || die "required command not found: $command_name"
    done
}

prepare_roots() {
    require_commands

    [[ $SOURCE_ROOT == /* ]] || die "SOURCE_ROOT must be an absolute path: $SOURCE_ROOT"
    [[ $DEST_ROOT == /* ]] || die "DEST_ROOT must be an absolute path: $DEST_ROOT"
    [[ -d $SOURCE_ROOT && ! -L $SOURCE_ROOT ]] || die "source is not a real directory: $SOURCE_ROOT"
    [[ ! -L $DEST_ROOT ]] || die "destination root must not be a symlink: $DEST_ROOT"

    mkdir -p -- "$DEST_ROOT"
    [[ -d $DEST_ROOT && ! -L $DEST_ROOT ]] || die "cannot create destination: $DEST_ROOT"

    SRC_REAL=$(realpath -e -- "$SOURCE_ROOT")
    DST_REAL=$(realpath -e -- "$DEST_ROOT")

    [[ $SRC_REAL != / ]] || die 'refusing to use / as the source root'
    [[ $DST_REAL != / ]] || die 'refusing to use / as the destination root'
    [[ $SRC_REAL != "$DST_REAL" ]] || die 'source and destination resolve to the same directory'
    [[ $DST_REAL != "$SRC_REAL"/* ]] || die 'destination must not be inside the source directory'
    [[ $SRC_REAL != "$DST_REAL"/* ]] || die 'source must not be inside the destination directory'
}

configure_project_filter() {
    local spec=${1-}
    local project
    local -a requested=()

    PROJECT_FILTER_ACTIVE=0
    PROJECT_SPEC=
    SELECTED_PROJECT_LIST=()
    SELECTED_PROJECT_SET=()

    [[ -n $spec ]] || return 0

    PROJECT_FILTER_ACTIVE=1
    PROJECT_SPEC=$spec
    IFS=, read -r -a requested <<< "$spec"
    (( ${#requested[@]} > 0 )) || die 'project list must not be empty'

    for project in "${requested[@]}"; do
        if [[ -z $project || $project == . || $project == .. || $project == .* ||
              $project == */* || $project == *\\* || $project == *$'\n'* ||
              $project == *$'\r'* || $project == *$'\t'* ]]; then
            die "invalid top-level project name: $project"
        fi
        [[ ${project,,} != ukb ]] || die 'the protected ukb project cannot be selected'
        [[ -d $SRC_REAL/$project && ! -L $SRC_REAL/$project ]] ||
            die "selected project is not a real source directory: $SOURCE_ROOT/$project"

        if [[ ! ${SELECTED_PROJECT_SET[$project]+_} ]]; then
            SELECTED_PROJECT_SET["$project"]=1
            SELECTED_PROJECT_LIST+=("$project")
        fi
    done
}

project_is_selected() {
    local project=$1

    [[ ${project,,} != ukb ]] || return 1
    if (( PROJECT_FILTER_ACTIVE )); then
        [[ ${SELECTED_PROJECT_SET[$project]+_} ]]
        return
    fi
    return 0
}

extension_of() {
    local base=${1##*/}
    if [[ $base == *.* ]]; then
        printf '%s' "${base##*.}" | tr '[:upper:]' '[:lower:]'
    fi
}

is_cache_or_temporary_component() {
    local lower=${1,,}

    [[ $lower == *cache* ]] && return 0
    case "$lower" in
        tmp|temp|temporary|*.tmp|*.temp|*.temporary|*.tmp.*|*.temp.*|*.temporary.*|\
        tmp.*|temp.*|temporary.*|*_tmp|*_temp|*_temporary|*_tmp_*|*_temp_*|*_temporary_*|\
        tmp_*|temp_*|temporary_*|*-tmp|*-temp|*-temporary|*-tmp-*|*-temp-*|*-temporary-*|\
        tmp-*|temp-*|temporary-*)
            return 0
            ;;
    esac
    return 1
}

matches_result_extension() {
    local ext=$1
    [[ $ext == csv || $ext == tsv || $ext == xlsx || $ext == png ]]
}


# 🚩 Publishing whitelist
# Adding a project should normally mean
# adding one narrow rule here, not weakening the safety checks below.
matches_publish_rule() {
    local rel=$1
    local lower=${rel,,}
    local ext part
    local -a parts
    IFS=/ read -r -a parts <<< "$rel"
    ext=$(extension_of "$rel")

    for part in "${parts[@]}"; do
        if [[ ${part,,} == mrlink2 || ${part,,} == dandelion_network ]]; then
            return 1
        fi
        if is_cache_or_temporary_component "$part"; then
            return 1
        fi
    done

    if [[ $lower == gu/normalize/* && ${#parts[@]} -ge 3 ]]; then
        return 0
    fi

    if [[ ${parts[0],,} == le8 &&
          ( ${#parts[@]} -eq 4 || ${#parts[@]} -eq 5 ) ]]; then
        matches_result_extension "$ext"
        return
    fi

    # Accept the actual source spelling and the conventional spelling.
    if [[ ( ${parts[0],,} == nexstrain || ${parts[0],,} == nextstrain ) && ${#parts[@]} -eq 2 ]]; then
        return 0
    fi

    # Other datasets may publish result files from the dataset root and from
    # first- or second-level result directories. le8 has one extra grouping
    # component and is handled by its explicit rule above.
    if [[ ${#parts[@]} -ge 2 && ${#parts[@]} -le 4 ]]; then
        matches_result_extension "$ext"
        return
    fi

    return 1
}


# 🚩 File safety checks
# Sets REJECT_REASON and FILE_SIZE. A nonzero return means the file is unsafe.
check_file_safety() {
    local rel=$1
    local source_path=$2
    local lower_base ext part current
    local -a parts

    REJECT_REASON=
    FILE_SIZE=0

    if [[ -z $rel || $rel == /* || $rel == *\\* || $rel == *$'\n'* || $rel == *$'\r'* || $rel == *$'\t'* ]]; then
        REJECT_REASON=unsafe_path
        return 1
    fi

    IFS=/ read -r -a parts <<< "$rel"
    if [[ ${#parts[@]} -lt 2 ]]; then
        REJECT_REASON=not_in_project_folder
        return 1
    fi

    for part in "${parts[@]}"; do
        if [[ -z $part || $part == . || $part == .. ]]; then
            REJECT_REASON=unsafe_path
            return 1
        fi
        if [[ ${part,,} == ukb ]]; then
            REJECT_REASON=ukb_directory
            return 1
        fi
        if [[ $part == .* ]]; then
            REJECT_REASON=hidden_path
            return 1
        fi
        if is_cache_or_temporary_component "$part"; then
            REJECT_REASON=cache_or_temporary_path
            return 1
        fi
    done

    current=$SRC_REAL
    for part in "${parts[@]}"; do
        current=$current/$part
        if [[ -L $current ]]; then
            REJECT_REASON=symlink
            return 1
        fi
    done
    if [[ ! -f $source_path ]]; then
        REJECT_REASON=not_a_regular_file
        return 1
    fi

    lower_base=${parts[${#parts[@]}-1],,}
    ext=$(extension_of "$lower_base")

    case "$lower_base" in
        .env|credentials|credentials.*|secrets|secrets.*|id_rsa|id_dsa|id_ecdsa|id_ed25519)
            REJECT_REASON=secret_filename
            return 1
            ;;
    esac

    case "$ext" in
        rds|rda|rdata|sav|dta|sas7bdat|parquet|feather|fst|sqlite|sqlite3|db|mdb|accdb|pem|key|p12|pfx|kdbx|gpg|age|zip|7z|rar|tar|tgz|bz2|xz|zst)
            REJECT_REASON=blocked_file_type
            return 1
            ;;
    esac

    FILE_SIZE=$(stat -c '%s' -- "$source_path") || {
        REJECT_REASON=cannot_read_size
        return 1
    }
    if (( FILE_SIZE > MAX_FILE_BYTES )); then
        REJECT_REASON=over_50_MiB
        return 1
    fi

    if [[ $ext == gz ]]; then
        if [[ ${rel,,} != gu/normalize/* ]]; then
            REJECT_REASON=gzip_outside_gu_normalize
            return 1
        fi
        if (( FILE_SIZE > MAX_GZIP_BYTES )); then
            REJECT_REASON=gzip_over_10_MiB
            return 1
        fi
    fi

    return 0
}

find_scan_candidates() {
    local root=$1
    local max_depth=$2
    local -a depth_args=()

    if (( max_depth > 0 )); then
        depth_args=(-maxdepth "$max_depth")
    fi

    find "$root" -mindepth 1 "${depth_args[@]}" \
        \( -type d \( -iname ukb -o -iname mrlink2 -o -iname dandelion_network -o \
            -iname '*cache*' -o -iname tmp -o -iname temp -o -iname temporary -o \
            -iname '*.tmp' -o -iname '*.temp' -o -iname '*.temporary' -o \
            -iname '*_tmp' -o -iname '*_temp' -o -iname '*_temporary' -o \
            -iname '*-tmp' -o -iname '*-temp' -o -iname '*-temporary' \) -prune \) -o \
        \( \( -type f -o -type l \) -print0 \)
}

scan_files() {
    local project_spec=$1
    local source_path rel reason
    local -a approved=()
    local -A approved_set=()
    local -A skipped=()

    prepare_roots
    configure_project_filter "$project_spec"

    # Do not descend into protected data, caches, temporary directories, or
    # excluded result-tool directories. Depth is capped at the publishable
    # result levels so a scan does not walk large raw/cache trees needlessly.
    while IFS= read -r -d '' source_path; do
        rel=${source_path#"$SRC_REAL"/}

        if ! matches_publish_rule "$rel"; then
            (( skipped[not_publish_result] += 1 ))
            continue
        fi
        if ! check_file_safety "$rel" "$source_path"; then
            reason=$REJECT_REASON
            (( skipped[$reason] += 1 ))
            continue
        fi
        if [[ ! ${approved_set[$rel]+_} ]]; then
            approved_set["$rel"]=1
            approved+=("$rel")
        fi
    done < <(
        local project lower_project normalize_root
        for project in "${SELECTED_PROJECT_LIST[@]}"; do
            lower_project=${project,,}
            case "$lower_project" in
                le8)
                    find_scan_candidates "$SRC_REAL/$project" 4
                    ;;
                gu)
                    find_scan_candidates "$SRC_REAL/$project" 3
                    normalize_root=$SRC_REAL/$project/normalize
                    if [[ -d $normalize_root && ! -L $normalize_root ]]; then
                        find_scan_candidates "$normalize_root" 0
                    fi
                    ;;
                nexstrain|nextstrain)
                    find_scan_candidates "$SRC_REAL/$project" 1
                    ;;
                *)
                    find_scan_candidates "$SRC_REAL/$project" 3
                    ;;
            esac
        done
    )

    TEMP_FILE=$(mktemp "$DST_REAL/.github_files.lst.tmp.XXXXXX")
    {
        printf '# Generated by github_copy.sh scan. Paths are relative to %s.\n' "$SOURCE_ROOT"
        printf '# Project scope: %s\n' "$PROJECT_SPEC"
        printf '# Review this allowlist before running github_copy.sh sync.\n'
        if (( ${#approved[@]} > 0 )); then
            printf '%s\n' "${approved[@]}" | LC_ALL=C sort -u
        fi
    } > "$TEMP_FILE"
    chmod 0644 "$TEMP_FILE"
    mv -f -- "$TEMP_FILE" "$DST_REAL/$MANIFEST_NAME"
    TEMP_FILE=

    printf 'SCAN OK: wrote %d approved paths to %s\n' \
        "${#approved[@]}" "$DST_REAL/$MANIFEST_NAME"
    printf 'Safety note: ukb, mrlink2, dandelion_network, cache, and temporary directories were pruned without inspection.\n'
    if (( ${#skipped[@]} > 0 )); then
        printf 'Skipped files by reason:\n'
        for reason in "${!skipped[@]}"; do
            printf '  %-28s %d\n' "$reason" "${skipped[$reason]}"
        done
    fi
}

check_destination_target() {
    local rel=$1
    local current=$DST_REAL
    local target=$DST_REAL/$rel
    local -a parts
    local i

    IFS=/ read -r -a parts <<< "$rel"
    for ((i = 0; i < ${#parts[@]} - 1; i++)); do
        current=$current/${parts[i]}
        [[ ! -L $current ]] || die "destination parent is a symlink: $current"
        [[ ! -e $current || -d $current ]] || die "destination parent is not a directory: $current"
    done

    [[ ! -L $target ]] || die "destination file is a symlink: $target"
    [[ ! -d $target ]] || die "destination file path is a directory: $target"
}

read_manifest_project_spec() {
    local manifest=$1
    local line value
    local found=0

    MANIFEST_PROJECT_SPEC=
    while IFS= read -r line || [[ -n $line ]]; do
        [[ $line == '# Project scope: '* ]] || continue
        value=${line#'# Project scope: '}
        [[ -n $value && $value != 'all top-level projects except ukb.' ]] ||
            die 'manifest has an obsolete or empty project scope; run scan with an explicit dataset list'
        ((found += 1))
        (( found == 1 )) || die 'manifest contains more than one project scope header'
        MANIFEST_PROJECT_SPEC=$value
    done < "$manifest"

    (( found == 1 )) || die 'manifest has no project scope header; run scan again'
}

manifest_scope_has_project() {
    local wanted_project=$1
    local project
    local -a manifest_projects=()

    IFS=, read -r -a manifest_projects <<< "$MANIFEST_PROJECT_SPEC"
    for project in "${manifest_projects[@]}"; do
        [[ $project == "$wanted_project" ]] && return 0
    done
    return 1
}

project_skips_ukb_content_check() {
    case "${1,,}" in
        gu|grid|ems120) return 0 ;;
        *) return 1 ;;
    esac
}


# 🚩 UKB content scan
# Scan copied destination files, not source files. The detector deliberately
# emits only a relative path and a reason; participant-like values are never
# printed to the terminal or written to the manifest.
run_ukb_content_scanner() {
    local input_file=$1
    local result_file=$2

    command -v python3 >/dev/null 2>&1 ||
        die 'python3 is required to inspect non-exempt destination files for UKB participant data'

    if ! python3 - "$DST_REAL" "$input_file" > "$result_file" <<'PY'
import csv
import io
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

root = Path(sys.argv[1])
input_path = Path(sys.argv[2])
SEVEN_DIGIT_ID = re.compile(r"^[0-9]{7}(?:[.]0+)?$")
CELL_REF = re.compile(r"^([A-Za-z]+)")
UKB_FIELD = re.compile(r"^(?:f_)?[0-9]{1,6}_[0-9]+_[0-9]+$")
STRONG_ID_HEADERS = {
    "eid", "e_id", "f_eid", "ukb_eid", "ukb_id", "participant_eid",
    "participant_id", "participantid", "eid_1", "eid_2",
}
GENERIC_ID_HEADERS = {
    "id", "id_1", "id_2", "iid", "fid", "sample_id", "sampleid",
    "subject_id", "subjectid", "person_id", "individual_id",
}
TEXT_EXTENSIONS = {
    ".csv", ".tsv", ".txt", ".dat", ".sample", ".fam", ".ped",
    ".map", ".json", ".jsonl", ".yaml", ".yml", ".html", ".htm",
    ".svg", ".md",
}
MAX_TEXT_BYTES = 4 * 1024 * 1024
MAX_ROWS = 250
MAX_COLUMNS = 256
MAX_XLSX_EXPANDED_BYTES = 250 * 1024 * 1024
MAX_XLSX_MEMBER_BYTES = 100 * 1024 * 1024


def clean_header(value):
    value = str(value or "").lstrip("\ufeff").strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def eid_like(value):
    return bool(SEVEN_DIGIT_ID.fullmatch(str(value or "").strip()))


def values_in_column(rows, column, start):
    values = set()
    for row in rows[start : start + 200]:
        if column < len(row) and eid_like(row[column]):
            values.add(str(row[column]).strip().split(".", 1)[0])
    return values


def detect_table(rows):
    rows = [list(row[:MAX_COLUMNS]) for row in rows if any(str(v).strip() for v in row)]
    if len(rows) < 2:
        return None

    for header_index in range(min(25, len(rows) - 1)):
        headers = [clean_header(value) for value in rows[header_index]]
        if not headers:
            continue

        for column, header in enumerate(headers):
            values = values_in_column(rows, column, header_index + 1)
            if header in STRONG_ID_HEADERS and values:
                return "eid_column_with_seven_digit_values"
            if header in GENERIC_ID_HEADERS and len(values) >= 3:
                return "id_column_with_multiple_seven_digit_values"

        has_ukb_fields = any(UKB_FIELD.fullmatch(header) for header in headers)
        if has_ukb_fields:
            for column in range(min(3, len(headers))):
                if len(values_in_column(rows, column, header_index + 1)) >= 2:
                    return "ukb_field_columns_with_eid_like_values"
    return None


def text_rows(path):
    with path.open("rb") as handle:
        raw = handle.read(MAX_TEXT_BYTES)
    if not raw or b"\x00" in raw[:65536]:
        return []

    text = raw.decode("utf-8-sig", errors="replace")
    lines = text.splitlines()[:MAX_ROWS]
    if not lines:
        return []

    sample = "\n".join(lines[:20])
    counts = {delimiter: sample.count(delimiter) for delimiter in ("\t", ",", ";", "|")}
    delimiter = max(counts, key=counts.get)
    if counts[delimiter] > 0:
        return list(csv.reader(lines, delimiter=delimiter))
    return [re.split(r"\s+", line.strip()) for line in lines if line.strip()]


def excel_column(reference):
    match = CELL_REF.match(reference or "")
    if not match:
        return 0
    value = 0
    for char in match.group(1).upper():
        value = value * 26 + ord(char) - 64
    return max(0, value - 1)


def shared_strings(workbook):
    name = "xl/sharedStrings.xml"
    if name not in workbook.namelist():
        return []
    strings = []
    with workbook.open(name) as handle:
        for _event, element in ET.iterparse(handle, events=("end",)):
            if element.tag.endswith("}si") or element.tag == "si":
                strings.append("".join(node.text or "" for node in element.iter()
                                       if node.tag.endswith("}t") or node.tag == "t"))
                element.clear()
    return strings


def cell_value(cell, strings):
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter()
                       if node.tag.endswith("}t") or node.tag == "t")
    value = next((node.text or "" for node in cell
                  if node.tag.endswith("}v") or node.tag == "v"), "")
    if cell_type == "s" and value.isdigit():
        index = int(value)
        return strings[index] if index < len(strings) else ""
    return value


def worksheet_rows(workbook, member, strings):
    rows = []
    with workbook.open(member) as handle:
        for _event, element in ET.iterparse(handle, events=("end",)):
            if not (element.tag.endswith("}row") or element.tag == "row"):
                continue
            values = {}
            for cell in element:
                if not (cell.tag.endswith("}c") or cell.tag == "c"):
                    continue
                column = excel_column(cell.attrib.get("r", ""))
                if column < MAX_COLUMNS:
                    values[column] = cell_value(cell, strings)
            if values:
                last_column = min(MAX_COLUMNS - 1, max(values))
                rows.append([values.get(column, "") for column in range(last_column + 1)])
            element.clear()
            if len(rows) >= MAX_ROWS:
                break
    return rows


def detect_xlsx(path):
    with zipfile.ZipFile(path) as workbook:
        infos = workbook.infolist()
        if (sum(info.file_size for info in infos) > MAX_XLSX_EXPANDED_BYTES or
                any(info.file_size > MAX_XLSX_MEMBER_BYTES for info in infos)):
            raise ValueError("expanded workbook is too large to inspect safely")
        strings = shared_strings(workbook)
        sheets = sorted(name for name in workbook.namelist()
                        if re.fullmatch(r"xl/worksheets/sheet[0-9]+[.]xml", name))
        for member in sheets[:50]:
            reason = detect_table(worksheet_rows(workbook, member, strings))
            if reason:
                return reason
    return None


def detect_file(path):
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        try:
            return detect_xlsx(path)
        except Exception:
            # A workbook that cannot be inspected must not be left ready to publish.
            return "xlsx_content_check_failed"
    if suffix in TEXT_EXTENSIONS:
        try:
            return detect_table(text_rows(path))
        except Exception:
            return "text_content_check_failed"
    return None


with input_path.open("r", encoding="utf-8") as entries:
    for raw_entry in entries:
        relative = raw_entry.rstrip("\n")
        if not relative:
            continue
        candidate = root.joinpath(*relative.split("/"))
        reason = detect_file(candidate)
        if reason:
            print(f"{relative}\t{reason}")
PY
    then
        die 'UKB destination content scanner failed'
    fi
}

rewrite_manifest_without_ukb_files() {
    local manifest=$1
    local line
    local -n rejected_ref=$2

    TEMP_FILE=$(mktemp "$DST_REAL/.github_files.lst.tmp.XXXXXX")
    while IFS= read -r line || [[ -n $line ]]; do
        if [[ -n $line && ${line:0:1} != '#' && ${rejected_ref[$line]+_} ]]; then
            continue
        fi
        printf '%s\n' "$line"
    done < "$manifest" > "$TEMP_FILE"
    chmod 0644 "$TEMP_FILE"
    mv -f -- "$TEMP_FILE" "$manifest"
    TEMP_FILE=
}


# 🚩 Synchronization
sync_files() {
    local project_spec=${1-}
    local manifest line source_path target parent candidate rel top reason project rejected_size
    local line_number=0 copied=0 unchanged=0 deleted=0 total_bytes=0
    local ukb_removed=0 ukb_checked_files=0 ukb_exempt_files=0 approved_count=0
    local -a files=()
    local -A wanted=()
    local -A file_sizes=()
    local -A copied_files=()
    local -A unchanged_files=()
    local -A ukb_rejected=()

    prepare_roots
    manifest=$DST_REAL/$MANIFEST_NAME
    [[ -f $manifest && ! -L $manifest ]] || die "manifest not found or unsafe; run scan first: $manifest"
    read_manifest_project_spec "$manifest"

    # Validate the complete recorded scope. With no explicit sync list, it is
    # also the effective scope. An explicit list may safely select a subset.
    configure_project_filter "$MANIFEST_PROJECT_SPEC"
    if [[ -n $project_spec ]]; then
        configure_project_filter "$project_spec"
        for project in "${SELECTED_PROJECT_LIST[@]}"; do
            manifest_scope_has_project "$project" ||
                die "sync dataset is outside the manifest scope ($MANIFEST_PROJECT_SPEC): $project"
        done
    fi

    # Complete preflight within the selected scope: no destination changes
    # happen until every selected manifest entry passes.
    while IFS= read -r line || [[ -n $line ]]; do
        ((line_number += 1))
        [[ -z $line || ${line:0:1} == '#' ]] && continue

        top=${line%%/*}
        manifest_scope_has_project "$top" ||
            die "manifest line $line_number is outside its recorded project scope: $line"
        project_is_selected "$top" || continue

        matches_publish_rule "$line" || die "manifest line $line_number is outside the publish whitelist: $line"
        source_path=$SRC_REAL/$line
        if ! check_file_safety "$line" "$source_path"; then
            die "manifest line $line_number failed safety check ($REJECT_REASON): $line"
        fi
        [[ ! ${wanted[$line]+_} ]] || die "duplicate manifest entry on line $line_number: $line"
        check_destination_target "$line"

        wanted["$line"]=1
        file_sizes["$line"]=$FILE_SIZE
        files+=("$line")
        total_bytes=$((total_bytes + FILE_SIZE))
    done < "$manifest"

    # Copy through a same-directory temporary file, then rename atomically.
    for rel in "${files[@]}"; do
        source_path=$SRC_REAL/$rel
        target=$DST_REAL/$rel
        parent=$(dirname -- "$target")
        mkdir -p -- "$parent"

        # Recheck the source immediately before use in case it changed after preflight.
        if ! check_file_safety "$rel" "$source_path"; then
            die "source changed or became unsafe during sync ($REJECT_REASON): $rel"
        fi

        if [[ -f $target && ! -L $target ]] && cmp -s -- "$source_path" "$target"; then
            ((unchanged += 1))
            unchanged_files["$rel"]=1
            continue
        fi

        TEMP_FILE=$(mktemp "$parent/.github_copy.tmp.XXXXXX")
        cp --preserve=mode,timestamps -- "$source_path" "$TEMP_FILE"
        mv -f -- "$TEMP_FILE" "$target"
        TEMP_FILE=
        ((copied += 1))
        copied_files["$rel"]=1
    done

    # UKB's public guidance describes a fast EID-pattern scan followed by a
    # contextual file review. Apply the same two signals to destination copies.
    # gu, grid, and ems120 are explicitly exempt because they do not use UKB.
    UKB_INPUT_FILE=$(mktemp "$DST_REAL/.github_copy.ukb_input.XXXXXX")
    for rel in "${files[@]}"; do
        top=${rel%%/*}
        if project_skips_ukb_content_check "$top"; then
            ((ukb_exempt_files += 1))
            continue
        fi
        printf '%s\n' "$rel" >> "$UKB_INPUT_FILE"
        ((ukb_checked_files += 1))
    done

    if [[ -s $UKB_INPUT_FILE ]]; then
        UKB_RESULT_FILE=$(mktemp "$DST_REAL/.github_copy.ukb_result.XXXXXX")
        run_ukb_content_scanner "$UKB_INPUT_FILE" "$UKB_RESULT_FILE"

        while IFS=$'\t' read -r rel reason || [[ -n $rel ]]; do
            [[ -n $rel && -n $reason ]] || die 'UKB scanner returned a malformed result'
            [[ ${wanted[$rel]+_} ]] || die "UKB scanner returned an unexpected path: $rel"

            target=$DST_REAL/$rel
            check_destination_target "$rel"
            [[ -f $target && ! -L $target ]] ||
                die "UKB scanner target disappeared or became unsafe: $rel"

            rm -f -- "$target"
            unset 'wanted[$rel]'
            ukb_rejected["$rel"]=$reason
            ((ukb_removed += 1))
            rejected_size=${file_sizes[$rel]}
            total_bytes=$((total_bytes - rejected_size))
            if [[ ${copied_files[$rel]+_} ]]; then
                ((copied -= 1))
            elif [[ ${unchanged_files[$rel]+_} ]]; then
                ((unchanged -= 1))
            fi
            printf 'UKB CHECK REMOVED: %s (%s)\n' "$rel" "$reason" >&2
        done < "$UKB_RESULT_FILE"
    fi

    if (( ukb_removed > 0 )); then
        rewrite_manifest_without_ukb_files "$manifest" ukb_rejected
        printf 'UKB CHECK: removed %d destination file(s) and the matching manifest entries.\n' \
            "$ukb_removed" >&2
    elif (( ukb_checked_files > 0 )); then
        printf 'UKB CHECK: no participant-level UKB pattern detected in %d checked file(s).\n' \
            "$ukb_checked_files"
    else
        printf 'UKB CHECK: no non-exempt files required content inspection.\n'
    fi
    if (( ukb_exempt_files > 0 )); then
        printf 'UKB CHECK: skipped %d file(s) in exempt datasets gu, grid, and ems120.\n' \
            "$ukb_exempt_files"
    fi

    rm -f -- "$UKB_INPUT_FILE"
    UKB_INPUT_FILE=
    if [[ -n $UKB_RESULT_FILE ]]; then
        rm -f -- "$UKB_RESULT_FILE"
        UKB_RESULT_FILE=
    fi

    # Keep root repository files and hidden top-level directories (.git, .github,
    # etc.). Below ordinary project directories, the manifest is authoritative.
    while IFS= read -r -d '' candidate; do
        rel=${candidate#"$DST_REAL"/}
        [[ $rel == */* ]] || continue
        top=${rel%%/*}
        [[ $top != .* ]] || continue
        project_is_selected "$top" || continue
        [[ ${wanted[$rel]+_} ]] && continue

        rm -f -- "$candidate"
        ((deleted += 1))
    done < <(
        find "$DST_REAL" \
            \( -path "$DST_REAL/.git" -o -path "$DST_REAL/.github" \) -prune -o \
            \( -type f -o -type l \) -print0
    )

    # Remove empty directories only within non-hidden project trees.
    for candidate in "$DST_REAL"/*; do
        [[ -d $candidate && ! -L $candidate ]] || continue
        top=${candidate##*/}
        project_is_selected "$top" || continue
        find "$candidate" -depth -type d -empty -delete
    done

    printf 'SYNC OK: %d copied/updated, %d unchanged, %d stale files removed.\n' \
        "$copied" "$unchanged" "$deleted"
    approved_count=$((${#files[@]} - ukb_removed))
    printf 'Approved payload: %d files, %d MiB.\n' \
        "$approved_count" "$((total_bytes / 1024 / 1024))"
}


# 🚩 CLI dispatch
main() {
    [[ $# -ge 1 && $# -le 2 ]] || {
        usage >&2
        exit 2
    }

    case $1 in
        scan)
            [[ $# -eq 2 && -n $2 ]] || die 'scan requires an explicit dataset list, for example: scan le8,maha'
            scan_files "$2"
            ;;
        sync) sync_files "${2-}" ;;
        -h|--help|help)
            [[ $# -eq 1 ]] || die 'help does not accept a project list'
            usage
            ;;
        *)
            usage >&2
            exit 2
            ;;
    esac
}

main "$@"
