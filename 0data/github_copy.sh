#!/usr/bin/env bash


# 🚩 Publish selected analysis results
#
# Usage:
#   ./github_copy.sh scan project,...    # inspect source files and rebuild the allowlist
#   ./github_copy.sh sync [project,...]  # inspect staged copies, then rebuild listed projects
# Delete each destination project represented by selected allowlist entries,
# then copy its listed files. Projects absent from the entries stay untouched.

set -Eeuo pipefail
IFS=$'\n\t'

# These can be overridden for testing, for example:
#   GITHUB_COPY_SOURCE_ROOT=/tmp/source GITHUB_COPY_DEST_ROOT=/tmp/dest ./github_copy.sh scan le8
SOURCE_ROOT=${GITHUB_COPY_SOURCE_ROOT:-/mnt/d/analysis}
DEST_ROOT=${GITHUB_COPY_DEST_ROOT:-/mnt/d/github/analysis}
MANIFEST_NAME=github_files.lst

# GitHub rejects files above 100 MiB. Use a deliberately conservative limit.
MAX_FILE_BYTES=$((50 * 1024 * 1024))
# Panome: admit non-participant results up to 100 MB (decimal), below GitHub's 100 MiB ceiling.
PANOME_MAX_FILE_BYTES=100000000
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
STAGING_ROOT=
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
    if [[ -n ${STAGING_ROOT:-} && -d $STAGING_ROOT ]]; then
        rm -rf -- "$STAGING_ROOT"
    fi
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

usage(){ cat <<'HELP'
cd /mnt/d/scripts/0data

./github_copy.sh scan le8,maha
./github_copy.sh sync le8,maha

# Panome: recursively inspect results; excludes participant IDs and files >100 MB.
./github_copy.sh scan panome
./github_copy.sh sync

# sync deletes and rebuilds destination projects with selected entries in github_files.lst.
# Old files within those projects are removed, even if not listed.
# Projects with no selected entries (for example grid/) are left untouched.
# sync rechecks staged file contents before changing any destination project.
HELP
}

require_commands() {
    local command_name
    for command_name in find sort stat realpath mktemp mv cp mkdir dirname rm chmod tr; do
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

    if [[ ${parts[0],,} == panome ]]; then
        case "$ext" in
            csv|tsv|xlsx|png|txt|md|json|jsonl|yaml|yml|html|htm|svg) return 0 ;;
            *) return 1 ;;
        esac
    fi

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
    local lower_base ext part current max_bytes size_reason
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

    if [[ ${parts[0],,} == panome ]]; then
        # Binary matrices/models cannot be inspected by this text/workbook scanner.
        case "$ext" in
            npy|npz|joblib|pkl|pickle|pt|pth|h5|hdf5|onnx|bin)
                REJECT_REASON=participant_matrix_or_model; return 1 ;;
        esac
        # Inspect readable outputs by content. Names, latent coordinates and
        # measurements alone are not participant identifiers.
    fi

    FILE_SIZE=$(stat -c '%s' -- "$source_path") || {
        REJECT_REASON=cannot_read_size
        return 1
    }
    max_bytes=$MAX_FILE_BYTES
    size_reason=over_50_MiB
    if [[ ${parts[0],,} == panome ]]; then
        max_bytes=$PANOME_MAX_FILE_BYTES
        size_reason=over_100_MB
    fi
    if (( FILE_SIZE > max_bytes )); then
        REJECT_REASON=$size_reason
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
    local ukb_checked_files=0 ukb_exempt_files=0
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

        if ! check_file_safety "$rel" "$source_path"; then
            reason=$REJECT_REASON
            (( skipped[$reason] += 1 ))
            continue
        fi
        if ! matches_publish_rule "$rel"; then
            (( skipped[not_publish_result] += 1 ))
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
                panome)
                    find_scan_candidates "$SRC_REAL/$project" 0
                    ;;
                le8)
                    # Includes le8/cvd_cad/prot/c2_cause/<result file>.
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

    # Inspect source content before writing the allowlist. Never copy or delete
    # analysis files during scan; a failed scanner leaves the old manifest intact.
    UKB_INPUT_FILE=$(mktemp "$DST_REAL/.github_copy.ukb_input.XXXXXX")
    UKB_RESULT_FILE=$(mktemp "$DST_REAL/.github_copy.ukb_result.XXXXXX")
    for rel in "${approved[@]}"; do
        if project_skips_ukb_content_check "${rel%%/*}"; then
            ((ukb_exempt_files += 1))
            continue
        fi
        printf '%s\n' "$rel" >> "$UKB_INPUT_FILE"
        ((ukb_checked_files += 1))
    done
    if [[ -s $UKB_INPUT_FILE ]]; then
        run_ukb_content_scanner "$UKB_INPUT_FILE" "$UKB_RESULT_FILE"
    fi
    while IFS=$'\t' read -r rel reason; do
        [[ -n $rel && -n $reason && ${approved_set[$rel]+_} ]] ||
            die 'UKB scanner returned an invalid scan result'
        unset 'approved_set[$rel]'
        (( skipped[$reason] += 1 ))
        printf 'UKB CHECK EXCLUDED: %s (%s)\n' "$rel" "$reason" >&2
    done < "$UKB_RESULT_FILE"
    approved=()
    for rel in "${!approved_set[@]}"; do
        approved+=("$rel")
    done
    printf 'UKB CHECK: inspected %d file(s); exempt gu/grid/ems120 files: %d.\n' \
        "$ukb_checked_files" "$ukb_exempt_files"

    TEMP_FILE=$(mktemp "$DST_REAL/.github_files.lst.tmp.XXXXXX")
    {
        printf '# Generated by github_copy.sh scan. Paths are relative to %s.\n' "$SOURCE_ROOT"
        printf '# Project scope: %s\n' "$PROJECT_SPEC"
        printf '# UKB content checks completed for non-exempt files during scan.\n'
        printf '# sync rechecks staged content before publishing non-exempt files.\n'
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

check_destination_project() {
    local project=$1
    local target=$DST_REAL/$project

    [[ -n $project && $project != .* && $project != */* && $project != *\\* &&
       ${project,,} != ukb ]] || die "unsafe destination project: $project"
    [[ ! -L $target ]] || die "destination project is a symlink: $target"
    [[ ! -e $target || -d $target ]] || die "destination project is not a directory: $target"
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
# Scan source files before admitting them to the allowlist. The detector deliberately
# emits only a relative path and a reason; participant-like values are never
# printed to the terminal or written to the manifest.
run_ukb_content_scanner() {
    local input_file=$1
    local result_file=$2
    local content_root=${3:-$SRC_REAL}

    command -v python3 >/dev/null 2>&1 ||
        die 'python3 is required to inspect non-exempt files for UKB participant data'

    if ! python3 - "$content_root" "$input_file" > "$result_file" <<'PY'
import csv
import json
import re
import sys
import zipfile
from html.parser import HTMLParser
from pathlib import Path
from xml.etree import ElementTree as ET

root = Path(sys.argv[1])
input_path = Path(sys.argv[2])
CELL_REF = re.compile(r"^([A-Za-z]+)")
PERSON_TOKENS = {
    "ukb", "participant", "sample", "subject", "patient", "person",
    "individual", "case", "control", "reference", "donor", "recipient",
}
EMPTY_ID_VALUES = {"", "na", "n/a", "nan", "none", "null", "."}
# Also inspect identifiers printed inside labels, YAML and other text fields.
INLINE_ID = re.compile(
    r'''["']?([A-Za-z][A-Za-z0-9_. -]*?)["']?\s*[:=]\s*["']?([0-9]{7})(?![0-9])'''
)
TEXT_EXTENSIONS = {
    ".csv", ".tsv", ".txt", ".dat", ".sample", ".fam", ".ped",
    ".map", ".json", ".jsonl", ".yaml", ".yml", ".html", ".htm",
    ".svg", ".md",
}
csv.field_size_limit(100000000)
MAX_XLSX_EXPANDED_BYTES = 250 * 1024 * 1024
MAX_XLSX_MEMBER_BYTES = 100 * 1024 * 1024


def clean_header(value):
    value = str(value if value is not None else "").lstrip("\ufeff").strip()
    value = re.sub(r"(?<=[a-z0-9])(?=[A-Z])", "_", value).lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def identifier_header(value):
    header = clean_header(value)
    tokens = header.split("_")
    if set(tokens) & {"eid", "iid", "fid"} or header == "e_id":
        return True
    if header in {"id", "sampleid", "subjectid", "participantid", "personid"}:
        return True
    if "id" in tokens and (set(tokens) & PERSON_TOKENS):
        return True
    # ID1/ID2/ID3 are real gene/protein names, not identifier column aliases.
    return bool(re.fullmatch(r"(?:id_[ab0-9]+|[ab0-9]+_id)", header))


def populated_id(value):
    if value is None or isinstance(value, bool):
        return False
    text = str(value).strip()
    return (text.lower() not in EMPTY_ID_VALUES and not identifier_header(text)
            and bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:@+-]*", text)))


def inline_identifier(value):
    return any(identifier_header(match.group(1)) for match in INLINE_ID.finditer(str(value)))


def detect_table(rows, panome=False):
    # Walk every row and column, including headers after introductory notes.
    # A single populated identifier is sufficient; never require three people.
    id_columns = set()
    for row in rows:
        if any(inline_identifier(value) for value in row):
            return "participant_identifier_value"
        if any(column < len(row) and populated_id(row[column]) for column in id_columns):
            return "participant_identifier_value"
        id_columns.update(column for column, value in enumerate(row) if identifier_header(value))
    return None


def text_rows(path):
    with path.open("rb") as handle:
        prefix = handle.read(4)
    encoding = "utf-16" if prefix.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8-sig"
    with path.open("r", encoding=encoding, errors="strict", newline="") as handle:
        sample = handle.read(65536)
        if "\x00" in sample:
            raise ValueError("unexpected binary text")
        handle.seek(0)
        # Respect quoted multiline CSV fields; inspect all rows, not a sample.
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters="\t,;|").delimiter
        except csv.Error:
            delimiter = {".csv": ",", ".tsv": "\t"}.get(path.suffix.lower())
        if delimiter:
            yield from csv.reader(handle, delimiter=delimiter, strict=True)
        else:
            for line in handle:
                if line.strip():
                    yield re.split(r"\s+", line.strip())


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
    if cell_type == "s":
        if not value.isdigit() or int(value) >= len(strings):
            raise ValueError("invalid shared-string reference")
        index = int(value)
        return strings[index]
    return value


def worksheet_rows(workbook, member, strings):
    with workbook.open(member) as handle:
        for _event, element in ET.iterparse(handle, events=("end",)):
            if not (element.tag.endswith("}row") or element.tag == "row"):
                continue
            values = {}
            for cell in element:
                if not (cell.tag.endswith("}c") or cell.tag == "c"):
                    continue
                column = excel_column(cell.attrib.get("r", ""))
                if column >= 16384:
                    raise ValueError("invalid Excel column")
                values[column] = cell_value(cell, strings)
            if values:
                yield [values.get(column, "") for column in range(max(values) + 1)]
            element.clear()


def detect_xlsx(path, panome=False):
    with zipfile.ZipFile(path) as workbook:
        infos = workbook.infolist()
        if (sum(info.file_size for info in infos) > MAX_XLSX_EXPANDED_BYTES or
                any(info.file_size > MAX_XLSX_MEMBER_BYTES for info in infos)):
            raise ValueError("expanded workbook is too large to inspect safely")
        strings = shared_strings(workbook)
        sheets = sorted(name for name in workbook.namelist()
                        if re.fullmatch(r"xl/worksheets/[^/]+[.]xml", name))
        if not sheets:
            raise ValueError("workbook has no inspectable sheets")
        for member in sheets:
            reason = detect_table(worksheet_rows(workbook, member, strings), panome)
            if reason:
                return reason
    return None


def detect_json(value):
    if isinstance(value, dict):
        for key, item in value.items():
            # Recognize both records and column-oriented identifier arrays.
            # Configuration such as {"id_col": "eid"} is not participant data.
            if identifier_header(key):
                values = item if isinstance(item, list) else [item]
                if any(not isinstance(v, (dict, list)) and populated_id(v) for v in values):
                    return "participant_identifier_value"
            if re.fullmatch(r"[0-9]{7}", str(key)) and isinstance(item, (dict, list)):
                return "participant_identifier_key"
            reason = detect_json(item)
            if reason:
                return reason
    elif isinstance(value, list):
        # Also support JSON tables expressed as a header row followed by rows.
        if value and all(isinstance(row, list) for row in value):
            reason = detect_table(value)
            if reason:
                return reason
        for item in value:
            reason = detect_json(item)
            if reason:
                return reason
    elif isinstance(value, str) and inline_identifier(value):
        return "participant_identifier_value"
    return None


class HTMLRows(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.cell = [], [], None

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.row = []
        elif tag in {"th", "td"}:
            self.cell = []

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in {"th", "td"} and self.cell is not None:
            self.row.append("".join(self.cell))
            self.cell = None
        elif tag == "tr":
            self.rows.append(self.row)


def detect_html(path):
    text = path.read_text(encoding="utf-8-sig")
    if inline_identifier(text):
        return "participant_identifier_value"
    parser = HTMLRows()
    parser.feed(text)
    return detect_table(parser.rows)


def detect_file(path, panome=False):
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        try:
            return detect_xlsx(path, panome)
        except Exception:
            # A workbook that cannot be inspected must not be left ready to publish.
            return "xlsx_content_check_failed"
    if suffix in {".json", ".jsonl"}:
        try:
            with path.open(encoding="utf-8-sig") as handle:
                if suffix == ".json":
                    return detect_json(json.load(handle))
                for line in handle:
                    if line.strip():
                        reason = detect_json(json.loads(line))
                        if reason:
                            return reason
            return None
        except Exception:
            return "json_content_check_failed"
    if suffix in TEXT_EXTENSIONS:
        try:
            if suffix in {".html", ".htm", ".svg"}:
                return detect_html(path)
            return detect_table(text_rows(path), panome)
        except Exception:
            return "text_content_check_failed"
    return None


with input_path.open("r", encoding="utf-8") as entries:
    for raw_entry in entries:
        relative = raw_entry.rstrip("\n")
        if not relative:
            continue
        candidate = root.joinpath(*relative.split("/"))
        reason = detect_file(candidate, relative.split("/", 1)[0].lower() == "panome")
        if reason:
            print(f"{relative}\t{reason}")
PY
    then
        die 'UKB content scanner failed'
    fi
}

# 🚩 Synchronization
sync_files() {
    local project_spec=${1-}
    local manifest line source_path target parent rel top project
    local line_number=0 copied=0 total_bytes=0
    local -a files=() projects=()
    local -A wanted=() projects_seen=()

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
        check_destination_project "$top"

        if [[ ! ${projects_seen[$top]+_} ]]; then
            projects_seen["$top"]=1
            projects+=("$top")
        fi

        wanted["$line"]=1
        files+=("$line")
        total_bytes=$((total_bytes + FILE_SIZE))
    done < "$manifest"

    # Snapshot outside the publishing tree, then inspect exactly those bytes.
    # A hand-edited/old allowlist and changes since scan must not bypass checks.
    # Never put an uninspected copy (even a temporary one) inside DEST_ROOT.
    STAGING_ROOT=$(mktemp -d /tmp/github_copy.stage.XXXXXX)
    UKB_INPUT_FILE=$(mktemp /tmp/github_copy.ukb_input.XXXXXX)
    UKB_RESULT_FILE=$(mktemp /tmp/github_copy.ukb_result.XXXXXX)
    for rel in "${files[@]}"; do
        source_path=$SRC_REAL/$rel
        if ! check_file_safety "$rel" "$source_path"; then
            die "source changed or became unsafe during sync ($REJECT_REASON): $rel"
        fi
        target=$STAGING_ROOT/$rel
        mkdir -p -- "$(dirname -- "$target")"
        cp --no-dereference --preserve=mode,timestamps -- "$source_path" "$target"
        if [[ -L $target ]] || ! check_file_safety "$rel" "$target"; then
            die "staged source failed safety check: $rel"
        fi
        if ! project_skips_ukb_content_check "${rel%%/*}"; then
            printf '%s\n' "$rel" >> "$UKB_INPUT_FILE"
        fi
    done
    if [[ -s $UKB_INPUT_FILE ]]; then
        run_ukb_content_scanner "$UKB_INPUT_FILE" "$UKB_RESULT_FILE" "$STAGING_ROOT"
    fi
    if [[ -s $UKB_RESULT_FILE ]]; then
        while IFS=$'\t' read -r rel reason; do
            printf 'UKB SYNC BLOCKED: %s (%s)\n' "$rel" "$reason" >&2
        done < "$UKB_RESULT_FILE"
        die 'staged content is unsafe; destination projects were not changed; run scan again'
    fi

    # Only projects with selected file entries are cleared. A scope header alone
    # (including an empty scan result) must not cause an unrelated tree deletion.
    for project in "${projects[@]}"; do
        check_destination_project "$project"
        rm -rf -- "$DST_REAL/$project"
    done

    # Copy through a same-directory temporary file, then rename atomically.
    for rel in "${files[@]}"; do
        source_path=$STAGING_ROOT/$rel
        target=$DST_REAL/$rel
        parent=$(dirname -- "$target")
        mkdir -p -- "$parent"

        # Only the inspected snapshot is published, never a fresh source copy.
        TEMP_FILE=$(mktemp "$parent/.github_copy.tmp.XXXXXX")
        mv -f -- "$source_path" "$TEMP_FILE"
        mv -f -- "$TEMP_FILE" "$target"
        TEMP_FILE=
        ((copied += 1))
    done

    printf 'SYNC OK: %d projects rebuilt, %d files copied.\n' "${#projects[@]}" "$copied"
    printf 'Allowlisted payload: %d files, %d MiB.\n' \
        "${#files[@]}" "$((total_bytes / 1024 / 1024))"
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
