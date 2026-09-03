#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C

WES_DIR_DEFAULT=/mnt/d/data.BIG/gwas/wes
SOURCE_MT_DEFAULT=gs://ukbb-exome-public/500k/results/results.mt
SINGLE_SOURCE_MT_DEFAULT=gs://ukbb-exome-public/500k/results/variant_results.mt
HAIL_VENV_DEFAULT="$HOME/venvs/hail"
# Connector 4.x is compiled for Java 17, while this Hail setup requires Java 11.
# Keep this on the Java 11-compatible 3.1 line.
GCS_CONNECTOR_DEFAULT=https://repo1.maven.org/maven2/com/google/cloud/bigdataoss/gcs-connector/3.1.15/gcs-connector-3.1.15-shaded.jar


# 🚩 Genebass WES download
usage() {
cat <<'USAGE'
Usage: ./wes_download.sh [options]

Export the complete GeneBass gene-based results (SKAT-O, Burden and SKAT),
one bgzipped TSV per phenotype.

  --dir-wes PATH        WES root [/mnt/d/data.BIG/gwas/wes]
  --dir-out PATH        Downloads root [<dir-wes>/Downloads]
  --tmpdir PATH         Hail temp [<dir-out>/tmp/hail]
  --hail-venv PATH      Hail virtual environment [~/venvs/hail]
  --source-mt URI       Override the official public MatrixTable
  --single-variant BOOL  Export single-variant GWAS after gene results [FALSE]
  --single-source-mt URI  Override the official variant MatrixTable
  --single-batch-size N  Variant phenotypes per export batch [same as --batch-size]
  --single-estimate-gib N
                        Estimated GiB per single-variant phenotype [0.65]
  --allow-low-space BOOL
                        Override the pre-export disk-space guard [FALSE]
  --gcs-connector SPEC  Maven coordinate or connector JAR URI [official 3.1.15]
  --proxy URL           HTTP(S) proxy; auto-detects WSL localhost:7897
  --gcp-project ID      Billing/quota project required by Requester Pays
  --pheno-regex REGEX   Filter trait_type/phenocode/category/description
  --batch-size N        Columns per Hail export batch [16]
  --spark-memory SIZE   Local Spark Java heap [16g]
  --spark-cores N       Concurrent local Spark tasks [8]
  --task-max-failures N  Attempts per Spark task before aborting [8]
  --gcs-connect-timeout DURATION
                        GCS HTTP connection timeout [30s]
  --gcs-read-timeout DURATION
                        GCS HTTP read timeout [60s]
  --gcs-max-retries N   Low-level GCS HTTP retries [20]
  --cleanup-interval N  Remove merged temp batches every N seconds; 0 disables [21600]
  --profile full|compact  All official fields, or compact analysis fields [full]
  --overwrite TRUE|FALSE  Discard resume state and restart from zero [FALSE]
  --check               Check dependencies only
  -h, --help            Show help

Default layout:
  wes/Downloads/genebass_500k/      project-level download staging
    files/*.tsv.bgz                 atomically committed phenotype results
    files/index.tsv                 merged completed-result index
    .resume_batch/                  current batch only; safe to recompute
    resume.state.json               progress and resume-compatibility state
  wes/Downloads/single_500k/
    files/*.tsv.bgz                 one single-variant GWAS file per phenotype
    files/index.tsv                 completed single-variant result index
    resume.state.json               independent variant export progress
  wes/Downloads/tmp/hail/           Hail temporary data
  wes/common/<trait>/raw/           split official single-variant data
  wes/common/<trait>/gwas/          cleaned common-variant output
  wes/rare/<trait>/raw/             split official rare-variant data
  wes/rare/<trait>/gwas/            cleaned rare-variant gene results
USAGE
}

die() { echo "ERROR: $*" >&2; exit 2; }
upper() { printf '%s' "$1" | tr '[:lower:]' '[:upper:]'; }
need_value() { [[ -n "${2-}" && "${2-}" != --* ]] || die "missing value for $1"; }

dir_wes="$WES_DIR_DEFAULT"; dir_out=""; tmpdir=""
source_mt="$SOURCE_MT_DEFAULT"; pheno_regex=""; batch_size=16
hail_venv="$HAIL_VENV_DEFAULT"
single_variant=FALSE
single_source_mt="$SINGLE_SOURCE_MT_DEFAULT"
single_batch_size=""
allow_low_space=FALSE
single_estimate_gib=0.65
gcs_connector="$GCS_CONNECTOR_DEFAULT"
profile=full; overwrite=FALSE; check_only=FALSE
spark_memory=16g
spark_cores=8
task_max_failures=8
gcs_connect_timeout=30s
gcs_read_timeout=60s
gcs_max_retries=20
cleanup_interval=21600
proxy_url="${HTTPS_PROXY:-${https_proxy:-}}"
gcp_project=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir-wes) need_value "$1" "${2-}"; dir_wes="$2"; shift 2;;
    --dir-out) need_value "$1" "${2-}"; dir_out="$2"; shift 2;;
    --tmpdir) need_value "$1" "${2-}"; tmpdir="$2"; shift 2;;
    --hail-venv) need_value "$1" "${2-}"; hail_venv="$2"; shift 2;;
    --source-mt) need_value "$1" "${2-}"; source_mt="$2"; shift 2;;
    --single-variant) need_value "$1" "${2-}"; single_variant="$(upper "$2")"; shift 2;;
    --single-source-mt) need_value "$1" "${2-}"; single_source_mt="$2"; shift 2;;
    --single-batch-size) need_value "$1" "${2-}"; single_batch_size="$2"; shift 2;;
    --single-estimate-gib) need_value "$1" "${2-}"; single_estimate_gib="$2"; shift 2;;
    --allow-low-space) need_value "$1" "${2-}"; allow_low_space="$(upper "$2")"; shift 2;;
    --gcs-connector) need_value "$1" "${2-}"; gcs_connector="$2"; shift 2;;
    --proxy) need_value "$1" "${2-}"; proxy_url="$2"; shift 2;;
    --gcp-project) need_value "$1" "${2-}"; gcp_project="$2"; shift 2;;
    --pheno-regex) need_value "$1" "${2-}"; pheno_regex="$2"; shift 2;;
    --batch-size) need_value "$1" "${2-}"; batch_size="$2"; shift 2;;
    --spark-memory) need_value "$1" "${2-}"; spark_memory="$2"; shift 2;;
    --spark-cores) need_value "$1" "${2-}"; spark_cores="$2"; shift 2;;
    --task-max-failures) need_value "$1" "${2-}"; task_max_failures="$2"; shift 2;;
    --gcs-connect-timeout) need_value "$1" "${2-}"; gcs_connect_timeout="$2"; shift 2;;
    --gcs-read-timeout) need_value "$1" "${2-}"; gcs_read_timeout="$2"; shift 2;;
    --gcs-max-retries) need_value "$1" "${2-}"; gcs_max_retries="$2"; shift 2;;
    --cleanup-interval) need_value "$1" "${2-}"; cleanup_interval="$2"; shift 2;;
    --profile) need_value "$1" "${2-}"; profile="$2"; shift 2;;
    --overwrite) need_value "$1" "${2-}"; overwrite="$(upper "$2")"; shift 2;;
    --check) check_only=TRUE; shift;;
    -h|--help) usage; exit 0;;
    *) die "unknown option: $1";;
  esac
done

[[ "$batch_size" =~ ^[1-9][0-9]*$ ]] || die "--batch-size must be a positive integer"
[[ -n "$single_batch_size" ]] || single_batch_size="$batch_size"
[[ "$single_batch_size" =~ ^[1-9][0-9]*$ ]] || die "--single-batch-size must be a positive integer"
[[ "$spark_memory" =~ ^[1-9][0-9]*[mMgG]$ ]] ||
  die "--spark-memory must look like 4096m or 16g"
[[ "$spark_cores" =~ ^[1-9][0-9]*$ ]] || die "--spark-cores must be a positive integer"
[[ "$task_max_failures" =~ ^[1-9][0-9]*$ ]] ||
  die "--task-max-failures must be a positive integer"
[[ "$gcs_connect_timeout" =~ ^(0|[1-9][0-9]*(ms|s|m|h))$ ]] ||
  die "--gcs-connect-timeout must look like 30s, 2m, or 0"
[[ "$gcs_read_timeout" =~ ^(0|[1-9][0-9]*(ms|s|m|h))$ ]] ||
  die "--gcs-read-timeout must look like 60s, 2m, or 0"
[[ "$gcs_max_retries" =~ ^[0-9]+$ ]] ||
  die "--gcs-max-retries must be a non-negative integer"
[[ "$cleanup_interval" =~ ^[0-9]+$ ]] || die "--cleanup-interval must be a non-negative integer"
[[ "$profile" == full || "$profile" == compact ]] || die "--profile must be full or compact"
[[ "$overwrite" == TRUE || "$overwrite" == FALSE ]] || die "--overwrite must be TRUE or FALSE"
single_variant="$(upper "$single_variant")"
allow_low_space="$(upper "$allow_low_space")"
[[ "$single_variant" == TRUE || "$single_variant" == FALSE ]] ||
  die "--single-variant must be TRUE or FALSE"
[[ "$allow_low_space" == TRUE || "$allow_low_space" == FALSE ]] ||
  die "--allow-low-space must be TRUE or FALSE"
[[ "$single_estimate_gib" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
  die "--single-estimate-gib must be a non-negative number"
[[ -n "$dir_out" ]] || dir_out="${dir_wes%/}/Downloads"
downloads_dir="${dir_out%/}"
dir_out="${downloads_dir%/}/genebass_500k"
single_dir_out="${downloads_dir%/}/single_500k"
[[ -n "$tmpdir" ]] || tmpdir="${downloads_dir%/}/tmp/hail"
if [[ -z "$gcp_project" ]] && command -v gcloud >/dev/null 2>&1; then
  gcp_project="$(gcloud config get-value project 2>/dev/null || true)"
  [[ "$gcp_project" == '(unset)' ]] && gcp_project=""
fi
if [[ "$source_mt" == gs://* ||
      ( "$single_variant" == TRUE && "$single_source_mt" == gs://* ) ]]; then
  [[ -n "$gcp_project" ]] || die "GeneBass is Requester Pays; supply --gcp-project ID"
  [[ -s "$HOME/.config/gcloud/application_default_credentials.json" ]] ||
    die "ADC missing; run: gcloud auth application-default login"
fi

[[ -r "$hail_venv/bin/activate" ]] || {
  echo "ERROR: Hail virtual environment not found: $hail_venv" >&2
  echo "Create it as documented in README.wes.md or use --hail-venv PATH." >&2
  exit 2
}
# shellcheck disable=SC1090
source "$hail_venv/bin/activate"

# On this workstation the Windows proxy is reachable through mirrored WSL
# networking. Use it only when no proxy was explicitly supplied/inherited.
if [[ -z "$proxy_url" ]] && command -v timeout >/dev/null 2>&1 &&
   timeout 2 bash -c '</dev/tcp/127.0.0.1/7897' >/dev/null 2>&1; then
  proxy_url=http://127.0.0.1:7897
fi

command -v python3 >/dev/null 2>&1 || die "python3 is required"
command -v java >/dev/null 2>&1 || die "Java 11 is required"
java_version="$(java -version 2>&1 | awk -F '[\".]' 'NR==1 {print $2; exit}')"
[[ "$java_version" == 11 ]] || die "Java 11 is required; detected: $(java -version 2>&1 | head -n 1)"
python3 -c 'import hail as hl; print(hl.__version__)' >/dev/null 2>&1 || {
  echo "ERROR: Python package 'hail' is unavailable. See README.wes.md." >&2
  exit 2
}

if [[ "$check_only" == TRUE ]]; then
  echo "Software check OK"
  echo "Python: $(python3 --version 2>&1)"
  echo "Java:   $(java -version 2>&1 | head -n 1)"
  echo "Hail:   $(python3 -c 'import hail as hl; print(hl.__version__)')"
  echo "Spark:  $spark_memory heap, $spark_cores cores"
  echo "Retry:  $task_max_failures task attempts, $gcs_max_retries GCS HTTP retries"
  echo "Timeout: GCS connect $gcs_connect_timeout, read $gcs_read_timeout"
  echo "Batch:  $batch_size phenotypes"
  echo "Single: $single_variant; batch $single_batch_size; source $single_source_mt"
  echo "Downloads root: $downloads_dir"
  echo "Gene output: $dir_out"
  echo "Variant output: $single_dir_out"
  echo "Temporary data: $tmpdir"
  echo "Cleanup: every ${cleanup_interval}s (0 disables)"
  exit 0
fi

files_dir="${dir_out%/}/files"
single_files_dir="${single_dir_out%/}/files"
if [[ "$overwrite" == TRUE ]]; then
  files_real="$(realpath -m -- "$files_dir")"
  out_real="$(realpath -m -- "$dir_out")"
  [[ "$out_real" != / ]] || die "refusing to overwrite output rooted at /"
  [[ "$(dirname -- "$files_real")" == "$out_real" && "$(basename -- "$files_real")" == files ]] ||
    die "refusing to remove unexpected path: $files_real"
  rm -rf -- "$files_real"
  rm -rf -- "${dir_out%/}/.resume_batch"
  rm -f -- "${dir_out%/}/resume.state.json" "${dir_out%/}/run.info.json"
  if [[ "$single_variant" == TRUE ]]; then
    single_files_real="$(realpath -m -- "$single_files_dir")"
    single_out_real="$(realpath -m -- "$single_dir_out")"
    [[ "$single_out_real" != / ]] || die "refusing to overwrite variant output rooted at /"
    [[ "$(dirname -- "$single_files_real")" == "$single_out_real" &&
       "$(basename -- "$single_files_real")" == files ]] ||
      die "refusing to remove unexpected variant path: $single_files_real"
    rm -rf -- "$single_files_real"
    rm -rf -- "${single_dir_out%/}/.resume_batch"
    rm -f -- "${single_dir_out%/}/resume.state.json" \
      "${single_dir_out%/}/run.info.json"
  fi
elif [[ -e "$files_dir" ]]; then
  [[ -d "$files_dir" && ! -L "$files_dir" ]] ||
    die "output is not a real directory: $files_dir"
  echo "[$(date '+%F %T')] resume : preserving existing output in $files_dir"
fi
if [[ "$single_variant" == TRUE && -e "$single_files_dir" ]]; then
  [[ -d "$single_files_dir" && ! -L "$single_files_dir" ]] ||
    die "variant output is not a real directory: $single_files_dir"
  echo "[$(date '+%F %T')] resume : preserving variant output in $single_files_dir"
fi

mkdir -p -- "$dir_out" "$tmpdir" "$downloads_dir/logs"
if [[ "$single_variant" == TRUE ]]; then
  mkdir -p -- "$single_dir_out"
fi
dir_out="$(realpath -m -- "$dir_out")"
single_dir_out="$(realpath -m -- "$single_dir_out")"
tmpdir="$(realpath -m -- "$tmpdir")"
downloads_dir="$(realpath -m -- "$downloads_dir")"

export WES_SOURCE_MT="$source_mt" WES_DIR_OUT="$dir_out" WES_TMPDIR="$tmpdir"
export WES_PHENO_REGEX="$pheno_regex" WES_BATCH_SIZE="$batch_size" WES_PROFILE="$profile"
export WES_SPARK_MEMORY="$spark_memory" WES_SPARK_CORES="$spark_cores"
export WES_TASK_MAX_FAILURES="$task_max_failures"
export WES_GCS_CONNECT_TIMEOUT="$gcs_connect_timeout"
export WES_GCS_READ_TIMEOUT="$gcs_read_timeout"
export WES_GCS_MAX_RETRIES="$gcs_max_retries"
export WES_GCS_CONNECTOR="$gcs_connector"
export WES_PROXY_URL="$proxy_url"
export WES_GCP_PROJECT="$gcp_project"
export WES_DOWNLOADS_DIR="$downloads_dir"
export WES_SINGLE_VARIANT="$single_variant"
export WES_SINGLE_SOURCE_MT="$single_source_mt"
export WES_SINGLE_DIR_OUT="$single_dir_out"
export WES_SINGLE_BATCH_SIZE="$single_batch_size"
export WES_ALLOW_LOW_SPACE="$allow_low_space"
export WES_SINGLE_ESTIMATE_GIB_PER_PHENO="$single_estimate_gib"

if [[ -n "$proxy_url" ]]; then
  proxy_hostport="${proxy_url#*://}"
  proxy_hostport="${proxy_hostport%%/*}"
  proxy_host="${proxy_hostport%:*}"
  proxy_port="${proxy_hostport##*:}"
  [[ -n "$proxy_host" && "$proxy_port" =~ ^[0-9]+$ ]] ||
    die "--proxy must look like http://HOST:PORT"
  export HTTP_PROXY="$proxy_url" HTTPS_PROXY="$proxy_url"
  export http_proxy="$proxy_url" https_proxy="$proxy_url"
  # Hail Python talks to its local Java backend over HTTP. Never send that RPC
  # through the external proxy, or Python receives an empty/non-Hail response.
  no_proxy_local="localhost,127.0.0.1,::1,10.255.255.254"
  export NO_PROXY="${NO_PROXY:+$NO_PROXY,}$no_proxy_local"
  export no_proxy="${no_proxy:+$no_proxy,}$no_proxy_local"
  java_proxy_opts="-Dhttp.proxyHost=$proxy_host -Dhttp.proxyPort=$proxy_port -Dhttps.proxyHost=$proxy_host -Dhttps.proxyPort=$proxy_port -Dhttp.nonProxyHosts=localhost|127.*|[::1]|10.255.255.254"
  export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:+$JAVA_TOOL_OPTIONS }$java_proxy_opts"
fi

echo "[$(date '+%F %T')] source : $source_mt"
echo "[$(date '+%F %T')] output : $files_dir"
echo "[$(date '+%F %T')] tmpdir : $tmpdir"
echo "[$(date '+%F %T')] profile: $profile"
echo "[$(date '+%F %T')] spark  : ${spark_memory} heap, ${spark_cores} cores"
echo "[$(date '+%F %T')] retry  : ${task_max_failures} task attempts; GCS ${gcs_max_retries} retries"
echo "[$(date '+%F %T')] timeout: GCS connect ${gcs_connect_timeout}, read ${gcs_read_timeout}"
echo "[$(date '+%F %T')] cleanup: every ${cleanup_interval}s (0 disables)"
echo "[$(date '+%F %T')] proxy  : ${proxy_url:-none}"
echo "[$(date '+%F %T')] project: ${gcp_project:-none}"
echo "[$(date '+%F %T')] single : $single_variant"
if [[ "$single_variant" == TRUE ]]; then
  echo "[$(date '+%F %T')] var src: $single_source_mt"
  echo "[$(date '+%F %T')] var out: $single_dir_out/files"
fi

# Remove old Hail fragments only after their merged root part exists.


# 🚩 Export cleanup
cleanup_completed_export_temp() {
  local requested_files_dir=$1
  local requested_tmp_dir="${requested_files_dir%/}/tmp"
  local files_real tmp_real path name stage max_stage=-1 merge_stage i
  local deleted_count=0
  local -a paths=() stages=() delete_paths=() batch=()
  local -A merged_stages=()

  [[ -d "$requested_tmp_dir" && ! -L "$requested_tmp_dir" ]] || return 0
  [[ -d "$requested_files_dir" && ! -L "$requested_files_dir" ]] || return 0
  files_real="$(realpath -e -- "$requested_files_dir")"
  tmp_real="$(realpath -e -- "$requested_tmp_dir")"
  [[ "$files_real" != / && "$tmp_real" != / ]] || return 2
  [[ "$(dirname -- "$tmp_real")" == "$files_real" && "$(basename -- "$tmp_real")" == tmp ]] ||
    return 2

  while IFS= read -r -d '' path; do
    name="${path##*/}"
    if [[ "$name" =~ ^part-[0-9]+-([0-9]+)- ]]; then
      merged_stages["${BASH_REMATCH[1]}"]=1
    fi
  done < <(find "$files_real" -mindepth 1 -maxdepth 1 -type f -name 'part-*' -print0)

  while IFS= read -r -d '' path; do
    name="${path##*/}"
    if [[ "$name" =~ ^part-[0-9]+-([0-9]+)- ]]; then
      stage="${BASH_REMATCH[1]}"
      paths+=("$path")
      stages+=("$stage")
      (( stage > max_stage )) && max_stage=$stage
    fi
  done < <(find "$tmp_real" -mindepth 1 -maxdepth 1 -type d -print0)

  (( max_stage >= 2 )) || return 0
  for ((i = 0; i < ${#paths[@]}; i++)); do
    stage="${stages[$i]}"
    merge_stage=$((stage + 1))
    if (( stage <= max_stage - 2 )) && [[ -n "${merged_stages[$merge_stage]+x}" ]]; then
      delete_paths+=("${paths[$i]}")
    fi
  done
  (( ${#delete_paths[@]} > 0 )) || return 0

  for path in "${delete_paths[@]}"; do
    [[ "$(dirname -- "$path")" == "$tmp_real" ]] || return 2
    [[ "${path##*/}" =~ ^part-[0-9]+-[0-9]+- ]] || return 2
    batch+=("$path")
    if (( ${#batch[@]} >= 128 )); then
      rm -rf -- "${batch[@]}"
      ((deleted_count += ${#batch[@]}))
      batch=()
    fi
  done
  if (( ${#batch[@]} > 0 )); then
    rm -rf -- "${batch[@]}"
    ((deleted_count += ${#batch[@]}))
  fi
  echo "[$(date '+%F %T')] cleanup: deleted $deleted_count completed-batch directories; retained stage $max_stage."
}

cleanup_watcher_pid=""
stop_cleanup_watcher() {
  if [[ -n "$cleanup_watcher_pid" ]]; then
    kill "$cleanup_watcher_pid" 2>/dev/null || true
    wait "$cleanup_watcher_pid" 2>/dev/null || true
    cleanup_watcher_pid=""
  fi
}
start_cleanup_watcher() {
  (( cleanup_interval > 0 )) || return 0
  (
    while true; do
      sleep "$cleanup_interval"
      cleanup_completed_export_temp "$files_dir" ||
        echo "[$(date '+%F %T')] WARNING: automatic temporary-file cleanup failed" >&2
    done
  ) &
  cleanup_watcher_pid=$!
}

start_cleanup_watcher
trap stop_cleanup_watcher EXIT

python3 - <<'PY'
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import threading
from datetime import datetime
from pathlib import Path

import hail as hl

source = os.environ['WES_SOURCE_MT']
out_dir = os.environ['WES_DIR_OUT']
tmp_dir = os.environ['WES_TMPDIR']
pheno_regex = os.environ['WES_PHENO_REGEX']
batch_size = int(os.environ['WES_BATCH_SIZE'])
spark_memory = os.environ['WES_SPARK_MEMORY']
spark_cores = int(os.environ['WES_SPARK_CORES'])
task_max_failures = int(os.environ['WES_TASK_MAX_FAILURES'])
gcs_connect_timeout = os.environ['WES_GCS_CONNECT_TIMEOUT']
gcs_read_timeout = os.environ['WES_GCS_READ_TIMEOUT']
gcs_max_retries = int(os.environ['WES_GCS_MAX_RETRIES'])
profile = os.environ['WES_PROFILE']
gcs_connector = os.environ['WES_GCS_CONNECTOR']
proxy_url = os.environ['WES_PROXY_URL']
gcp_project = os.environ['WES_GCP_PROJECT']
downloads_dir = Path(os.environ['WES_DOWNLOADS_DIR'])

out_path = Path(out_dir)
files_path = out_path / 'files'
work_path = out_path / '.resume_batch'
state_path = out_path / 'resume.state.json'
run_info_path = out_path / 'run.info.json'
data_suffix = '.tsv.bgz'
safe_file_re = re.compile(r'^[A-Za-z0-9._-]+$')
cleanup_paths = []


def atomic_write_text(path, text):
    """Write state and indexes atomically on the same filesystem."""
    path = Path(path)
    tmp = path.with_name(f'.{path.name}.tmp-{os.getpid()}')
    with tmp.open('w', encoding='utf-8', newline='\n') as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def atomic_write_json(path, value):
    atomic_write_text(
        path, json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def read_header_json(path):
    """Read and validate the phenotype JSON stored in a BGZF first line."""
    path = Path(path)
    try:
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            line = handle.readline()
    except (OSError, EOFError) as exc:
        raise RuntimeError(f'invalid BGZF result {path}: {exc}') from exc
    if not line.startswith('#'):
        raise RuntimeError(f'missing phenotype JSON header: {path}')
    try:
        metadata = json.loads(line[1:])
    except json.JSONDecodeError as exc:
        raise RuntimeError(f'invalid phenotype JSON in {path}: {exc}') from exc
    file_name = metadata.get('_file')
    if not isinstance(file_name, str) or not safe_file_re.fullmatch(file_name):
        raise RuntimeError(f'unsafe or missing _file in {path}: {file_name!r}')
    return metadata


def crc_path(path):
    path = Path(path)
    return path.with_name(f'.{path.name}.crc')


def remove_tree(path, expected_parent):
    """Remove only a fixed, non-symlink child of the expected directory."""
    path = Path(path)
    expected_parent = Path(expected_parent).resolve()
    if not path.exists():
        return
    if path.is_symlink() or path.parent.resolve() != expected_parent:
        raise RuntimeError(f'refusing to remove unexpected directory: {path}')
    shutil.rmtree(path)


def quarantine_directory(path, prefix):
    """Atomically detach a large incomplete tree for background deletion."""
    path = Path(path)
    if not path.exists():
        return None
    if path.is_symlink() or path.parent.resolve() not in {
            files_path.resolve(), out_path.resolve()}:
        raise RuntimeError(f'refusing to quarantine unexpected path: {path}')
    destination = out_path / (
        f'.{prefix}-{datetime.now():%Y%m%d-%H%M%S}-{os.getpid()}')
    if destination.exists():
        raise RuntimeError(f'cleanup quarantine already exists: {destination}')
    os.replace(path, destination)
    cleanup_paths.append(destination)
    print(f'Quarantined temporary tree: {destination.name}', flush=True)
    return destination


def cleanup_quarantined_tree(path):
    """Delete a quarantined tree incrementally without blocking resume."""
    path = Path(path)
    try:
        if path.is_symlink() or path.parent.resolve() != out_path.resolve():
            raise RuntimeError(f'unsafe cleanup quarantine: {path}')
        removed = 0
        for child in list(path.iterdir()):
            try:
                if child.is_dir() and not child.is_symlink():
                    shutil.rmtree(child)
                else:
                    child.unlink(missing_ok=True)
            except FileNotFoundError:
                pass
            removed += 1
            if removed % 2000 == 0:
                print(
                    f'Background cleanup {path.name}: '
                    f'{removed:,} top-level entries removed.', flush=True)
        path.rmdir()
        print(f'Background cleanup completed: {path.name}', flush=True)
    except Exception as exc:  # Leave the remainder for the next restart.
        print(
            f'WARNING: background cleanup paused for {path}: {exc}',
            file=sys.stderr, flush=True)


def start_background_cleanup():
    threads = []
    unique_paths = list(dict.fromkeys(cleanup_paths))
    for path in unique_paths:
        if not path.exists():
            continue
        thread = threading.Thread(
            target=cleanup_quarantined_tree, args=(path,), daemon=True,
            name=f'cleanup-{path.name}')
        thread.start()
        threads.append(thread)
    return threads


def move_result(source_path, destination_dir, expected_names):
    """Atomically commit one completed result and its optional Hadoop CRC."""
    source_path = Path(source_path)
    metadata = read_header_json(source_path)
    file_name = metadata['_file']
    if file_name not in expected_names:
        raise RuntimeError(
            f'result {source_path} is not part of the selected phenotype set')
    destination = Path(destination_dir) / f'{file_name}{data_suffix}'
    source_crc = crc_path(source_path)
    destination_crc = crc_path(destination)

    if destination.exists():
        existing = read_header_json(destination)
        if existing['_file'] != file_name:
            raise RuntimeError(f'conflicting resume result: {destination}')
        source_path.unlink()
    else:
        os.replace(source_path, destination)

    if source_crc.exists():
        if destination_crc.exists():
            source_crc.unlink()
        else:
            os.replace(source_crc, destination_crc)
    return metadata


def scan_stable_results(expected_names):
    """Return validated, phenotype-named results already committed to files/."""
    found = {}
    files_path.mkdir(parents=True, exist_ok=True)
    for path in files_path.glob(f'*{data_suffix}'):
        metadata = read_header_json(path)
        file_name = metadata['_file']
        expected_path = files_path / f'{file_name}{data_suffix}'
        if path.name != expected_path.name:
            raise RuntimeError(
                f'resume filename/header mismatch: {path.name} != {expected_path.name}')
        if file_name not in expected_names:
            raise RuntimeError(
                f'existing result is not part of this run: {path}')
        if file_name in found:
            raise RuntimeError(f'duplicate completed phenotype: {file_name}')
        found[file_name] = metadata
    return found


def recover_finished_work(expected_names, expected_batch=None):
    """Commit a batch that finished export but was interrupted during rename."""
    if not work_path.exists():
        return 0
    if work_path.is_symlink() or work_path.parent.resolve() != out_path.resolve():
        raise RuntimeError(f'unsafe resume work directory: {work_path}')
    if not (work_path / 'index.tsv').is_file():
        print('Discarding one incomplete resume batch.', flush=True)
        quarantine_directory(work_path, 'resume-batch-cleanup')
        return 0

    sources = list(work_path.glob(f'*{data_suffix}'))
    source_names = {read_header_json(path)['_file'] for path in sources}
    if expected_batch is not None and source_names != set(expected_batch):
        raise RuntimeError(
            'completed resume batch does not match its expected phenotype set')
    for path in sources:
        move_result(path, files_path, expected_names)
    remove_tree(work_path, out_path)
    if sources:
        print(f'Recovered completed staging batch: {len(sources)} phenotypes.',
              flush=True)
    return len(sources)


def write_index(all_names, completed):
    lines = []
    for file_name in all_names:
        metadata = completed.get(file_name)
        if metadata is None:
            continue
        path = (files_path / f'{file_name}{data_suffix}').resolve().as_uri()
        metadata_json = json.dumps(
            metadata, ensure_ascii=False, separators=(',', ':'))
        lines.append(f'{path}\t{metadata_json}\n')
    atomic_write_text(files_path / 'index.tsv', ''.join(lines))


def write_resume_state(signature, completed_count, total_count, status, hail_log):
    atomic_write_json(state_path, {
        **signature,
        'status': status,
        'completed_phenotypes': completed_count,
        'remaining_phenotypes': total_count - completed_count,
        'batch_size': batch_size,
        'spark_memory': spark_memory,
        'spark_cores': spark_cores,
        'hail_log': str(hail_log),
        'updated_at': datetime.now().astimezone().isoformat(timespec='seconds'),
    })


# Never let Hail place its default log in the source-code working directory.
log_dir = downloads_dir / 'logs'
log_dir.mkdir(parents=True, exist_ok=True)
hail_log = log_dir / (
    f'hail-genebass-{datetime.now():%Y%m%d-%H%M%S}-{os.getpid()}.log')

if pheno_regex:
    try:
        re.compile(pheno_regex, re.IGNORECASE)
    except re.error as exc:
        sys.exit(f'Invalid --pheno-regex: {exc}')

spark_conf = {
    'spark.driver.memory': spark_memory,
    'spark.executor.memory': spark_memory,
    'spark.task.maxFailures': str(task_max_failures),
    'spark.hadoop.fs.gs.impl':
        'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem',
    'spark.hadoop.fs.AbstractFileSystem.gs.impl':
        'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS',
    'spark.hadoop.fs.gs.auth.type': 'APPLICATION_DEFAULT',
    'spark.hadoop.fs.gs.requester.pays.mode': 'ENABLED',
    'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
    'spark.hadoop.fs.gs.http.connect-timeout': gcs_connect_timeout,
    'spark.hadoop.fs.gs.http.read-timeout': gcs_read_timeout,
    'spark.hadoop.fs.gs.http.max.retry': str(gcs_max_retries),
}
if '://' in gcs_connector or gcs_connector.endswith('.jar'):
    spark_conf['spark.jars'] = gcs_connector
else:
    spark_conf['spark.jars.packages'] = gcs_connector
if proxy_url:
    from urllib.parse import urlparse
    proxy = urlparse(proxy_url)
    spark_conf['spark.hadoop.fs.gs.proxy.address'] = (
        f'{proxy.hostname}:{proxy.port}')

hl.init(master=f'local[{spark_cores},{task_max_failures}]',
        app_name='genebass-gene-results-export',
        tmp_dir='file://' + tmp_dir, local_tmpdir=tmp_dir,
        log=str(hail_log), spark_conf=spark_conf)

try:
    hl.default_reference('GRCh38')
    mt = hl.read_matrix_table(source)

    if pheno_regex:
        searchable = hl.delimit([
            hl.or_else(mt.trait_type, ''), hl.or_else(mt.phenocode, ''),
            hl.or_else(mt.category, ''), hl.or_else(mt.description, '')
        ], delimiter='|')
        mt = mt.filter_cols(searchable.matches('(?i).*' + pheno_regex + '.*'))

    # Turn the official five-field column key into a safe, unique file name.
    raw_name = hl.delimit([
        hl.or_else(mt.trait_type, 'NA'), hl.or_else(mt.phenocode, 'NA'),
        hl.or_else(mt.pheno_sex, 'NA'), hl.or_else(mt.coding, 'NA'),
        hl.or_else(mt.modifier, 'NA')
    ], delimiter='__')
    mt = mt.annotate_cols(_file=raw_name.replace(r'[^A-Za-z0-9._-]+', '_'))
    mt = mt.key_cols_by(mt._file)

    if profile == 'compact':
        compact_rows = [
            'interval', 'total_variants', 'CAF', 'mean_coverage',
            'keep_gene_burden', 'keep_gene_coverage',
            'keep_gene_expected_ac', 'keep_gene_n_var']
        compact_entries = [
            'Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT',
            'BETA_Burden', 'SE_Burden',
            'Pvalue.NA', 'Pvalue_Burden.NA', 'Pvalue_SKAT.NA',
            'BETA_Burden.NA', 'SE_Burden.NA',
            'total_variants_pheno', 'expected_AC', 'keep_entry_expected_ac']
        row_fields = set(mt.row_value.dtype.fields)
        entry_fields = set(mt.entry.dtype.fields)
        selected_rows = [name for name in compact_rows if name in row_fields]
        selected_entries = [
            name for name in compact_entries if name in entry_fields]
        missing_rows = [name for name in compact_rows if name not in row_fields]
        missing_entries = [
            name for name in compact_entries if name not in entry_fields]
        if missing_rows:
            print('Compact profile: unavailable row fields: ' +
                  ', '.join(missing_rows), flush=True)
        if missing_entries:
            print('Compact profile: unavailable entry fields: ' +
                  ', '.join(missing_entries), flush=True)
        mt = mt.select_rows(
            **{name: mt.row_value[name] for name in selected_rows})
        mt = mt.select_entries(
            **{name: mt.entry[name] for name in selected_entries})

    # Table.select() retains key fields automatically.  Passing '_file'
    # explicitly would be treated by Hail as an attempt to overwrite the key.
    column_rows = mt.cols().select().collect()
    all_names = [row._file for row in column_rows]
    n_cols = len(all_names)
    if n_cols == 0:
        sys.exit('No phenotype matched --pheno-regex.')
    if len(set(all_names)) != n_cols:
        raise RuntimeError('generated phenotype filenames are not unique')
    expected_names = set(all_names)
    names_digest = hashlib.sha256(
        ('\n'.join(all_names) + '\n').encode('utf-8')).hexdigest()
    signature = {
        'source_matrix_table': source,
        'profile': profile,
        'phenotype_regex': pheno_regex or None,
        'phenotype_count': n_cols,
        'phenotype_names_sha256': names_digest,
        'reference_genome': 'GRCh38',
    }

    files_path.mkdir(parents=True, exist_ok=True)
    for pattern in ('.resume-batch-cleanup-*',):
        for path in out_path.glob(pattern):
            if path.is_dir() and not path.is_symlink():
                cleanup_paths.append(path)
    recover_finished_work(expected_names)
    completed = scan_stable_results(expected_names)
    write_index(all_names, completed)
    write_resume_state(
        signature, len(completed), n_cols, 'running', hail_log)
    cleanup_threads = start_background_cleanup()

    print(
        f'Resume inventory: {len(completed):,}/{n_cols:,} phenotypes complete; '
        f'{n_cols - len(completed):,} remaining.', flush=True)
    missing_names = [name for name in all_names if name not in completed]
    total_resume_batches = (
        (len(missing_names) + batch_size - 1) // batch_size)
    for batch_number, start in enumerate(
            range(0, len(missing_names), batch_size), start=1):
        batch_names = missing_names[start:start + batch_size]
        print(
            f'Resume batch {batch_number}/{total_resume_batches}: '
            f'{len(batch_names)} phenotypes; '
            f'{len(completed)}/{n_cols} already committed.', flush=True)

        batch_name_set = set(batch_names)
        batch_mt = mt.filter_cols(
            hl.literal(batch_name_set).contains(mt._file))
        hl.experimental.export_entries_by_col(
            batch_mt, work_path.resolve().as_uri(),
            batch_size=len(batch_names), bgzip=True,
            header_json_in_file=True, use_string_key_as_file_name=True)
        recover_finished_work(
            expected_names, expected_batch=batch_name_set)
        completed = scan_stable_results(expected_names)
        write_index(all_names, completed)
        write_resume_state(
            signature, len(completed), n_cols, 'running', hail_log)
        print(
            f'Committed resume batch {batch_number}: '
            f'{len(completed):,}/{n_cols:,} complete.', flush=True)

    completed = scan_stable_results(expected_names)
    if len(completed) != n_cols:
        raise RuntimeError(
            f'export ended with {len(completed)}/{n_cols} completed phenotypes')
    write_index(all_names, completed)
    write_resume_state(signature, n_cols, n_cols, 'complete', hail_log)
except KeyboardInterrupt:
    print(
        'Interrupted by user; committed result files, index, and resume state '
        'remain safe.',
        file=sys.stderr, flush=True)
    raise SystemExit(130)
finally:
    hl.stop()

atomic_write_json(run_info_path, {
    **signature,
    'status': 'complete',
    'completed_phenotypes': n_cols,
    'batch_size': batch_size,
    'spark_memory': spark_memory,
    'spark_cores': spark_cores,
    'gcs_connector': gcs_connector,
    'proxy': proxy_url or None,
    'gcp_project': gcp_project,
    'hail_log': str(hail_log),
    'resumable': True,
})
print('Export completed.', flush=True)
PY

stop_cleanup_watcher
trap - EXIT

echo "[$(date '+%F %T')] gene-based export completed"
echo "Gene data:     $files_dir/*.tsv.bgz"
echo "Gene manifest: $files_dir/index.tsv"
echo "Gene run info: $dir_out/run.info.json"


# 🚩 Optional single-variant export
# Keep the multi-terabyte single-variant export opt-in.
single_variant="${single_variant:-FALSE}"
if [[ "$single_variant" == TRUE ]]; then
  echo "[$(date '+%F %T')] starting single-variant GWAS export"
  python3 - <<'PY_SINGLE_VARIANT'
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import threading
from datetime import datetime
from pathlib import Path

import hail as hl

source = os.environ['WES_SINGLE_SOURCE_MT']
out_path = Path(os.environ['WES_SINGLE_DIR_OUT'])
tmp_dir = os.environ['WES_TMPDIR']
downloads_dir = Path(os.environ['WES_DOWNLOADS_DIR'])
pheno_regex = os.environ['WES_PHENO_REGEX']
batch_size = int(os.environ['WES_SINGLE_BATCH_SIZE'])
spark_memory = os.environ['WES_SPARK_MEMORY']
spark_cores = int(os.environ['WES_SPARK_CORES'])
task_max_failures = int(os.environ['WES_TASK_MAX_FAILURES'])
gcs_connect_timeout = os.environ['WES_GCS_CONNECT_TIMEOUT']
gcs_read_timeout = os.environ['WES_GCS_READ_TIMEOUT']
gcs_max_retries = int(os.environ['WES_GCS_MAX_RETRIES'])
profile = os.environ['WES_PROFILE']
gcs_connector = os.environ['WES_GCS_CONNECTOR']
proxy_url = os.environ['WES_PROXY_URL']
gcp_project = os.environ['WES_GCP_PROJECT']
allow_low_space = os.environ['WES_ALLOW_LOW_SPACE'] == 'TRUE'
configured_gib_per_pheno = float(
    os.environ['WES_SINGLE_ESTIMATE_GIB_PER_PHENO'])

files_path = out_path / 'files'
work_path = out_path / '.resume_batch'
state_path = out_path / 'resume.state.json'
run_info_path = out_path / 'run.info.json'
data_suffix = '.tsv.bgz'
safe_file_re = re.compile(r'^[A-Za-z0-9._-]+$')
cleanup_paths = []


def atomic_write_text(path, text):
    path = Path(path)
    tmp = path.with_name(f'.{path.name}.tmp-{os.getpid()}')
    with tmp.open('w', encoding='utf-8', newline='\n') as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def atomic_write_json(path, value):
    atomic_write_text(
        path, json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def read_header_json(path):
    path = Path(path)
    try:
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            line = handle.readline()
    except (OSError, EOFError) as exc:
        raise RuntimeError(f'invalid BGZF result {path}: {exc}') from exc
    if not line.startswith('#'):
        raise RuntimeError(f'missing phenotype JSON header: {path}')
    try:
        metadata = json.loads(line[1:])
    except json.JSONDecodeError as exc:
        raise RuntimeError(f'invalid phenotype JSON in {path}: {exc}') from exc
    file_name = metadata.get('_file')
    if not isinstance(file_name, str) or not safe_file_re.fullmatch(file_name):
        raise RuntimeError(f'unsafe or missing _file in {path}: {file_name!r}')
    return metadata


def crc_path(path):
    path = Path(path)
    return path.with_name(f'.{path.name}.crc')


def quarantine_work(path):
    path = Path(path)
    if not path.exists():
        return None
    if path.is_symlink() or path.parent.resolve() != out_path.resolve():
        raise RuntimeError(f'refusing to quarantine unexpected path: {path}')
    destination = out_path / (
        f'.resume-batch-cleanup-{datetime.now():%Y%m%d-%H%M%S}-'
        f'{os.getpid()}')
    os.replace(path, destination)
    cleanup_paths.append(destination)
    print(f'Quarantined incomplete variant batch: {destination.name}',
          flush=True)
    return destination


def cleanup_quarantined_tree(path):
    path = Path(path)
    try:
        if path.is_symlink() or path.parent.resolve() != out_path.resolve():
            raise RuntimeError(f'unsafe cleanup quarantine: {path}')
        shutil.rmtree(path)
        print(f'Background cleanup completed: {path.name}', flush=True)
    except Exception as exc:
        print(f'WARNING: background cleanup paused for {path}: {exc}',
              file=sys.stderr, flush=True)


def start_background_cleanup():
    threads = []
    for path in list(dict.fromkeys(cleanup_paths)):
        if not path.exists():
            continue
        thread = threading.Thread(
            target=cleanup_quarantined_tree, args=(path,), daemon=True,
            name=f'cleanup-{path.name}')
        thread.start()
        threads.append(thread)
    return threads


def move_result(source_path, expected_names):
    source_path = Path(source_path)
    metadata = read_header_json(source_path)
    file_name = metadata['_file']
    if file_name not in expected_names:
        raise RuntimeError(
            f'variant result is outside selected phenotypes: {source_path}')
    destination = files_path / f'{file_name}{data_suffix}'
    source_crc = crc_path(source_path)
    destination_crc = crc_path(destination)
    if destination.exists():
        existing = read_header_json(destination)
        if existing['_file'] != file_name:
            raise RuntimeError(f'conflicting resume result: {destination}')
        source_path.unlink()
    else:
        os.replace(source_path, destination)
    if source_crc.exists():
        if destination_crc.exists():
            source_crc.unlink()
        else:
            os.replace(source_crc, destination_crc)
    return metadata


def scan_stable_results(expected_names):
    found = {}
    files_path.mkdir(parents=True, exist_ok=True)
    for path in files_path.glob(f'*{data_suffix}'):
        metadata = read_header_json(path)
        file_name = metadata['_file']
        if path.name != f'{file_name}{data_suffix}':
            raise RuntimeError(
                f'resume filename/header mismatch: {path.name}')
        if file_name not in expected_names:
            raise RuntimeError(f'unexpected existing variant result: {path}')
        if file_name in found:
            raise RuntimeError(f'duplicate variant phenotype: {file_name}')
        found[file_name] = metadata
    return found


def recover_finished_work(expected_names, expected_batch=None):
    if not work_path.exists():
        return 0
    if work_path.is_symlink() or work_path.parent.resolve() != out_path.resolve():
        raise RuntimeError(f'unsafe resume work directory: {work_path}')
    if not (work_path / 'index.tsv').is_file():
        print('Discarding one incomplete single-variant batch.', flush=True)
        quarantine_work(work_path)
        return 0
    sources = list(work_path.glob(f'*{data_suffix}'))
    source_names = {read_header_json(path)['_file'] for path in sources}
    if expected_batch is not None and source_names != set(expected_batch):
        raise RuntimeError(
            'completed variant batch differs from its expected phenotypes')
    for path in sources:
        move_result(path, expected_names)
    shutil.rmtree(work_path)
    if sources:
        print(f'Recovered completed variant batch: {len(sources)} phenotypes.',
              flush=True)
    return len(sources)


def write_index(all_names, completed):
    lines = []
    for file_name in all_names:
        metadata = completed.get(file_name)
        if metadata is None:
            continue
        uri = (files_path / f'{file_name}{data_suffix}').resolve().as_uri()
        metadata_json = json.dumps(
            metadata, ensure_ascii=False, separators=(',', ':'))
        lines.append(f'{uri}\t{metadata_json}\n')
    atomic_write_text(files_path / 'index.tsv', ''.join(lines))


def write_resume_state(signature, completed_count, total_count, status,
                       hail_log, space=None):
    value = {
        **signature,
        'status': status,
        'completed_phenotypes': completed_count,
        'remaining_phenotypes': total_count - completed_count,
        'batch_size': batch_size,
        'spark_memory': spark_memory,
        'spark_cores': spark_cores,
        'hail_log': str(hail_log),
        'updated_at': datetime.now().astimezone().isoformat(
            timespec='seconds'),
    }
    if space:
        value['space_estimate'] = space
    atomic_write_json(state_path, value)


log_dir = downloads_dir / 'logs'
log_dir.mkdir(parents=True, exist_ok=True)
hail_log = log_dir / (
    f'hail-genebass-single-variant-'
    f'{datetime.now():%Y%m%d-%H%M%S}-{os.getpid()}.log')

if pheno_regex:
    try:
        re.compile(pheno_regex, re.IGNORECASE)
    except re.error as exc:
        sys.exit(f'Invalid --pheno-regex: {exc}')

spark_conf = {
    'spark.driver.memory': spark_memory,
    'spark.executor.memory': spark_memory,
    'spark.task.maxFailures': str(task_max_failures),
    'spark.hadoop.fs.gs.impl':
        'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem',
    'spark.hadoop.fs.AbstractFileSystem.gs.impl':
        'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS',
    'spark.hadoop.fs.gs.auth.type': 'APPLICATION_DEFAULT',
    'spark.hadoop.fs.gs.requester.pays.mode': 'ENABLED',
    'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
    'spark.hadoop.fs.gs.http.connect-timeout': gcs_connect_timeout,
    'spark.hadoop.fs.gs.http.read-timeout': gcs_read_timeout,
    'spark.hadoop.fs.gs.http.max.retry': str(gcs_max_retries),
}
if '://' in gcs_connector or gcs_connector.endswith('.jar'):
    spark_conf['spark.jars'] = gcs_connector
else:
    spark_conf['spark.jars.packages'] = gcs_connector
if proxy_url:
    from urllib.parse import urlparse
    proxy = urlparse(proxy_url)
    spark_conf['spark.hadoop.fs.gs.proxy.address'] = (
        f'{proxy.hostname}:{proxy.port}')

hl.init(master=f'local[{spark_cores},{task_max_failures}]',
        app_name='genebass-single-variant-results-export',
        tmp_dir='file://' + tmp_dir, local_tmpdir=tmp_dir,
        log=str(hail_log), spark_conf=spark_conf)

try:
    hl.default_reference('GRCh38')
    mt = hl.read_matrix_table(source)
    if pheno_regex:
        searchable = hl.delimit([
            hl.or_else(mt.trait_type, ''), hl.or_else(mt.phenocode, ''),
            hl.or_else(mt.category, ''), hl.or_else(mt.description, '')
        ], delimiter='|')
        mt = mt.filter_cols(searchable.matches('(?i).*' + pheno_regex + '.*'))

    raw_name = hl.delimit([
        hl.or_else(mt.trait_type, 'NA'), hl.or_else(mt.phenocode, 'NA'),
        hl.or_else(mt.pheno_sex, 'NA'), hl.or_else(mt.coding, 'NA'),
        hl.or_else(mt.modifier, 'NA')
    ], delimiter='__')
    mt = mt.annotate_cols(
        _file=raw_name.replace(r'[^A-Za-z0-9._-]+', '_'))
    mt = mt.key_cols_by(mt._file)

    if profile == 'compact':
        compact_rows = [
            'markerID', 'gene', 'annotation', 'call_stats',
            'keep_var_expected_ac', 'keep_var_annt']
        compact_entries = [
            'AC', 'AF', 'BETA', 'SE', 'AF.Cases', 'AF.Controls', 'Pvalue',
            'expected_AC', 'keep_entry_expected_ac']
        row_fields = set(mt.row_value.dtype.fields)
        entry_fields = set(mt.entry.dtype.fields)
        selected_rows = [name for name in compact_rows if name in row_fields]
        selected_entries = [
            name for name in compact_entries if name in entry_fields]
        missing_rows = [name for name in compact_rows if name not in row_fields]
        missing_entries = [
            name for name in compact_entries if name not in entry_fields]
        if missing_rows:
            print('Variant compact profile unavailable row fields: ' +
                  ', '.join(missing_rows), flush=True)
        if missing_entries:
            print('Variant compact profile unavailable entry fields: ' +
                  ', '.join(missing_entries), flush=True)
        mt = mt.select_rows(
            **{name: mt.row_value[name] for name in selected_rows})
        mt = mt.select_entries(
            **{name: mt.entry[name] for name in selected_entries})

    column_rows = mt.cols().select().collect()
    all_names = [row._file for row in column_rows]
    n_cols = len(all_names)
    if n_cols == 0:
        sys.exit('No single-variant phenotype matched --pheno-regex.')
    if len(set(all_names)) != n_cols:
        raise RuntimeError('generated variant phenotype names are not unique')
    expected_names = set(all_names)
    names_digest = hashlib.sha256(
        ('\n'.join(all_names) + '\n').encode('utf-8')).hexdigest()
    signature = {
        'result_type': 'single_variant',
        'source_matrix_table': source,
        'profile': profile,
        'phenotype_regex': pheno_regex or None,
        'phenotype_count': n_cols,
        'phenotype_names_sha256': names_digest,
        'reference_genome': 'GRCh38',
    }

    files_path.mkdir(parents=True, exist_ok=True)
    for path in out_path.glob('.resume-batch-cleanup-*'):
        if path.is_dir() and not path.is_symlink():
            cleanup_paths.append(path)
    recover_finished_work(expected_names)
    completed = scan_stable_results(expected_names)
    write_index(all_names, completed)
    cleanup_threads = start_background_cleanup()

    missing_names = [name for name in all_names if name not in completed]
    existing_bytes = sum(
        (files_path / f'{name}{data_suffix}').stat().st_size
        for name in completed)
    configured_bytes = configured_gib_per_pheno * 1024 ** 3
    if completed:
        observed_bytes = existing_bytes / len(completed)
        bytes_per_pheno = max(configured_bytes, observed_bytes)
        estimate_basis = 'max(configured, observed completed-file average)'
    else:
        bytes_per_pheno = configured_bytes
        estimate_basis = 'configured initial estimate'
    estimated_remaining = int(bytes_per_pheno * len(missing_names))
    free_bytes = shutil.disk_usage(out_path).free
    reserve_bytes = int(max(50 * 1024 ** 3,
                            bytes_per_pheno * batch_size * 2))
    required_bytes = estimated_remaining + reserve_bytes
    space = {
        'basis': estimate_basis,
        'estimated_gib_per_phenotype': round(bytes_per_pheno / 1024 ** 3, 3),
        'estimated_remaining_gib': round(estimated_remaining / 1024 ** 3, 1),
        'temporary_and_safety_reserve_gib': round(reserve_bytes / 1024 ** 3, 1),
        'free_gib': round(free_bytes / 1024 ** 3, 1),
        'low_space_override': allow_low_space,
    }
    write_resume_state(
        signature, len(completed), n_cols, 'running', hail_log, space)
    print(
        f'Variant resume inventory: {len(completed):,}/{n_cols:,} complete; '
        f'{len(missing_names):,} remaining.', flush=True)
    print(
        f'Variant disk estimate: {space["estimated_remaining_gib"]:,.1f} GiB '
        f'remaining plus {space["temporary_and_safety_reserve_gib"]:,.1f} GiB '
        f'reserve; {space["free_gib"]:,.1f} GiB free.', flush=True)
    if missing_names and free_bytes < required_bytes and not allow_low_space:
        write_resume_state(
            signature, len(completed), n_cols, 'insufficient_space',
            hail_log, space)
        raise SystemExit(
            'ERROR: insufficient disk space for the single-variant export. '
            'Choose a larger --dir-out, restrict --pheno-regex, or use '
            '--allow-low-space TRUE to override the safety check.')

    total_batches = (
        (len(missing_names) + batch_size - 1) // batch_size)
    for batch_number, start in enumerate(
            range(0, len(missing_names), batch_size), start=1):
        batch_names = missing_names[start:start + batch_size]
        print(
            f'Variant resume batch {batch_number}/{total_batches}: '
            f'{len(batch_names)} phenotypes; '
            f'{len(completed)}/{n_cols} already committed.', flush=True)
        batch_name_set = set(batch_names)
        batch_mt = mt.filter_cols(
            hl.literal(batch_name_set).contains(mt._file))
        hl.experimental.export_entries_by_col(
            batch_mt, work_path.resolve().as_uri(),
            batch_size=len(batch_names), bgzip=True,
            header_json_in_file=True, use_string_key_as_file_name=True)
        recover_finished_work(
            expected_names, expected_batch=batch_name_set)
        completed = scan_stable_results(expected_names)
        write_index(all_names, completed)
        write_resume_state(
            signature, len(completed), n_cols, 'running', hail_log, space)
        print(
            f'Committed variant batch {batch_number}: '
            f'{len(completed):,}/{n_cols:,} complete.', flush=True)

    completed = scan_stable_results(expected_names)
    if len(completed) != n_cols:
        raise RuntimeError(
            f'variant export ended with {len(completed)}/{n_cols} complete')
    write_index(all_names, completed)
    write_resume_state(
        signature, n_cols, n_cols, 'complete', hail_log, space)
except KeyboardInterrupt:
    print(
        'Single-variant export interrupted; committed files and resume state '
        'remain safe.', file=sys.stderr, flush=True)
    raise SystemExit(130)
finally:
    hl.stop()

atomic_write_json(run_info_path, {
    **signature,
    'status': 'complete',
    'completed_phenotypes': n_cols,
    'batch_size': batch_size,
    'spark_memory': spark_memory,
    'spark_cores': spark_cores,
    'gcs_connector': gcs_connector,
    'proxy': proxy_url or None,
    'gcp_project': gcp_project,
    'hail_log': str(hail_log),
    'space_estimate': space,
    'resumable': True,
})
print('Single-variant export completed.', flush=True)
PY_SINGLE_VARIANT
else
  echo "[$(date '+%F %T')] single-variant export skipped; enable with --single-variant TRUE"
fi

echo "[$(date '+%F %T')] completed"
if [[ "$single_variant" == TRUE ]]; then
  echo "Variant data:     $single_dir_out/files/*.tsv.bgz"
  echo "Variant manifest: $single_dir_out/files/index.tsv"
  echo "Variant run info: $single_dir_out/run.info.json"
fi
