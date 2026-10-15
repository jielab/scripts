#!/usr/bin/env bash
# Reproducible Linux/WSL environment and opt-in reference preparation.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=f/conda_env.sh
source "$HERE/f/conda_env.sh"

usage() {
  cat <<'EOF'
用法：./install.sh [options]
  默认                 安装/检查固定版本的 metag 环境，不下载参考或运行全量分析
  --prepare            另行下载/建立共享参考、ICTV 映射并审计
  --datasets CSV       与 --prepare 配合，整理指定已下载数据集的 metadata/raw
  --data-root DIR      数据集父目录，默认 /mnt/f/metag
  --reference-root DIR 共享参考目录，默认 DATA_ROOT/ref_kraken2_ictv
  --refresh-metadata   路径迁移后备份并重建清单；人工注释需从备份合并
  --threads N          建库/审计线程，默认 16
  --download-threads N 参考下载线程，默认 4
  --max-db-gb N        hash 表上限 GiB，默认 24；不限制下载和临时文件
  --check              只检查现有环境，不安装、不联网、不准备数据
  --dry-run            只打印安装/准备计划，不写文件或联网
  --direct             此次安装/准备使用直连，保留父 shell 的代理设置
  -h, --help           显示帮助

新电脑（Linux x86_64 或 Windows 的 WSL2）：
  ./install.sh
  ./install.sh --prepare --data-root /mnt/f/metag --datasets emp,wbe,长江,厦门污水,全球污水,纳米布沙漠

没有 Conda 时，安装固定版本 Miniforge 到 ~/miniforge3；已有 Conda 自动复用。
可用 METAG_CONDA_ROOT=/其他位置/miniforge3 显式指定；之后运行 metag.sh 时也须设置。
无需 source、conda activate 或 conda init。已有环境不自动升级或删除。
环境下载失败可重跑；参考准备很大，失败时停止，排查后重跑同一命令。
--prepare 不下载原始 FASTQ，不替你运行 pilot/full；完成后按 metag.sh --help 操作。
完整复现旧结果须迁移旧参考、ICTV/accession 缓存、metadata 与原始 FASTQ；
重新下载 NCBI 参考会随上游更新，不能保证与旧库逐字节一致。
EOF
}

PREPARE=0; CHECK=0; DRY_RUN=0; DIRECT=0; REFRESH=0
DATA_ROOT=/mnt/f/metag; REF_ROOT=""; DATASETS=""
THREADS=16; DOWNLOAD_THREADS=4; MAX_GB=24
fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
need_value() { [[ $# -ge 2 && -n "$2" && "$2" != --* ]] || fail "$1 需要参数"; }
while (($#)); do
  case "$1" in
    --prepare) PREPARE=1; shift;;
    --check) CHECK=1; shift;;
    --dry-run) DRY_RUN=1; shift;;
    --direct) DIRECT=1; shift;;
    --refresh-metadata) REFRESH=1; shift;;
    --datasets) need_value "$@"; DATASETS="$2"; shift 2;;
    --data-root) need_value "$@"; DATA_ROOT="$2"; shift 2;;
    --reference-root) need_value "$@"; REF_ROOT="$2"; shift 2;;
    --threads) need_value "$@"; THREADS="$2"; shift 2;;
    --download-threads) need_value "$@"; DOWNLOAD_THREADS="$2"; shift 2;;
    --max-db-gb) need_value "$@"; MAX_GB="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) fail "未知参数：$1";;
  esac
done
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || fail '锁定环境仅支持 Linux x86_64 / WSL2 x86_64。'
[[ "$THREADS" =~ ^[1-9][0-9]*$ && "$DOWNLOAD_THREADS" =~ ^[1-9][0-9]*$ && "$MAX_GB" =~ ^[1-9][0-9]*$ ]] || fail '线程和容量须为正整数。'
if ((CHECK && PREPARE)); then fail '--check 不能和 --prepare 一起使用。'; fi
if [[ -n "$DATASETS" ]] && (( ! PREPARE )); then fail '--datasets 需要 --prepare。'; fi
if ((REFRESH)) && [[ -z "$DATASETS" ]]; then fail '--refresh-metadata 需要 --prepare --datasets。'; fi
DATA_ROOT="$(realpath -m -- "$DATA_ROOT")"
REF_ROOT="$(realpath -m -- "${REF_ROOT:-$DATA_ROOT/ref_kraken2_ictv}")"
DB="$REF_ROOT/kraken2_broad_v1"
MAP="$REF_ROOT/virus_map_v1/virus_genome_types.sqlite"
LOCKFILE="$HERE/f/conda-linux-64.lock"
[[ -s "$LOCKFILE" ]] || fail "缺少环境锁文件：$LOCKFILE"
CONDA_ROOT="$(metag_find_conda || true)"
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
ENV_PREFIX="$CONDA_ROOT/envs/metag"
export METAG_CONDA_ROOT="$CONDA_ROOT"

print_cmd() { printf '%q ' "$@"; printf '\n'; }
if ((DRY_RUN)); then
  printf 'PLAN ONLY：Conda=%s；环境=%s；锁文件=%s\n' "$CONDA_ROOT" "$ENV_PREFIX" "$LOCKFILE"
  if ((CHECK)); then echo '检查包版本及命令，不联网。'; exit 0; fi
  echo '无 Conda 时下载并校验固定版本 Miniforge；缺少环境时按锁文件安装，已有环境只检查。'
  if [[ -n "$DATASETS" ]]; then
    plan=("$HERE/metag.sh" prepare "$DATASETS" --data-root "$DATA_ROOT" --execute)
    if ((REFRESH)); then plan+=(--refresh); fi
    print_cmd "${plan[@]}"
  fi
  if ((PREPARE)); then
    print_cmd "$HERE/metag.sh" build-db --db "$DB" --threads "$THREADS" --download-threads "$DOWNLOAD_THREADS" --max-db-gb "$MAX_GB" --execute
    echo "获取固定 ICTV VMR，续用 accession 缓存并建立/核验映射：$MAP"
    print_cmd "$HERE/metag.sh" db-audit --db "$DB" --virus-map "$MAP" --threads "$THREADS"
  fi
  exit 0
fi

# Script-local only: users need not unset anything in their shell.
if ((DIRECT)); then
  unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY
fi
retry() {
  local attempt
  for attempt in 1 2 3; do
    if "$@"; then return 0; fi
    if ((attempt < 3)); then
      printf '步骤失败，%s 秒后重试（%s/3）…\n' "$((attempt * 5))" "$attempt" >&2
      sleep "$((attempt * 5))"
    fi
  done
  return 1
}

if (( ! CHECK )); then
  CACHE="${XDG_CACHE_HOME:-$HOME/.cache}/metag"
  mkdir -p "$CACHE"
  exec 8>"$CACHE/install.lock"
  flock -n 8 || fail '另一个 install.sh 正在运行。'
  if [[ ! -r "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
    [[ ! -e "$CONDA_ROOT" ]] || fail "Conda 目录存在但不完整，已保留：$CONDA_ROOT"
    command -v curl >/dev/null || fail '需要系统 curl（Ubuntu/WSL 可安装 curl ca-certificates）。'
    INSTALLER="$CACHE/Miniforge3-26.7.2-0-Linux-x86_64.sh"
    INSTALLER_SHA=281b0ac7d550802efc81af633225a5e6116d29ae72f3ab4eae7168c3931a4c05
    check_installer() { [[ -f "$INSTALLER" && "$(sha256sum "$INSTALLER" | cut -d' ' -f1)" == "$INSTALLER_SHA" ]]; }
    if ! check_installer; then
      curl --fail --location --show-error --retry 5 --retry-all-errors --connect-timeout 30 \
        --speed-limit 1024 --speed-time 60 \
        'https://github.com/conda-forge/miniforge/releases/download/26.7.2-0/Miniforge3-26.7.2-0-Linux-x86_64.sh' \
        --output "$INSTALLER.partial"
      [[ "$(sha256sum "$INSTALLER.partial" | cut -d' ' -f1)" == "$INSTALLER_SHA" ]] || fail 'Miniforge SHA256 校验失败，未执行安装包。'
      mv -- "$INSTALLER.partial" "$INSTALLER"
    fi
    bash "$INSTALLER" -b -p "$CONDA_ROOT"
  fi
  if [[ ! -f "$ENV_PREFIX/conda-meta/history" ]]; then
    retry "$CONDA_ROOT/bin/conda" create --yes --prefix "$ENV_PREFIX" --no-default-packages \
      --override-channels -c conda-forge -c bioconda --file "$LOCKFILE" || fail '环境安装失败；网络恢复后重跑，包缓存会保留。'
  fi
fi
[[ -x "$CONDA_ROOT/bin/conda" && -f "$ENV_PREFIX/conda-meta/history" ]] || fail 'metag 环境未安装，先运行 ./install.sh。'
metag_activate
# Installation checks the managed interpreter, independent of an analysis override.
"$CONDA_PREFIX/bin/python3" - "$LOCKFILE" "$CONDA_PREFIX" <<'PY'
import json, sys, urllib.parse
from pathlib import Path
lock, prefix = map(Path, sys.argv[1:])
expected = {}
for line in lock.read_text().splitlines():
    if not line.startswith('https://'):
        continue
    url, checksum = line.rsplit('#', 1)
    filename = urllib.parse.unquote(urllib.parse.urlsplit(url).path.rsplit('/', 1)[1])
    stem = filename.removesuffix('.conda').removesuffix('.tar.bz2')
    name, version, build = stem.rsplit('-', 2)
    expected[name] = (version, build, checksum)
actual = {}
for path in (prefix / 'conda-meta').glob('*.json'):
    record = json.loads(path.read_text())
    actual[record['name']] = (record['version'], record['build'], record.get('md5'))
different = [name for name, spec in expected.items() if actual.get(name) != spec]
extra = sorted(set(actual) - set(expected))
if different or extra:
    raise SystemExit('现有环境与锁文件不同；不自动修改正在使用的环境。差异：' + ', '.join(different + extra)
                     + '\n请用 METAG_CONDA_ROOT 指定另一个 Conda 安装位置。')
import openpyxl
print(f'PASS：{len(expected)} 个 Conda 包版本/build/MD5 与锁文件一致；Python {sys.version.split()[0]}，openpyxl {openpyxl.__version__}')
PY
for program in kraken2 kraken2-build kraken2-inspect k2 fastp seqkit pigz bowtie2 bowtie2-build samtools multiqc blastn; do
  command -v "$program" >/dev/null || fail "缺少程序：$program"
done
if ! "$CONDA_PREFIX/bin/python3" "$HERE/f/k2_compat.py" --help >/dev/null 2>&1; then
  fail 'k2 下载兼容层与当前安装的软件不匹配。'
fi
kraken2 --version
fastp --version
"$HERE/metag.sh" --version
echo 'PASS：环境和入口检查完成。'
if ((CHECK || ! PREPARE)); then exit 0; fi

mkdir -p "$REF_ROOT"
exec 7>"$REF_ROOT/.metag_prepare.lock"
flock -n 7 || fail '该共享参考已有 install.sh --prepare 在运行。'
if [[ -e "$DB/.metag_build.lock" ]]; then
  flock -n "$DB/.metag_build.lock" true || fail '该数据库正在构建；等待现有任务，不重复启动。'
fi
if [[ -n "$DATASETS" ]]; then
  "$CONDA_PREFIX/bin/python3" - "$HERE" "$DATA_ROOT" "$DATASETS" "$REFRESH" <<'PY'
import csv, sys
from pathlib import Path
here, root = map(Path, sys.argv[1:3])
names = sys.argv[3].split(',')
with (here / 'f/datasets.example.tsv').open() as stream:
    known = {row['dataset'] for row in csv.DictReader(stream, delimiter='\t')}
if len(names) != len(set(names)) or any(name not in known for name in names):
    raise SystemExit('数据集名称为空、重复或未登记；先查看 ./metag.sh datasets。')
for name in names:
    folder = root / name
    if not folder.is_dir():
        raise SystemExit(f'先放置已下载的原始 FASTQ：{folder}')
    manifest = folder / 'metadata/files.tsv'
    if manifest.exists() and sys.argv[4] != '1':
        with manifest.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if not Path(row['raw_path']).is_relative_to(folder):
                    raise SystemExit(f'{name} 清单仍指向旧路径；迁移原始 FASTQ 后用 --refresh-metadata 重建，并从备份合并人工注释。')
    if sys.argv[4] == '1' and (folder / 'raw').exists():
        for path in (folder / 'raw').rglob('*'):
            if path.is_symlink() and not path.resolve().is_relative_to(folder):
                raise SystemExit(f'旧 raw 链接仍指向其他位置：{path}。迁移原始 FASTQ，保留旧链接备份后再重建；不要删除原始数据。')
PY
  args=(prepare "$DATASETS" --data-root "$DATA_ROOT" --execute)
  if ((REFRESH)); then args+=(--refresh); fi
  "$HERE/metag.sh" "${args[@]}"
fi
"$HERE/metag.sh" build-db --db "$DB" --threads "$THREADS" --download-threads "$DOWNLOAD_THREADS" --max-db-gb "$MAX_GB" --execute

VMR_NAME=VMR_MSL41.v1.20260729.xlsx
VMR_SHA=b79b5d82a1b3b8e9dd5e19afe8fe1a8f441267474918a7cefa8ae4913adf45bb
VMR_URL="https://ictv.global/sites/default/files/VMR/$VMR_NAME"
VMR_RECEIPT="$REF_ROOT/ictv_vmr_v1/${VMR_NAME%.xlsx}.source.json"
if [[ -f "$VMR_RECEIPT" ]]; then
  # Older local fetches used /vmr/current; fetch validates their saved receipt.
  VMR_URL="$("$CONDA_PREFIX/bin/python3" -c 'import json,sys; print(json.load(open(sys.argv[1]))["source_url"])' "$VMR_RECEIPT")"
fi
VMR="$("$HERE/metag.sh" virus fetch --url "$VMR_URL" --out-dir "$REF_ROOT/ictv_vmr_v1")"
[[ "$(sha256sum "$VMR" | cut -d' ' -f1)" == "$VMR_SHA" ]] || fail 'ICTV 文件与固定版本校验值不同；已停止，请使用新的参考目录。'
if [[ ! -e "$REF_ROOT/virus_map_v1" ]]; then
  "$HERE/metag.sh" virus cache-accessions --vmr "$VMR" --cache-dir "$REF_ROOT/ictv_accessions_v1"
  "$HERE/metag.sh" virus build --db "$DB" --vmr "$VMR" \
    --accession-taxid "$REF_ROOT/ictv_accessions_v1/accession_taxid.tsv" --out "$REF_ROOT/virus_map_v1"
fi
"$HERE/metag.sh" db-audit --db "$DB" --virus-map "$MAP" --threads "$THREADS"
"$CONDA_PREFIX/bin/python3" - "$DB/metag.audit.json" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    audit = json.load(stream)
if (audit.get('schema') != 'v2' or audit.get('missing_groups') != []
        or audit.get('missing_virus_classes') != [] or not audit.get('virus_types_tested')):
    raise SystemExit('参考覆盖未通过：要求七大基础类以及病毒 DNA/RNA 均有实际索引支持。')
PY
echo 'PASS：参考准备和覆盖审计完成。请先运行 pilot，成功后再 full。'
printf 'DB=%s\nMAP=%s\n' "$DB" "$MAP"
