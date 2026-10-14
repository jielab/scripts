#!/usr/bin/env bash
# Public entry point; runtime modules live in f/.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
METAG_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'HELP'
用法：./metag.sh COMMAND [数据集1,数据集2,...] [options]
换电脑先运行 ./install.sh；环境/参考准备选项见 ./install.sh --help。

命令：
  datasets          只读查看数据目录、下载标记、非 ENA 文件名和已有清单
  prepare           建立/复用清单并整理 raw 软链接；默认计划，--execute 执行
  pilot             每数据集选至多 2 个有效 run，校验并试跑；--execute 执行
  full              按数据集顺序审计、MD5、QC、分类、汇总；--execute 执行
  build-db          下载参考序列并建库；默认计划，--execute 执行
  virus             获取 ICTV 数据、构建病毒 DNA/RNA/unknown 映射
  db-audit          核查数据库基础类别和病毒 DNA/RNA 覆盖
  manifest          单数据集的 ENA 文件清单与下载状态检查
  organize / verify 单数据集 raw 整理 / MD5、gzip 校验
  run / summary     单数据集分类分流 / 汇总
  metadata-upgrade  保留人工值并补充实验字段
  align             候选序列的独立竞争比对

Usage examples：共享参考准备完成后，复制以下命令
metag.sh 自动激活已安装的 metag 环境；六个数据集的 metadata/raw 已准备好。
当前共享参考仍在下载/建库，正式库审计和真实数据分类试跑尚未完成。
准备完成后，完整复制第 1–3 步的整个括号代码块；任一步失败就停止。
第一步会检查参考，未完成时停止；不再重复启动正在运行的建库任务。

(
set -e

# 1. 选择六个已下载的数据集，并检查共享参考是否准备完成。
DATASETS=emp,wbe,长江,厦门污水,全球污水,纳米布沙漠
./metag.sh datasets "$DATASETS"
REF=/mnt/f/metag/reference
DB="$REF/kraken2_broad_v1"
MAP="$REF/virus_map_v1/virus_genome_types.sqlite"
if [[ ! -f "$DB/metag_build_state/all.complete" || ! -s "$DB/hash.k2d" || ! -s "$MAP" ]]; then
  echo '共享参考尚未准备完成，请等待；不要重复启动 build-db。' >&2
  exit 1
fi
./metag.sh db-audit --db "$DB" --virus-map "$MAP" --threads 16

# 2. 六个数据集逐个试跑；每个最多 2 个 run，每输入流最多 100000 records。
./metag.sh pilot "$DATASETS" --db "$DB" --virus-map "$MAP" --threads 16 --execute

# 3. 试跑全部成功后，全量处理这六个数据集。
./metag.sh full "$DATASETS" --db "$DB" --virus-map "$MAP" --threads 16 --execute
)

以下单独留到 ocean 和 sub 均下载完成后执行；先完成上面的共享参考准备。
(
set -e

# 4. 两个目录均下载完成后，重建清单并整理 raw；旧 metadata 会先备份。
./metag.sh prepare ocean,sub --refresh --execute
./metag.sh pilot ocean,sub --execute
./metag.sh full ocean,sub --execute
)

目录默认在 /mnt/f/metag；可用 --data-root 修改（共享数据库路径仍由 --db/--virus-map 指定）。
数据库默认 /mnt/f/metag/reference/kraken2_broad_v1。
病毒映射默认 /mnt/f/metag/reference/virus_map_v1/virus_genome_types.sqlite。
结果分别写到每个数据集的 results_v2/；试跑写到各自的 pilot_virus_v2/。
批量命令中某个数据集失败，会继续其余已选项，最后列出失败项并返回非零状态。
ocean 和 sub 仍在下载，前三步只处理已列出的六个数据集。
neon 有非 ENA 文件名，需核对文件归属与配对关系后再加入处理列表。
参考下载和映射构建成功后复用，不必重复构建。
若按 Ctrl+Z 暂停了建库，它仍占锁；在原终端用 jobs -l 查看，fg 恢复，不要重复启动或删除锁文件。

manifest = 建立 FASTQ 清单：样本/run/单双端/ENA 大小与 MD5/下载状态。
prepare = manifest（已有清单默认保留）+ organize（raw 软链接）。
valid=yes 仅表示清单检查通过；pilot/full 还会完整读取文件做 MD5 校验。
full = 审计、清单检查、MD5、QC、分类、分流、汇总；不带 --execute 只打印计划。
它处理本地清单，不会自动下载尚未获取的样本。

仍可使用单目录形式：./metag.sh full --root DIR --db DIR --virus-map FILE。
用 ./metag.sh prepare --help、./metag.sh pilot --help 查看批量参数。
批量 full 参数可用 ./metag.sh full emp --help 查看。
METAG_PYTHON / PYTHON_BIN 指定 Python；默认使用 metag 环境的 python3。
详细说明见 README.md。
HELP
}

case "${1:-}" in
  ''|-h|--help|help) usage; exit 0;;
esac

# Activate inside this process, including PATH for fastp, Kraken2 and shell modules.
# shellcheck source=f/conda_env.sh
source "$METAG_ROOT/f/conda_env.sh"
metag_activate

case "$1" in
  build-db) shift; exec bash "$METAG_ROOT/f/build_db.sh" "$@";;
  datasets|prepare|pilot) module=datasets.py;;
  full)
    if [[ -n "${2:-}" && "${2:0:1}" != '-' ]]; then
      module=datasets.py
    else
      shift; exec bash "$METAG_ROOT/f/run_full.sh" "$@"
    fi
    ;;
  align) shift; exec bash "$METAG_ROOT/f/align_candidates.sh" "$@";;
  virus)
    shift
    if (($# == 0)); then set -- --help; fi
    module=virus_metadata.py
    ;;
  manifest|metadata-upgrade|organize|verify|db-audit|run|summary|--version)
    module=metag.py
    ;;
  *) printf 'Unknown command: %s\n' "$1" >&2; usage >&2; exit 2;;
esac

command -v "$METAG_PYTHON" >/dev/null 2>&1 || {
  printf 'Python not found: %s. Check the metag environment or METAG_PYTHON setting.\n' "$METAG_PYTHON" >&2
  exit 127
}
exec "$METAG_PYTHON" "$METAG_ROOT/f/$module" "$@"
