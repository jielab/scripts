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
  pilot             每数据集选至多 2 个有效 run 试跑；--execute 执行
  full              按数据集顺序审计、QC、分类、汇总；--execute 执行
  prep_ref          下载后恢复建库、审计及六组试跑；不会启动 full
  build-db          下载参考序列并建库；默认计划，--execute 执行
  reference-urls    从已有参考清单导出多份 aria2 .url 文件，不启动下载
  virus             获取 ICTV 数据、构建病毒 DNA/RNA/unknown 映射
  db-audit          核查数据库基础类别和病毒 DNA/RNA 覆盖
  manifest          单数据集的 ENA 文件清单与下载状态检查
  organize / verify 单数据集 raw 整理 / MD5、gzip 校验
  run / summary     单数据集分类分流 / 汇总
  metadata-upgrade  保留人工值并补充实验字段
  align             候选序列的独立竞争比对

Usage examples：下载 → prep_ref → 正式分析
先完成 README 中四份下载：ref_microbes.url、ref_plants.url、ref_vertebrates.url、ref_invertebrates.url。
共享参考目录为 /mnt/f/metag/ref_kraken2_ictv。本机环境已安装，无需再运行 install.sh。
prep_ref 内部完成建库、审计和六组试跑，不会启动 full。
当前微生物精简模式按候选物种各保留一个基因组；下载清单与建库选择一致。
下面两条按顺序分别执行；第一条失败时先排查，不执行第二条。
metag.sh 自动激活环境；pilot/full 不重复校验已由用户确认通过 MD5 的原始 FASTQ。

# 1. 四份 aria2 下载全部成功后：恢复建库、审计并完成六组真实试跑。
/mnt/d/scripts/metag/metag.sh prep_ref

# 2. 上一条成功后，由你启动六组正式全量分析。
systemd-run --user --wait --pipe --collect --unit="metag-full-$(date +%s%N)" \
  --property=MemoryHigh=30G --property=MemoryMax=36G \
  --property=MemorySwapMax=4G --property=OOMPolicy=stop --property=KillMode=control-group \
  /mnt/d/scripts/metag/metag.sh full emp,wbe,长江,厦门污水,全球污水,纳米布沙漠 \
  --threads 16 --execute

prep_ref 与全量命令限制为 36 GiB 内存和 4 GiB swap；其他程序不受此限制。
正式参考索引、覆盖审计及真实试跑仍待完成，自动检查保持暂停。
建库进度：
  systemctl --user status metag-reference-build.service --no-pager
  tail -n 20 /mnt/f/metag/ref_kraken2_ictv/kraken2_broad_v1/metag_build_logs/supervised-build.log

目录默认在 /mnt/f/metag；可用 --data-root 修改（共享数据库路径仍由 --db/--virus-map 指定）。
数据库默认 /mnt/f/metag/ref_kraken2_ictv/kraken2_broad_v1。
病毒映射默认 /mnt/f/metag/ref_kraken2_ictv/virus_map_v1/virus_genome_types.sqlite。
结果分别写到每个数据集的 results_v2/；试跑写到各自的 pilot_virus_v2/。
批量命令中某个数据集失败，会继续其余已选项，最后列出失败项并返回非零状态。
ocean 和 sub 仍在下载，上述命令只处理已列出的六个数据集。
neon 有非 ENA 文件名，需核对文件归属与配对关系后再加入处理列表。
参考下载和映射构建成功后复用，不必重复构建。
若按 Ctrl+Z 暂停了建库，它仍占锁；在原终端用 jobs -l 查看，fg 恢复，不要重复启动或删除锁文件。

manifest = 建立 FASTQ 清单：样本/run/单双端/ENA 大小与 MD5/下载状态。
prepare = manifest（已有清单默认保留）+ organize（raw 软链接）。
valid=yes 仅表示清单检查通过，不代表 MD5 校验通过。
pilot/full 默认不重复检查 FASTQ MD5；如需重查，添加 --check-md5。
full = 审计、清单检查、QC、分类、分流、汇总；不带 --execute 只打印计划。
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
  prep_ref) shift; exec bash "$METAG_ROOT/f/prep_ref.sh" "$@";;
  build-db) shift; exec bash "$METAG_ROOT/f/build_db.sh" "$@";;
  datasets|prepare|pilot) module=datasets.py;;
  reference-urls) shift; module=reference_urls.py;;
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
