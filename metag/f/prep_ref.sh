#!/usr/bin/env bash
# Current-machine post-download preparation; full analysis is a separate user command.
case "${1:-}" in
  -h|--help)
    echo 'Usage: ./metag.sh prep_ref'
    echo '先完成 README 中四份 aria2 下载。本脚本恢复建库、审计并对六组数据试跑，不运行 full。'
    echo '每项任务限制为 36 GiB 内存和 4 GiB swap；失败立即停止。'
    exit 0;;
  '') ;;
  *) echo '不接受参数；运行 --help 查看说明。' >&2; exit 2;;
esac
set -euo pipefail

# 1. 选择六个已下载的数据集，设置路径和本机内存上限。
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
METAG="$HERE/../metag.sh"
DATASETS=emp,wbe,长江,厦门污水,全球污水,纳米布沙漠
REF=/mnt/f/metag/ref_kraken2_ictv
DB="$REF/kraken2_broad_v1"
MAP="$REF/virus_map_v1/virus_genome_types.sqlite"
UNIT=metag-reference-build.scope
LOG="$DB/metag_build_logs/supervised-build.log"
LIMITS=(--property=MemoryHigh=30G --property=MemoryMax=36G
        --property=MemorySwapMax=4G --property=OOMPolicy=stop
        --property=KillMode=control-group)
[[ -d "$DB" && -s "$MAP" ]] || { echo '参考目录或已准备的病毒映射不存在，请先检查共享参考路径。' >&2; exit 1; }
exec 8>"$REF/.prepare_reference.lock"
flock -n 8 || { echo '另一个 metag.sh prep_ref 正在执行，请等待它完成。' >&2; exit 1; }
"$METAG" datasets "$DATASETS"
systemctl --user show --property=Version --value >/dev/null

run_metag() {
  systemd-run --user --scope --quiet --collect \
    --unit="metag-$1-$(date +%s%N)" "${LIMITS[@]}" "$METAG" "$@"
}

# 2. 前台恢复建库；不能接管历史上独立启动的后台服务。
case "$(systemctl --user show metag-reference-build.service --property=ActiveState --value)" in
  active|activating|deactivating)
    echo '历史后台建库服务仍在运行；请先明确停止该服务，再以前台方式运行 prep_ref。' >&2
    exit 1 ;;
esac
if [[ ! -f "$DB/metag_build_state/all.complete" ]]; then
  mkdir -p "$DB/metag_build_logs"
  systemctl --user reset-failed "$UNIT" 2>/dev/null || true
  printf '前台恢复参考建库；日志：%s\n' "$LOG"
  systemd-run --user --scope --quiet --collect --unit="$UNIT" "${LIMITS[@]}" \
    "$METAG" build-db --db "$DB" --threads 16 --download-threads 8 \
    --max-db-gb 24 --execute 2>&1 | tee -a "$LOG"
fi
if [[ ! -f "$DB/metag_build_state/all.complete" || ! -s "$DB/hash.k2d" ||
      ! -s "$DB/opts.k2d" || ! -s "$DB/taxo.k2d" || ! -s "$MAP" ]]; then
  printf '参考尚未完成，已停止。请查看日志：%s\n' "$LOG" >&2
  systemctl --user show "$UNIT" --property=Result,ExecMainStatus >&2
  exit 1
fi

# 3. 检查正式索引及病毒 DNA/RNA 覆盖。
run_metag db-audit --db "$DB" --virus-map "$MAP" --threads 16

# 4. 六组真实试跑：每组最多 2 个 run，每输入流最多 100000 records；不重复检查 FASTQ MD5。
run_metag pilot "$DATASETS" --db "$DB" --virus-map "$MAP" --threads 16 \
  --limit 2 --max-records 100000 --execute

printf '\n参考建库、覆盖审计和六组试跑全部成功。接下来按 README 运行 metag.sh full。\n'
