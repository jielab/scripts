#!/usr/bin/env bash
# No disk-wide loops: one explicit dataset, audited database, frozen virus map.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT=/mnt/f/metag/ocean; DB=""; VIRUS_MAP=""; OUT=""; THREADS=16; EXECUTE=0
usage() { cat <<'EOF'
批量形式：./metag.sh full emp,wbe [--execute] [--db DIR] [--virus-map FILE]
单目录形式：./metag.sh full --db DIR --virus-map FILE [--root DIR] [--out DIR] [--threads 16] [--execute]
默认 ROOT=/mnt/f/metag/ocean，OUT=ROOT/results_v2。
不带 --execute 只打印计划；带 --execute 依次审计、检查清单、完整 MD5 校验、分类、汇总。
执行前须准备好环境、数据库、病毒映射，并完成 manifest 和 organize。
清单有 valid=no（缺文件、.aria2 等）时先完成下载并更新清单。
首次使用请运行 ./metag.sh --help，按其中步骤准备共享参考并试跑。

EOF
}
while (($#)); do
  case "$1" in
    --root) ROOT="$2"; shift 2;; --db) DB="$2"; shift 2;; --virus-map) VIRUS_MAP="$2"; shift 2;;
    --out) OUT="$2"; shift 2;; --threads) THREADS="$2"; shift 2;; --execute) EXECUTE=1; shift;;
    -h|--help) usage; exit 0;; *) usage >&2; exit 2;;
  esac
done
[[ -n "$DB" && -n "$VIRUS_MAP" && "$THREADS" =~ ^[1-9][0-9]*$ ]] || { usage >&2; exit 2; }
ROOT="$(realpath -m -- "$ROOT")"; DB="$(realpath -m -- "$DB")"; VIRUS_MAP="$(realpath -m -- "$VIRUS_MAP")"
OUT="${OUT:-$ROOT/results_v2}"
step() {
  printf '%q ' "$@"; printf '\n'
  if ((EXECUTE)); then "$@"; fi
}
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" db-audit --db "$DB" --virus-map "$VIRUS_MAP" --threads "$THREADS"
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" run --root "$ROOT" --db "$DB" --virus-map "$VIRUS_MAP" --out "$OUT" --threads "$THREADS" --dry-run
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" verify --root "$ROOT" --md5
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" run --root "$ROOT" --db "$DB" --virus-map "$VIRUS_MAP" --out "$OUT" --threads "$THREADS"
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" summary --out "$OUT"
