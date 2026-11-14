#!/usr/bin/env bash
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case ${1:-} in
        -h|--help|help)
            cat <<'ALLELE_HELP'
Usage: bash f/scan_annotations.sh [OUT_DIR]
下载并校验 UCSC hg38 rmskJoinedCurrent.txt.gz；已有完整文件则复用。
默认目录：$DIR_ANNOT/allele3/repeats/hg38；未配置 DIR_ANNOT 时必须指定 OUT_DIR。
示例：bash f/scan_annotations.sh /work/sph-huangj/data/annot/allele3/repeats/hg38
需要 curl、flock、gzip、sha256sum 和网络；输出下载文件的绝对路径。
-h / --help 仅显示帮助，不创建文件或下载。
ALLELE_HELP
            exit 0 ;;
    esac
fi
# Shared all-repeat joined annotation. Atomic publication; safe concurrent reuse.
set -euo pipefail
OUT=${1:-${DIR_ANNOT}/allele3/repeats/hg38}
mkdir -p "$OUT"
exec 9>"$OUT/.joined-download.lock"
flock -x 9
FILE="$OUT/rmskJoinedCurrent.txt.gz"
URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmskJoinedCurrent.txt.gz
if [[ -s $FILE ]]; then gzip -t "$FILE"; printf '%s\n' "$FILE"; exit 0; fi
# Each attempt re-detects the current length, retaining progress after truncation.
for attempt in {1..9}; do
    if curl --fail --location --connect-timeout 30 --speed-time 60 --speed-limit 1024 --continue-at - -o "$FILE.partial" "$URL"; then break; fi
    if (( attempt == 9 )); then exit 1; fi
    sleep 2
done
gzip -t "$FILE.partial"
mv "$FILE.partial" "$FILE"
sha256sum "$FILE" > "$FILE.sha256"
printf '%s\n' "$FILE"
