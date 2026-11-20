#!/usr/bin/env bash
# One-pass chromosome-sharded summary of the existing merged GRCh38 Atlas VCFs.
if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    cat <<'HELP'
./4.avi_scan.sh --grch 38  # 完整 AVI 扫描
./4.avi_scan.sh --grch 38 MODULE  # scan|merge|plan
./4.avi_scan.sh --grch 38 --jobs 4 --replace FALSE

Usage: ./4.avi_scan.sh --grch 38 [OPTIONS]
No action: scan existing GRCh38 AlphaGenome VCFs into AVI windows/hotspots.
  --jobs N [4] --vcf-dir DIR --out DIR
  --replace FALSE|TRUE [FALSE]
  --plan-only           Inspect without submitting
  --queue NAME --job-hours N
  --help-full           All engine options (including merge/plan)
HPC is detected automatically. Job name: allele3.avi_scan
HELP
    exit 0
fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=avi-scan
source "$ROOT/f/0.run.sh"
action=scan
if (($#)) && [[ $1 != -* ]]; then action=$1;shift;fi
out=${ALLELE_AVI_STORE:-$ALLELE_ANALYSIS_ROOT/avi_scan}
python=${ALLELE_PYTHON:-};memory=${ALLELE_MEMORY_LIMIT_GB:-16};interval=300
extra=();have_dir=0
while (($#)); do
    case "$1" in
        --cores) shift 2;;
        --out) out=$2;shift 2;;
        --python) python=$2;shift 2;;
        --memory-limit-gb) memory=$2;shift 2;;
        --screen-log-interval) interval=$2;shift 2;;
        --vcf|--vcf-dir) have_dir=1;extra+=("$1" "$2");shift 2;;
        --help-full) exec "${python:-python3}" "$ROOT/f/avi_scan.py" --help;;
        *) extra+=("$1");shift;;
    esac
done
((have_dir)) || extra+=(--vcf-dir "$ALLELE_ATLAS_DIR")
python=$(workflow_python "$python" "$HOME/anaconda3/envs/allele3-evolution/bin/python")
export WORKFLOW_LOG_DIR="$out/logs"
workflow_launch avi-scan "$python" "$ROOT/f/avi_scan.py" "$out" "$memory" "$interval" 0 0 1 \
    "$action" --out "$out" --replace "${ALLELE_REPLACE:-FALSE}" "${extra[@]}"
