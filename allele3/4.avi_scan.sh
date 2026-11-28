#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
ALLELE3: summarize existing AlphaGenome GRCh38 VCFs.

Usage:
  ./4.avi_scan.sh run|MODULE [options]

Modules:
  run | scan                       Scan VCFs into AVI windows and hotspots.
  merge | plan                     Merge completed shards or inspect the plan.

Runtime: analysis/allele3/4.avi_scan/{jobs,logs,qc,tmp}.
  Permanent products: analysis/allele3/4.avi_scan/.
  Preserve shards, results, reports, figures and metadata together.
  Explicit --out takes precedence.

Options:
  --dir-out-avi_scan DIR           AVI store [${ALLELE_ANALYSIS_ROOT}/4.avi_scan].
  --vcf-dir DIR                    Existing Atlas VCF directory.
  --out DIR                        Output directory.
  --jobs N                         Concurrent workers [4].
  --grch 37|38                     Genome build [38].
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME --cores N           Queue and CPU budget from shared 0f settings.
  LSF controls HPC job duration and memory.
  --replace TRUE|FALSE             Recompute or reuse completed outputs [FALSE].
  --plan-only                      Write commands without running or submitting.
  --help-full                      Show all engine options.
  -h, --help                       Show this help.

Permanent paths above are defaults under the configured data root (HPC:
  /work/sph-huangj/data). CLI directory options override environment defaults;
  explicit --out overrides the selected task output. Relative CLI paths use cwd.
  --replace FALSE (default) checks completion receipts and required files:
  complete outputs are reused; missing/incomplete scan products are generated
  at the configured destination. --replace TRUE explicitly requests rebuilding.

Examples:
  ./4.avi_scan.sh run --grch 38 --plan-only
  ./4.avi_scan.sh run --grch 38 --jobs 4
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- scan "$@"; fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=avi-scan
source "$ROOT/f/0.run.sh"
action=scan
if (($#)) && [[ $1 != -* ]]; then action=$1;shift;fi
out=${ALLELE_AVI_STORE:-${ALLELE_ANALYSIS_ROOT}/4.avi_scan}
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
