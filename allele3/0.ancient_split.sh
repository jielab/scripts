#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    cat <<'HELP'
./0.ancient_split.sh --grch 38  # pfile 拆分和 GRCh37 -> GRCh38
./0.ancient_split.sh --grch 38 --chromosomes 1,2  # 默认 auto，选择部分染色体
./0.ancient_split.sh --grch 38 --jobs 4 --threads 4 --replace FALSE

Usage: ./0.ancient_split.sh --grch 38 [--submit hpc] [--plan-only] [OPTIONS]
Split ancient/37/ancient.p* into 37/chr*.pgen/.pvar/.psam,
then lift to 38/chr*.p* using the UKB row-preserving policy.
Only POS changes; rejected rows keep POS=-1. PGEN/PSAM use relative links.
Also provides 38/ancient.p* in original global row order for existing PCA.

--input-prefix PREFIX  --out37 DIR  --out38 DIR
--runtime DIR   logs/jobs only [ANALYSIS_ROOT/allele3/ancient_split]
--scratch-root DIR [/tmp on compute node]; intermediates are cleaned on success
--chain FILE --fasta FILE --chromosomes auto|1,2,...
--jobs N [4] --threads N [4] --memory-mb N [8192 per PLINK process]
--queue NAME [33] --job-hours N [120]
--replace TRUE|FALSE [FALSE] --plan-only
Receipts and rejection tables: ancient/{37,38}/.ancient_split/.
One LSF job runs up to --jobs chromosomes in parallel. Requires LSF allocation.
HELP
    exit 0
fi
export submit=hpc SUBMIT=hpc
source "${ALLELE_OS_CONFIG:-$ROOT/../0f/0os.sh}"
export PATH="${DIR_WORK}/.conda/envs/allele3-evolution/bin:${DIR_SOFTWARE}/bin:$PATH"
exec "${ALLELE_PYTHON:-${DIR_WORK}/.conda/envs/allele3-evolution/bin/python3}" "$ROOT/f/ancient_split.py" "$@"
