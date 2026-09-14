#!/usr/bin/env bash
set -euo pipefail
F=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ACTION=${1:-run}
[[ ${GU_BUILD:-37} == 37 ]] || { echo 'ERROR: GWAS PhyML sequences must use GRCh37' >&2; exit 2; }
args=(--lead-table "${GU_PHYML_LEAD_TABLE:?original COJO lead table required}"
      --loci "${GU_LOCI_CORE_FILE:?}" --vcf-dir "${GU_TARGET_VCF_DIR:?}"
      --archaic-root "${GU_ARCHAIC_ROOT:?}" --sample-file "${GU_SAMPLE_PANEL:?}"
      --out "${PHYML_OUT:?}" --dataset "${GU_TARGET:-1kg}" --action "$ACTION"
      --plot-phy "${PHYML_PLOT_PHY:-TRUE}" --timeout "${PHYML_TREE_TIMEOUT:-86400}" --cpus "${PHYML_TREE_CPUS:-4}")
[[ ${GU_CHRX_MALE_ONLY:-0} != 1 ]] || args+=(--x-male-only)
[[ ${GU_CHRX_PAR_DIPLOID:-0} != 1 ]] || args+=(--x-par-diploid)
exec python3 "$F/phyml_gwas.py" "${args[@]}"
