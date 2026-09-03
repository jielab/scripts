#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
F=$ROOT/f
ACTION=${1:-run}
case "$ACTION" in run|check|match) ;; *) echo "ERROR: invalid phyml action: $ACTION" >&2; exit 2;; esac

plot_phy=0
if [[ $ACTION == run ]]; then
  case "${PHYML_PLOT_PHY:-FALSE}" in
    TRUE) plot_phy=1;;
    FALSE) ;;
    *) echo "ERROR: PHYML_PLOT_PHY must be TRUE or FALSE" >&2; exit 2;;
  esac
fi
case "${PHYML_REPLACE:-FALSE}" in TRUE|FALSE) ;; *) echo "ERROR: internal PHYML_REPLACE must be TRUE or FALSE" >&2; exit 2;; esac

for cmd in python3 bcftools bash awk; do command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing command: $cmd" >&2; exit 1; }; done
if (( plot_phy )) && ! command -v phyml >/dev/null 2>&1; then
  echo "ERROR: phyml is required for tree inference; use 'phyml --plot-phy FALSE' for direct haplotype comparison only" >&2
  exit 1
fi
if (( plot_phy )); then
  command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript is required to render the PhyML tree PNG" >&2; exit 1; }
  Rscript -e 'quit(status=if(requireNamespace("ape", quietly=TRUE) && capabilities("png")) 0 else 1)' >/dev/null 2>&1 || {
    echo "ERROR: R package 'ape' and PNG graphics support are required to render the PhyML tree" >&2
    exit 1
  }
fi
PHE_F=${PHE_F:-$ROOT/../0f/0phe.f.sh}
[[ -s $PHE_F ]] || { echo "ERROR: shared helper missing: $PHE_F" >&2; exit 1; }
# shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
source "$PHE_F"
declare -F match_HAP >/dev/null || { echo "ERROR: $PHE_F does not define match_HAP" >&2; exit 1; }

OUT=${PHYML_OUT:-${GU_ANALYSIS_ROOT:-/mnt/d/analysis/gu}/phyml/${GU_SCOPE_ID:-default}}
phyml_input_paths(){
  local flat=$OUT/loci/haplotypes.phy
  if [[ -e $flat ]]; then
    printf '%s\n' "$flat"
  elif [[ -d $OUT/loci ]]; then
    find "$OUT/loci" -mindepth 2 -maxdepth 2 -type f -name 'haplotypes.phy' | sort
  fi
}
args=(--vcf-dir "${GU_TARGET_VCF_DIR:?}" --archaic-root "${GU_ARCHAIC_ROOT:?}" --out "$OUT" --match-helper "$PHE_F"
      --build "${GU_BUILD:-37}"
      --chr-window-bp "${PHYML_CHR_WINDOW_BP:-500000}"
      --jobs "${PHYML_JOBS:-1}"
      --refs "${PHYML_REFS:-Altai Chagyr Vindija Denisova}" --min-mac "${PHYML_MIN_MAC:-2}"
      --min-compared "${PHYML_MIN_COMPARED:-2}" --min-prop "${PHYML_MIN_PROP:-0.8}"
      --ld-r2 "${PHYML_LD_R2:-0.98}" --min-ld-sites "${PHYML_MIN_LD_SITES:-2}"
      --region-mode "${PHYML_REGION_MODE:-core}")
[[ -n ${GU_SAMPLE_PANEL:-} && -s ${GU_SAMPLE_PANEL:-} ]] && args+=(--sample-file "$GU_SAMPLE_PANEL")
[[ ${GU_CHRX_MALE_ONLY:-0} == 1 ]] && args+=(--x-male-only)
[[ ${GU_CHRX_PAR_DIPLOID:-0} == 1 ]] && args+=(--x-par-diploid)
[[ ${PHYML_REQUIRE_ALL_ARCHAIC:-0} == 1 ]] && args+=(--require-all-archaic)
[[ ${GU_OUTPUT_LAYOUT:-} == per_locus ]] && args+=(--flat-locus-dir)
if [[ -n ${GU_LOCI_FILE:-} ]]; then
  args+=(--loci "$GU_LOCI_FILE")
  [[ -n ${GU_LOCI_MAP_FILE:-} && -s ${GU_LOCI_MAP_FILE:-} ]] && args+=(--loci-map "$GU_LOCI_MAP_FILE")
else
  args+=(--chrs "${GU_CHRS:?}")
fi

if [[ $ACTION == check ]]; then exec python3 "$F/phyml.py" "${args[@]}" --check; fi

comparison_record(){
  printf 'argv'
  printf '\t%q' "${args[@]}"
  printf '\n'
  stat -c 'script\t%n:%s:%Y' "$F/phyml.py" "$F/phyml_tree_summary.py" "$F/phyml_tree_plot.R" "$PHE_F"
  [[ -z ${GU_LOCI_FILE:-} ]] || stat -c 'input\t%n:%s:%Y' "$GU_LOCI_FILE"
  [[ -z ${GU_LOCI_MAP_FILE:-} ]] || stat -c 'input\t%n:%s:%Y' "$GU_LOCI_MAP_FILE"
  [[ -z ${GU_SAMPLE_PANEL:-} ]] || stat -c 'input\t%n:%s:%Y' "$GU_SAMPLE_PANEL"
  find "${GU_TARGET_VCF_DIR:?}" "${GU_ARCHAIC_ROOT:?}" -type f \( -name '*.vcf.gz' -o -name '*.bcf' \) -printf 'input\t%p:%s:%T@\n' 2>/dev/null | sort
}
comparison_meta=$OUT/final/comparison.input.meta.tsv
expected_meta=$(mktemp "${TMPDIR:-/tmp}/gu.phyml.meta.XXXXXX")
trap 'rm -f -- "$expected_meta"' EXIT
comparison_record > "$expected_meta"
comparison_ready=0
if [[ ${PHYML_REPLACE:-FALSE} == TRUE ]]; then
  echo "[GU PHYML] comparison=OVERWRITE reason=--replace-phyml output=$OUT"
elif [[ -s $OUT/final/haplotypes.tsv && -s $OUT/final/loci.tsv && -s $OUT/final/haplotype_samples.tsv && -s $OUT/final/skipped_loci.tsv && -s $comparison_meta ]] && cmp -s "$expected_meta" "$comparison_meta"; then
  comparison_ready=1
fi
if (( comparison_ready )); then
  echo "[GU PHYML] comparison=SKIP reason=outputs_exist output=$OUT"
else
  echo "[GU PHYML] comparison=RUN reason=missing_or_stale_outputs output=$OUT"
  python3 "$F/phyml.py" "${args[@]}"
  mkdir -p "$OUT/final"
  mv -f -- "$expected_meta" "$comparison_meta"
  expected_meta=''
  find "$OUT/loci" -type f \( -name '*_phyml_tree.txt' -o -name '*_phyml_stats.txt' -o -name '*_phyml_tree.png' -o -name '*.phyml.log' \) -delete 2>/dev/null || true
fi

if (( plot_phy )); then
  echo "[GU PHYML] plot_phy=TRUE tree_inference=RUN bootstrap=${PHYML_BOOT:-100}"
  status=0
  while IFS= read -r phy; do
    [[ -s ${phy}_phyml_tree.txt ]] && continue
    phyml -i "$phy" -m HKY85 -c 4 -a e -v e -b "${PHYML_BOOT:-100}" > "${phy}.phyml.log" 2>&1 || status=1
  done < <(phyml_input_paths)
  (( status == 0 )) || { echo "ERROR: one or more PhyML jobs failed" >&2; exit 1; }
  python3 "$F/phyml_tree_summary.py" --out "$OUT"
elif [[ $ACTION == run ]]; then
  echo "[GU PHYML] plot_phy=FALSE tree_inference=SKIP; use --plot-phy TRUE to build trees"
  python3 "$F/phyml_tree_summary.py" --out "$OUT" --default-status not_requested
fi
