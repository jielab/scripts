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
case "${PHYML_RUN_ALL_UNIQUE_TREE:-FALSE}" in TRUE|FALSE) ;; *) echo "ERROR: PHYML_RUN_ALL_UNIQUE_TREE must be TRUE or FALSE" >&2; exit 2;; esac
case "${PHYML_TREE_FAILURE_POLICY:-warn}" in warn|error) ;; *) echo "ERROR: PHYML_TREE_FAILURE_POLICY must be warn or error" >&2; exit 2;; esac

for cmd in python3 bcftools bash awk find; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing command: $cmd" >&2; exit 1; }
done
for required in phyml.py phyml/core.py phyml/tree_summary.py phyml/region_scan.py phyml/evidence.py phyml/anchor.py phyml/layered_plot.py phyml/tree_plot.R; do
  [[ -s $F/$required ]] || { echo "ERROR: PhyML module missing: $F/$required" >&2; exit 1; }
done
if (( plot_phy )) && ! command -v phyml >/dev/null 2>&1; then
  echo "ERROR: phyml is required for --plot-phy TRUE" >&2
  exit 1
fi
if (( plot_phy )); then
  python3 - <<'PY' >/dev/null 2>&1 || {
import matplotlib
PY
    echo "ERROR: matplotlib is required to render layered PhyML QC plots" >&2
    exit 1
  }
fi

PHE_F=${PHE_F:-$ROOT/../0f/0phe.f.sh}
[[ -s $PHE_F ]] || { echo "ERROR: shared helper missing: $PHE_F" >&2; exit 1; }
# shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
source "$PHE_F"
declare -F match_HAP >/dev/null || { echo "ERROR: $PHE_F does not define match_HAP" >&2; exit 1; }

OUT=${PHYML_OUT:-${GU_ANALYSIS_ROOT:-/mnt/d/analysis/gu}/phyml/${GU_SCOPE_ID:-default}}

# Use every locally available high-coverage archaic reference unless the caller
# explicitly supplies PHYML_REFS.  The earlier GU default omitted Denisova25,
# although the original baldness analysis and manuscript used two Denisovans.
# Auto-detection is directory-based and the downstream phyml.py preflight still
# verifies that the requested chromosome VCF exists for every selected ref.
resolve_phyml_refs(){
  local configured=${PHYML_REFS:-} root=${GU_ARCHAIC_ROOT:?} ref d first
  local -a selected=() known=(Altai Chagyr Vindija Denisova Denisova25)
  if [[ -n $configured ]]; then
    printf '%s\n' "$configured"
    return 0
  fi
  local chr chr_list=${GU_CHRS:-}
  chr_list=${chr_list//,/ }
  for ref in "${known[@]}"; do
    for d in "$root/$ref" "$root/avcf/$ref"; do
      [[ -d $d ]] || continue
      first=$(find "$d" -maxdepth 1 -type f \( -name '*.vcf.gz' -o -name '*.bcf' \) -size +0c -print -quit 2>/dev/null || true)
      [[ -n $first ]] || continue
      # Do not auto-select a partially installed reference panel.  Explicit
      # PHYML_REFS remains available when the caller intentionally wants one.
      local complete=1
      for chr in $chr_list; do
        chr=${chr#chr}; chr=${chr#CHR}; chr=${chr^^}; [[ $chr == 23 ]] && chr=X
        first=$(find "$d" -maxdepth 1 -type f -size +0c \
          \( -name "*chr${chr}_*.vcf.gz" -o -name "*chr${chr}.*.vcf.gz" \
             -o -name "*chr${chr}.vcf.gz" -o -name "*${chr}.vcf.gz" \
             -o -name "*chr${chr}_*.bcf" -o -name "*chr${chr}.*.bcf" \
             -o -name "*chr${chr}.bcf" -o -name "*${chr}.bcf" \) \
          -print -quit 2>/dev/null || true)
        if [[ -z $first ]]; then
          complete=0
          break
        fi
      done
      if (( complete )); then
        selected+=("$ref")
        break
      fi
    done
  done
  if (( ${#selected[@]} == 0 )); then
    # Preserve the historical four-reference contract so phyml.py emits the
    # usual precise missing-file diagnostic rather than an empty-ref error.
    printf '%s\n' 'Altai Chagyr Vindija Denisova'
  else
    printf '%s\n' "${selected[*]}"
  fi
}
PHYML_REFS_EFFECTIVE=$(resolve_phyml_refs)
if [[ -n ${PHYML_REFS:-} ]]; then
  PHYML_REFS_SOURCE=environment_override
else
  PHYML_REFS_SOURCE=auto_detect_local_reference_directories
fi
echo "[GU PHYML] archaic_refs=$PHYML_REFS_EFFECTIVE source=$PHYML_REFS_SOURCE"

# BED loci are always scanned across the complete core before any anchor-based
# zoom.  A requested LD mode is retained as provenance but is not allowed to
# erase regional matches when BED column 4 is absent or incorrect.
requested_region_mode=${PHYML_REGION_MODE:-core}
base_region_mode=$requested_region_mode
if [[ -n ${GU_LOCI_FILE:-} ]]; then
  base_region_mode=core
  if [[ $requested_region_mode != core ]]; then
    echo "[GU PHYML] region_mode=core_for_stage1 requested=$requested_region_mode reason=anchor_independent_complete_BED_scan"
  fi
fi

raw_phy_inputs(){
  local flat=$OUT/loci/haplotypes.phy
  if [[ -e $flat ]]; then
    printf '%s\n' "$flat"
  elif [[ -d $OUT/loci ]]; then
    find "$OUT/loci" -mindepth 2 -maxdepth 2 -type f -name 'haplotypes.phy' -size +0c | sort
  fi
}
region_common_phy_inputs(){
  [[ -d $OUT/loci ]] || return 0
  find "$OUT/loci" -type f -name 'haplotypes.region_common*.phy' -size +0c | sort
}
evidence_phy_inputs(){
  [[ -d $OUT/loci ]] || return 0
  find "$OUT/loci" -type f -name 'haplotypes.evidence.*.phy' -size +0c | sort
}

args=(
  --vcf-dir "${GU_TARGET_VCF_DIR:?}"
  --archaic-root "${GU_ARCHAIC_ROOT:?}"
  --out "$OUT"
  --match-helper "$PHE_F"
  --build "${GU_BUILD:-37}"
  --chr-window-bp "${PHYML_CHR_WINDOW_BP:-500000}"
  --jobs "${PHYML_JOBS:-1}"
  --refs "$PHYML_REFS_EFFECTIVE"
  --min-mac "${PHYML_MIN_MAC:-2}"
  --min-compared "${PHYML_MIN_COMPARED:-2}"
  --min-prop "${PHYML_MIN_PROP:-0.8}"
  --ld-r2 "${PHYML_LD_R2:-0.98}"
  --min-ld-sites "${PHYML_MIN_LD_SITES:-2}"
  --region-mode "$base_region_mode"
)
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

region_scan_args=(
  --out "$OUT"
  --scan-min-callable-sites "${PHYML_REGION_SCAN_MIN_SITES:-12}"
  --scan-seed-prop "${PHYML_REGION_SCAN_MIN_PROP:-0.90}"
  --scan-max-gap-bp "${PHYML_REGION_SCAN_MAX_GAP_BP:-50000}"
  --full-min-callable-sites "${PHYML_REGION_FULL_MIN_SITES:-10}"
  --full-min-prop "${PHYML_REGION_FULL_MIN_PROP:-0.80}"
  --family-min-overlap "${PHYML_REGION_FAMILY_MIN_OVERLAP:-0.50}"
  --common-tree-min-copies "${PHYML_COMMON_TREE_MIN_COPIES:-11}"
  --common-tree-min-sites "${PHYML_COMMON_TREE_MIN_SITES:-10}"
  --min-tree-bootstrap "${PHYML_REGION_TREE_BOOTSTRAP_MIN:-70}"
  --min-tree-modern-types "${PHYML_REGION_TREE_MIN_MODERN_TYPES:-2}"
)

evidence_args=(
  --out "$OUT"
  --outgroup "${PHYML_OUTGROUP:-YRI}"
  --outgroup-fallback "${PHYML_OUTGROUP_FALLBACK:-AFR}"
  --max-outgroup-allele-freq "${PHYML_MAX_OUTGROUP_ALLELE_FREQ:-0.05}"
  --max-outgroup-haplotype-freq "${PHYML_MAX_OUTGROUP_HAPLOTYPE_FREQ:-0.05}"
  --min-outgroup-copies "${PHYML_MIN_OUTGROUP_COPIES:-20}"
  --min-diagnostic-sites "${PHYML_MIN_DIAGNOSTIC_SITES:-3}"
  --min-diagnostic-match-prop "${PHYML_MIN_DIAGNOSTIC_MATCH_PROP:-0.80}"
  --min-candidate-copies "${PHYML_MIN_CANDIDATE_COPIES:-2}"
  --max-diagnostic-gap-bp "${PHYML_MAX_DIAGNOSTIC_GAP_BP:-50000}"
  --min-neanderthal-refs "${PHYML_MIN_NEANDERTHAL_REFS:-2}"
  --min-other-lineage-refs "${PHYML_MIN_OTHER_LINEAGE_REFS:-1}"
  --evidence-controls "${PHYML_EVIDENCE_CONTROLS:-24}"
  --max-supporting-contexts "${PHYML_EVIDENCE_MAX_SUPPORTING_CONTEXTS:-24}"
  --min-control-copies "${PHYML_EVIDENCE_MIN_CONTROL_COPIES:-2}"
  --max-candidate-haplotype-types "${PHYML_MAX_CANDIDATE_HAPLOTYPE_TYPES:-1000}"
  --min-tree-sites "${PHYML_MIN_EVIDENCE_TREE_SITES:-10}"
  --tree-flank-bp "${PHYML_EVIDENCE_TREE_FLANK_BP:-5000}"
  --min-tree-bootstrap "${PHYML_EVIDENCE_BOOTSTRAP_MIN:-70}"
  --strong-tree-bootstrap "${PHYML_EVIDENCE_BOOTSTRAP_STRONG:-90}"
  --strong-min-diagnostic-sites "${PHYML_STRONG_MIN_DIAGNOSTIC_SITES:-10}"
  --strong-min-candidate-copies "${PHYML_STRONG_MIN_CANDIDATE_COPIES:-5}"
  --strong-min-archaic-tips "${PHYML_STRONG_MIN_ARCHAIC_TIPS:-2}"
  --min-candidate-purity "${PHYML_EVIDENCE_MIN_PURITY:-0.80}"
  --min-candidate-sensitivity "${PHYML_EVIDENCE_MIN_SENSITIVITY:-0.50}"
  --min-tree-candidate-types "${PHYML_EVIDENCE_MIN_CANDIDATE_TYPES:-1}"
)
EVIDENCE_SAMPLE_PANEL=${PHYML_EVIDENCE_SAMPLE_PANEL:-${GU_SAMPLE_PANEL:-}}
[[ -n $EVIDENCE_SAMPLE_PANEL && -s $EVIDENCE_SAMPLE_PANEL ]] && evidence_args+=(--sample-file "$EVIDENCE_SAMPLE_PANEL")

evidence_enabled=0
if [[ -n ${GU_LOCI_FILE:-} || ${PHYML_EVIDENCE_WHOLE_CHROMOSOME:-FALSE} == TRUE ]]; then
  evidence_enabled=1
fi

anchor_args=()
if [[ -n ${GU_LOCI_FILE:-} ]]; then
  anchor_args=(
    --loci "$GU_LOCI_FILE"
    --vcf-dir "${GU_TARGET_VCF_DIR:?}"
    --archaic-root "${GU_ARCHAIC_ROOT:?}"
    --out "$OUT"
    --refs "$PHYML_REFS_EFFECTIVE"
    --common-min-copies "${PHYML_COMMON_TREE_MIN_COPIES:-11}"
  )
  [[ -n ${GU_LOCI_MAP_FILE:-} && -s ${GU_LOCI_MAP_FILE:-} ]] && anchor_args+=(--loci-map "$GU_LOCI_MAP_FILE")
  [[ ${GU_CHRX_MALE_ONLY:-0} == 1 ]] && anchor_args+=(--x-male-only)
  [[ ${GU_CHRX_PAR_DIPLOID:-0} == 1 ]] && anchor_args+=(--x-par-diploid)
  [[ ${GU_OUTPUT_LAYOUT:-} == per_locus ]] && anchor_args+=(--flat-locus-dir)
fi

if [[ $ACTION == check ]]; then
  python3 "$F/phyml.py" region prepare --help >/dev/null
  python3 "$F/phyml.py" evidence prepare --help >/dev/null
  exec python3 "$F/phyml.py" compare "${args[@]}" --check
fi

mkdir -p "$OUT/final"
ref_selection_tmp=$OUT/final/.archaic_references.selected.tsv.tmp.$$
{
  printf 'reference_order\treference\tlineage\tselection_source\n'
  ref_order=0
  for ref_name in $PHYML_REFS_EFFECTIVE; do
    ref_order=$((ref_order+1))
    case "${ref_name,,}" in
      *altai*|*vindija*|*chagyr*|*neand*) ref_lineage=Neanderthal ;;
      *denis*) ref_lineage=Denisovan ;;
      *) ref_lineage=Archaic ;;
    esac
    printf '%s\t%s\t%s\t%s\n' "$ref_order" "$ref_name" "$ref_lineage" "$PHYML_REFS_SOURCE"
  done
} > "$ref_selection_tmp"
mv -f -- "$ref_selection_tmp" "$OUT/final/archaic_references.selected.tsv"

comparison_record(){
  local i c source_meta
  local -a record_args=("${args[@]}")
  for (( i=0; i<${#record_args[@]}; i++ )); do
    if [[ ${record_args[$i]} == --vcf-dir && $((i+1)) -lt ${#record_args[@]} ]]; then
      record_args[$((i+1))]='PREPARED_TARGET_VCF_DIR'
    fi
  done
  printf 'argv'
  printf '\t%q' "${record_args[@]}"
  printf '\n'
  printf 'requested_region_mode\t%s\n' "$requested_region_mode"
  printf 'archaic_refs_source\t%s\n' "$PHYML_REFS_SOURCE"
  stat -c 'script\t%n:%s:%Y' "$F/phyml.py" "$F/phyml/"*.py "$F/phyml/tree_plot.R" "$PHE_F"
  [[ -z ${GU_LOCI_FILE:-} ]] || stat -c 'input\t%n:%s:%Y' "$GU_LOCI_FILE"
  [[ -z ${GU_LOCI_MAP_FILE:-} ]] || stat -c 'input\t%n:%s:%Y' "$GU_LOCI_MAP_FILE"
  [[ -z ${GU_SAMPLE_PANEL:-} ]] || stat -c 'input\t%n:%s:%Y' "$GU_SAMPLE_PANEL"
  if [[ -n ${GU_TARGET_TMP_DIR:-} ]]; then
    for c in ${GU_CHRS:-}; do
      source_meta=$GU_TARGET_TMP_DIR/chr${c}/source.tsv
      [[ -s $source_meta ]] || { echo "ERROR: missing stable target provenance: $source_meta" >&2; return 1; }
      awk -v c="$c" 'BEGIN{OFS="\t"}{print "target","chr" c,$0}' "$source_meta"
    done
  else
    find "${GU_TARGET_VCF_DIR:?}" -type f \( -name '*.vcf.gz' -o -name '*.bcf' \) -printf 'target\t%p:%s:%T@\n' 2>/dev/null | sort
  fi
  find "${GU_ARCHAIC_ROOT:?}" -type f \( -name '*.vcf.gz' -o -name '*.bcf' \) -printf 'input\t%p:%s:%T@\n' 2>/dev/null | sort
}

comparison_record_without_ephemeral_vcf_dir(){
  awk -F'\t' 'BEGIN{OFS="\t"}$1=="argv"{for(i=1;i<NF;i++)if($i=="--vcf-dir")$(i+1)="PREPARED_TARGET_VCF_DIR"}1' "$1"
}

legacy_target_sources_not_newer_than(){
  local old_meta=$1 cutoff c source_meta line record stamp
  cutoff=$(stat -c '%Y' "$old_meta")
  for c in ${GU_CHRS:-}; do
    source_meta=$GU_TARGET_TMP_DIR/chr${c}/source.tsv
    [[ -s $source_meta ]] || return 1
    while IFS= read -r line; do
      case "$line" in
        pgen$'\t'*|pvar$'\t'*|psam$'\t'*|vcf$'\t'*|index$'\t'*)
          record=${line#*$'\t'}
          stamp=${record##*:}; stamp=${stamp%%.*}
          [[ $stamp =~ ^[0-9]+$ && $stamp -le $cutoff ]] || return 1
          ;;
      esac
    done < "$source_meta"
  done
}

comparison_meta=$OUT/final/comparison.input.meta.tsv
expected_meta=$(mktemp "${TMPDIR:-/tmp}/gu.phyml.meta.XXXXXX")
legacy_expected=$(mktemp "${TMPDIR:-/tmp}/gu.phyml.legacy.XXXXXX")
legacy_actual=$(mktemp "${TMPDIR:-/tmp}/gu.phyml.actual.XXXXXX")
trap 'rm -f -- "$expected_meta" "$legacy_expected" "$legacy_actual"' EXIT
comparison_record > "$expected_meta"
comparison_ready=0
if [[ ${PHYML_REPLACE:-FALSE} == TRUE ]]; then
  echo "[GU PHYML] comparison=OVERWRITE reason=--replace-phyml output=$OUT"
elif [[ -s $OUT/final/haplotypes.tsv && -s $OUT/final/loci.tsv && -s $OUT/final/haplotype_samples.tsv && -s $OUT/final/skipped_loci.tsv && -s $comparison_meta ]]; then
  if cmp -s "$expected_meta" "$comparison_meta"; then
    comparison_ready=1
  else
    grep -v $'^target\t' "$expected_meta" > "$legacy_expected"
    comparison_record_without_ephemeral_vcf_dir "$comparison_meta" > "$legacy_actual"
    if cmp -s "$legacy_expected" "$legacy_actual" && legacy_target_sources_not_newer_than "$comparison_meta"; then
      comparison_ready=1
      cp "$expected_meta" "$comparison_meta.tmp.$$"
      mv -f -- "$comparison_meta.tmp.$$" "$comparison_meta"
      echo "[GU PHYML] comparison_metadata=MIGRATED stable_target_provenance=TRUE output=$OUT"
    fi
  fi
fi

if (( comparison_ready )); then
  echo "[GU PHYML] comparison=SKIP reason=outputs_exist output=$OUT"
else
  echo "[GU PHYML] comparison=RUN reason=missing_or_stale_outputs output=$OUT"
  rm -f -- "$OUT/final/trees.tsv" "$OUT/final/evidence_trees.tsv" "$OUT/final/region_common_trees.tsv"
  python3 "$F/phyml.py" compare "${args[@]}"
  mkdir -p "$OUT/final"
  mv -f -- "$expected_meta" "$comparison_meta"
  expected_meta=''
  find "$OUT/loci" -type f \( -name '*_phyml_tree.txt' -o -name '*_phyml_stats.txt' -o -name '*_phyml_tree.png' -o -name '*.phyml.log' -o -name '*.phyml.run.status.tsv' \) -delete 2>/dev/null || true
fi

# Stage 1: complete anchor-independent inventory. This runs before exact-anchor
# extraction and before YRI/AFR/ancestral filters by design.
if (( evidence_enabled )); then
  python3 "$F/phyml.py" region prepare "${region_scan_args[@]}"

  # Stage 2a: optional exact BED-column-4 marker. Missing or incorrect markers
  # are recorded but never remove Stage-1 regional matches.
  if (( ${#anchor_args[@]} )); then
    python3 "$F/phyml.py" anchor "${anchor_args[@]}"
  fi

  # Stage 2b: derived-state, outgroup-frequency and recurrence filters; Stage 3
  # candidate-focused trees are prepared from these rows.
  python3 "$F/phyml.py" evidence prepare "${evidence_args[@]}"
else
  echo "[GU PHYML] layered_evidence=SKIP reason=whole_chromosome_scan; set PHYML_EVIDENCE_WHOLE_CHROMOSOME=TRUE to enable"
fi

run_phyml_one(){
  local phy=$1 scope=$2 tree=${phy}_phyml_tree.txt log=${phy}.phyml.log status_file=${phy}.phyml.run.status.tsv
  local rc=0 start end timeout_value=${PHYML_TREE_TIMEOUT:-0}
  if [[ -s $tree ]]; then
    printf 'scope\tstatus\trc\tphy\ttree\tlog\n%s\tSKIP_EXISTING\t0\t%s\t%s\t%s\n' "$scope" "$phy" "$tree" "$log" > "$status_file"
    return 0
  fi
  start=$(date +%s)
  echo "[GU PHYML] tree=RUN scope=$scope input=$phy bootstrap=${PHYML_BOOT:-100}"
  if [[ $timeout_value != 0 && -n $timeout_value ]] && command -v timeout >/dev/null 2>&1; then
    timeout "$timeout_value" phyml -i "$phy" -m HKY85 -c 4 -a e -v e -b "${PHYML_BOOT:-100}" > "$log" 2>&1 || rc=$?
  else
    phyml -i "$phy" -m HKY85 -c 4 -a e -v e -b "${PHYML_BOOT:-100}" > "$log" 2>&1 || rc=$?
  fi
  end=$(date +%s)
  if (( rc == 0 )) && [[ -s $tree ]]; then
    printf 'scope\tstatus\trc\telapsed_seconds\tphy\ttree\tlog\n%s\tCOMPLETE\t0\t%s\t%s\t%s\t%s\n' "$scope" "$((end-start))" "$phy" "$tree" "$log" > "$status_file"
    return 0
  fi
  printf 'scope\tstatus\trc\telapsed_seconds\tphy\ttree\tlog\n%s\tFAILED\t%s\t%s\t%s\t%s\t%s\n' "$scope" "$rc" "$((end-start))" "$phy" "$tree" "$log" > "$status_file"
  echo "WARNING: PhyML tree failed scope=$scope rc=$rc input=$phy log=$log" >&2
  return 1
}

status=0
if (( plot_phy )); then
  if (( evidence_enabled )); then
    while IFS= read -r phy; do
      [[ -n $phy ]] || continue
      run_phyml_one "$phy" region_common_unfiltered || status=1
    done < <(region_common_phy_inputs)
    while IFS= read -r phy; do
      [[ -n $phy ]] || continue
      run_phyml_one "$phy" candidate_filtered || status=1
    done < <(evidence_phy_inputs)
  fi
  if [[ ${PHYML_RUN_ALL_UNIQUE_TREE:-FALSE} == TRUE ]]; then
    while IFS= read -r phy; do
      [[ -n $phy ]] || continue
      run_phyml_one "$phy" all_unique_exploratory || status=1
    done < <(raw_phy_inputs)
  else
    echo "[GU PHYML] all_unique_tree=SKIP default=FALSE reason=large_tree_is_exploratory; set PHYML_RUN_ALL_UNIQUE_TREE=TRUE to request"
  fi
fi

if (( plot_phy )); then
  python3 "$F/phyml.py" plot --out "$OUT" \
    --dpi "${PHYML_TREE_PLOT_DPI:-220}" \
    --max-labels "${PHYML_TREE_LABEL_LIMIT:-180}" || status=1
fi

# Always write summaries, including explicit not-requested/not-run states.
if (( evidence_enabled )); then
  tree_default=not_requested
  (( plot_phy )) && tree_default=not_run
  python3 "$F/phyml.py" region summarize "${region_scan_args[@]}" --default-tree-status "$tree_default"

  # Establish one compatibility row per locus even when the huge all-unique tree
  # is intentionally not run. phyml/evidence.py then replaces its interpretation
  # with the candidate-focused authoritative fields.
  python3 "$F/phyml.py" tree --out "$OUT" --skip-render --default-status not_requested
  python3 "$F/phyml.py" evidence summarize "${evidence_args[@]}" --default-tree-status "$tree_default"
else
  python3 "$F/phyml.py" tree --out "$OUT" --skip-render --default-status not_requested
fi

if (( status != 0 )); then
  echo "WARNING: one or more trees failed; all completed comparisons and partial tree summaries were preserved" >&2
  [[ ${PHYML_TREE_FAILURE_POLICY:-warn} == warn ]] || exit 1
fi
