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

restore_completed_tree_from_summary(){
  local phy=$1 tree=${phy}_phyml_tree.txt locus tmp
  [[ ! -s $tree && -s $OUT/final/trees.tsv ]] || return 0
  if [[ $(dirname -- "$phy") == "$OUT/loci" ]]; then
    locus=$(awk -F'\t' 'NR>1&&$1!=""{print $1; exit}' "$OUT/final/loci.tsv")
  else
    locus=$(basename -- "$(dirname -- "$phy")")
  fi
  [[ -n $locus ]] || return 0
  tmp=$tree.restore.$$
  awk -F'\t' -v locus="$locus" '
    NR==1{for(i=1;i<=NF;i++){if($i=="locus_id")id=i;else if($i=="tree_status")status=i;else if($i=="tree_newick")newick=i};next}
    $id==locus && $status=="complete" && $newick!=""{print $newick; found=1; exit}
    END{exit !found}' "$OUT/final/trees.tsv" > "$tmp" || { rm -f -- "$tmp"; return 0; }
  [[ -s $tmp ]] || { rm -f -- "$tmp"; return 0; }
  mv -f -- "$tmp" "$tree"
  echo "[GU PHYML] tree=RESTORE reason=completed_final_summary file=$tree"
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
  local i c source_meta
  local -a record_args=("${args[@]}")
  # GU_TARGET_VCF_DIR is a disposable run.XXXXXX directory.  Its path must not
  # invalidate an otherwise identical scientific request.
  for (( i=0; i<${#record_args[@]}; i++ )); do
    if [[ ${record_args[$i]} == --vcf-dir && $((i+1)) -lt ${#record_args[@]} ]]; then
      record_args[$((i+1))]='PREPARED_TARGET_VCF_DIR'
    fi
  done
  printf 'argv'
  printf '\t%q' "${record_args[@]}"
  printf '\n'
  stat -c 'script\t%n:%s:%Y' "$F/phyml.py" "$F/phyml_tree_summary.py" "$F/phyml_tree_plot.R" "$PHE_F"
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
        pgen\\t*|pvar\\t*|psam\\t*|vcf\\t*|index\\t*)
          record=${line#*\\t}
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
    # Migrate metadata written before target provenance became stable.  Accept
    # it only when every current target source predates that completed record.
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
  rm -f -- "$OUT/final/trees.tsv"
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
    (( comparison_ready )) && restore_completed_tree_from_summary "$phy"
    [[ -s ${phy}_phyml_tree.txt ]] && continue
    phyml -i "$phy" -m HKY85 -c 4 -a e -v e -b "${PHYML_BOOT:-100}" > "${phy}.phyml.log" 2>&1 || status=1
  done < <(phyml_input_paths)
  (( status == 0 )) || { echo "ERROR: one or more PhyML jobs failed" >&2; exit 1; }
  python3 "$F/phyml_tree_summary.py" --out "$OUT"
elif [[ $ACTION == run ]]; then
  echo "[GU PHYML] plot_phy=FALSE tree_inference=SKIP; use --plot-phy TRUE to build trees"
  python3 "$F/phyml_tree_summary.py" --out "$OUT" --default-status not_requested
fi
