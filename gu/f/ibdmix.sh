#!/usr/bin/env bash
# Lean IBDmix caller: computation only. Reporting, summaries and visualization live in Shiny.
set -Eeuo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
F=$ROOT/f
PHE_F=${PHE_F:-$ROOT/../0f/0phe.f.sh}
[[ -s $PHE_F ]] || { echo "ERROR: missing $PHE_F" >&2; exit 1; }
# shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
source "$PHE_F"

action=${GU_ACTION:-ibdmix_run}
case "$action" in ibdmix_run|ibdmix_check) ;; *) echo "ERROR: unsupported IBDmix action: $action" >&2; exit 2;; esac
for cmd in awk bcftools gzip python3 sort comm find stat cmp grep samtools sha256sum flock; do command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing command: $cmd" >&2; exit 1; }; done

dir0=${dir0:-/mnt/d}
dir_ref=${dir_ref:-/mnt/e/refGen}
dirmod=${dirmod:-${GU_TARGET_ROOT:-$dir_ref/1kg/${GRCH:-37}}}
dirarch=${dirarch:-${GU_ARCHAIC_ROOT:-$dir_ref/archaic/${GRCH:-37}/vcf}}
sample_file=${sample_file:-$dirmod/samples.txt}
target_vcf_dir=${GU_TARGET_VCF_DIR:-$dirmod/vcf}
dirsoft=${dirsoft:-$dir0/software/gu/IBDmix}
dirout=${dirout:-${GU_ANALYSIS_ROOT:-$dir0/analysis/gu}/ibdmix/${GU_SCOPE_ID:-genome}}
genome_build=${genome_build:-b${GRCH:-37}}
loci_file=${GU_LOCI_FILE:-}
chrs_arg=${chrs:-${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}}
profile=${IBDMIX_PROFILE:-multi_reference}
default_refs="Altai Chagyr Vindija Denisova Denisova25"
[[ $profile == multi_reference ]] || default_refs="Altai Denisova"
refs_arg=${refs:-$default_refs}
read -r -a refs <<< "$refs_arg"
for i in "${!refs[@]}"; do
  case "${refs[$i],,}" in
    altai) refs[$i]=Altai;; denisova|denisovan) refs[$i]=Denisova;; denisova25|den25) refs[$i]=Denisova25;;
    vindija) refs[$i]=Vindija;; chagyr|chagyrskaya) refs[$i]=Chagyr;;
    *) echo "ERROR: unknown IBDmix reference: ${refs[$i]}" >&2; exit 2;;
  esac
done
refs_arg="${refs[*]}"
lod_cut=${lod_cut:-4}
len_cut=${len_cut:-50000}
emit_lod_cut=${emit_lod_cut:-$lod_cut}
minor_allele_count=${minor_allele_count:-1}
archaic_error=${archaic_error:-0.01}
modern_error_max=${modern_error_max:-0.002}
modern_error_proportion=${modern_error_proportion:-2}
# Keep analysis units serial by default to limit temporary I/O and memory.
job_of_chr=${job_of_chr:-1}
job_in_chr=${job_in_chr:-1} # populations processed independently; bounded workers per chromosome.
IBDMIX_LOCUS_FLANK_BP=${IBDMIX_LOCUS_FLANK_BP:-100000}
IBDMIX_REPLACE=${IBDMIX_REPLACE:-0}
GU_CHRX_PAR_DIPLOID=${GU_CHRX_PAR_DIPLOID:-0}

generate_gt=$dirsoft/build/src/generate_gt
ibdmix_bin=$dirsoft/build/src/ibdmix
helper=$F/ibdmix_workflow.py
# Semantic reuse contract: bump for changes to genotype preparation, calling,
# filtering or output meaning; path/log/comment edits do not change it.
pipeline_version=2026-09-14.1
# Published resources only; generated masks live under $dirout/mask/derived.
[[ -z ${IBDMIX_MASK_CACHE:-} ]] || { echo "ERROR: IBDMIX_MASK_CACHE is retired; use IBDMIX_MASK_ROOT" >&2; exit 2; }
mask_root=${IBDMIX_MASK_ROOT:-$(dirname -- "$dirarch")/mask}
custom_masks=${IBDMIX_MASK_DIR:-}
ref_fasta=${IBDMIX_REFERENCE_FASTA:-$dir_ref/fasta/GRCH37.fasta}
background_filter=${IBDMIX_AFR_DENISOVAN_FILTER:-1}
case "$background_filter" in 0|1) ;; *) echo "ERROR: IBDMIX_AFR_DENISOVAN_FILTER must be 0 or 1" >&2; exit 2;; esac
case "$profile" in
  cell2020)
    [[ $genome_build == b37 && $refs_arg == 'Altai Denisova' && $background_filter == 1 ]] || { echo "ERROR: cell2020 requires GRCh37, refs='Altai Denisova' and African Denisovan filtering" >&2; exit 2; }
    [[ -z $custom_masks ]] || { echo "ERROR: supplied masks require IBDMIX_PROFILE=custom; the default profile builds the documented Cell masks" >&2; exit 2; }
    ;;
  custom)
    [[ -n $custom_masks ]] || { echo "ERROR: custom IBDmix requires IBDMIX_MASK_DIR with excluded-site BED files <ref>/chrN.bed" >&2; exit 2; }
    ;;
  multi_reference)
    [[ $genome_build == b37 ]] || { echo "ERROR: multi_reference masks require GRCh37" >&2; exit 2; }
    ;;
  *) echo "ERROR: IBDMIX_PROFILE must be cell2020, multi_reference or custom" >&2; exit 2;;
esac
export_denisovan=0
[[ $profile != multi_reference ]] || export_denisovan=1
denisovan_args=(); [[ $export_denisovan == 0 ]] || denisovan_args=(--export-denisovan)
if [[ $background_filter == 1 && " $refs_arg " != *' Denisova '* ]]; then echo "ERROR: background filtering requires Denisova calls" >&2; exit 2; fi
background_args=(); [[ $background_filter == 1 ]] || background_args=(--no-background)
[[ $genome_build == b37 || $genome_build == b38 ]] || { echo "ERROR: genome_build must be b37 or b38" >&2; exit 1; }
[[ $IBDMIX_LOCUS_FLANK_BP =~ ^[0-9]+$ ]] || { echo "ERROR: IBDMIX_LOCUS_FLANK_BP must be a non-negative integer" >&2; exit 1; }
[[ $IBDMIX_REPLACE == 0 || $IBDMIX_REPLACE == 1 ]] || { echo "ERROR: internal replace-ibdmix state must be 0 or 1" >&2; exit 1; }
[[ $job_of_chr =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: job_of_chr must be a positive integer" >&2; exit 1; }
[[ $job_in_chr =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: job_in_chr must be a positive integer" >&2; exit 1; }

final_output=$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz

for file in "$generate_gt" "$ibdmix_bin" "$helper" "$F/ibdmix_provenance.py" "$F/vcf_gt_fix.py" "$sample_file"; do [[ -s $file ]] || { echo "ERROR: missing required file: $file" >&2; exit 1; }; done
[[ -d $target_vcf_dir ]] || { echo "ERROR: missing target VCF directory: $target_vcf_dir" >&2; exit 1; }
mkdir -p "$dirout"/{samples,genotype,raw,segments,final,log,tmp}
exec 9> "$dirout/.run.lock"
flock -n 9 || { echo "ERROR: another IBDmix run is using $dirout" >&2; exit 1; }

exec > >(tee "$dirout/ibdmix.log") 2>&1
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
log "CONFIG profile=$profile refs=$refs_arg export_denisovan=$export_denisovan background_filter=$background_filter"
gzip_ok(){ [[ -s $1 ]] && gzip -t "$1" >/dev/null 2>&1; }

declare -a analysis_units=()
declare -A unit_chr=() unit_start=() unit_end=() unit_core_start=() unit_core_end=() unit_sex=() unit_locus=()
add_unit(){
  local key=$1 chr=$2 start=$3 end=$4 sex=$5 locus=$6 core_start=${7:-0} core_end=${8:-0}
  (( end >= start )) || return 0
  analysis_units+=("$key"); unit_chr[$key]=$chr; unit_start[$key]=$start; unit_end[$key]=$end
  unit_core_start[$key]=$core_start; unit_core_end[$key]=$core_end; unit_sex[$key]=$sex; unit_locus[$key]=$locus
}
if [[ -n $loci_file ]]; then
  [[ -s $loci_file ]] || { echo "ERROR: GU_LOCI_FILE missing: $loci_file" >&2; exit 1; }
  index=0
  while IFS=$'\t' read -r chr start end name; do
    [[ -n $chr ]] || continue
    index=$((index + 1)); stem=$(printf 'L%05d' "$index")
    if [[ $chr == X ]]; then
      scan_start=$((start + 1 - IBDMIX_LOCUS_FLANK_BP)); (( scan_start < 1 )) && scan_start=1
      scan_end=$((end + IBDMIX_LOCUS_FLANK_BP))
      if [[ $GU_CHRX_PAR_DIPLOID == 1 ]]; then
        add_unit "${stem}_X_PAR" X "$scan_start" "$scan_end" all "$name" "$((start+1))" "$end"
      else
        add_unit "${stem}_X_MALE" X "$scan_start" "$scan_end" male "$name" "$((start+1))" "$end"
      fi
    else
      scan_start=$((start + 1 - IBDMIX_LOCUS_FLANK_BP)); (( scan_start < 1 )) && scan_start=1
      scan_end=$((end + IBDMIX_LOCUS_FLANK_BP))
      add_unit "$stem" "$chr" "$scan_start" "$scan_end" all "$name" "$((start+1))" "$end"
    fi
  done < "$loci_file"
else
  read -r -a requested <<< "${chrs_arg//,/ }"
  for chr in "${requested[@]}"; do
    chr=${chr#chr}; [[ $chr == 23 ]] && chr=X
    if [[ $chr == X ]]; then
      add_unit X_MALE X 0 0 male chrX
    else add_unit "C$chr" "$chr" 0 0 all "chr$chr"
    fi
  done
fi
if [[ $profile == cell2020 ]]; then
  for u in "${analysis_units[@]}"; do
    [[ ${unit_chr[$u]} =~ ^([1-9]|1[0-9]|2[0-2]|X)$ ]] || { echo "ERROR: unsupported IBDmix chromosome: ${unit_chr[$u]}" >&2; exit 2; }
  done
fi
(( ${#analysis_units[@]} )) || { echo "ERROR: no IBDmix analysis units" >&2; exit 1; }

unit_region(){ local u=$1; (( unit_start[$u] > 0 )) && printf '%s:%s-%s\n' "${unit_chr[$u]}" "${unit_start[$u]}" "${unit_end[$u]}" || printf '\n'; }
modern_vcf(){ printf '%s/chr%s.vcf.gz\n' "$target_vcf_dir" "$1"; }
ref_dir(){ case "${1,,}" in altai) printf '%s/Altai\n' "$dirarch";; denisova|denisovan) printf '%s/Denisova\n' "$dirarch";; *) printf '%s/%s\n' "$dirarch" "$1";; esac; }
archaic_vcf(){
  local ref=$1 chr=$2 d f
  for d in "$(ref_dir "$ref")" "$dirarch/avcf/$(basename "$(ref_dir "$ref")")"; do
    [[ -d $d ]] || continue
    f=$(find "$d" -maxdepth 1 -type f \( -name "*chr${chr}_*.vcf.gz" -o -name "*chr${chr}.*.vcf.gz" -o -name "*chr${chr}.vcf.gz" \) | sort | head -1)
    [[ -n $f ]] && { printf '%s\n' "$f"; return 0; }
  done
  return 1
}
vcf_contig(){
  local vcf=$1 wanted=$2
  bcftools view -h "$vcf" | awk -F'[=,>]' -v c="$wanted" '
    /^##contig=<ID=/{id=$3; if(id==c || id=="chr" c || (c=="X" && id=="23") || (c=="23" && id=="X")){print id; found=1; exit}}
    END{if(!found)exit 1}'
}
adapt_region(){
  local vcf=$1 region=$2 contig
  [[ -n $region ]] || { printf '\n'; return; }
  contig=$(vcf_contig "$vcf" "$unit_chr_current") || { echo "ERROR: contig for chr$unit_chr_current not declared in $vcf" >&2; return 1; }
  printf '%s\n' "${region/${unit_chr_current}:/$contig:}"
}

prepare_samples(){
  local u vcf d
  for u in "${analysis_units[@]}"; do
    vcf=$(modern_vcf "${unit_chr[$u]}"); d=$dirout/samples/$u
    mkdir -p "$d"
    bcftools query -l "$vcf" > "$d/actual.txt"
    python3 "$helper" samples --panel "$sample_file" --actual "$d/actual.txt" --output "$d" "${background_args[@]}"
  done
  # Compatibility roster for downstream read-only population summaries.
  cp "$dirout/samples/${analysis_units[0]}/ALL.txt" "$dirout/samples/ALL.txt"
  if [[ ${unit_sex[${analysis_units[0]}]} == male ]]; then cp "$dirout/samples/ALL.txt" "$dirout/samples/male.txt"; fi
}
validate_inputs(){
  local u chr ref vcf mask
  local -a mask_chroms=() mask_args=()
  for u in "${analysis_units[@]}"; do mask_chroms+=("${unit_chr[$u]}"); done
  mask_args=(--root "$mask_root" --archaic-root "$dirarch" --chroms "${mask_chroms[@]}" --refs "${refs[@]}")
  [[ -z $custom_masks ]] || mask_args+=(--custom-masks "$custom_masks")
  python3 "$helper" check-masks "${mask_args[@]}" || return $?
  for u in "${analysis_units[@]}"; do
    chr=${unit_chr[$u]}; vcf=$(modern_vcf "$chr")
    [[ -s $vcf ]] || { echo "ERROR: missing modern VCF: $vcf" >&2; return 1; }
    for ref in "${refs[@]}"; do
      local avcf
      avcf=$(archaic_vcf "$ref" "$chr") || { echo "ERROR: missing archaic VCF ref=$ref chr=$chr" >&2; return 1; }
      [[ $(bcftools query -l "$avcf" | wc -l) == 1 ]] || { echo "ERROR: expected exactly one archaic sample in $avcf" >&2; return 1; }
      if [[ -n $custom_masks ]]; then
        mask=$custom_masks/$ref/chr$chr.bed
        [[ -f $mask ]] || { echo "ERROR: excluded-site BED missing: $mask" >&2; return 1; }
        python3 "$helper" validate-mask --path "$mask" --chrom "$chr" --build "$genome_build"
      fi
    done
  done
  if [[ -z $custom_masks ]]; then
    for u in "${analysis_units[@]}"; do
      [[ ${unit_chr[$u]} == X ]] && continue
      [[ -s $ref_fasta && -s $ref_fasta.fai ]] || { echo "ERROR: indexed GRCh37 FASTA required: $ref_fasta" >&2; return 1; }
    done
  fi
  prepare_samples
  log "INPUT CHECK PASSED: per-population calls; profile=$profile units=${analysis_units[*]}"
}
run_record(){
  local u chr ref vcf index meta line x_mode=not_requested
  for u in "${analysis_units[@]}"; do
    [[ ${unit_chr[$u]} != X ]] || x_mode=$([[ ${unit_sex[$u]} == male ]] && printf male_haploid_nonpar || printf par_diploid_chrxy)
  done
  printf 'pipeline_version\t%s\nprofile\t%s\ncoordinate_system\t0-based-half-open\npopulation_mode\tpopulation\nbackground_filter\t%s\nmask_root\t%s\ncustom_masks\t%s\n' "$pipeline_version" "$profile" "$background_filter" "$mask_root" "$custom_masks"
  printf 'export_denisovan\t%s\ndenisovan_call_definition\tnative_reference_matching_no_African_Denisova_subtraction\n' "$export_denisovan"
  printf 'genome_build\t%s\nx_mode\t%s\nloci_flank_bp\t%s\nrefs\t%s\nlod_cut\t%s\nlen_cut\t%s\nemit_lod_cut\t%s\nminor_allele_count\t%s\narchaic_error\t%s\nmodern_error_max\t%s\nmodern_error_proportion\t%s\n' \
    "$genome_build" "$x_mode" "$IBDMIX_LOCUS_FLANK_BP" "$refs_arg" "$lod_cut" "$len_cut" "$emit_lod_cut" "$minor_allele_count" "$archaic_error" "$modern_error_max" "$modern_error_proportion"
  stat -c 'sample_file\t%n:%s:%Y' "$sample_file"
  [[ -z $loci_file ]] || stat -c 'loci_file\t%n:%s:%Y' "$loci_file"
  stat -c 'software\t%n:%s:%Y' "$generate_gt" "$ibdmix_bin" "$helper" "$F/vcf_gt_fix.py" "$F/ibdmix.sh"
  for u in "${analysis_units[@]}"; do
    [[ ! -f $dirout/mask/$u/manifest.json ]] || printf 'mask_manifest_sha256\t%s\t%s\n' "$u" "$(sha256sum "$dirout/mask/$u/manifest.json" | cut -d' ' -f1)"
    for ref in "${refs[@]}"; do
      printf 'excluded_mask_sha256\t%s\t%s\t%s\n' "$u" "$ref" "$(sha256sum "$dirout/mask/$u/$ref.exclude.bed" | cut -d' ' -f1)"
    done
  done
  [[ ! -f $ref_fasta ]] || stat -Lc 'reference_fasta\t%n:%s:%Y' "$ref_fasta"
  declare -A seen=()
  for u in "${analysis_units[@]}"; do
    chr=${unit_chr[$u]}; [[ -z ${seen[$chr]:-} ]] || continue; seen[$chr]=1
    vcf=$(modern_vcf "$chr"); meta=${GU_TARGET_TMP_DIR:-$dirout/tmp}/chr${chr}/source.tsv
    if [[ -s $meta ]]; then
      # The VCF is a disposable pfile export. Track the source pgen/pvar/psam
      # conversion contract, not the regenerated VCF/index mtimes.
      while IFS= read -r line; do printf 'modern_source\tchr%s\t%s\n' "$chr" "$line"; done < "$meta"
    else
      # A native VCF has no conversion record and remains the source input.
      stat -c 'modern_vcf\t%n:%s:%Y' "$vcf"
      for index in "$vcf.tbi" "$vcf.csi"; do [[ -s $index ]] && stat -c 'modern_index\t%n:%s:%Y' "$index"; done
    fi
    for ref in "${refs[@]}"; do
      vcf=$(archaic_vcf "$ref" "$chr"); stat -c 'archaic_vcf\t%n:%s:%Y' "$vcf"
      [[ -z $custom_masks ]] || stat -Lc 'excluded_mask\t%n:%s:%Y' "$custom_masks/$ref/chr$chr.bed"
    done
  done
}
check_run_provenance(){
  local current=$dirout/run.meta.tsv candidate=$dirout/tmp/run.meta.current.tsv
  run_record > "$candidate"
  if [[ -s $current ]] && ! python3 "$F/ibdmix_provenance.py" "$current" "$candidate"; then
    if [[ $IBDMIX_REPLACE != 1 ]] && [[ -e $dirout/.complete || -n $(find "$dirout/genotype" "$dirout/raw" "$dirout/segments" "$dirout/final" -type f ! -name '*.part*' -print -quit) ]]; then
      echo "ERROR: IBDmix inputs, profile or implementation changed; old results cannot be reused. Rerun with --replace-ibdmix TRUE." >&2
      return 1
    fi
  fi
  if [[ $IBDMIX_REPLACE == 1 ]]; then
    rm -f "$dirout/.complete"
    rm -rf -- "$dirout/genotype" "$dirout/raw" "$dirout/segments" "$dirout/final"
    mkdir -p "$dirout"/{genotype,raw,segments,final,mask}
  elif [[ ! -s $current ]] && [[ -e $dirout/.complete || -n $(find "$dirout/genotype" "$dirout/raw" "$dirout/segments" "$dirout/final" -type f -print -quit) ]]; then
    echo "ERROR: cached IBDmix files have no provenance; rerun with --replace-ibdmix TRUE" >&2; return 1
  fi
  mv -f "$candidate" "$current"
}
completed_output_ok(){
  [[ -e $dirout/.complete ]] && gzip_ok "$final_output"
}

prepare_modern(){
  local u=$1 out=$2 chr vcf region adapted
  chr=${unit_chr[$u]}; vcf=$(modern_vcf "$chr"); region=$(unit_region "$u")
  unit_chr_current=$chr; adapted=$(adapt_region "$vcf" "$region")
  local -a view=(view -Ou); [[ -z $adapted ]] || view+=(-r "$adapted");view+=("$vcf")
  # Keep non-SNP records until generate_gt sees them: deleting a modern indel
  # beforehand can incorrectly turn that site into a homozygous-reference call.
  local -a fix=()
  [[ $chr != X ]] || fix=(--chrom 23 --duplicate-haploid slash)
  bcftools "${view[@]}" | bcftools +setGT -Ou -- -t . -n . |
    bcftools +fixploidy -Ou -- -f 2 |
    bcftools annotate -x INFO,^FORMAT/GT -Ov |
    if [[ $chr == X ]]; then python3 "$F/vcf_gt_fix.py" "${fix[@]}"; else awk 'BEGIN{OFS="\t"} /^#/{print;next}{sub(/^chr/,"",$1);print}'; fi > "$out"
}
prepare_archaic(){
  local u=$1 ref=$2 out=$3 chr vcf region adapted
  chr=${unit_chr[$u]}; vcf=$(archaic_vcf "$ref" "$chr");region=$(unit_region "$u")
  unit_chr_current=$chr; adapted=$(adapt_region "$vcf" "$region")
  local -a view=(view -Ou); [[ -z $adapted ]] || view+=(-r "$adapted");view+=("$vcf")
  # ALT='.' / 0/0 sites carry essential discordant-homozygote evidence.
  # Upstream generate_gt itself handles biallelic/SNV compatibility.
  bcftools "${view[@]}" | bcftools +setGT -Ou -- -t . -n . |
    bcftools +fixploidy -Ou -- -f 2 |
    bcftools annotate -x INFO,^FORMAT/GT -Ov |
    if [[ $chr == X ]]; then python3 "$F/vcf_gt_fix.py" --chrom 23 --duplicate-haploid slash; else awk 'BEGIN{OFS="\t"} /^#/{print;next}{sub(/^chr/,"",$1);print}'; fi > "$out"
}
prepare_unit_masks(){
  local u=$1 ref chr=${unit_chr[$1]} out=$dirout/mask/$1
  mkdir -p "$out"
  if [[ -n $custom_masks ]]; then
    for ref in "${refs[@]}"; do
      if [[ $chr == X ]]; then awk 'BEGIN{OFS="\t"}{$1=23;print}' "$custom_masks/$ref/chr$chr.bed" > "$out/$ref.exclude.bed"
      else cp "$custom_masks/$ref/chr$chr.bed" "$out/$ref.exclude.bed"; fi
    done
  else
    python3 "$helper" mask --root "$mask_root" --chrom "$chr" --modern "$(modern_vcf "$chr")" --fasta "$ref_fasta" --upstream "$dirsoft" --output "$out" --archaic-root "$dirarch" --refs "${refs[@]}"
  fi
}

run_population(){
  local u=$1 ref=$2 pop=$3 sample_list=$4 gt=$5 tmp=$6 raw
  raw=$dirout/raw/$u/$ref.$pop.raw.txt.gz
  if gzip_ok "$raw"; then return 0; fi
  log "CALL unit=$u ref=$ref population=$pop"
  "$ibdmix_bin" --genotype <(gzip -dc "$gt") --output "$tmp/$ref.$pop.raw.txt" --sample "$sample_list" \
    --mask "$dirout/mask/$u/$ref.exclude.bed" --LOD-threshold "$emit_lod_cut" \
    --minor-allele-count-threshold "$minor_allele_count" --archaic-error "$archaic_error" --modern-error-max "$modern_error_max" \
    --modern-error-proportion "$modern_error_proportion" --more-stats > "$dirout/log/$ref.$u.$pop.ibdmix.log" 2>&1 || return 1
  gzip -c "$tmp/$ref.$pop.raw.txt" > "$raw.part" && mv "$raw.part" "$raw" || return 1
  rm -f "$tmp/$ref.$pop.raw.txt"
}
run_ref_unit(){
  local u=$1 ref=$2 tmp=$3 gt producer archaic_fd modern_fd modern_pid pop _super _n sample_list status=0 producer_status=0 modern_status=0 pid
  local -a pids=()
  gt=$dirout/genotype/$ref.$u.gt.txt.gz
  if ref_calls_complete "$u" "$ref"; then
    rm -f "$gt" "$gt.part.gz"
    log "REUSE unit=$u ref=$ref population_calls=complete genotype=SKIP"
    return 0
  fi
  if ! gzip_ok "$gt"; then
    log "GENOTYPE unit=$u ref=$ref stage=start"
    # Stream directly from the source; do not cache another modern VCF.
    exec {modern_fd}< <(prepare_modern "$u" /dev/stdout)
    modern_pid=$!
    exec {archaic_fd}< <(exec {modern_fd}<&-; prepare_archaic "$u" "$ref" /dev/stdout)
    producer=$!
    if ! "$generate_gt" -a "/dev/fd/$archaic_fd" -m "/dev/fd/$modern_fd" -o - 2> "$dirout/log/$ref.$u.generate_gt.log" |
        gzip -1c > "$gt.part.gz"; then
      exec {archaic_fd}<&-; exec {modern_fd}<&-
      wait "$producer" 2>/dev/null || true; wait "$modern_pid" 2>/dev/null || true
      rm -f "$gt.part.gz"; return 1
    fi
    exec {archaic_fd}<&-; exec {modern_fd}<&-
    # Modern EOF is required. The trailing all-sites archaic stream may
    # legitimately receive SIGPIPE once generate_gt has consumed modern EOF.
    wait "$producer" || producer_status=$?
    wait "$modern_pid" || modern_status=$?
    if [[ $modern_status != 0 || ( $producer_status != 0 && $producer_status != 141 ) ]]; then
      rm -f "$gt.part.gz"; return 1
    fi
    python3 "$helper" validate-genotypes --path "$gt.part.gz" --output "$dirout/genotype/$ref.$u.qc.json"
    mv "$gt.part.gz" "$gt"
    log "GENOTYPE unit=$u ref=$ref stage=complete"
  fi
  mkdir -p "$dirout/raw/$u"
  while IFS=$'\t' read -r pop _super _n sample_list; do
    [[ $pop != population ]] || continue
    [[ $export_denisovan == 1 || $background_filter == 0 || $ref != Denisova* || $pop =~ ^(ESN|GWD|LWK|MSL|YRI)$ ]] || continue
    run_population "$u" "$ref" "$pop" "$sample_list" "$gt" "$tmp" & pids+=("$!")
    if (( ${#pids[@]} >= job_in_chr )); then wait "${pids[0]}" || status=1; pids=("${pids[@]:1}"); fi
  done < "$dirout/samples/$u/populations.tsv"
  for pid in "${pids[@]}"; do wait "$pid" || status=1; done
  rm -f "$gt" "$gt.part.gz"
  (( status == 0 ))
}

ref_calls_complete(){
  local u=$1 ref=$2 pop _super _n sample_list count=0
  while IFS=$'\t' read -r pop _super _n sample_list; do
    [[ $pop != population ]] || continue
    [[ $export_denisovan == 1 || $background_filter == 0 || $ref != Denisova* || $pop =~ ^(ESN|GWD|LWK|MSL|YRI)$ ]] || continue
    gzip_ok "$dirout/raw/$u/$ref.$pop.raw.txt.gz" || return 1
    count=$((count+1))
  done < "$dirout/samples/$u/populations.tsv"
  (( count > 0 ))
}

cleanup_unit(){
  local ref
  rm -rf -- "$tmp"
  for ref in "${refs[@]}"; do
    rm -f -- "$dirout/genotype/$ref.$u.gt.txt.gz" "$dirout/genotype/$ref.$u.gt.txt.gz.part.gz"
  done
}

run_unit(){
  local u=$1 tmp ref rc=0
  trap 'rc=$?; log "ERROR unit=$u rc=$rc command=$BASH_COMMAND"; exit "$rc"' ERR
  tmp=$dirout/tmp/$u
  trap cleanup_unit EXIT
  trap 'exit 129' HUP
  trap 'exit 130' INT
  trap 'exit 143' TERM
  log "START unit=$u chr=${unit_chr[$u]} region=$(unit_region "$u") sex=${unit_sex[$u]}"
  rm -rf "$tmp"; mkdir -p "$tmp"
  log "PREPARE unit=$u modern_vcf=streamed genotype=temporary_compressed"
  for ref in "${refs[@]}"; do run_ref_unit "$u" "$ref" "$tmp"; done
  local -a core_args=()
  if (( unit_core_end[$u] > 0 )); then core_args=(--core-start "$((unit_core_start[$u]-1))" --core-end "${unit_core_end[$u]}"); fi
  python3 "$helper" finalize --raw-dir "$dirout/raw/$u" --populations "$dirout/samples/$u/populations.tsv" \
    --output "$dirout/segments/$u.segments.tsv.gz" --chrom "${unit_chr[$u]}" --build "$genome_build" --locus "${unit_locus[$u]}" \
    --min-bp "$len_cut" --lod "$lod_cut" "${core_args[@]}" "${background_args[@]}" "${denisovan_args[@]}" --refs "${refs[@]}"
  cleanup_unit
  trap - ERR EXIT HUP INT TERM
}
run_scan(){
  local status=0 u final pid
  local -a pids=()
  validate_inputs
  # Resolve and validate local masks before replacing any old calls.
  # Their contents participate in provenance, so changed masks invalidate raw
  # calls as well as genotypes. A missing input leaves prior results intact.
  for u in "${analysis_units[@]}"; do prepare_unit_masks "$u"; done
  check_run_provenance
  final=$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz
  if completed_output_ok; then
    log "SKIP reason=output_exists final=$final"
    return 0
  fi
  rm -f "$dirout/.complete"
  for u in "${analysis_units[@]}"; do
    run_unit "$u" & pids+=("$!")
    if (( ${#pids[@]} >= job_of_chr )); then wait "${pids[0]}" || status=1; pids=("${pids[@]:1}"); fi
  done
  for pid in "${pids[@]}"; do wait "$pid" || status=1; done
  (( status == 0 )) || { echo "ERROR: one or more IBDmix units failed" >&2; exit 1; }
  mapfile -t files < <(find "$dirout/segments" -type f -name '*.segments.tsv.gz' | sort)
  (( ${#files[@]} )) || { echo "ERROR: no IBDmix segment files" >&2; exit 1; }
  gzip -dc "${files[@]}" | awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} $1!="ID"{print}' | gzip -c > "$final.part"
  mv "$final.part" "$final"
  touch "$dirout/.complete"
  log "ALL DONE: final=$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz"
}

if [[ $action == ibdmix_check ]]; then validate_inputs; exit 0; fi
run_scan
