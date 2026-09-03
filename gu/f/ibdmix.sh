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
for cmd in awk bcftools gzip python3 sort comm find stat cmp grep; do command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing command: $cmd" >&2; exit 1; }; done

dir0=${dir0:-/mnt/d}
dir_ref=${dir_ref:-$dir0/data.BIG/refGen}
dirmod=${dirmod:-${GU_TARGET_ROOT:-$dir_ref/1kg/${GRCH:-37}}}
dirarch=${dirarch:-${GU_ARCHAIC_ROOT:-$dir_ref/archaic/${GRCH:-37}/vcf}}
sample_file=${sample_file:-$dirmod/samples.txt}
target_vcf_dir=${GU_TARGET_VCF_DIR:-$dirmod/vcf}
dirsoft=${dirsoft:-$dir0/software/gu/IBDmix}
dirout=${dirout:-${GU_ANALYSIS_ROOT:-$dir0/analysis/gu}/ibdmix/${GU_SCOPE_ID:-genome}}
genome_build=${genome_build:-b${GRCH:-37}}
loci_file=${GU_LOCI_FILE:-}
chrs_arg=${chrs:-${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}}
refs_arg=${refs:-"Altai Chagyr Vindija Denisova"}
read -r -a refs <<< "$refs_arg"
lod_cut=${lod_cut:-4}
len_cut=${len_cut:-50000}
emit_lod_cut=${emit_lod_cut:-$lod_cut}
minor_allele_count=${minor_allele_count:-1}
archaic_error=${archaic_error:-0.01}
modern_error_max=${modern_error_max:-0.002}
modern_error_proportion=${modern_error_proportion:-2}
# Keep analysis units serial by default to limit temporary I/O and memory.
job_of_chr=${job_of_chr:-1}
job_in_chr=${job_in_chr:-1} # retained as a configuration field; there is one ALL-sample call per unit/ref.
IBDMIX_LOCUS_FLANK_BP=${IBDMIX_LOCUS_FLANK_BP:-100000}
IBDMIX_REPLACE=${IBDMIX_REPLACE:-0}
GU_CHRX_PAR_DIPLOID=${GU_CHRX_PAR_DIPLOID:-0}

generate_gt=$dirsoft/build/src/generate_gt
ibdmix_bin=$dirsoft/build/src/ibdmix
summary_sh=$dirsoft/src/summary.sh
[[ $genome_build == b37 || $genome_build == b38 ]] || { echo "ERROR: genome_build must be b37 or b38" >&2; exit 1; }
[[ $IBDMIX_LOCUS_FLANK_BP =~ ^[0-9]+$ ]] || { echo "ERROR: IBDMIX_LOCUS_FLANK_BP must be a non-negative integer" >&2; exit 1; }
[[ $IBDMIX_REPLACE == 0 || $IBDMIX_REPLACE == 1 ]] || { echo "ERROR: internal replace-ibdmix state must be 0 or 1" >&2; exit 1; }
[[ $job_of_chr =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: job_of_chr must be a positive integer" >&2; exit 1; }

final_output=$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz
if [[ $action == ibdmix_run && $IBDMIX_REPLACE == 0 && -s $final_output ]]; then
  printf '[IBDmix] SKIP reason=output_exists final=%s\n' "$final_output"
  exit 0
fi

for file in "$generate_gt" "$ibdmix_bin" "$summary_sh" "$F/vcf_gt_fix.py" "$sample_file"; do [[ -s $file ]] || { echo "ERROR: missing required file: $file" >&2; exit 1; }; done
[[ -d $target_vcf_dir ]] || { echo "ERROR: missing target VCF directory: $target_vcf_dir" >&2; exit 1; }
mkdir -p "$dirout"/{samples,genotype,raw,segments,final,log,tmp}
if [[ $action == ibdmix_run && $IBDMIX_REPLACE == 1 ]]; then
  rm -rf -- "$dirout/genotype" "$dirout/raw" "$dirout/segments" "$dirout/final" "$dirout/tmp"
  mkdir -p "$dirout"/{genotype,raw,segments,final,tmp}
fi
exec > >(tee "$dirout/ibdmix.log") 2>&1
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
gzip_ok(){ [[ -s $1 ]] && gzip -t "$1" >/dev/null 2>&1; }
data_rows(){ gzip -dc "$1" | awk 'NR>1{n++} END{print n+0}'; }

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
  local needs_x=0 u
  for u in "${analysis_units[@]}"; do [[ ${unit_chr[$u]} != X || ${unit_sex[$u]} != male ]] || needs_x=1; done
  : > "$dirout/samples/ALL.txt"; : > "$dirout/samples/male.txt"
  awk -v needs_x="$needs_x" 'BEGIN{FS="[[:space:]]+"}
    NR==1{for(i=1;i<=NF;i++){h[tolower($i)]=i}; s=h["sample"]; x=h["sex"];
      if(!s || (needs_x && !x)){print "ERROR: samples.txt needs sample, plus sex for chrX" > "/dev/stderr"; exit 2}; next}
    {id=$s; if(id=="")next; print id > all; if(!x)next; sex=tolower($x);
      if(sex=="1"||sex=="m"||sex=="male") print id > male;
      else if(sex=="2"||sex=="f"||sex=="female") next;
      else if(needs_x){print "ERROR: invalid sex for " id > "/dev/stderr"; bad=1}}
    END{exit bad}' all="$dirout/samples/ALL.txt" male="$dirout/samples/male.txt" "$sample_file"
  [[ -s $dirout/samples/ALL.txt ]] || { echo "ERROR: failed to create sample list" >&2; exit 1; }
  if (( needs_x )); then
    [[ -s $dirout/samples/male.txt ]] || { echo "ERROR: chrX.male analysis needs at least one sample marked male" >&2; exit 1; }
  fi
}
unit_has_samples(){
  case ${unit_sex[$1]} in
    male) [[ -s $dirout/samples/male.txt ]] ;;
    *) [[ -s $dirout/samples/ALL.txt ]] ;;
  esac
}
validate_vcf_samples(){
  local vcf=$1 expected=$2 work=$dirout/tmp/sample_check
  mkdir -p "$work"; bcftools query -l "$vcf" | sort -u > "$work/vcf.txt"; sort -u "$expected" > "$work/expected.txt"
  comm -23 "$work/expected.txt" "$work/vcf.txt" > "$work/missing.txt"
  comm -13 "$work/expected.txt" "$work/vcf.txt" > "$work/unexpected.txt"
  [[ ! -s $work/missing.txt ]] || { echo "ERROR: expected sample IDs absent from $vcf (first examples):" >&2; head -10 "$work/missing.txt" >&2; return 1; }
  [[ ! -s $work/unexpected.txt ]] || { echo "ERROR: VCF sample IDs absent from the selected sample metadata (first examples):" >&2; head -10 "$work/unexpected.txt" >&2; return 1; }
}
validate_inputs(){
  local u chr ref vcf expected
  declare -A seen=()
  prepare_samples
  for u in "${analysis_units[@]}"; do
    chr=${unit_chr[$u]}; [[ -z ${seen[$chr]:-} ]] || continue; seen[$chr]=1
    vcf=$(modern_vcf "$chr"); [[ -s $vcf ]] || { echo "ERROR: missing modern VCF: $vcf" >&2; return 1; }
    expected=$dirout/samples/ALL.txt
    for u in "${analysis_units[@]}"; do [[ ${unit_chr[$u]} != "$chr" || ${unit_sex[$u]} != male ]] || expected=$dirout/samples/male.txt; done
    validate_vcf_samples "$vcf" "$expected"
    for ref in "${refs[@]}"; do archaic_vcf "$ref" "$chr" >/dev/null || { echo "ERROR: missing archaic VCF ref=$ref chr=$chr" >&2; return 1; }; done
  done
  log "INPUT CHECK PASSED: metadata_samples=$(awk 'END{print NR}' "$dirout/samples/ALL.txt") units=${analysis_units[*]}"
}
run_record(){
  local u chr ref vcf index meta line x_mode=not_requested
  for u in "${analysis_units[@]}"; do
    [[ ${unit_chr[$u]} != X ]] || x_mode=$([[ ${unit_sex[$u]} == male ]] && printf male_haploid_nonpar || printf par_diploid_chrxy)
  done
  printf 'genome_build\t%s\nx_mode\t%s\nloci_flank_bp\t%s\nrefs\t%s\nlod_cut\t%s\nlen_cut\t%s\nemit_lod_cut\t%s\nminor_allele_count\t%s\narchaic_error\t%s\nmodern_error_max\t%s\nmodern_error_proportion\t%s\n' \
    "$genome_build" "$x_mode" "$IBDMIX_LOCUS_FLANK_BP" "$refs_arg" "$lod_cut" "$len_cut" "$emit_lod_cut" "$minor_allele_count" "$archaic_error" "$modern_error_max" "$modern_error_proportion"
  stat -c 'sample_file\t%n:%s:%Y' "$sample_file"
  [[ -z $loci_file ]] || stat -c 'loci_file\t%n:%s:%Y' "$loci_file"
  stat -c 'software\t%n:%s:%Y' "$generate_gt" "$ibdmix_bin" "$summary_sh" "$F/vcf_gt_fix.py"
  declare -A seen=()
  for u in "${analysis_units[@]}"; do
    chr=${unit_chr[$u]}; [[ -z ${seen[$chr]:-} ]] || continue; seen[$chr]=1
    vcf=$(modern_vcf "$chr"); meta=${GU_TARGET_TMP_DIR:?}/chr${chr}/source.tsv
    if [[ -s $meta ]]; then
      # The VCF is a disposable pfile export. Track the source pgen/pvar/psam
      # conversion contract, not the regenerated VCF/index mtimes.
      while IFS= read -r line; do printf 'modern_source\tchr%s\t%s\n' "$chr" "$line"; done < "$meta"
    else
      # A native VCF has no conversion record and remains the source input.
      stat -c 'modern_vcf\t%n:%s:%Y' "$vcf"
      for index in "$vcf.tbi" "$vcf.csi"; do [[ -s $index ]] && stat -c 'modern_index\t%n:%s:%Y' "$index"; done
    fi
    for ref in "${refs[@]}"; do vcf=$(archaic_vcf "$ref" "$chr"); stat -c 'archaic_vcf\t%n:%s:%Y' "$vcf"; done
  done
}
check_run_provenance(){
  local current=$dirout/run.meta.tsv
  [[ -s $current ]] || run_record > "$current"
}

completed_output_ok(){
  local final=$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz u ref header
  [[ -e $dirout/.complete ]] || return 1
  gzip_ok "$final" || return 1
  header=$(gzip -dc "$final" | head -1)
  [[ $header == $'ID\tchrom\tstart\tend\tlength\tslod\tsites\tpositive_lods\tnegative_lods\tsample_set\tsuper_pop\tanc\tlocus_id\tgenome_build' ]] || return 1
  for u in "${analysis_units[@]}"; do
    for ref in "${refs[@]}"; do
      gzip_ok "$dirout/genotype/$ref.$u.gt.txt.gz" || return 1
      gzip_ok "$dirout/raw/$ref.$u.raw.txt.gz" || return 1
      gzip_ok "$dirout/segments/$ref.$u.segments.tsv.gz" || return 1
    done
  done
}

prepare_modern(){
  local u=$1 out=$2 chr sex vcf region adapted
  chr=${unit_chr[$u]}; sex=${unit_sex[$u]}
  local -a sample_args=() fix=()
  vcf=$(modern_vcf "$chr"); region=$(unit_region "$u"); unit_chr_current=$chr; adapted=$(adapt_region "$vcf" "$region")
  [[ $sex == male ]] && sample_args=(-S "$dirout/samples/male.txt")
  if [[ $chr == X ]]; then
    fix=(--chrom 23)
    [[ $sex != male ]] || fix+=(--duplicate-haploid slash)
  fi
  local -a view=(view -m2 -M2 -v snps); [[ -n $adapted ]] && view+=(-r "$adapted"); view+=("${sample_args[@]}" "$vcf")
  if [[ $chr == X ]]; then bcftools "${view[@]}" | python3 "$F/vcf_gt_fix.py" "${fix[@]}" > "$out"
  else bcftools "${view[@]}" | awk 'BEGIN{OFS="\t"} /^#/{print;next}{sub(/^chr/,"",$1);print}' > "$out"
  fi
  [[ -s $out ]]
}
prepare_archaic(){
  local u=$1 ref=$2 out=$3 chr vcf region adapted
  chr=${unit_chr[$u]}
  vcf=$(archaic_vcf "$ref" "$chr"); region=$(unit_region "$u"); unit_chr_current=$chr; adapted=$(adapt_region "$vcf" "$region")
  local -a view=(view -m2 -M2 -v snps); [[ -n $adapted ]] && view+=(-r "$adapted"); view+=("$vcf")
  if [[ $chr == X ]]; then bcftools "${view[@]}" | python3 "$F/vcf_gt_fix.py" --chrom 23 --duplicate-haploid slash > "$out"
  else bcftools "${view[@]}" | awk 'BEGIN{OFS="\t"} /^#/{print;next}{sub(/^chr/,"",$1);print}' > "$out"
  fi
  [[ -s $out ]]
}

run_ref_unit(){
  local u=$1 ref=$2 tmp=$3 gt raw raw_txt seg modern archaic sample_list outchr
  gt=$dirout/genotype/$ref.$u.gt.txt.gz; raw=$dirout/raw/$ref.$u.raw.txt.gz; seg=$dirout/segments/$ref.$u.segments.tsv.gz
  modern=$tmp/modern.vcf; archaic=$tmp/$ref.vcf
  if ! gzip_ok "$gt"; then
    prepare_archaic "$u" "$ref" "$archaic"
    "$generate_gt" -a "$archaic" -m "$modern" -o "$tmp/$ref.gt.txt" > "$dirout/log/$ref.$u.generate_gt.log" 2>&1
    gzip -c "$tmp/$ref.gt.txt" > "$gt"; rm -f "$tmp/$ref.gt.txt" "$archaic"
  fi
  sample_list=$dirout/samples/ALL.txt; [[ ${unit_sex[$u]} == male ]] && sample_list=$dirout/samples/male.txt
  if ! gzip_ok "$raw"; then
    raw_txt=$tmp/$ref.raw.txt
    "$ibdmix_bin" --genotype <(gzip -dc "$gt") --output "$raw_txt" --sample "$sample_list" --LOD-threshold "$emit_lod_cut" \
      --minor-allele-count-threshold "$minor_allele_count" --archaic-error "$archaic_error" --modern-error-max "$modern_error_max" \
      --modern-error-proportion "$modern_error_proportion" --more-stats > "$dirout/log/$ref.$u.ibdmix.log" 2>&1
    gzip -c "$raw_txt" > "$raw"; rm -f "$raw_txt"
  fi
  outchr=${unit_chr[$u]}
  gzip -dc "$raw" | bash "$summary_sh" "$len_cut" "$lod_cut" ALL - - 2> "$dirout/log/$ref.$u.segment.log" |
    awk -v ref="$ref" -v outchr="$outchr" -v locus="${unit_locus[$u]}" -v core_start="${unit_core_start[$u]}" -v core_end="${unit_core_end[$u]}" -v build="$genome_build" 'BEGIN{FS=OFS="\t"}
      NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "ID","chrom","start","end","length","slod","sites","positive_lods","negative_lods","sample_set","super_pop","anc","locus_id","genome_build";next}
      $1!="ID"{
        id=(h["ID"]?$(h["ID"]):$1); st=(h["start"]?$(h["start"]):$2); en=(h["end"]?$(h["end"]):$3)
        if(core_start>0 && st<core_start)st=core_start; if(core_end>0 && en>core_end)en=core_end; if(en<=st)next
        len=en-st; slod=(h["slod"]?$(h["slod"]):"")
        sites=(h["sites"]?$(h["sites"]):""); pos=(h["positive_lods"]?$(h["positive_lods"]):""); neg=(h["negative_lods"]?$(h["negative_lods"]):"")
        print id,outchr,st,en,len,slod,sites,pos,neg,"ALL","ALL",ref,locus,build
      }' | gzip -c > "$seg"
  log "unit=$u chr=$outchr region=$(unit_region "$u") sex=${unit_sex[$u]} ref=$ref segments=$(data_rows "$seg")"
}

unit_needs_modern_vcf(){
  local u=$1 ref gt
  for ref in "${refs[@]}"; do
    gt=$dirout/genotype/$ref.$u.gt.txt.gz
    gzip_ok "$gt" || return 0
  done
  return 1
}

run_unit(){
  local u=$1 tmp ref
  trap 'rc=$?; log "ERROR unit=$u rc=$rc command=$BASH_COMMAND"; exit "$rc"' ERR
  tmp=$dirout/tmp/$u
  log "START unit=$u chr=${unit_chr[$u]} region=$(unit_region "$u") sex=${unit_sex[$u]}"
  rm -rf "$tmp"; mkdir -p "$tmp"
  if unit_needs_modern_vcf "$u"; then
    prepare_modern "$u" "$tmp/modern.vcf"
  else
    log "REUSE unit=$u genotype=complete modern_vcf=SKIP"
  fi
  for ref in "${refs[@]}"; do run_ref_unit "$u" "$ref" "$tmp"; done
  rm -rf "$tmp"
  trap - ERR
}
run_scan(){
  local running=0 status=0 u final
  validate_inputs
  check_run_provenance
  final=$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz
  if completed_output_ok; then
    log "SKIP reason=output_exists final=$final"
    return 0
  fi
  rm -f "$dirout/.complete"
  for u in "${analysis_units[@]}"; do
    if ! unit_has_samples "$u"; then log "SKIP unit=$u sex=${unit_sex[$u]}: no samples"; continue; fi
    run_unit "$u" & running=$((running+1))
    if (( running >= job_of_chr )); then wait -n || status=1; running=$((running-1)); fi
  done
  while (( running )); do wait -n || status=1; running=$((running-1)); done
  (( status == 0 )) || { echo "ERROR: one or more IBDmix units failed" >&2; exit 1; }
  mapfile -t files < <(find "$dirout/segments" -type f -name '*.segments.tsv.gz' | sort)
  (( ${#files[@]} )) || { echo "ERROR: no IBDmix segment files" >&2; exit 1; }
  gzip -dc "${files[@]}" | awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} $1!="ID"{print}' | gzip -c > "$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz"
  touch "$dirout/.complete"
  log "ALL DONE: final=$dirout/final/all_archaic_refs.lod${lod_cut}.len${len_cut}.segments.tsv.gz"
}

if [[ $action == ibdmix_check ]]; then validate_inputs; exit 0; fi
run_scan
