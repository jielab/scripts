#!/usr/bin/env bash
# One-time repair for chrX rows produced with an autosome-only HapMap3 list.
# Rebuilds chrX from the original GRCh38 GWAS, then adds chrX MAGMA and
# clump/COJO results without recalculating or replacing chromosomes 1-22.

set -euo pipefail
export LC_ALL=C

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
gwas_post="$script_dir/gwas_post.sh"
hm3=/mnt/d/data.BIG/refGen/hm3/hapmap3_r3.snp
hm3_chr36=/mnt/d/data.BIG/refGen/hm3/hapmap3_r3_grch36.snplist
hm3_pos=/mnt/d/data.BIG/refGen/hm3/hapmap3_r3_grch38.snplist
raw_base=/mnt/g/gwas
project_base=/mnt/d/data.BIG/gwas
work_root=/mnt/d/data.BIG/gwas/.fix_X_work
magma_ref=/mnt/d/data.BIG/refLD/magma/g1000_eur
fix_version=fix_X_v3_magma
labels=all
gwas_filter=
jobs=4
apply=FALSE

usage() {
  cat <<'USAGE'
Usage: ./fix_X.sh [options]

  --label prot|met|all  project(s) to repair [all]
  --gwas LIST           comma-separated traits [all]
  --jobs INT            parallel raw scans / MAGMA / lead jobs [4]
  --apply TRUE|FALSE    TRUE atomically replaces real chrX, MAGMA, and lead tables [FALSE]
  --work-root PATH      resumable candidates, audit, and chrX staging

Examples:
  ./fix_X.sh --label prot --gwas A1BG --apply FALSE --work-root /mnt/d/scripts/0data/_fix_X_test
  ./fix_X.sh --label all --jobs 8 --apply TRUE
USAGE
}

need_value() {
  [[ -n "${2-}" && "${2-}" != --* ]] || { echo "ERROR: missing value for $1" >&2; exit 2; }
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --label) need_value "$1" "${2-}"; labels=${2,,}; shift 2 ;;
    --gwas) need_value "$1" "${2-}"; gwas_filter=$2; shift 2 ;;
    --jobs) need_value "$1" "${2-}"; jobs=$2; shift 2 ;;
    --apply) need_value "$1" "${2-}"; apply=${2^^}; shift 2 ;;
    --work-root) need_value "$1" "${2-}"; work_root=${2%/}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ "$labels" == all || "$labels" == prot || "$labels" == met ]] || {
  echo "ERROR: --label must be prot, met, or all" >&2; exit 2;
}
[[ "$apply" == TRUE || "$apply" == FALSE ]] || {
  echo "ERROR: --apply must be TRUE or FALSE" >&2; exit 2;
}
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --jobs must be a positive integer" >&2; exit 2; }
[[ -s "$gwas_post" ]] || { echo "ERROR: missing gwas_post.sh: $gwas_post" >&2; exit 1; }
[[ -s "$hm3" ]] || { echo "ERROR: missing HapMap3 SNP list: $hm3" >&2; exit 1; }
[[ -s "$hm3_chr36" ]] || { echo "ERROR: missing original GRCh36 HapMap3 list: $hm3_chr36" >&2; exit 1; }
[[ -s "$hm3_pos" ]] || { echo "ERROR: missing GRCh38 HapMap3 coordinate list: $hm3_pos" >&2; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip is required" >&2; exit 127; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix is required" >&2; exit 127; }
[[ -s "$magma_ref.bed" && -s "$magma_ref.bim" && -s "$magma_ref.fam" ]] || {
  echo "ERROR: missing MAGMA LD reference: $magma_ref.{bed,bim,fam}" >&2; exit 1;
}

if [[ "$labels" == all ]]; then
  label_list=(prot met)
else
  label_list=("$labels")
fi

mkdir -p "$work_root/audit" "$work_root/candidates" "$work_root/stage_magma" "$work_root/stage_lead"
printf 'fix_X: apply=%s labels=%s jobs=%s hm3=%s hm3-pos=%s work=%s\n' \
  "$apply" "${label_list[*]}" "$jobs" "$hm3" "$hm3_pos" "$work_root" >&2

# Reading two 1.46-million-row reference files once per trait would dominate a
# 3,000-trait repair.  Build one resumable chrX-only ID/POS cache up front.
hm3_x="$work_root/reference/hapmap3_r3_grch38.X.tsv"
hm3_x_version=fix_X_hm3_cache_v3
mkdir -p "$(dirname "$hm3_x")"
if [[ ! -s "$hm3_x" || "$hm3" -nt "$hm3_x" || "$hm3_chr36" -nt "$hm3_x" || "$hm3_pos" -nt "$hm3_x" ]] || \
   [[ "$(head -n 1 "$hm3_x" 2>/dev/null || true)" != "$hm3_x_version" ]]; then
  hm3_x_tmp="$hm3_x.tmp.$$"
  printf '%s\n' "$hm3_x_version" > "$hm3_x_tmp"
  awk -v OFS='\t' '
    function normchr(x){x=toupper(x);gsub(/^CHR/,"",x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="XY")x=25;return x}
    FNR==NR{
      gsub(/\r/,"");split($0,a,/[ \t]+/)
      if(a[1]!=""&&a[1]!="SNP")valid_id[a[1]]=1
      next
    }
    {
      gsub(/\r/,"");split($0,a,/[ \t]+/)
      if((normchr(a[1])==23||normchr(a[1])==25)&&(a[2] in valid_id))print "ID",a[2]
    }
  ' "$hm3" "$hm3_chr36" >> "$hm3_x_tmp"
  awk -v OFS='\t' '
    function normchr(x){x=toupper(x);gsub(/^CHR/,"",x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="XY")x=25;return x}
    {
      gsub(/\r/,"");split($0,a,/[ \t]+/)
      if((normchr(a[1])==23||normchr(a[1])==25)&&a[2]!=""&&a[4]~/^[0-9]+$/&&a[4]+0>0)print "POS",a[4]+0,a[2]
    }
  ' "$hm3_pos" >> "$hm3_x_tmp"
  [[ -s "$hm3_x_tmp" ]] || { echo "ERROR: no chrX rows in HapMap3 references" >&2; exit 1; }
  mv -f -- "$hm3_x_tmp" "$hm3_x"
fi
printf 'fix_X: chrX HM3 cache IDs/POS=%s file=%s\n' \
  "$(awk -F '\t' '{n[$1]++}END{print n["ID"]+0 "/" n["POS"]+0}' "$hm3_x")" "$hm3_x" >&2

selected_trait() {
  local trait="$1"
  [[ -z "$gwas_filter" || ",$gwas_filter," == *",$trait,"* ]]
}

# find(1) -path globs allow '*' to cross '/', and the project root itself is
# named "gwas".  Validate each exact three-level layout explicitly so qc/pgs
# gzip files can never be mistaken for <trait>/gwas/<trait>.gz.
list_final_gwas() {
  local project="$1" path parent subdir trait_parent trait base
  while IFS= read -r -d '' path; do
    parent=${path%/*}; subdir=${parent##*/}
    trait_parent=${parent%/*}; trait=${trait_parent##*/}; base=${path##*/}
    [[ "$subdir" == gwas && "$base" == "$trait.gz" ]] || continue
    printf '%s\0' "$path"
  done < <(find "$project/common" -mindepth 3 -maxdepth 3 -type f -name '*.gz' -print0)
}

has_final_gwas() {
  IFS= read -r -d '' < <(list_final_gwas "$1")
}

build_x_candidate() {
  set -euo pipefail
  local label="$1" final="$2" trait raw candidate audit done_dir done_file final_done=FALSE
  local tmp_dir raw_x sorted_x candidate_tmp old_x=NA sig_x kept_x magma_stage_gwas lead_stage_gwas
  trait=$(basename "$final" .gz)
  selected_trait "$trait" || return 0
  raw="$raw_base/$label/raw/$trait.gz"
  candidate="$work_root/candidates/$label/$trait.gz"
  audit="$work_root/audit/$label/$trait.tsv"
  done_dir="$work_root/audit/$label"
  done_file="$done_dir/$trait.final.done"
  mkdir -p "$(dirname "$candidate")" "$done_dir"

  if [[ "$apply" == TRUE && -s "$done_file" ]] && awk -F '\t' -v v="$fix_version" 'NR==2&&$7==v{ok=1}END{exit !ok}' "$done_file"; then
    printf 'fix_X: REUSE repaired final %s/%s\n' "$label" "$trait" >&2
    final_done=TRUE
  fi
  [[ -s "$raw" ]] || { echo "ERROR: missing raw GWAS: $raw" >&2; return 1; }
  [[ -s "$final" ]] || { echo "ERROR: missing formatted GWAS: $final" >&2; return 1; }

  tmp_dir="$work_root/tmp/$label/$trait.$$"
  mkdir -p "$tmp_dir"
  raw_x="$tmp_dir/raw_x.tsv"
  sorted_x="$tmp_dir/sorted_x.tsv"

  if [[ ! -s "$candidate" || ! -s "$audit" || "$raw" -nt "$candidate" || "$hm3_x" -nt "$candidate" ]]; then
    printf 'fix_X: extract %s/%s from %s\n' "$label" "$trait" "$raw" >&2
    gzip -cd -- "$raw" | awk -v FS='\t' -v OFS='\t' \
      -v hm3_x_file="$hm3_x" -v p_small='1e-4' -v p_lead='5e-8' -v audit="$audit.tmp.$$" '
      function clean(x){gsub(/\r/,"",x);return x}
      function get(c,x){x=(c>0?$c:"");return clean(x)}
      function value(c,x){x=get(c);return x==""?"NA":x}
      function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
      function normchr(x){x=toupper(clean(x));gsub(/^CHR/,"",x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="XY")x=25;return x}
      BEGIN{
        while((getline x<hm3_x_file)>0){
          x=clean(x);split(x,a,/[ \t]+/)
          if(a[1]=="ID"&&a[2]!="")hm3[a[2]]=1
          else if(a[1]=="POS"&&a[2]~/^[0-9]+$/&&a[2]+0>0&&a[3]!=""&&!((a[2]+0) in hm3_pos))hm3_pos[a[2]+0]=a[3]
        }
        close(hm3_x_file)
      }
      NR==1{
        for(i=1;i<=NF;i++){
          h=tolower(clean($i));gsub(/^#/,"",h)
          if(h=="snp")snp_col=i;else if(h=="rsid"&&snp_col==0)snp_col=i
          if(h=="chr"||h=="chrom"||h=="chromosome")chr_col=i
          if(h=="pos"||h=="bp"||h=="position"||h=="base_pair_location")pos_col=i
          if(h=="ea"||h=="a1"||h=="allele1"||h=="effect_allele")ea_col=i
          if(h=="nea"||h=="a2"||h=="allele0"||h=="other_allele")nea_col=i
          if(h=="eaf"||h=="a1freq"||h=="effect_allele_frequency")eaf_col=i
          if(h=="n"||h=="neff")n_col=i
          if(h=="beta"||h=="effect")beta_col=i
          if(h=="se"||h=="standard_error")se_col=i
          if(h=="p"||h=="pval"||h=="p_value")p_col=i
          if(h=="log10p"||h=="neg_log_10_p_value")logp_col=i
        }
        if(!snp_col||!chr_col||!pos_col||!ea_col||!nea_col||!beta_col||!se_col||(!p_col&&!logp_col)){
          print "ERROR: unsupported raw header" > "/dev/stderr"; exit 2
        }
        print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P"
        next
      }
      {
        chr=normchr(get(chr_col));if(chr!=23)next
        raw_n++
        pos=get(pos_col);if(pos!~/^[0-9]+$/||pos+0<1){bad_coord++;next}
        raw_snp=get(snp_col);snp=raw_snp;ea=get(ea_col);nea=get(nea_col)
        in_hm3_id=(raw_snp in hm3)
        in_hm3_pos=((pos+0) in hm3_pos)
        if(snp==""||snp=="NA"||snp=="."){
          if(in_hm3_pos)snp=hm3_pos[pos+0]
          else{snp=chr ":" pos;if(ea!="")snp=snp ":" ea;if(nea!="")snp=snp ":" nea}
        }
        p=(p_col?get(p_col):"");lp=(logp_col?get(logp_col):"")
        if(p==""&&isnum(lp))p=10^(-lp);if(lp==""&&isnum(p)&&p>0)lp=-log(p)/log(10)
        in_hm3=(in_hm3_id||in_hm3_pos)
        in_p=(isnum(p)&&p+0<=p_small+0)
        if(!in_hm3&&!in_p)next
        if(ea=="")ea="NA";if(nea=="")nea="NA";if(p=="")p="NA";if(lp=="")lp="NA"
        print snp,23,pos+0,ea,nea,value(eaf_col),value(n_col),value(beta_col),value(se_col),p,lp
        kept++;hm3_id_n+=in_hm3_id;hm3_pos_n+=in_hm3_pos;hm3_n+=in_hm3;p_n+=in_p;if(isnum(p)&&p+0<=p_lead+0)sig_n++
      }
      END{
        print "RAW_X\tKEPT_X\tHM3_ID_X\tHM3_POS_X\tHM3_X\tP_SMALL_X\tP_LEAD_X\tBAD_COORD_X" > audit
        print raw_n+0 "\t" kept+0 "\t" hm3_id_n+0 "\t" hm3_pos_n+0 "\t" hm3_n+0 "\t" p_n+0 "\t" sig_n+0 "\t" bad_coord+0 >> audit
      }
    ' > "$raw_x"
    {
      IFS= read -r header
      printf '%s\n' "$header"
      sort -T "$tmp_dir" -S 512M -t $'\t' -k3,3n -k1,1 -k4,4 -k5,5 -u
    } < "$raw_x" > "$sorted_x"
    candidate_tmp="$candidate.tmp.$$"
    bgzip -@ 1 -c "$sorted_x" > "$candidate_tmp"
    bgzip -t "$candidate_tmp"
    tabix -f -s 2 -b 3 -e 3 -S 1 "$candidate_tmp"
    mv -f -- "$candidate_tmp" "$candidate"
    mv -f -- "$candidate_tmp.tbi" "$candidate.tbi"
    mv -f -- "$audit.tmp.$$" "$audit"
  fi

  kept_x=$(awk -F '\t' 'NR==2{print $2+0}' "$audit")
  sig_x=$(awk -F '\t' 'NR==2{print $7+0}' "$audit")
  (( kept_x > 0 )) || { echo "ERROR: no replacement chrX rows for $label/$trait" >&2; return 1; }
  if [[ -s "$final.tbi" || -s "$final.csi" ]]; then
    old_x=$(tabix "$final" 23 2>/dev/null | awk 'END{print NR+0}')
  fi

  if [[ "$apply" == TRUE && "$final_done" != TRUE ]]; then
    local merged_tmp="$final.fix_X.tmp.$$"
    awk -v FS='\t' -v OFS='\t' '
      NR==FNR{if(FNR>1)x[++nx]=$0;next}
      FNR==1{print;next}
      {
        ch=toupper($2);gsub(/^CHR/,"",ch);if(ch=="X")ch=23;else if(ch=="Y")ch=24;else if(ch=="MT"||ch=="M")ch=25
        if(!inserted&&ch+0>=23){for(i=1;i<=nx;i++)print x[i];inserted=1}
        if(ch+0!=23)print
      }
      END{if(!inserted)for(i=1;i<=nx;i++)print x[i]}
    ' <(bgzip -cd "$candidate") <(gzip -cd "$final") | bgzip -@ 1 -c > "$merged_tmp"
    bgzip -t "$merged_tmp"
    tabix -f -s 2 -b 3 -e 3 -S 1 "$merged_tmp"
    chmod --reference="$final" "$merged_tmp" 2>/dev/null || true
    rm -f -- "$final.tbi" "$final.csi"
    mv -f -- "$merged_tmp" "$final"
    mv -f -- "$merged_tmp.tbi" "$final.tbi"
    printf 'GWAS\tOLD_X\tNEW_X\tP_LEAD_X\tHM3\tHM3_POS\tFIX_VERSION\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$trait" "$old_x" "$kept_x" "$sig_x" "$hm3" "$hm3_pos" "$fix_version" > "$done_file.tmp.$$"
  fi

  magma_stage_gwas="$work_root/stage_magma/$label/common/$trait/gwas"
  mkdir -p "$magma_stage_gwas"
  cp -f -- "$candidate" "$magma_stage_gwas/$trait.gz"
  cp -f -- "$candidate.tbi" "$magma_stage_gwas/$trait.gz.tbi"
  if (( sig_x > 0 )); then
    lead_stage_gwas="$work_root/stage_lead/$label/common/$trait/gwas"
    mkdir -p "$lead_stage_gwas"
    cp -f -- "$candidate" "$lead_stage_gwas/$trait.gz"
    cp -f -- "$candidate.tbi" "$lead_stage_gwas/$trait.gz.tbi"
  fi
  if [[ "$apply" == TRUE && "$final_done" != TRUE ]]; then
    mv -f -- "$done_file.tmp.$$" "$done_file"
  fi
  rm -f -- "$raw_x" "$sorted_x"
  rmdir -- "$tmp_dir"
  printf 'fix_X: %s/%s old_X=%s new_X=%s lead_X=%s apply=%s\n' \
    "$label" "$trait" "$old_x" "$kept_x" "$sig_x" "$apply" >&2
}

export -f selected_trait build_x_candidate
export hm3 hm3_chr36 hm3_pos hm3_x raw_base project_base work_root magma_ref fix_version gwas_filter apply

for label in "${label_list[@]}"; do
  project="$project_base/$label"
  raw_dir="$raw_base/$label/raw"
  [[ -d "$project/common" ]] || { echo "ERROR: missing project common directory: $project/common" >&2; exit 1; }
  [[ -d "$raw_dir" ]] || { echo "ERROR: missing mounted raw directory: $raw_dir" >&2; exit 1; }
  printf 'fix_X: scan project %s\n' "$label" >&2
  while IFS= read -r -d '' final; do
    trait=$(basename "$final" .gz)
    if selected_trait "$trait"; then printf '%s\0' "$final"; fi
  done < <(list_final_gwas "$project" | sort -z) | \
    xargs -0 -r -n 1 -P "$jobs" bash -c "build_x_candidate \"\$1\" \"\$2\"" _ "$label"
done

if [[ "$apply" != TRUE ]]; then
  printf 'fix_X: dry run complete; candidates and audits retained under %s\n' "$work_root" >&2
  exit 0
fi

filter_args=()
[[ -z "$gwas_filter" ]] || filter_args=(--gwas "$gwas_filter")

# The staged inputs contain chrX only. MAGMA therefore computes only chrX gene
# rows even though the shared LD reference/annotation also contains autosomes.
for label in "${label_list[@]}"; do
  magma_project="$work_root/stage_magma/$label"
  if has_final_gwas "$magma_project"; then
    printf 'fix_X: run staged chrX MAGMA for %s\n' "$label" >&2
    "$gwas_post" magma --dir-out "$magma_project" --label "$label" --category common \
      --grch 38 --magma-ref "$magma_ref" --magma-annot-cache "$work_root/reference/magma_annotation" \
      --replace TRUE --run-cmd TRUE --foreground TRUE --jobs "$jobs" "${filter_args[@]}"
  else
    echo "ERROR: no staged chrX MAGMA inputs for $label" >&2
    exit 1
  fi
done

# Run the existing, reference-matching lead implementation only on traits with
# a genome-wide significant chrX association.  Their staged GWAS files contain
# chrX only, so --chr X cannot overwrite or recalculate autosomal results.
for label in "${label_list[@]}"; do
  stage_project="$work_root/stage_lead/$label"
  if has_final_gwas "$stage_project"; then
    printf 'fix_X: run staged chrX clump/COJO for %s\n' "$label" >&2
    "$gwas_post" lead --dir-out "$stage_project" --label "$label" --category common \
      --grch 38 --chr X --replace TRUE --run-cmd TRUE --foreground TRUE --jobs "$jobs" "${filter_args[@]}"
  else
    printf 'fix_X: no genome-wide significant chrX traits for %s\n' "$label" >&2
  fi
done

merge_x_table() {
  local old="$1" new="$2" tmp header chr_col
  [[ -s "$old" || -s "$new" ]] || return 0
  if [[ -s "$old" ]]; then header=$(head -n 1 "$old"); else header=$(head -n 1 "$new"); fi
  tmp="$old.fix_X.tmp.$$"
  {
    printf '%s\n' "$header"
    if [[ -s "$old" ]]; then
      chr_col=$(awk -v FS='\t' 'NR==1{for(i=1;i<=NF;i++){x=toupper($i);gsub(/^#/,"",x);if(x=="CHR"||x=="CHROM"){print i;exit}}}' "$old")
      [[ -n "$chr_col" ]] || { echo "ERROR: no CHR/CHROM column in lead table: $old" >&2; return 1; }
      awk -v FS='\t' -v c="$chr_col" 'NR>1{x=toupper($c);gsub(/^CHR/,"",x);if(x!="X"&&x!="23")print}' "$old"
    fi
    if [[ -s "$new" ]]; then tail -n +2 "$new"; fi
  } > "$tmp"
  mv -f -- "$tmp" "$old"
}

merge_x_magma() {
  local old="$1" new="$2" kind="$3" tmp n_x
  [[ -s "$new" ]] || { echo "ERROR: missing staged chrX MAGMA output: $new" >&2; return 1; }
  case "$kind" in
    out)
      n_x=$(awk 'NR>1{x=toupper($2);gsub(/^CHR/,"",x);if(x=="X"||x=="23")n++}END{print n+0}' "$new")
      (( n_x > 0 )) || { echo "ERROR: staged MAGMA genes.out has no chrX rows: $new" >&2; return 1; }
      tmp="$old.fix_X.tmp.$$"
      {
        if [[ -s "$old" ]]; then head -n 1 "$old"; else head -n 1 "$new"; fi
        if [[ -s "$old" ]]; then
          awk 'NR>1{x=toupper($2);gsub(/^CHR/,"",x);if(x!="X"&&x!="23")print}' "$old"
        fi
        tail -n +2 "$new"
      } > "$tmp"
      ;;
    raw)
      n_x=$(awk '!/^#/{x=toupper($2);gsub(/^CHR/,"",x);if(x=="X"||x=="23")n++}END{print n+0}' "$new")
      (( n_x > 0 )) || { echo "ERROR: staged MAGMA genes.raw has no chrX rows: $new" >&2; return 1; }
      tmp="$old.fix_X.tmp.$$"
      {
        if [[ -s "$old" ]]; then awk '/^#/{print}' "$old"; else awk '/^#/{print}' "$new"; fi
        if [[ -s "$old" ]]; then
          awk '!/^#/{x=toupper($2);gsub(/^CHR/,"",x);if(x!="X"&&x!="23")print}' "$old"
        fi
        awk '!/^#/{print}' "$new"
      } > "$tmp"
      ;;
    *) echo "ERROR: invalid MAGMA merge type: $kind" >&2; return 2 ;;
  esac
  mv -f -- "$tmp" "$old"
}

for label in "${label_list[@]}"; do
  project="$project_base/$label"
  while IFS= read -r -d '' final; do
    trait=$(basename "$final" .gz)
    selected_trait "$trait" || continue
    real_dir=$(dirname "$final")
    lead_stage_dir="$work_root/stage_lead/$label/common/$trait/gwas"
    magma_stage_dir="$work_root/stage_magma/$label/common/$trait/magma"
    merge_x_magma "$project/common/$trait/magma/$trait.genes.out" "$magma_stage_dir/$trait.genes.out" out
    merge_x_magma "$project/common/$trait/magma/$trait.genes.raw" "$magma_stage_dir/$trait.genes.raw" raw
    merge_x_table "$real_dir/$trait.awk.snp" "$lead_stage_dir/$trait.awk.snp"
    merge_x_table "$real_dir/$trait.clumps" "$lead_stage_dir/$trait.clumps"
    merge_x_table "$real_dir/$trait.jma.cojo" "$lead_stage_dir/$trait.jma.cojo"
    printf 'GWAS\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n%s\t5e-8\t1000000\tautosome,X\tcomplete_fix_X\n' "$trait" \
      > "$real_dir/$trait.lead.done.tmp.$$"
    mv -f -- "$real_dir/$trait.lead.done.tmp.$$" "$real_dir/$trait.lead.done"
    if [[ -s "$real_dir/$trait.clumps" ]]; then
      printf 'GWAS\tPHASE\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n%s\tclump\t5e-8\t1000000\tautosome,X\tcomplete_fix_X\n' "$trait" \
        > "$real_dir/$trait.clump.done"
    fi
    if [[ -s "$real_dir/$trait.jma.cojo" || -s "$real_dir/$trait.ldr.cojo" ]]; then
      printf 'GWAS\tPHASE\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n%s\tcojo\t5e-8\t1000000\tautosome,X\tcomplete_fix_X\n' "$trait" \
        > "$real_dir/$trait.cojo.done"
    fi
  done < <(list_final_gwas "$project" | sort -z)
done

printf 'fix_X: repair complete; audit and staging retained under %s\n' "$work_root" >&2
