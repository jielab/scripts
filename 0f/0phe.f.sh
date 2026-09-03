# This file is a source-only library. Loading it more than once is a no-op, so
# later sources cannot replace functions that the caller defined in between.
if [[ "${PHE_F_SH_LOADED:-0}" == 1 ]]; then
  return 0 2>/dev/null || exit 0
fi
PHE_F_SH_LOADED=1


# 🚩 Shared genomic coordinates
# Defaults preserve values supplied by the caller.
: "${X_LEN_b37:=155270560}"
: "${X_LEN_b38:=156040895}"

: "${X_PAR1_START_b37:=60001}"
: "${X_PAR1_END_b37:=2699520}"
: "${X_PAR2_START_b37:=154931044}"
: "${X_PAR2_END_b37:=155260560}"
: "${X_NONPAR1_START_b37:=1}"
: "${X_NONPAR1_END_b37:=60000}"
: "${X_NONPAR2_START_b37:=2699521}"
: "${X_NONPAR2_END_b37:=154931043}"
: "${X_NONPAR3_START_b37:=155260561}"
: "${X_NONPAR3_END_b37:=155270560}"

: "${X_PAR1_START_b38:=10001}"
: "${X_PAR1_END_b38:=2781479}"
: "${X_PAR2_START_b38:=155701383}"
: "${X_PAR2_END_b38:=156030895}"
: "${X_NONPAR1_START_b38:=1}"
: "${X_NONPAR1_END_b38:=10000}"
: "${X_NONPAR2_START_b38:=2781480}"
: "${X_NONPAR2_END_b38:=155701382}"
: "${X_NONPAR3_START_b38:=156030896}"
: "${X_NONPAR3_END_b38:=156040895}"

: "${MHC_START_b37:=28477897}"
: "${MHC_END_b37:=33448354}"
: "${MHC_START_b38:=28510120}"
: "${MHC_END_b38:=33480577}"

: "${X_PAR_b37:=X:${X_PAR1_START_b37}-${X_PAR1_END_b37},X:${X_PAR2_START_b37}-${X_PAR2_END_b37}}"
: "${X_NONPAR_b37:=X:${X_NONPAR1_START_b37}-${X_NONPAR1_END_b37},X:${X_NONPAR2_START_b37}-${X_NONPAR2_END_b37},X:${X_NONPAR3_START_b37}-${X_NONPAR3_END_b37}}"

: "${X_PAR_b38:=X:${X_PAR1_START_b38}-${X_PAR1_END_b38},X:${X_PAR2_START_b38}-${X_PAR2_END_b38}}"
: "${X_NONPAR_b38:=X:${X_NONPAR1_START_b38}-${X_NONPAR1_END_b38},X:${X_NONPAR2_START_b38}-${X_NONPAR2_END_b38},X:${X_NONPAR3_START_b38}-${X_NONPAR3_END_b38}}"

: "${MHC_b37:=6:${MHC_START_b37}-${MHC_END_b37}}"
: "${MHC_b38:=6:${MHC_START_b38}-${MHC_END_b38}}"


# 🚩 Logging and tabular-file helpers
phe_log() { echo "[$(date '+%F %T')] $*"; }

phe_zcat() {
  local file="$1"

  if declare -F gwas_clean_zcat >/dev/null 2>&1; then
    gwas_clean_zcat "$file"
    return
  fi

  case "$file" in
    *.zst)
      command -v zstdcat >/dev/null 2>&1 || {
        echo "ERROR: zstdcat is required to read: $file" >&2
        return 2
      }
      zstdcat -- "$file"
      ;;
    *.gz|*.bgz)
      if command -v pigz >/dev/null 2>&1; then
        pigz -dc -p "${GWAS_CLEAN_DECOMP_THREADS:-1}" "$file"
      else
        gzip -dc "$file"
      fi
      ;;
    *)
      cat "$file"
      ;;
  esac
}

# Return success only when FILE has a non-empty, current tabix index that
# tabix itself can read.  Both .tbi and .csi are accepted.
phe_tabix_index_valid() {
  local file="${1:?phe_tabix_index_valid requires FILE}" idx=""
  command -v tabix >/dev/null 2>&1 || return 1
  [[ -s "$file" ]] || return 1
  if [[ -s "$file.tbi" && ! "$file" -nt "$file.tbi" ]]; then
    idx="$file.tbi"
  elif [[ -s "$file.csi" && ! "$file" -nt "$file.csi" ]]; then
    idx="$file.csi"
  fi
  [[ -n "$idx" ]] || return 1
  tabix -l "$file" >/dev/null 2>&1
}

# Extract only GRCh sentinel positions from a coordinate-indexed standardized
# GWAS. Query both GRCh37 and GRCh38 coordinates; the caller then compares the
# returned positions with the gold table.  This turns a whole-file scan into at
# most 78 tiny tabix queries while preserving the legacy scan as a fallback.
phe_check_grch_tabix_rows() {
  local file="${1:?phe_check_grch_tabix_rows requires FILE}"
  local gold="${2:?phe_check_grch_tabix_rows requires GOLD}" out="${3:?phe_check_grch_tabix_rows requires OUT}"
  local work contigs regions header rc=0
  local -a query_regions=()

  phe_tabix_index_valid "$file" || return 1
  work=$(mktemp -d) || return 1
  contigs="$work/contigs"; regions="$work/regions"
  if ! tabix -l "$file" > "$contigs" 2>/dev/null || [[ ! -s "$contigs" ]]; then
    rm -rf -- "$work"
    return 1
  fi
  awk 'BEGIN{FS=OFS="\t"}
    function norm(x,y){y=toupper(x);sub(/^CHR/,"",y);if(y=="X")y="23";else if(y=="Y")y="24";else if(y=="M"||y=="MT")y="25";sub(/^0+/,"",y);return y}
    NR==FNR{actual[norm($1)]=$1;next}
    FNR>1{c=actual[norm($3)];if(c!=""){if($4~/^[0-9]+$/)print c ":" $4 "-" $4;if($5~/^[0-9]+$/)print c ":" $5 "-" $5}}
  ' "$contigs" "$gold" | sort -u > "$regions"
  mapfile -t query_regions < "$regions"
  if (( ${#query_regions[@]} == 0 )); then
    rm -rf -- "$work"
    return 1
  fi

  header=$( (set +o pipefail; phe_zcat "$file" | head -n 1) 2>/dev/null ) || header=""
  if [[ -z "$header" ]]; then
    rm -rf -- "$work"
    return 1
  fi
  if (
    set -o pipefail
    {
      printf '%s\n' "$header"
      tabix "$file" "${query_regions[@]}"
    } 2>/dev/null | awk 'BEGIN{FS=OFS="\t"}
    NR==1{
      for(i=1;i<=NF;i++){
        x=toupper($i);gsub(/^#/,"",x);gsub(/\r/,"",x)
        p=0;if(x=="SNP")p=1;else if(x=="RSID")p=2;else if(x=="ID")p=3;else if(x=="VARIANT_ID")p=4;else if(x=="VARIANT_IDS")p=5;if(p&&(!idp||p<idp)){id=i;idp=p}
        p=0;if(x=="CHR")p=1;else if(x=="CHROM")p=2;else if(x=="CHROMOSOME")p=3;else if(x=="CHR_NAME")p=4;if(p&&(!chrp||p<chrp)){chr=i;chrp=p}
        p=0;if(x=="POS")p=1;else if(x=="BP")p=2;else if(x=="BASE_PAIR")p=3;else if(x=="BASE_PAIR_LOCATION")p=4;else if(x=="POSITION")p=5;else if(x=="GENPOS")p=6;else if(x=="HM_POS")p=7;if(p&&(!posp||p<posp)){pos=i;posp=p}
      }
      if(!id||!chr||!pos)exit 2
      next
    }
    {v=$id;gsub(/\r/,"",v);print v,$chr,$pos}
    '
  ) > "$out"; then
    rc=0
  else
    rc=$?
  fi
  if (( rc == 0 )); then
    awk 'BEGIN{FS="\t"} NR==FNR{if(FNR>1)want[$1]=1;next} ($1 in want){found=1;exit} END{exit(found?0:1)}' "$gold" "$out" || rc=1
  fi
  rm -rf -- "$work"
  (( rc == 0 )) || : > "$out"
  return "$rc"
}

phe_header_names() {
  local file="$1" header_line header_fields matches match_count column value i pattern
  local -a names=(SNP CHR POS EA NEA EAF N BETA SE P LOG10P)
  local -a preferred_patterns=()
  local -a patterns=(
    '^snp$|^rsid$|^id$|^variant_id$|^variant_ids$'
    '^chr$|^chrom$|^chromosome$|^chr_name$'
    '^pos$|^bp$|^base_pair$|^base_pair_location$|^genpos$|^hm_pos$'
    '^a1$|^ea$|^eff.allele$|^effect_allele$|^allele1$'
    '^OMITTED$|^nea$|^non_effect_allele$|^other_allele$|^allele0$|^allele2$|^ref.allele$|^reference_allele$|^ref$'
    '^eaf$|^a1freq$|^a1_freq$|^effect_allele_frequency'
    '^n$|^obs_ct$|^Neff$'
    '^beta$|^effect_weight$'
    '^se$|^standard_error'
    '^p$|^pval$|^p_value$|^p_bolt_lmm$'
    '^log10p$|^neg.log.10.p.value$|^neg.log10.p.value$|^negative.log.10.p.value$|^negative.log10.p.value$|^minus.log10.p$|^mlog10p$'
  )

  header_line=$(set +o pipefail; phe_zcat "$file" | head -n 1 | tr '\t' ' ' | sed 's/\r$//')
  header_fields=$(tr ' ' '\n' <<< "$header_line")

  for i in "${!names[@]}"; do
    matches=""
    IFS='|' read -r -a preferred_patterns <<< "${patterns[$i]}"
    for pattern in "${preferred_patterns[@]}"; do
      matches=$(grep -Einw "$pattern" <<< "$header_fields" | cut -d: -f1 || true)
      match_count=$(grep -c -v '^$' <<< "$matches" || true)
      (( match_count == 0 )) || break
    done

    match_count=$(grep -c -v '^$' <<< "$matches" || true)
    if (( match_count > 0 )); then
      column=$(head -n 1 <<< "$matches")
      value=$(cut -d' ' -f "$column" <<< "$header_line")
    else
      column=""
      value=""
    fi
    printf -v "${names[$i]}" '%s' "$value"
    printf -v "${names[$i]}_col" '%s' "$column"
  done

  echo "dat $file, snp $SNP $SNP_col, chr $CHR $CHR_col, pos $POS $POS_col, ea $EA $EA_col, nea $NEA $NEA_col, eaf $EAF $EAF_col, n $N $N_col, beta $BETA $BETA_col, se $SE $SE_col, p $P $P_col, log10p $LOG10P $LOG10P_col"
}


# Calculate one PGS input against chromosome-split PLINK 2 pfiles.
#
# Required input headers (case-insensitive aliases are accepted):
#   CHR, SNP/ID, EA/A1/refA, BETA/bJ
# The function owns per-chromosome scoring and IID-checked merging only.  Batch
# discovery, scheduling, completion markers, and project metadata stay in the
# caller.  STATUS_FILE reports the values wrappers need after a successful run.
pgs_plink_calc() (
  set +e
  set -o pipefail
  local input="" pfile_dir="" output="" label="" work_dir="" status_file="" threads=1
  local tmp_dir="" normalized="" output_tmp="" status_tmp="" chr score_chr score_n pfile out_prefix norm_file rc=0
  local matched_status=scored intermediate_status=removed_after_success chr_csv=""
  local -a score_chrs=() norm_files=() work_prefixes=() pfile_modifier=()

  while (( $# )); do
    case "$1" in
      --input) [[ $# -ge 2 ]] || { echo "ERROR: pgs_plink_calc --input requires a value" >&2; return 2; }; input=$2; shift 2 ;;
      --pfile-dir) [[ $# -ge 2 ]] || { echo "ERROR: pgs_plink_calc --pfile-dir requires a value" >&2; return 2; }; pfile_dir=$2; shift 2 ;;
      --output) [[ $# -ge 2 ]] || { echo "ERROR: pgs_plink_calc --output requires a value" >&2; return 2; }; output=$2; shift 2 ;;
      --label) [[ $# -ge 2 ]] || { echo "ERROR: pgs_plink_calc --label requires a value" >&2; return 2; }; label=$2; shift 2 ;;
      --work-dir) [[ $# -ge 2 ]] || { echo "ERROR: pgs_plink_calc --work-dir requires a value" >&2; return 2; }; work_dir=$2; shift 2 ;;
      --status-file) [[ $# -ge 2 ]] || { echo "ERROR: pgs_plink_calc --status-file requires a value" >&2; return 2; }; status_file=$2; shift 2 ;;
      --threads) [[ $# -ge 2 ]] || { echo "ERROR: pgs_plink_calc --threads requires a value" >&2; return 2; }; threads=$2; shift 2 ;;
      -h|--help)
        cat <<'PGS_USAGE'
Usage: pgs_plink_calc --input FILE --pfile-dir DIR --output FILE --label LABEL
       [--work-dir DIR] [--status-file FILE] [--threads INT]
PGS_USAGE
        return 0
        ;;
      *) echo "ERROR: unknown pgs_plink_calc option: $1" >&2; return 2 ;;
    esac
  done

  [[ -n "$input" && -s "$input" ]] || { echo "ERROR: PGS input is missing/empty: ${input:-<none>}" >&2; return 2; }
  [[ -n "$pfile_dir" && -d "$pfile_dir" ]] || { echo "ERROR: PGS pfile directory is missing: ${pfile_dir:-<none>}" >&2; return 2; }
  [[ -n "$output" ]] || { echo "ERROR: pgs_plink_calc requires --output" >&2; return 2; }
  [[ -n "$label" && "$label" != */* ]] || { echo "ERROR: pgs_plink_calc requires a path-free --label" >&2; return 2; }
  [[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: PGS threads must be a positive integer: $threads" >&2; return 2; }
  command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 is required by pgs_plink_calc" >&2; return 127; }
  if ! { command -v awk >/dev/null 2>&1 && command -v gzip >/dev/null 2>&1 && command -v paste >/dev/null 2>&1; }; then
    echo "ERROR: pgs_plink_calc requires awk, gzip, and paste" >&2
    return 127
  fi

  [[ -n "$work_dir" ]] || work_dir=$(dirname -- "$output")
  mkdir -p -- "$work_dir" "$(dirname -- "$output")" || return 2
  if [[ -n "$status_file" ]]; then mkdir -p -- "$(dirname -- "$status_file")" || return 2; fi
  tmp_dir=$(mktemp -d "$work_dir/.pgs_plink_calc.${label}.XXXXXX") || return 2
  normalized="$tmp_dir/score.normalized.tsv"
  output_tmp="${output}.tmp.$$"
  [[ -z "$status_file" ]] || status_tmp="${status_file}.tmp.$$"
  cleanup_pgs_plink_calc(){
    [[ -z "$tmp_dir" ]] || rm -rf -- "$tmp_dir"
    [[ -z "$output_tmp" ]] || rm -f -- "$output_tmp"
    [[ -z "$status_tmp" ]] || rm -f -- "$status_tmp"
  }
  trap cleanup_pgs_plink_calc EXIT
  trap 'exit 129' HUP
  trap 'exit 130' INT
  trap 'exit 143' TERM

  phe_zcat "$input" | awk -v FS='[ \t]+' -v OFS='\t' '
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x==23)x="X";else if(x==24)x="Y";return x}
    NR==1{
      for(i=1;i<=NF;i++){
        h=toupper($i);sub(/^#/,"",h);gsub(/\r/,"",h)
        if(!chr&&(h=="CHR"||h=="CHROM"||h=="CHROMOSOME"))chr=i
        else if(!snp&&(h=="SNP"||h=="RSID"||h=="ID"||h=="VARIANT_ID"||h=="VARIANT_IDS"))snp=i
        else if(!ea&&(h=="EA"||h=="A1"||h=="EFFECT_ALLELE"||h=="REFA"))ea=i
        else if(!beta&&(h=="BETA"||h=="EFFECT_WEIGHT"||h=="BJ"))beta=i
      }
      if(!chr||!snp||!ea||!beta){
        print "ERROR: PGS input requires CHR, SNP/ID, EA/A1/refA, and BETA/bJ headers" > "/dev/stderr"
        bad=1;exit 2
      }
      print "CHR","SNP","EA","BETA";next
    }
    NF{
      c=normchr($chr);id=$snp;a=toupper($ea);w=$beta;gsub(/\r/,"",w)
      if(c!~/^([1-9]|1[0-9]|2[0-2]|X|Y)$/||id==""||id=="."||toupper(id)=="NA"||a!~/^[ACGT]+$/||!isnum(w)){
        print "ERROR: invalid PGS row in " FILENAME " at line " NR ": " $0 > "/dev/stderr";bad=1;next
      }
      key=c SUBSEP id
      if(seen[key]++){print "ERROR: duplicate PGS SNP within chr" c ": " id > "/dev/stderr";bad=1;next}
      print c,id,a,w;n++
    }
    END{if(n==0){print "ERROR: no valid PGS rows in " FILENAME > "/dev/stderr";bad=1}exit bad?2:0}
  ' > "$normalized"
  rc=$?
  (( rc == 0 )) || return "$rc"

  mapfile -t score_chrs < <(awk -F '\t' 'NR>1{print $1}' "$normalized" | sort -uV)
  (( ${#score_chrs[@]} > 0 )) || { echo "ERROR: no PGS chromosomes found in $input" >&2; return 2; }
  chr_csv=$(IFS=,; echo "${score_chrs[*]}")

  for chr in "${score_chrs[@]}"; do
    score_chr="$work_dir/$label.chr$chr.score.tsv"
    awk -v FS='\t' -v OFS='\t' -v chr="$chr" 'NR==1{print "SNP","EA","BETA";next}$1==chr{print $2,$3,$4}' "$normalized" > "$score_chr"
    score_n=$(($(wc -l < "$score_chr")-1))
    (( score_n > 0 )) || continue

    pfile="${pfile_dir%/}/chr$chr"
    [[ -s "${pfile}.pgen" && -s "${pfile}.psam" ]] || { echo "ERROR: incomplete PGS pfile for chr$chr: $pfile" >&2; return 2; }
    [[ -s "${pfile}.pvar" || -s "${pfile}.pvar.zst" ]] || { echo "ERROR: missing PGS pvar for chr$chr: $pfile" >&2; return 2; }
    pfile_modifier=()
    [[ ! -s "${pfile}.pvar.zst" ]] || pfile_modifier=(vzs)

    out_prefix="$work_dir/$label.chr$chr"
    work_prefixes+=("$out_prefix")
    rm -f -- "${out_prefix}.sscore" "${out_prefix}.sscore.vars" "${out_prefix}.log"
    plink2 --threads "$threads" --pfile "$pfile" "${pfile_modifier[@]}" --chr "$chr" \
      --score "$score_chr" 1 2 3 header-read no-mean-imputation ignore-dup-ids \
      list-variants cols=+scoresums --out "$out_prefix"
    rc=$?
    if (( rc != 0 )); then
      if [[ -s "${out_prefix}.log" ]] && grep -Eqi 'no ([[:alnum:]_-]+ )?variants|0 variants|no variants (remaining|loaded|processed)' "${out_prefix}.log"; then
        phe_log "No chr$chr variants were scored; retained PLINK2 log: ${out_prefix}.log" >&2
        continue
      fi
      echo "ERROR: plink2 PGS failed for $label chr$chr (exit $rc); see ${out_prefix}.log" >&2
      return "$rc"
    fi
    if [[ ! -s "${out_prefix}.sscore" ]]; then
      phe_log "No chr$chr score output; retained PLINK2 log: ${out_prefix}.log" >&2
      continue
    fi

    norm_file="$tmp_dir/$label.chr$chr.norm.tsv"
    awk -v FS='[ \t]+' -v OFS='\t' '
      NR==1{
        for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);if(h=="IID")iid=i;else if(h=="ALLELE_CT")ct=i;else if(h=="SCORE1_SUM")sum=i;else if(!sum&&h~/_SUM$/)sum=i}
        if(!iid||!ct||!sum){print "ERROR: .sscore lacks IID, ALLELE_CT, or *_SUM in " FILENAME > "/dev/stderr";exit 2}
        print "#IID","ALLELE_CT","SCORE_SUM";next
      }
      {if($iid==""||$sum=="NA"||$sum=="."){print "ERROR: invalid .sscore row in " FILENAME " at line " NR > "/dev/stderr";exit 2}print $iid,$ct,$sum}
    ' "${out_prefix}.sscore" > "$norm_file"
    rc=$?
    (( rc == 0 )) || return "$rc"
    norm_files+=("$norm_file")
  done

  if (( ${#norm_files[@]} == 0 )); then
    printf '#IID\t%s.ALLELE_CT\t%s.SCORE_SUM\n' "$label" "$label" > "$tmp_dir/output.txt"
    matched_status=none
    intermediate_status=plink_logs_retained_no_matched_variants
  else
    paste "${norm_files[@]}" | awk -v FS='\t' -v OFS='\t' -v n="${#norm_files[@]}" -v label="$label" '
      NR==1{print "#IID",label ".ALLELE_CT",label ".SCORE_SUM";next}
      {
        iid=$1;allele_ct=0;score_sum=0
        for(i=0;i<n;i++){
          base=i*3
          if($(base+1)!=iid){print "ERROR: IID order mismatch while merging chromosome PGS at row " NR > "/dev/stderr";exit 2}
          allele_ct+=$(base+2);score_sum+=$(base+3)
        }
        print iid,allele_ct,sprintf("%.15g",score_sum)
      }
    ' > "$tmp_dir/output.txt"
    rc=$?
    (( rc == 0 )) || return "$rc"
  fi

  if [[ "$output" == *.gz || "$output" == *.bgz ]]; then
    gzip -c "$tmp_dir/output.txt" > "$output_tmp"
    rc=$?
    (( rc == 0 )) || return "$rc"
    gzip -t "$output_tmp" || return $?
  else
    cp -f -- "$tmp_dir/output.txt" "$output_tmp" || return $?
  fi
  mv -f -- "$output_tmp" "$output" || return $?
  output_tmp=""

  for out_prefix in "${work_prefixes[@]}"; do
    if [[ "$matched_status" == none ]]; then
      find "$work_dir" -maxdepth 1 -type f -name "$(basename -- "$out_prefix").*" ! -name '*.log' -delete
    else
      rm -f -- "${out_prefix}".*
    fi
  done
  rm -f -- "$work_dir/$label.chr"*.score.tsv

  if [[ -n "$status_file" ]]; then
    {
      printf 'key\tvalue\n'
      printf 'matched_variants\t%s\n' "$matched_status"
      printf 'chromosomes\t%s\n' "$chr_csv"
      printf 'chromosome_intermediates\t%s\n' "$intermediate_status"
      printf 'output\t%s\n' "$output"
    } > "$status_tmp"
    mv -f -- "$status_tmp" "$status_file" || return $?
    status_tmp=""
  fi

  trap - EXIT HUP INT TERM
  rm -rf -- "$tmp_dir"
  tmp_dir=""
  return 0
)


# 🚩 Genome-build detection
# Detect GRCh37/38 from sentinel variants; export CHECK_GRCH_RESULT.
check_GRCH() {
  local input="${1:?check_GRCH requires INPUT}" expected="${2:-}" gold="${CHECK_GRCH_SNP_LIST:-/mnt/d/data/ukb/phe/common/snp.lst}"
  local src="$input" kind=text tmp n m37 m38 best mismatch build expected_n chr candidate pv
  local -a ref_bims=() ref_pvars=()
  [[ -s "$gold" ]] || { echo "ERROR: check_GRCH gold standard is missing/empty: $gold" >&2; return 2; }

  if [[ -d "$src" ]]; then
    kind=refdir
  elif [[ ! -s "$src" ]]; then
    for src in "$input.pvar" "$input.pvar.gz" "$input.pvar.bgz" "$input.pvar.zst" "$input.bim"; do [[ -s "$src" ]] && break; done
  fi
  [[ -d "$src" || -s "$src" ]] || { echo "ERROR: check_GRCH input is missing/empty: $input" >&2; return 2; }
  if [[ "$kind" != refdir ]]; then
    case "$src" in *.bim) kind=bim;; *.pvar|*.pvar.gz|*.pvar.bgz|*.pvar.zst) kind=pvar;; *.vcf|*.vcf.gz|*.vcf.bgz) kind=vcf;; esac
  fi
  tmp=$(mktemp)

  case "$kind" in
    refdir)
      # One sentinel-rich chromosome is sufficient for split references.
      while read -r _ chr; do
        [[ -n "$chr" ]] || continue
        candidate="$src/chr${chr}.bim"
        if [[ -s "$candidate" ]]; then
          ref_bims+=("$candidate")
          break
        fi
        for candidate in "$src/chr${chr}.pvar" "$src/chr${chr}.pvar.gz" "$src/chr${chr}.pvar.bgz" "$src/chr${chr}.pvar.zst"; do
          if [[ -s "$candidate" ]]; then ref_pvars+=("$candidate"); break; fi
        done
        (( ${#ref_pvars[@]} == 0 )) || break
      done < <(awk -F '\t' 'NR>1&&$3!=""{n[$3]++}END{for(c in n)print n[c],c}' "$gold" | sort -k1,1nr -k2,2V)
      (( ${#ref_bims[@]} + ${#ref_pvars[@]} > 0 )) || {
        echo "ERROR: check_GRCH found no chromosome-split .bim/.pvar files under: $src" >&2
        rm -f "$tmp"
        return 2
      }
      if (( ${#ref_bims[@]} > 0 )); then
        awk 'BEGIN{FS=OFS="\t"}
          NR==FNR{if(FNR>1){want[$1]=1;p37[$3 SUBSEP $4]=$1;p38[$3 SUBSEP $5]=$1};next}
          NF>=4{id=$2;chr=$1;pos=$4;key=chr SUBSEP pos
            if(id in want)print id,chr,pos;else if(key in p37)print p37[key],chr,pos;else if(key in p38)print p38[key],chr,pos}' \
          "$gold" "${ref_bims[@]}" >> "$tmp"
      fi
      for pv in "${ref_pvars[@]}"; do
        phe_zcat "$pv" | awk -v gold="$gold" 'BEGIN{FS=OFS="\t";while((getline line<gold)>0){n=split(line,a,"\t");if(n&&a[1]!="SNP"){want[a[1]]=1;p37[a[3] SUBSEP a[4]]=a[1];p38[a[3] SUBSEP a[5]]=a[1]}}close(gold)}
          $1=="#CHROM"{for(i=1;i<=NF;i++){x=$i;sub(/^#/,"",x);c[x]=i};next} /^##/{next}
          ("ID" in c)&&("POS" in c){id=$(c["ID"]);chr=$(c["CHROM"]);pos=$(c["POS"]);key=chr SUBSEP pos
            if(id in want)print id,chr,pos;else if(key in p37)print p37[key],chr,pos;else if(key in p38)print p38[key],chr,pos}' >> "$tmp"
      done
      ;;
    bim)
      awk 'BEGIN{OFS="\t"} NF>=4{print $2,$1,$4}' "$src" > "$tmp" ;;
    pvar)
      phe_zcat "$src" | awk 'BEGIN{FS=OFS="\t"} $1=="#CHROM"{for(i=1;i<=NF;i++){x=$i;sub(/^#/,"",x);c[x]=i};next} /^##/{next} ("ID" in c)&&("POS" in c){print $(c["ID"]),$(c["CHROM"]),$(c["POS"])}' > "$tmp" ;;
    vcf)
      phe_zcat "$src" | awk 'BEGIN{FS=OFS="\t"} !/^#/{print $3,$1,$2}' > "$tmp" ;;
    text)
      if phe_check_grch_tabix_rows "$src" "$gold" "$tmp"; then
        phe_log "check_GRCH: used tabix sentinel lookup for $src"
      else
        phe_zcat "$src" | awk 'BEGIN{FS=OFS="\t"}
          NR==1{
            for(i=1;i<=NF;i++){
              x=toupper($i);gsub(/^#/,"",x);gsub(/\r/,"",x)
              p=0;if(x=="SNP")p=1;else if(x=="RSID")p=2;else if(x=="ID")p=3;else if(x=="VARIANT_ID")p=4;else if(x=="VARIANT_IDS")p=5;if(p&&(!idp||p<idp)){id=i;idp=p}
              p=0;if(x=="CHR")p=1;else if(x=="CHROM")p=2;else if(x=="CHROMOSOME")p=3;else if(x=="CHR_NAME")p=4;if(p&&(!chrp||p<chrp)){chr=i;chrp=p}
              p=0;if(x=="POS")p=1;else if(x=="BP")p=2;else if(x=="BASE_PAIR")p=3;else if(x=="BASE_PAIR_LOCATION")p=4;else if(x=="POSITION")p=5;else if(x=="GENPOS")p=6;else if(x=="HM_POS")p=7;if(p&&(!posp||p<posp)){pos=i;posp=p}
            }
            next
          }
          id&&chr&&pos{v=$id;gsub(/\r/,"",v);print v,$chr,$pos}' > "$tmp"
      fi ;;
  esac

  read -r n m37 m38 < <(awk 'BEGIN{FS="\t"}
    NR==FNR{if(FNR>1){gchr[$1]=$3;g37[$1]=$4;g38[$1]=$5}next}
    {id=$1;chr=$2;pos=$3;sub(/^chr/,"",chr);if(!(id in gchr)||seen[id]++)next;if(chr!=gchr[id])next;n++;if(pos==g37[id])m37++;if(pos==g38[id])m38++}
    END{print n+0,m37+0,m38+0}' "$gold" "$tmp")
  (( n > 0 )) || { echo "ERROR: check_GRCH found 0 of the 39 gold-standard SNPs in $src; GRCh cannot be verified." >&2; rm -f "$tmp"; return 3; }
  if (( m37 > m38 )); then build=37; best=$m37; elif (( m38 > m37 )); then build=38; best=$m38; else
    echo "ERROR: $n SNPs exist in snp.lst, GRCh37 matches=$m37 and GRCh38 matches=$m38; build is ambiguous, quit." >&2; rm -f "$tmp"; return 4
  fi
  mismatch=$((n-best))
  if (( mismatch > 1 )); then
    echo "ERROR: $n SNPs exist in snp.lst, $best match GRCH.$build, $mismatch mismatch, quit." >&2; rm -f "$tmp"; return 5
  elif (( mismatch == 1 )); then
    phe_log "$n SNPs exist in snp.lst, $best match GRCH.$build, 1 mismatch is allowed, so GRCH is set to $build"
  else
    phe_log "$n SNPs exist in snp.lst, all match GRCH.$build, so GRCH is set to $build"
  fi
  expected_n=$(printf '%s' "$expected" | tr '[:upper:]' '[:lower:]')
  case "$expected_n" in 37|b37|grch37|hg19) expected_n=37;; 38|b38|grch38|hg38) expected_n=38;; "") expected_n="";; *) echo "ERROR: invalid expected GRCh: $expected" >&2; rm -f "$tmp"; return 2;; esac
  [[ -z "$expected_n" || "$expected_n" == "$build" ]] || { echo "ERROR: configured GRCH=$expected_n but check_GRCH inferred GRCH=$build from $src, quit." >&2; rm -f "$tmp"; return 6; }
  CHECK_GRCH_RESULT="$build"; export CHECK_GRCH_RESULT
  rm -f "$tmp"
  return 0
}


# 🚩 Background-process control
phe_background_job_pattern() {
  local x b
  for x in "$@"; do
    b="$(basename -- "$x")"
    case "$b" in
      *.sh|*.py) printf '%s\n' "$b"; return 0 ;;
    esac
  done
  b="$(basename -- "${1:-background}")"
  printf '%s\n' "$b"
}

phe_running_tasks() {
  local pattern="$1" self="${2:-$$}" bash_self="${3:-${BASHPID:-$$}}" parent="${4:-${PPID:-}}"
  ps -eo pid=,ppid=,args= 2>/dev/null | PHE_TASK_PATTERN="$pattern" awk -v self="$self" -v bash_self="$bash_self" -v parent="$parent" '
    BEGIN { pattern = ENVIRON["PHE_TASK_PATTERN"] }
    {
      pid = $1
      ppid = $2
      $1 = ""
      $2 = ""
      sub(/^[[:space:]]+/, "")
      cmd = $0
    }
    pid == self { next }
    pid == bash_self { next }
    pid == parent { next }
    ppid == self { next }
    ppid == bash_self { next }
    index(cmd, pattern) > 0 { print pid "\t" cmd }
  '
}

phe_check_existing_tasks() {
  local pattern="$1" old pids self bash_self parent
  self="$$"
  bash_self="${BASHPID:-$$}"
  parent="${PPID:-}"
  old="$(phe_running_tasks "$pattern" "$self" "$bash_self" "$parent" || true)"
  [[ -z "$old" ]] && return 0
  pids="$(printf '%s\n' "$old" | awk '{print $1}' | paste -sd ' ' -)"
  {
    echo "ERROR: existing $pattern task(s) are still running:"
    printf '%s\n' "$old" | sed 's/^/  /'
    echo
    echo "Stop them first, then rerun this command:"
    echo "  kill $pids"
    echo
    echo "If they do not exit after a few seconds, use:"
    echo "  kill -9 $pids"
  } >&2
  return 1
}

phe_run_background() {
  local title="background job" job_pattern
  if [[ "${1-}" == "--title" ]]; then
    title="$2"
    shift 2
  fi
  local log_file="$1"
  shift || {
    echo "ERROR: phe_run_background requires: LOG_FILE COMMAND [ARGS...]" >&2
    return 2
  }
  [[ $# -gt 0 ]] || {
    echo "ERROR: phe_run_background missing command" >&2
    return 2
  }
  mkdir -p -- "$(dirname -- "$log_file")"
  echo "Starting $title. Log: $log_file"
  nohup "$@" > "$log_file" 2>&1 < /dev/null &
  echo "PID: $!"
  echo "Progress: tail -f $log_file"
  job_pattern="$(phe_background_job_pattern "$@")"
  echo "Job: pgrep -af '$job_pattern'"
}


# 🚩 Variant and haplotype matching
# Compare query/reference haplotype strings one character at a time.
# Input: query_id, reference_id, query_haplotype, reference_haplotype (tab-separated).
# Any symbol outside A/C/G/T is ambiguous/missing and excluded from the denominator.
match_HAP() {
  local has_header=0 input=-
  while (( $# )); do
    case "$1" in
      --header) has_header=1; shift ;;
      --no-header) has_header=0; shift ;;
      -|--stdin) input=-; shift ;;
      --) shift; [[ $# -le 1 ]] || { echo "ERROR: match_HAP accepts one input file" >&2; return 2; }; input=${1:--}; shift || true ;;
      -*) echo "ERROR: unknown match_HAP option: $1" >&2; return 2 ;;
      *) [[ $input == - ]] || { echo "ERROR: match_HAP accepts one input file" >&2; return 2; }; input=$1; shift ;;
    esac
  done
  [[ $input == - || -s $input ]] || { echo "ERROR: haplotype input is missing or empty: $input" >&2; return 2; }
  awk -v header="$has_header" 'BEGIN{FS=OFS="\t"; print "query_id","reference_id","query_haplotype","reference_haplotype","n_compared","n_match","n_mismatch","prop_match"}
    header && NR==1 {next}
    NF<4 {print "ERROR: match_HAP expects four tab-separated columns at line " NR > "/dev/stderr"; bad=1; next}
    {q=toupper($3);r=toupper($4);if(length(q)!=length(r)){print "ERROR: unequal haplotype lengths at line " NR > "/dev/stderr";bad=1;next};n=0;m=0
     for(i=1;i<=length(q);i++){a=substr(q,i,1);b=substr(r,i,1);if(a!~/^[ACGT]$/||b!~/^[ACGT]$/)continue;n++;if(a==b)m++}
     prop=n?sprintf("%.12g",m/n):"NA";print $1,$2,$3,$4,n,m,n-m,prop}
    END{if(bad)exit 2}' "$input"
}

# Disk-backed SNP matching.  QUERY rows and lookup keys are normalized to
# temporary files, externally sorted with a bounded buffer, and merge-joined to
# a streaming normalization of REFERENCE.  No AWK process retains QUERY-sized
# arrays.
match_SNP() (
  local query="" ref="" output="-" audit="" unmatched="" position_only=0 id_any_position=0
  local memory_mb="${MATCH_SNP_MEMORY_MB:-4096}" sort_memory_mb="${MATCH_SNP_SORT_MEMORY_MB:-512}"
  local sort_threads="${MATCH_SNP_SORT_THREADS:-1}" tmp_root="${MATCH_SNP_TMP_DIR:-/mnt/d/tmp}"
  local warn_query_rows="${MATCH_SNP_WARN_QUERY_ROWS:-500000}" max_query_rows="${MATCH_SNP_MAX_QUERY_ROWS:-5000000}"
  local allow_large_query=0 ref_file="" ref_kind="" tmp_dir="" rc=0 ext
  local qn rn snp_col chr_col pos_col ea_col nea_col ncols join_fields i
  local need_IDF=0 need_IDA=0 need_ID1=0 need_ID1A=0 need_ID0=0 need_ID0A=0
  local need_A=0 need_S=0 need_O=0 need_T=0 need_P=0

  while (( $# )); do
    case "$1" in
      -r|--reference) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; ref=$2; shift 2 ;;
      -o|--output) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; output=$2; shift 2 ;;
      --audit) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; audit=$2; shift 2 ;;
      --unmatched) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; unmatched=$2; shift 2 ;;
      --position-only) position_only=1; shift ;;
      --id-any-position) id_any_position=1; shift ;;
      --memory-mb) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires MiB" >&2; return 2; }; memory_mb=$2; shift 2 ;;
      --sort-memory-mb) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires MiB" >&2; return 2; }; sort_memory_mb=$2; shift 2 ;;
      --sort-threads) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an integer" >&2; return 2; }; sort_threads=$2; shift 2 ;;
      --tmp-dir) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a directory" >&2; return 2; }; tmp_root=$2; shift 2 ;;
      --warn-query-rows) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an integer" >&2; return 2; }; warn_query_rows=$2; shift 2 ;;
      --max-query-rows) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an integer" >&2; return 2; }; max_query_rows=$2; shift 2 ;;
      --allow-large-query) allow_large_query=1; shift ;;
      --) shift; [[ $# -eq 1 ]] || { echo "ERROR: match_SNP requires one QUERY file" >&2; return 2; }; query=$1; shift ;;
      -*) echo "ERROR: unknown match_SNP option: $1" >&2; return 2 ;;
      *) [[ -z "$query" ]] || { echo "ERROR: match_SNP accepts one QUERY file" >&2; return 2; }; query=$1; shift ;;
    esac
  done

  [[ -n "$query" && -s "$query" ]] || { echo "ERROR: match_SNP QUERY is missing/empty: ${query:-<none>}" >&2; return 2; }
  [[ -n "$ref" ]] || { echo "ERROR: match_SNP requires --reference REF" >&2; return 2; }
  for i in memory_mb sort_memory_mb sort_threads warn_query_rows max_query_rows; do
    [[ "${!i}" =~ ^[0-9]+$ ]] || { echo "ERROR: --${i//_/-} must be a non-negative integer: ${!i}" >&2; return 2; }
  done
  (( sort_memory_mb > 0 && sort_threads > 0 )) || { echo "ERROR: sort memory and threads must be positive" >&2; return 2; }
  command -v sort >/dev/null 2>&1 && command -v join >/dev/null 2>&1 && command -v awk >/dev/null 2>&1 || {
    echo "ERROR: match_SNP requires GNU sort, join, and awk" >&2; return 2;
  }

  case "$ref" in
    *.pvar|*.pvar.gz|*.pvar.bgz|*.pvar.zst) ref_file=$ref; ref_kind=pvar ;;
    *.bim|*.bim.gz|*.bim.bgz) ref_file=$ref; ref_kind=bim ;;
    *)
      for ext in .pvar .pvar.gz .pvar.bgz .pvar.zst .bim .bim.gz .bim.bgz; do
        if [[ -s "${ref}${ext}" ]]; then ref_file="${ref}${ext}"; [[ "$ext" == .pvar* ]] && ref_kind=pvar || ref_kind=bim; break; fi
      done
      ;;
  esac
  [[ -n "$ref_file" && -s "$ref_file" ]] || { echo "ERROR: no .pvar/.bim reference found from: $ref" >&2; return 2; }
  mkdir -p -- "$tmp_root" || return 2
  [[ -d "$tmp_root" && -w "$tmp_root" ]] || { echo "ERROR: temporary directory is not writable: $tmp_root" >&2; return 2; }
  tmp_dir=$(mktemp -d "${tmp_root%/}/match_SNP.XXXXXX") || return 2
  trap 'rm -rf -- "$tmp_dir"' EXIT HUP INT TERM

  # A function subshell prevents this limit from leaking into the caller's
  # terminal.  ulimit is per process; systemd-run MemoryMax is stronger when a
  # strict aggregate cap for the whole pipeline is required.
  if (( memory_mb > 0 )); then
    ulimit -Sv "$((memory_mb * 1024))" || { echo "ERROR: cannot set match_SNP memory cap to ${memory_mb} MiB" >&2; return 2; }
  fi

  : > "$tmp_dir/query.rows.tsv"; : > "$tmp_dir/query.keys.tsv"; : > "$tmp_dir/classes"
  if (set -o pipefail; phe_zcat "$query" | awk -v FS='\t' -v OFS='\t' \
      -v rows="$tmp_dir/query.rows.tsv" -v keys="$tmp_dir/query.keys.tsv" \
      -v header_file="$tmp_dir/header.tsv" -v cols="$tmp_dir/cols.tsv" \
      -v count_file="$tmp_dir/query.count" -v classes="$tmp_dir/classes" \
      -v pos_only="$position_only" -v id_any="$id_any_position" \
      -v max_rows="$max_query_rows" -v allow_large="$allow_large_query" '
    function hname(x){x=toupper(x);sub(/^#/,"",x);gsub(/\r/,"",x);return x}
    function normchr(x){x=toupper(x);sub(/^CHR/,"",x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="M"||x=="MT")x=25;sub(/^0+/,"",x);return x==""?"0":x}
    function valid(x){x=toupper(x);return x!=""&&x!="."&&x!="NA"&&x~/^[ACGT]+$/}
    function rc(x, i,c,y){x=toupper(x);y="";for(i=length(x);i>0;i--){c=substr(x,i,1);if(c=="A")c="T";else if(c=="T")c="A";else if(c=="C")c="G";else if(c=="G")c="C";else return "";y=y c}return y}
    function pair(a,b){a=toupper(a);b=toupper(b);return a<=b?a"~"b:b"~"a}
    function cleanid(x){x=toupper(x);gsub(/\r/,"",x);gsub(/%/,"%25",x);gsub(/[|]/,"%7C",x);return x}
    function emit(k,p,t,c){print k,rowkey,p,t >> keys;if(!(c in used)){print c >> classes;used[c]=1}}
    NR==1{
      print $0 > header_file;ncol=NF
      for(j=1;j<=NF;j++){
        x=hname($j)
        if(!snp&&(x=="SNP"||x=="RSID"||x=="ID"||x=="VARIANT_ID"))snp=j
        else if(!chr&&(x=="CHR"||x=="CHROM"||x=="CHROMOSOME"))chr=j
        else if(!pos&&(x=="POS"||x=="BP"||x=="POSITION"||x=="BASE_PAIR_LOCATION"))pos=j
        else if(!ea&&(x=="EA"||x=="A1"||x=="EFFECT_ALLELE"||x=="REFA"))ea=j
        else if(!nea&&(x=="NEA"||x=="A2"||x=="OTHER_ALLELE"||x=="NON_EFFECT_ALLELE"))nea=j
      }
      if(!snp||!chr||!pos){print "ERROR: match_SNP QUERY requires SNP, CHR, and POS columns" > "/dev/stderr";exit 2}
      print ncol,snp,chr,pos,ea+0,nea+0 > cols;next
    }
    {
      n++;if(max_rows>0&&!allow_large&&n>max_rows){print "ERROR: QUERY exceeds --max-query-rows " max_rows "; verify QUERY/REFERENCE or use --allow-large-query" > "/dev/stderr";exit 3}
      rowkey=sprintf("%012d",n);print rowkey,$0 >> rows
      id=cleanid($snp);ch=normchr($chr);bp=$pos;a=ea?toupper($ea):"";b=nea?toupper($nea):""
      full=valid(a)&&valid(b);one=!full&&(valid(a)||valid(b));o=valid(a)?a:b
      goodid=(id!=""&&id!="."&&id!="NA");goodpos=(ch!="0"&&bp~/^[0-9]+$/)
      if(goodid){
        if(full){
          emit("IDF|"id"|"ch"|"bp"|"pair(a,b),1,"ID_ALLELES","IDF")
          emit("IDF|"id"|"ch"|"bp"|"pair(rc(a),rc(b)),1,"ID_STRAND","IDF")
          if(id_any){emit("IDA|"id"|"pair(a,b),1,"ID_ANY_POSITION_ALLELES","IDA");emit("IDA|"id"|"pair(rc(a),rc(b)),1,"ID_ANY_POSITION_STRAND","IDA")}
        }else if(one){
          emit("ID1|"id"|"ch"|"bp"|"o,2,"ID_ALLELE1","ID1");emit("ID1|"id"|"ch"|"bp"|"rc(o),2,"ID_ALLELE1_STRAND","ID1")
          if(id_any){emit("ID1A|"id"|"o,2,"ID_ANY_POSITION_ALLELE1","ID1A");emit("ID1A|"id"|"rc(o),2,"ID_ANY_POSITION_ALLELE1_STRAND","ID1A")}
        }else{
          emit("ID0|"id"|"ch"|"bp,2,"ID_POSITION","ID0");if(id_any)emit("ID0A|"id,2,"ID_ANY_POSITION","ID0A")
        }
      }
      if(goodpos){
        if(full){emit("A|"ch"|"bp"|"pair(a,b),3,"ALLELES","A");emit("S|"ch"|"bp"|"pair(rc(a),rc(b)),4,"STRAND","S")}
        else if(one){emit("O|"ch"|"bp"|"o,5,"ALLELE1","O");emit("T|"ch"|"bp"|"rc(o),6,"ALLELE1_STRAND","T")}
        if(pos_only)emit("P|"ch"|"bp,7,"POSITION_ONLY","P")
      }
    }
    END{print n+0 > count_file}
  '); then :; else rc=$?; return "$rc"; fi

  read -r ncols snp_col chr_col pos_col ea_col nea_col < "$tmp_dir/cols.tsv"
  qn=$(<"$tmp_dir/query.count")
  if (( warn_query_rows > 0 && qn >= warn_query_rows )); then
    echo "WARNING: match_SNP QUERY has $qn rows. Confirm QUERY and REFERENCE were not reversed." >&2
  fi
  while IFS= read -r i; do [[ -n "$i" ]] && printf -v "need_$i" 1; done < "$tmp_dir/classes"

  : > "$tmp_dir/reference.keys.tsv"
  if (set -o pipefail; phe_zcat "$ref_file" | awk -v FS='\t' -v OFS='\t' -v kind="$ref_kind" \
      -v out="$tmp_dir/reference.keys.tsv" -v count_file="$tmp_dir/reference.count" \
      -v IDF="$need_IDF" -v IDA="$need_IDA" -v ID1="$need_ID1" -v ID1A="$need_ID1A" -v ID0="$need_ID0" -v ID0A="$need_ID0A" \
      -v A="$need_A" -v S="$need_S" -v O="$need_O" -v T="$need_T" -v P="$need_P" '
    function hname(x){x=toupper(x);sub(/^#/,"",x);gsub(/\r/,"",x);return x}
    function normchr(x){x=toupper(x);sub(/^CHR/,"",x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="M"||x=="MT")x=25;sub(/^0+/,"",x);return x==""?"0":x}
    function valid(x){x=toupper(x);return x!=""&&x!="."&&x!="NA"&&x~/^[ACGT]+$/}
    function pair(a,b){a=toupper(a);b=toupper(b);return a<=b?a"~"b:b"~"a}
    function cleanid(x){x=toupper(x);gsub(/\r/,"",x);gsub(/%/,"%25",x);gsub(/[|]/,"%7C",x);return x}
    function emit(k){print k,rid,ch,bp,r1,r2,rid"|"ch"|"bp"|"pair(r1,r2) >> out}
    function add(id,c,p,a,b, pp){
      rid=id;ch=normchr(c);bp=p;r1=toupper(a);r2=toupper(b);iid=cleanid(rid)
      if(iid==""||iid=="."||bp!~/^[0-9]+$/||!valid(r1)||!valid(r2))return
      pp=pair(r1,r2);n++
      if(IDF)emit("IDF|"iid"|"ch"|"bp"|"pp);if(IDA)emit("IDA|"iid"|"pp)
      if(ID1){emit("ID1|"iid"|"ch"|"bp"|"r1);emit("ID1|"iid"|"ch"|"bp"|"r2)}
      if(ID1A){emit("ID1A|"iid"|"r1);emit("ID1A|"iid"|"r2)}
      if(ID0)emit("ID0|"iid"|"ch"|"bp);if(ID0A)emit("ID0A|"iid)
      if(A)emit("A|"ch"|"bp"|"pp);if(S)emit("S|"ch"|"bp"|"pp)
      if(O){emit("O|"ch"|"bp"|"r1);emit("O|"ch"|"bp"|"r2)}
      if(T){emit("T|"ch"|"bp"|"r1);emit("T|"ch"|"bp"|"r2)}
      if(P)emit("P|"ch"|"bp)
    }
    kind=="pvar"{
      if(/^##/)next;nr=split($0,z,/[ \t]+/)
      if(/^#CHROM/){for(j=1;j<=nr;j++){x=hname(z[j]);if(x=="CHROM")cc=j;else if(x=="POS")pc=j;else if(x=="ID")ic=j;else if(x=="REF")ac=j;else if(x=="ALT")bc=j}seen_header=1;next}
      if(!seen_header||!cc||!pc||!ic||!ac||!bc){bad=1;next}
      na=split(z[bc],alts,",");for(k=1;k<=na;k++)add(z[ic],z[cc],z[pc],z[ac],alts[k]);next
    }
    kind=="bim"{nr=split($0,z,/[ \t]+/);if(nr>=6)add(z[2],z[1],z[4],z[5],z[6]);next}
    END{print n+0 > count_file;if(kind=="pvar"&&!seen_header){print "ERROR: invalid .pvar header" > "/dev/stderr";exit 2}else if(bad)exit 2}
  '); then :; else rc=$?; return "$rc"; fi
  rn=$(<"$tmp_dir/reference.count")
  if (( rn > 0 && qn > rn )); then
    echo "WARNING: QUERY rows ($qn) exceed REFERENCE variants ($rn); QUERY/REFERENCE may be reversed." >&2
  fi

  LC_ALL=C sort -t $'\t' -k1,1 -k2,2 -k3,3n -k4,4 --parallel="$sort_threads" -S "${sort_memory_mb}M" -T "$tmp_dir" \
    "$tmp_dir/query.keys.tsv" > "$tmp_dir/query.keys.sorted.tsv" || return $?
  LC_ALL=C sort -t $'\t' -k1,1 -k7,7 -u --parallel="$sort_threads" -S "${sort_memory_mb}M" -T "$tmp_dir" \
    "$tmp_dir/reference.keys.tsv" > "$tmp_dir/reference.keys.sorted.tsv" || return $?
  LC_ALL=C join -t $'\t' -1 1 -2 1 -o 1.2,1.3,1.4,2.2,2.3,2.4,2.5,2.6,2.7 \
    "$tmp_dir/query.keys.sorted.tsv" "$tmp_dir/reference.keys.sorted.tsv" > "$tmp_dir/candidates.tsv" || return $?
  LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k9,9 -k3,3 --parallel="$sort_threads" -S "${sort_memory_mb}M" -T "$tmp_dir" \
    "$tmp_dir/candidates.tsv" > "$tmp_dir/candidates.sorted.tsv" || return $?
  awk -v FS='\t' -v OFS='\t' '
    function flush(){if(cur=="")return;if(nsig==1)print cur,"MATCHED",rid,type;else print cur,"AMBIGUOUS",".","."}
    $1!=cur{flush();cur=$1;minp=$2;lastsig="";nsig=0;rid=".";type="."}
    $2!=minp{next}
    $9!=lastsig{nsig++;lastsig=$9;if(nsig==1){rid=$4;type=$3}}
    END{flush()}
  ' "$tmp_dir/candidates.sorted.tsv" > "$tmp_dir/resolution.tsv" || return $?

  join_fields="1.1";for ((i=2;i<=ncols+1;i++));do join_fields+=",1.$i";done;join_fields+=",2.2,2.3,2.4"
  printf '%s\n' "$(<"$tmp_dir/header.tsv")" > "$tmp_dir/output.tsv"
  [[ -z "$audit" ]] || printf 'ROW\tSNP_GWAS\tSNP_REFERENCE\tCHR\tPOS\tEA\tNEA\tSTATUS\tMATCH_TYPE\n' > "$tmp_dir/audit.tsv"
  [[ -z "$unmatched" ]] || printf '%s\tMATCH_STATUS\n' "$(<"$tmp_dir/header.tsv")" > "$tmp_dir/unmatched.tsv"
  LC_ALL=C join -t $'\t' -a 1 -e '' -o "$join_fields" "$tmp_dir/query.rows.tsv" "$tmp_dir/resolution.tsv" |
    awk -v FS='\t' -v OFS='\t' -v n="$ncols" -v sc="$snp_col" -v cc="$chr_col" -v pc="$pos_col" -v ec="$ea_col" -v nc="$nea_col" \
        -v out="$tmp_dir/output.tsv" -v audit="$tmp_dir/audit.tsv" -v unmatched="$tmp_dir/unmatched.tsv" -v summary="$tmp_dir/summary.tsv" '
      function rc(x, i,c,y){x=toupper(x);y="";for(i=length(x);i>0;i--){c=substr(x,i,1);if(c=="A")c="T";else if(c=="T")c="A";else if(c=="C")c="G";else if(c=="G")c="C";else return x;y=y c}return y}
      {
        row=$1+0;status=$(n+2);refid=$(n+3);type=$(n+4);if(status=="")status="UNMATCHED"
        qid=$(sc+1);qchr=$(cc+1);qpos=$(pc+1);qea=ec?$(ec+1):"";qnea=nc?$(nc+1):""
        if(status=="MATCHED"){
          $(sc+1)=refid;if(type~/STRAND$/){if(ec)$(ec+1)=rc($(ec+1));if(nc)$(nc+1)=rc($(nc+1))}
          line=$2;for(j=3;j<=n+1;j++)line=line OFS $j;print line >> out;matched++
        }else{if(status=="AMBIGUOUS")ambiguous++;else unmatched_n++}
        if(audit!="")print row,qid,(refid==""?".":refid),qchr,qpos,qea,qnea,status,(type==""?".":type) >> audit
        if(status!="MATCHED"&&unmatched!=""){line=$2;for(j=3;j<=n+1;j++)line=line OFS $j;print line,status >> unmatched}
      }
      END{print matched+0,unmatched_n+0,ambiguous+0 > summary}
    ' || return $?
  read -r matched_n unmatched_n ambiguous_n < "$tmp_dir/summary.tsv"

  if [[ "$output" == - ]]; then cat "$tmp_dir/output.tsv"; else mkdir -p -- "$(dirname -- "$output")"; mv -f -- "$tmp_dir/output.tsv" "$output"; fi
  if [[ -n "$audit" ]]; then mkdir -p -- "$(dirname -- "$audit")"; mv -f -- "$tmp_dir/audit.tsv" "$audit"; fi
  if [[ -n "$unmatched" ]]; then mkdir -p -- "$(dirname -- "$unmatched")"; mv -f -- "$tmp_dir/unmatched.tsv" "$unmatched"; fi
  echo "match_SNP: engine=external-sort-merge, reference=$ref_kind, total=$qn, matched=$matched_n, unmatched=$unmatched_n, ambiguous=$ambiguous_n, sort_memory=${sort_memory_mb}MiB" >&2
)

# Fill missing EAF values from a chromosome-split PLINK 2 reference.
match_EAF() {
  local query="" ref="" keep="" output="" audit="" tmp_parent tmp_dir map tmp_out header
  local chr label bed_unsorted bed prefix freq_prefix afreq plink_log rc=0
  local -a pfile_modifier=() keep_args=()

  while (( $# )); do
    case "$1" in
      -r|--reference) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; ref=$2; shift 2 ;;
      -o|--output) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; output=$2; shift 2 ;;
      --keep) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; keep=$2; shift 2 ;;
      --audit) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; audit=$2; shift 2 ;;
      --) shift; [[ $# -eq 1 ]] || { echo "ERROR: match_EAF requires one QUERY file" >&2; return 2; }; query=$1; shift ;;
      -*) echo "ERROR: unknown match_EAF option: $1" >&2; return 2 ;;
      *) [[ -z "$query" ]] || { echo "ERROR: match_EAF accepts one QUERY file" >&2; return 2; }; query=$1; shift ;;
    esac
  done

  [[ -n "$query" && -s "$query" ]] || { echo "ERROR: match_EAF QUERY is missing/empty: ${query:-<none>}" >&2; return 2; }
  [[ -n "$ref" ]] || { echo "ERROR: match_EAF requires --reference" >&2; return 2; }
  [[ -n "$output" ]] || { echo "ERROR: match_EAF requires --output" >&2; return 2; }
  [[ -z "$keep" || -s "$keep" ]] || { echo "ERROR: match_EAF keep file is missing/empty: $keep" >&2; return 2; }
  [[ -z "$keep" ]] || keep_args=(--keep "$keep")
  command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 is required by match_EAF" >&2; return 2; }
  command -v join >/dev/null 2>&1 || { echo "ERROR: join is required by match_EAF" >&2; return 2; }

  tmp_parent=$(dirname -- "$output")
  mkdir -p -- "$tmp_parent"
  tmp_dir=$(mktemp -d "$tmp_parent/.match_EAF.XXXXXX") || return 2
  map="$tmp_dir/reference.eaf.map.tsv"
  : > "$map"
  [[ -n "$audit" ]] || audit="$tmp_dir/audit.tsv"

  # Collect only coordinates whose EAF actually needs filling. The first pass
  # also validates the standardized header; at most 23 files stay open.
  if ! phe_zcat "$query" | awk -v FS='\t' -v OFS='\t' -v d="$tmp_dir" '
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="M"||x=="MT")x=25;return x}
    NR==1{
      for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}
      req[1]="SNP";req[2]="CHR";req[3]="POS";req[4]="EA";req[5]="NEA";req[6]="EAF"
      for(i=1;i<=6;i++)if(!(req[i] in c)){print "ERROR: match_EAF QUERY lacks " req[i] > "/dev/stderr";bad=1}
      if(bad)exit 2
      next
    }
    {
      f=$(c["EAF"]);if(isnum(f)&&f+0>=0&&f+0<=1)next
      chr=normchr($(c["CHR"]));pos=$(c["POS"])
      if(chr~/^([1-9]|1[0-9]|2[0-3])$/&&pos~/^[0-9]+$/&&pos+0>0){lab=(chr==23?"X":chr);print lab,pos,pos >> (d "/chr" lab ".bed.unsorted")}
    }
  '; then
    rc=$?
    rm -rf -- "$tmp_dir"
    return "$rc"
  fi

  for chr in {1..22} 23; do
    [[ "$chr" == 23 ]] && label=X || label=$chr
    bed_unsorted="$tmp_dir/chr${label}.bed.unsorted"
    [[ -s "$bed_unsorted" ]] || continue
    bed="$tmp_dir/chr${label}.bed"
    LC_ALL=C sort -t $'\t' -k2,2n -u "$bed_unsorted" > "$bed"
    rm -f -- "$bed_unsorted"

    if [[ -d "$ref" || "$ref" == */ ]]; then
      prefix="${ref%/}/chr${label}"
    elif [[ -s "${ref}chr${label}.pgen" ]]; then
      prefix="${ref}chr${label}"
    else
      prefix="${ref}${label}"
    fi
    pfile_modifier=()
    if [[ -s "${prefix}.pvar.zst" ]]; then
      pfile_modifier=(vzs)
    elif [[ ! -s "${prefix}.pvar" ]]; then
      echo "ERROR: match_EAF reference has no .pvar or .pvar.zst: $prefix" >&2
      rm -rf -- "$tmp_dir"
      return 2
    fi
    [[ -s "${prefix}.pgen" && -s "${prefix}.psam" ]] || {
      echo "ERROR: match_EAF reference pfile is incomplete: $prefix" >&2
      rm -rf -- "$tmp_dir"
      return 2
    }

    freq_prefix="$tmp_dir/chr${label}.freq"
    plink_log="$tmp_dir/chr${label}.plink.stdout"
    if ! plink2 --pfile "$prefix" "${pfile_modifier[@]}" --extract bed1 "$bed" \
        "${keep_args[@]}" --freq cols=chrom,pos,ref,alt,eq,nobs --out "$freq_prefix" > "$plink_log" 2>&1; then
      cat "$plink_log" >&2
      rm -rf -- "$tmp_dir"
      return 1
    fi
    afreq="${freq_prefix}.afreq"
    [[ -s "$afreq" ]] || { echo "ERROR: match_EAF did not create $afreq" >&2; rm -rf -- "$tmp_dir"; return 1; }
    awk -v FS='\t' -v OFS='\t' '
      function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="M"||x=="MT")x=25;return x}
      NR==1{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}next}
      {chr=normchr($(c["CHROM"]));pos=$(c["POS"]);if(chr~/^[0-9]+$/&&pos~/^[0-9]+$/)printf "%02d:%012d\t%s\t%s\t%s\n",chr,pos,$(c["REF"]),$(c["ALT"]),$(c["FREQS"])}
    ' "$afreq" >> "$map"
  done

  header=$(set +o pipefail; phe_zcat "$query" | head -n 1 | tr -d '\r' || true)
  tmp_out="${output}.tmp.$$"

  # Stream the coordinate-sorted join into allele resolution and gzip.
  if (
    set -o pipefail
    {
      printf '%s\n' "$header"
      LC_ALL=C join -t $'\t' -a 1 -e '' \
        -o '1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4' \
        <(phe_zcat "$query" | awk -v FS='\t' -v OFS='\t' '
          function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="M"||x=="MT")x=25;return x}
          NR==1{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}next}
          {chr=normchr($(c["CHR"]));pos=$(c["POS"]);if(chr!~/^[0-9]+$/)chr=0;if(pos!~/^[0-9]+$/)pos=0;printf "%02d:%012d\t%d\t%s\n",chr,pos,NR-1,$0}
        ') "$map" | awk -v FS='\t' -v OFS='\t' -v audit="$audit" '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
      function rc(x, i,c,y){x=toupper(x);y="";for(i=length(x);i>0;i--){c=substr(x,i,1);if(c=="A")c="T";else if(c=="T")c="A";else if(c=="C")c="G";else if(c=="G")c="C";else return "";y=y c}return y}
      function addfreq(x, k){if(!isnum(x)||x+0<0||x+0>1)return;k=sprintf("%.12g",x+0);if(!(k in cand)){cand[k]=x;nc++}}
      function consider(freqs, ea,nea, n,a,b,i,direct,cea,cnea){
        delete af;n=split(freqs,a,",");for(i=1;i<=n;i++){split(a[i],b,"=");if(length(b[1])&&isnum(b[2]))af[toupper(b[1])]=b[2]}
        ea=toupper(ea);nea=toupper(nea);direct=((ea in af)&&(nea in af))
        if(direct){addfreq(af[ea]);return}
        cea=rc(ea);cnea=rc(nea);if(cea!=""&&(cea in af)&&(cnea in af))addfreq(af[cea])
      }
      function flush( i,k){
        if(cur=="")return
        f=q[6]
        if(isnum(f)&&f+0>=0&&f+0<=1)existing++
        else if(nc==1){for(k in cand)q[6]=cand[k];filled++}
        else if(nc>1)ambiguous++
        else unmatched++
        print q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11]
      }
      {
        if($1!=cur){flush();cur=$1;delete q;delete cand;nc=0;for(i=1;i<=11;i++)q[i]=$(i+1)}
        if($15!=""&&!(isnum(q[6])&&q[6]+0>=0&&q[6]+0<=1))consider($15,q[4],q[5])
      }
      END{
        flush()
        print "STATUS","N" > audit
        print "existing",existing+0 >> audit
        print "filled_from_1KG",filled+0 >> audit
        print "unmatched",unmatched+0 >> audit
        print "ambiguous",ambiguous+0 >> audit
      }
      '
    } | gzip -c > "$tmp_out"
  ); then
    rc=0
  else
    rc=$?
    rm -f -- "$audit"
  fi

  if (( rc == 0 )); then
    gzip -t "$tmp_out" || rc=$?
  fi
  if (( rc == 0 )); then
    mv -f -- "$tmp_out" "$output"
    [[ "$audit" == "$tmp_dir/audit.tsv" ]] || mkdir -p -- "$(dirname -- "$audit")"
    echo "match_EAF: reference=$ref output=$output summary=$audit" >&2
  else
    rm -f -- "$tmp_out"
  fi
  rm -rf -- "$tmp_dir"
  return "$rc"
}

# Replace query SNP IDs from a reference table; coordinate fallback requires one build.
match_GRCH() (
  local query="" ref="" output="-" audit="" unmatched="" position_only=0
  local ref_file="" tmp_dir="" ref_pvar="" rc=0
  local memory_mb="${MATCH_SNP_MEMORY_MB:-4096}" sort_memory_mb="${MATCH_SNP_SORT_MEMORY_MB:-512}"
  local sort_threads="${MATCH_SNP_SORT_THREADS:-1}" match_tmp_dir="${MATCH_SNP_TMP_DIR:-/mnt/d/tmp}"
  local warn_query_rows="${MATCH_SNP_WARN_QUERY_ROWS:-500000}" max_query_rows="${MATCH_SNP_MAX_QUERY_ROWS:-5000000}"
  local allow_large_query=0

  while (( $# )); do
    case "$1" in
      -r|--reference) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; ref=$2; shift 2 ;;
      -o|--output) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; output=$2; shift 2 ;;
      --audit) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; audit=$2; shift 2 ;;
      --unmatched) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; return 2; }; unmatched=$2; shift 2 ;;
      --position-only) position_only=1; shift ;;
      --memory-mb) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires MiB" >&2; return 2; }; memory_mb=$2; shift 2 ;;
      --sort-memory-mb) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires MiB" >&2; return 2; }; sort_memory_mb=$2; shift 2 ;;
      --sort-threads) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an integer" >&2; return 2; }; sort_threads=$2; shift 2 ;;
      --tmp-dir) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a directory" >&2; return 2; }; match_tmp_dir=$2; shift 2 ;;
      --warn-query-rows) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an integer" >&2; return 2; }; warn_query_rows=$2; shift 2 ;;
      --max-query-rows) [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an integer" >&2; return 2; }; max_query_rows=$2; shift 2 ;;
      --allow-large-query) allow_large_query=1; shift ;;
      --) shift; [[ $# -eq 1 ]] || { echo "ERROR: match_GRCH requires one QUERY file" >&2; return 2; }; query=$1; shift ;;
      -*) echo "ERROR: unknown match_GRCH option: $1" >&2; return 2 ;;
      *) [[ -z "$query" ]] || { echo "ERROR: match_GRCH accepts one QUERY file" >&2; return 2; }; query=$1; shift ;;
    esac
  done

  [[ -n "$query" && -s "$query" ]] || { echo "ERROR: match_GRCH QUERY is missing/empty: ${query:-<none>}" >&2; return 2; }
  [[ -n "$ref" ]] || { echo "ERROR: match_GRCH requires --reference REF" >&2; return 2; }
  [[ "$memory_mb" =~ ^[0-9]+$ ]] || { echo "ERROR: --memory-mb must be a non-negative integer: $memory_mb" >&2; return 2; }
  if (( memory_mb > 0 )); then
    ulimit -Sv "$((memory_mb * 1024))" || { echo "ERROR: cannot set match_GRCH memory cap to ${memory_mb} MiB" >&2; return 2; }
  fi

  case "$ref" in
    *.pvar|*.pvar.gz|*.pvar.bgz|*.pvar.zst|*.bim|*.bim.gz|*.bim.bgz)
      ref_file=$ref
      ;;
    *)
      if [[ ! -s "$ref" ]]; then
        for ref_file in "$ref.pvar" "$ref.pvar.gz" "$ref.pvar.bgz" "$ref.pvar.zst" "$ref.bim" "$ref.bim.gz" "$ref.bim.bgz"; do
          [[ -s "$ref_file" ]] && break
        done
        [[ -s "$ref_file" ]] || { echo "ERROR: match_GRCH reference is missing/empty: $ref" >&2; return 2; }
      else
        mkdir -p -- "$match_tmp_dir" || return 2
        tmp_dir=$(mktemp -d "${match_tmp_dir%/}/match_GRCH.XXXXXX") || return 2
        trap '[[ -z "$tmp_dir" ]] || rm -rf -- "$tmp_dir"' EXIT HUP INT TERM
        ref_pvar="$tmp_dir/reference.pvar"
        if phe_zcat "$ref" | awk -v FS='\t' -v OFS='\t' '
          function hname(x){x=toupper(x);sub(/^#/,"",x);gsub(/\r/,"",x);return x}
          function valid_allele(x){return x!=""&&x!="."&&x!="NA"&&x~/^[ACGT]+$/}
          NR==1{
            for(i=1;i<=NF;i++){
              x=hname($i)
              if(x=="SNP"||x=="RSID"||x=="ID"||x=="VARIANT_ID")snp=i
              else if(x=="CHR"||x=="CHROM"||x=="CHROMOSOME")chr=i
              else if(x=="POS"||x=="BP"||x=="POSITION"||x=="BASE_PAIR_LOCATION")pos=i
              else if(x=="EA"||x=="A1"||x=="ALT"||x=="EFFECT_ALLELE")ea=i
              else if(x=="NEA"||x=="A2"||x=="REF"||x=="OTHER_ALLELE"||x=="NON_EFFECT_ALLELE"||x=="REFERENCE_ALLELE")nea=i
            }
            if(!snp||!chr||!pos||!ea||!nea){print "ERROR: match_GRCH REFERENCE requires SNP, CHR, POS, EA, and NEA columns" > "/dev/stderr";fatal=1;exit 2}
            print "#CHROM","POS","ID","REF","ALT";next
          }
          {
            id=$snp;ch=$chr;bp=$pos;a1=toupper($ea);a2=toupper($nea);sub(/^chr/,"",ch)
            if(id!=""&&id!="."&&bp~/^[0-9]+$/&&valid_allele(a1)&&valid_allele(a2))print ch,bp,id,a2,a1
          }
        ' > "$ref_pvar"; then
          :
        else
          rc=$?
          rm -rf -- "$tmp_dir"
          return "$rc"
        fi
        ref_file=$ref_pvar
      fi
      ;;
  esac

  local -a args=(--reference "$ref_file" --output "$output" --id-any-position
    --memory-mb "$memory_mb" --sort-memory-mb "$sort_memory_mb" --sort-threads "$sort_threads"
    --tmp-dir "$match_tmp_dir" --warn-query-rows "$warn_query_rows" --max-query-rows "$max_query_rows")
  [[ -z "$audit" ]] || args+=(--audit "$audit")
  [[ -z "$unmatched" ]] || args+=(--unmatched "$unmatched")
  (( position_only == 0 )) || args+=(--position-only)
  (( allow_large_query == 0 )) || args+=(--allow-large-query)
  match_SNP "${args[@]}" "$query" || rc=$?
  [[ -z "$tmp_dir" ]] || rm -rf -- "$tmp_dir"
  return "$rc"
)
