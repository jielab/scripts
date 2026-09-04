#!/usr/bin/env bash
# Performance overrides for gwas_post.sh at multi-thousand-GWAS scale.
# This file is sourced by each generated per-GWAS command after 0data.f.sh.

# Preserve the original fill function so integrated format-time N filling can
# skip only its redundant rewrite while retaining match_EAF behavior.
if declare -F gwas_post_fill_missing_fields >/dev/null 2>&1 &&
   ! declare -F gwas_post_fill_missing_fields_legacy >/dev/null 2>&1; then
  eval "$(declare -f gwas_post_fill_missing_fields | sed '1s/gwas_post_fill_missing_fields/gwas_post_fill_missing_fields_legacy/')"
fi
if declare -F gwas_clean_liftover_small_to_final >/dev/null 2>&1 &&
   ! declare -F gwas_clean_liftover_small_to_final_legacy >/dev/null 2>&1; then
  eval "$(declare -f gwas_clean_liftover_small_to_final | sed '1s/gwas_clean_liftover_small_to_final/gwas_clean_liftover_small_to_final_legacy/')"
fi

: "${GWAS_POST_SORT_MEMORY:=512M}"
: "${GWAS_CLEAN_COMP_THREADS:=1}"
GWAS_POST_N_FILLED_DURING_FORMAT=FALSE
GWAS_POST_VIEWS_READY=FALSE
GWAS_POST_MAGMA_ROWS=""
GWAS_POST_MPLOT_INPUT=""
GWAS_POST_LEAD_VIEW=""

gwas_clean_compress() {
  if command -v bgzip >/dev/null 2>&1; then
    bgzip -@ "$GWAS_CLEAN_COMP_THREADS" -c
  elif command -v pigz >/dev/null 2>&1; then
    pigz -c -p "$GWAS_CLEAN_COMP_THREADS"
  else
    gzip -c
  fi
}

gwas_clean_gzip_ok() {
  local file="$1"
  [[ -s "$file" ]] || return 1
  if { [[ -s "${file}.tbi" ]] || [[ -s "${file}.csi" ]]; } &&
     command -v tabix >/dev/null 2>&1; then
    tabix -l "$file" >/dev/null 2>&1
  else
    gzip -t "$file" >/dev/null 2>&1
  fi
}

gwas_post_sort_bgzf() {
  local output="$1" chr_col="$2" pos_col="$3" tmp
  tmp="${output}.tmp.$$"
  {
    IFS= read -r header
    printf '%s\n' "$header"
    sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' \
      -k"${chr_col}","${chr_col}"n -k"${pos_col}","${pos_col}"n -k1,1
  } | gwas_clean_compress > "$tmp"
  [[ -s "$tmp" ]] || { echo "ERROR: empty BGZF output: $output" >&2; return 1; }
  mv -f -- "$tmp" "$output"
}

# Full-format mode: one raw scan does schema conversion, coordinate filtering,
# optional N filling, and QC counting. No post-format counting scan is needed.
std_format() {
  local src="$1" out="$2" header_out="$3" input_fs tmp audit
  gwas_post_need_file "$src"
  gwas_post_log "Format full GWAS (single scan): $src -> $out"
  gwas_clean_header_names "$src" "$header_out"
  SNP_col=$(gwas_clean_col "${SNP_col:-}"); CHR_col=$(gwas_clean_col "${CHR_col:-}")
  POS_col=$(gwas_clean_col "${POS_col:-}"); EA_col=$(gwas_clean_col "${EA_col:-}")
  NEA_col=$(gwas_clean_col "${NEA_col:-}"); EAF_col=$(gwas_clean_col "${EAF_col:-}")
  N_col=$(gwas_clean_col "${N_col:-}"); BETA_col=$(gwas_clean_col "${BETA_col:-}")
  SE_col=$(gwas_clean_col "${SE_col:-}"); P_col=$(gwas_clean_col "${P_col:-}")
  LOG10P_col=$(gwas_clean_col "${LOG10P_col:-}")
  [[ "$SNP_col" -gt 0 && "$CHR_col" -gt 0 && "$POS_col" -gt 0 && "$EA_col" -gt 0 &&
     "$NEA_col" -gt 0 && "$BETA_col" -gt 0 && "$SE_col" -gt 0 && "$P_col" -gt 0 ]] || {
    echo "ERROR: required GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: $src" >&2
    return 1
  }
  input_fs=$(gwas_clean_detect_fs "$src")
  tmp="${out}.tmp.$$"
  audit="${QC_PREFIX}.format.audit.tsv"
  gwas_post_zcat "$src" | awk -v FS="$input_fs" -v OFS='\t' \
    -v snp_col="$SNP_col" -v chr_col="$CHR_col" -v pos_col="$POS_col" \
    -v ea_col="$EA_col" -v nea_col="$NEA_col" -v eaf_col="$EAF_col" -v n_col="$N_col" \
    -v beta_col="$BETA_col" -v se_col="$SE_col" -v p_col="$P_col" -v logp_col="$LOG10P_col" \
    -v fill_n="$FILL_N" -v audit="$audit" -v gwas="$GWAS" '
    function get(c, x){x=(c>0 ? $c : "");gsub(/\r/,"",x);return x}
    function val(c, x){x=get(c);return x=="" ? "NA" : x}
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function validn(x){return x ~ /^[0-9]+([.][0-9]+)?$/ && x+0>0}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P";next}
    {
      chr=normchr(get(chr_col));pos=get(pos_col)
      if(chr!~/^[0-9]+$/ || chr+0<1 || chr+0>25 || pos!~/^[0-9]+$/ || pos+0<1){badcoord++;next}
      snp=get(snp_col);ea=get(ea_col);nea=get(nea_col)
      if(snp==""||snp=="NA"||snp=="."){snp=chr ":" pos;if(ea!="")snp=snp ":" ea;if(nea!="")snp=snp ":" nea}
      p=(p_col>0 ? get(p_col) : "");lp=(logp_col>0 ? get(logp_col) : "")
      if(p=="" && isnum(lp))p=10^(-lp);if(lp=="" && isnum(p) && p>0)lp=-log(p)/log(10)
      if(ea=="")ea="NA";if(nea=="")nea="NA";if(p=="")p="NA";if(lp=="")lp="NA"
      n=val(n_col);if(fill_n!="" && !validn(n)){n=fill_n;nfilled++}else if(validn(n))nkept++
      print snp,chr+0,pos+0,ea,nea,val(eaf_col),n,val(beta_col),val(se_col),p,lp
      kept++
    }
    END{
      print "GWAS\tN_OUTPUT\tN_DROPPED_BAD_COORD\tN_EXISTING\tN_FILLED" > audit
      print gwas "\t" kept+0 "\t" badcoord+0 "\t" nkept+0 "\t" nfilled+0 >> audit
    }' | {
      IFS= read -r format_header
      printf '%s\n' "$format_header"
      sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k2,2n -k3,3n -k1,1
    } | gwas_clean_compress > "$tmp"
  [[ -s "$tmp" && -s "$audit" ]] || { echo "ERROR: format failed: $src" >&2; rm -f -- "$tmp"; return 1; }
  mv -f -- "$tmp" "$out"
  rm -f -- "${out}.tbi" "${out}.csi"
  awk -F '\t' 'NR==2{print $1 "\t" $2}' "$audit" > "${QC_PREFIX}.format.nrow.tsv"
  [[ -z "$FILL_N" ]] || GWAS_POST_N_FILLED_DURING_FORMAT=TRUE
}

# Skip auxiliary panels such as chrXY, chrX.female, and chrX.male while
# enumerating a reference directory.  lead/MAGMA use only the canonical
# chr1..chr22/chrX/chrY panels; otherwise several prefixes can collapse to the
# same normalized chromosome and overwrite the same per-chromosome files.
ref_clump_pfiles() {
  local ref="$1" f prefix label chr found=FALSE
  if [[ -f "${ref}.pgen" && -f "${ref}.psam" && ( -f "${ref}.pvar" || -f "${ref}.pvar.zst" ) ]]; then
    label=$(basename "$ref"); chr=$(ref_chr "$label") || return 1
    printf '%s %s %s\n' "$ref" "$label" "$chr"; return 0
  fi
  if [[ -d "$ref" ]]; then
    while IFS= read -r f; do
      prefix=${f%.pgen}
      [[ -f "${prefix}.psam" && ( -f "${prefix}.pvar" || -f "${prefix}.pvar.zst" ) ]] || continue
      label=$(basename "$prefix"); chr=$(ref_chr "$label") || continue
      [[ "$label" == "$(chr_label "$chr")" ]] || continue
      printf '%s %s %s\n' "$prefix" "$label" "$chr"; found=TRUE
    done < <(find "$ref" -maxdepth 1 -type f -name 'chr*.pgen' | sort -V)
  else
    while IFS= read -r f; do
      prefix=${f%.pgen}
      [[ -f "${prefix}.psam" && ( -f "${prefix}.pvar" || -f "${prefix}.pvar.zst" ) ]] || continue
      label=$(basename "$prefix"); chr=$(ref_chr "$label") || continue
      [[ "$label" == "$(chr_label "$chr")" ]] || continue
      printf '%s %s %s\n' "$prefix" "$label" "$chr"; found=TRUE
    done < <(compgen -G "${ref}chr*.pgen" | sort -V)
  fi
  [[ "$found" == TRUE ]] || { echo "ERROR: no chromosome pfile reference found: $ref" >&2; return 1; }
}

# --small TRUE: filter and standardize in one raw scan. Keep exactly the union
# of HapMap3 IDs/coordinates and variants meeting --p-small, then deduplicate
# exact rows.  Coordinates are the fallback for files whose rsID column is NA.
std_small() {
  local src="$1" out="$2" header_out="$3" input_fs tmp_tsv tmp audit n_out
  gwas_clean_need_file "$src"; gwas_clean_need_file "$HM3"
  [[ -z "${HM3_POS:-}" ]] || gwas_clean_need_file "$HM3_POS"
  gwas_post_log "Small GWAS (single scan): $src -> $out"
  gwas_clean_header_names "$src" "$header_out"
  SNP_col=$(gwas_clean_col "${SNP_col:-}"); CHR_col=$(gwas_clean_col "${CHR_col:-}")
  POS_col=$(gwas_clean_col "${POS_col:-}"); EA_col=$(gwas_clean_col "${EA_col:-}")
  NEA_col=$(gwas_clean_col "${NEA_col:-}"); EAF_col=$(gwas_clean_col "${EAF_col:-}")
  N_col=$(gwas_clean_col "${N_col:-}"); BETA_col=$(gwas_clean_col "${BETA_col:-}")
  SE_col=$(gwas_clean_col "${SE_col:-}"); P_col=$(gwas_clean_col "${P_col:-}")
  LOG10P_col=$(gwas_clean_col "${LOG10P_col:-}")
  [[ "$CHR_col" -gt 0 && "$POS_col" -gt 0 && ( "$P_col" -gt 0 || "$LOG10P_col" -gt 0 ) ]] || {
    echo "ERROR: CHR/POS and P/LOG10P are required: $src" >&2; return 1;
  }
  input_fs=$(gwas_clean_detect_fs "$src")
  tmp_tsv="$GWAS_POST_TMP/${GWAS}.small.sorted.tsv"
  tmp="${out}.tmp.$$"
  audit="${QC_PREFIX}.small.audit.tsv"
  gwas_clean_zcat "$src" | awk -v FS="$input_fs" -v OFS='\t' \
    -v hm3_file="$HM3" -v hm3_pos_file="${HM3_POS:-}" -v pthr="$P_SMALL" \
    -v snp_col="$SNP_col" -v chr_col="$CHR_col" -v pos_col="$POS_col" \
    -v ea_col="$EA_col" -v nea_col="$NEA_col" -v eaf_col="$EAF_col" -v n_col="$N_col" \
    -v beta_col="$BETA_col" -v se_col="$SE_col" -v p_col="$P_col" -v logp_col="$LOG10P_col" \
    -v fill_n="$FILL_N" -v audit="$audit" -v gwas="$GWAS" '
    function get(c,x){x=(c>0?$c:"");gsub(/\r/,"",x);return x}
    function val(c,x){x=get(c);return x==""?"NA":x}
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function validn(x){return x~/^[0-9]+([.][0-9]+)?$/&&x+0>0}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    BEGIN{
      while((getline x<hm3_file)>0){gsub(/\r/,"",x);split(x,a,/[ \t]+/);if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1}close(hm3_file)
      if(hm3_pos_file!="")while((getline x<hm3_pos_file)>0){
        gsub(/\r/,"",x);split(x,a,/[ \t]+/);ch=normchr(a[1]);bp=a[4]
        key=ch SUBSEP (bp+0)
        if(ch~/^[0-9]+$/&&bp~/^[0-9]+$/&&bp+0>0&&!(key in hm3pos))hm3pos[key]=a[2]
      }
      if(hm3_pos_file!="")close(hm3_pos_file)
    }
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P";next}
    {
      chr=normchr(get(chr_col));pos=get(pos_col);if(chr!~/^[0-9]+$/||chr+0<1||chr+0>25||pos!~/^[0-9]+$/||pos+0<1){badcoord++;next}
      key=chr SUBSEP (pos+0);in_hm3_pos=(key in hm3pos)
      snp=get(snp_col);ea=get(ea_col);nea=get(nea_col);if(snp==""||snp=="NA"||snp=="."){
        if(in_hm3_pos&&hm3pos[key]!="")snp=hm3pos[key]
        else{snp=chr ":" pos;if(ea!="")snp=snp ":" ea;if(nea!="")snp=snp ":" nea}
      }
      p=(p_col>0?get(p_col):"");lp=(logp_col>0?get(logp_col):"");if(p==""&&isnum(lp))p=10^(-lp);if(lp==""&&isnum(p)&&p>0)lp=-log(p)/log(10)
      keep=(snp in hm3)||in_hm3_pos||(isnum(p)&&p+0<=pthr+0);if(!keep)next
      if(ea=="")ea="NA";if(nea=="")nea="NA";if(p=="")p="NA";if(lp=="")lp="NA"
      n=val(n_col);if(fill_n!=""&&!validn(n)){n=fill_n;nfilled++}else if(validn(n))nkept++
      print snp,chr+0,pos+0,ea,nea,val(eaf_col),n,val(beta_col),val(se_col),p,lp;kept++
    }
    END{print "GWAS\tN_OUTPUT\tN_DROPPED_BAD_COORD\tN_EXISTING\tN_FILLED">audit;print gwas "\t" kept+0 "\t" badcoord+0 "\t" nkept+0 "\t" nfilled+0>>audit}
  ' | {
    IFS= read -r header; printf '%s\n' "$header"
    sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' \
      -k2,2n -k3,3n -k1,1 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 -k11,11 -u
  } > "$tmp_tsv"
  n_out=$(( $(wc -l < "$tmp_tsv") - 1 )); (( n_out < 0 )) && n_out=0
  gwas_clean_compress < "$tmp_tsv" > "$tmp"
  [[ -s "$tmp" ]] || { echo "ERROR: small format failed: $src" >&2; return 1; }
  mv -f -- "$tmp" "$out"
  rm -f -- "${out}.tbi" "${out}.csi"
  printf '%s\t%s\n' "$GWAS" "$n_out" > "${QC_PREFIX}.small.nrow.tsv"
  [[ -z "$FILL_N" ]] || GWAS_POST_N_FILLED_DURING_FORMAT=TRUE
}

if [[ "${SMALL_MODE:-FALSE}" == FALSE ]]; then
  std_small() { std_format "$@"; }
fi

gwas_post_fill_missing_fields() {
  local target="$1" saved_fill_n="$FILL_N"
  if [[ "$GWAS_POST_N_FILLED_DURING_FORMAT" == TRUE ]]; then FILL_N=""; fi
  gwas_post_fill_missing_fields_legacy "$target"
  FILL_N="$saved_fill_n"
  if [[ -n "$saved_fill_n" || "$FILL_EAF" == TRUE ]]; then
    rm -f -- "${target}.tbi" "${target}.csi"
  fi
}

gwas_clean_liftover_small_to_final() {
  gwas_clean_liftover_small_to_final_legacy
  rm -f -- "${FINAL}.tbi" "${FINAL}.csi"
}

gwas_post_ensure_index() {
  local file="$1"
  [[ -s "$file" ]] || return 0
  if ! gwas_index_valid "$file"; then gwas_index_file "$file"; fi
}

gwas_post_cis_regions() {
  local output="$1"
  awk -v name="$GWAS" -v flank="$CIS_FLANK" 'BEGIN{OFS="\t"}
    function chr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    /^[[:space:]]*($|#)/{next}
    $4==name{s=$2-flank;if(s<1)s=1;e=$3+flank;print chr($1),s,e}
  ' "$CIS_BED" | sort -k1,1V -k2,2n -k3,3n | awk 'BEGIN{OFS="\t"}
    NR==1{c=$1;s=$2;e=$3;next}
    $1==c&&$2<=e+1{if($3>e)e=$3;next}
    {print c,s,e;c=$1;s=$2;e=$3}
    END{if(NR)print c,s,e}' > "$output"
}

# Indexed cis extraction is proportional to the requested loci, not to the
# whole GWAS. A full-scan fallback remains for legacy/non-indexable inputs.
gwas_clean_make_cis() {
  local regions="$GWAS_POST_TMP/${GWAS}.cis.regions.tsv" body="$GWAS_POST_TMP/${GWAS}.cis.body.tsv" tmp header
  [[ -n "${CIS_BED:-}" ]] || { echo "ERROR: --cis-bed is required for cis output" >&2; return 2; }
  if [[ "$REPLACE" != TRUE && -s "$CIS_OUT" ]]; then
    if gwas_index_valid "$CIS_OUT"; then
      gwas_clean_log "indexed cis file exists: $CIS_OUT"; return 0
    fi
    gwas_clean_log "repair missing/stale cis BGZF index: $CIS_OUT"
    if gwas_index_file "$CIS_OUT" && gwas_index_valid "$CIS_OUT"; then
      gwas_clean_log "indexed cis file ready: $CIS_OUT"; return 0
    fi
    echo "ERROR: failed to repair cis BGZF/index: $CIS_OUT" >&2
    return 1
  fi
  gwas_clean_need_file "$FINAL"; gwas_clean_need_file "$CIS_BED"
  gwas_post_cis_regions "$regions"
  tmp="${CIS_OUT}.tmp.$$"
  header=$(gwas_index_header "$FINAL")
  [[ -n "$header" ]] || { echo "ERROR: unreadable GWAS header: $FINAL" >&2; return 1; }
  : > "$body"
  if [[ -s "$regions" ]]; then
    gwas_post_ensure_index "$FINAL"
    while IFS=$'\t' read -r chr start end; do
      tabix "$FINAL" "${chr}:${start}-${end}"
    done < "$regions" > "$body"
  fi
  { printf '%s\n' "$header"; sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k2,2n -k3,3n -k1,1 -u "$body"; } |
    gwas_clean_compress > "$tmp"
  mv -f -- "$tmp" "$CIS_OUT"
  rm -f -- "${CIS_OUT}.tbi" "${CIS_OUT}.csi"
  gwas_index_file "$CIS_OUT"
  printf '%s\t%s\n' "$GWAS" "$(wc -l < "$body" | tr -d ' ')" > "${QC_PREFIX}.cis.nrow.tsv"
  gwas_clean_log "cis subset via tabix: $FINAL -> $CIS_OUT"
}

gwas_post_prepare_views() {
  local want_magma=FALSE want_mplot=FALSE want_lead=FALSE raw_plot plot_tsv
  [[ "$GWAS_POST_VIEWS_READY" == TRUE ]] && return 0
  [[ ",${DO_STEP}," == *,magma,* ]] && want_magma=TRUE
  [[ "$DO_STEP" == all || ",${DO_STEP}," == *,mplot,* ]] && want_mplot=TRUE
  [[ "$DO_STEP" == all || ",${DO_STEP}," == *,lead,* ]] && want_lead=TRUE
  [[ "$want_magma" == TRUE || "$want_mplot" == TRUE || "$want_lead" == TRUE ]] || { GWAS_POST_VIEWS_READY=TRUE; return 0; }
  gwas_post_need_file "$FINAL"
  if [[ "$want_mplot" == TRUE ]]; then
    gwas_post_need_file "$HM3"
    [[ -z "${HM3_POS:-}" ]] || gwas_post_need_file "$HM3_POS"
  fi
  GWAS_POST_MAGMA_ROWS="$GWAS_POST_TMP/${GWAS}.magma.rows.tsv"
  GWAS_POST_LEAD_VIEW="$GWAS_POST_TMP/${GWAS}.lead.tsv"
  raw_plot="$GWAS_POST_TMP/${GWAS}.mplot.tsv"
  plot_tsv="${QC_PREFIX}.mplot.small.gz"
  : > "$GWAS_POST_MAGMA_ROWS"; : > "$GWAS_POST_LEAD_VIEW"; : > "$raw_plot"
  gwas_post_log "single full-GWAS scan for requested module views: magma=$want_magma mplot=$want_mplot lead=$want_lead"
  gwas_post_zcat "$FINAL" | awk -v FS='\t' -v OFS='\t' \
    -v want_magma="$want_magma" -v want_mplot="$want_mplot" -v want_lead="$want_lead" \
    -v magma_out="$GWAS_POST_MAGMA_ROWS" -v plot_out="$raw_plot" -v lead_out="$GWAS_POST_LEAD_VIEW" \
    -v hm3_file="$HM3" -v hm3_pos_file="${HM3_POS:-}" -v plot_p='0.001' -v lead_p="$P_LEAD" '
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    BEGIN{
      if(want_mplot=="TRUE"){
        while((getline x<hm3_file)>0){split(x,a,/[ \t]+/);if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1}close(hm3_file)
        if(hm3_pos_file!="")while((getline x<hm3_pos_file)>0){
          gsub(/\r/,"",x);split(x,a,/[ \t]+/);ch=normchr(a[1]);bp=a[4]
          if(ch~/^[0-9]+$/&&bp~/^[0-9]+$/&&bp+0>0)hm3pos[ch SUBSEP (bp+0)]=1
        }
        if(hm3_pos_file!="")close(hm3_pos_file)
      }
    }
    NR==1{
      for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}
      if(want_magma=="TRUE"&&(!("SNP" in c)||!("P" in c))){print "ERROR: MAGMA view requires SNP/P" > "/dev/stderr";exit 2}
      if(want_mplot=="TRUE")print > plot_out
      if(want_lead=="TRUE")print > lead_out
      next
    }
    {
      snp=("SNP" in c?$(c["SNP"]):"");p=("P" in c?$(c["P"]):"");n=("N" in c?$(c["N"]):"NA")
      chr=("CHR" in c?normchr($(c["CHR"])):"");pos=("POS" in c?$(c["POS"]):"")
      if(want_magma=="TRUE"&&snp!=""&&snp!="NA"&&snp!="."&&isnum(p)&&p+0>0&&p+0<=1)print snp,p,n > magma_out
      if(want_mplot=="TRUE"&&((snp in hm3)||((chr SUBSEP (pos+0)) in hm3pos)||(isnum(p)&&p+0<plot_p)))print > plot_out
      if(want_lead=="TRUE"&&isnum(p)&&p+0<=lead_p)print > lead_out
    }'
  if [[ "$want_mplot" == TRUE ]]; then
    gwas_clean_compress < "$raw_plot" > "$plot_tsv"
    GWAS_POST_MPLOT_INPUT="$plot_tsv"
  fi
  GWAS_POST_VIEWS_READY=TRUE
}

gwas_post_magma_annotation() {
  local ref_tag window_tag cache_key cache_dir annot meta lock_file lock_fd tmp_prefix
  local snploc_raw snploc pvar prefix nloc hash meta_tmp
  ref_tag=$(printf '%s' "$REFGEN_CLUMP" | sha256sum | awk '{print substr($1,1,16)}')
  window_tag=$(printf '%s' "$MAGMA_WINDOW" | tr -c 'A-Za-z0-9._-' '_')
  cache_key="v2.GRCh${GRCH}.window_${window_tag}.ref_${ref_tag}"
  cache_dir="$MAGMA_ANNOT_CACHE/v2/GRCh$GRCH/window_$window_tag/ref_$ref_tag"
  annot="$cache_dir/genes.annot"; meta="$cache_dir/annotation.meta.tsv"
  mkdir -p "$(dirname "$cache_dir")"
  lock_file="${cache_dir}.lock"; exec {lock_fd}> "$lock_file"; flock "$lock_fd"
  if [[ ! -s "$annot" || ! -s "$meta" ]]; then
    mkdir -p "$cache_dir"
    snploc_raw="$GWAS_POST_TMP/reference.snp.loc.raw"; snploc="$GWAS_POST_TMP/reference.snp.loc"
    : > "$snploc_raw"
    while read -r prefix _rest; do
      if [[ -s "${prefix}.pvar.zst" ]]; then pvar="${prefix}.pvar.zst"; else pvar="${prefix}.pvar"; fi
      gwas_clean_zcat "$pvar" | awk -v FS='[ \t]+' -v OFS='\t' '
        /^##/{next}
        /^#CHROM/{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}next}
        {id=$(c["ID"]);ch=$(c["CHROM"]);pos=$(c["POS"]);gsub(/^chr/,"",ch);if(id!=""&&id!="."&&ch~/^([1-9]|1[0-9]|2[0-5])$/&&pos~/^[0-9]+$/)print id,ch,pos}' >> "$snploc_raw"
    done < <(ref_clump_pfiles "$REFGEN_CLUMP")
    sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2n -k3,3n "$snploc_raw" | awk '$1!=last{print;last=$1}' > "$snploc"
    nloc=$(wc -l < "$snploc" | tr -d ' '); (( nloc > 1000 )) || { echo "ERROR: too few shared MAGMA reference SNPs: $nloc" >&2; return 1; }
    hash=$(sha256sum "$snploc" | awk '{print $1}')
    tmp_prefix="$GWAS_POST_TMP/annotation.$$"
    gwas_post_log "Build shared MAGMA annotation from GRCh$GRCH reference coordinates: $annot"
    if [[ "$MAGMA_WINDOW" == "0,0" || "$MAGMA_WINDOW" == 0 ]]; then
      magma --annotate --snp-loc "$snploc" --gene-loc "$GENE_LOC" --out "$tmp_prefix"
    else
      magma --annotate window="$MAGMA_WINDOW" --snp-loc "$snploc" --gene-loc "$GENE_LOC" --out "$tmp_prefix"
    fi
    gwas_post_need_file "$tmp_prefix.genes.annot"
    mv -f -- "$tmp_prefix.genes.annot" "${annot}.tmp.$$"; mv -f -- "${annot}.tmp.$$" "$annot"
    meta_tmp="${meta}.tmp.$$"
    printf 'key\tvalue\ngrch\t%s\nwindow_kb\t%s\ngene_loc\t%s\nrefgen_clump\t%s\nsnp_loc_n\t%s\nsnp_loc_sha256\t%s\n' \
      "$GRCH" "$MAGMA_WINDOW" "$GENE_LOC" "$REFGEN_CLUMP" "$nloc" "$hash" > "$meta_tmp"
    mv -f -- "$meta_tmp" "$meta"
  else
    gwas_post_log "Reuse shared MAGMA reference annotation: $annot"
  fi
  MAGMA_ANNOT="$annot"; MAGMA_ANNOT_KEY="$cache_key"
  MAGMA_SNPLOC_N=$(awk -F '\t' '$1=="snp_loc_n"{print $2;exit}' "$meta")
  MAGMA_SNPLOC_HASH=$(awk -F '\t' '$1=="snp_loc_sha256"{print $2;exit}' "$meta")
  flock -u "$lock_fd"; exec {lock_fd}>&-
}

gwas_post_magma() {
  local pval header narg has_usable_n nloc npval meta_tmp
  [[ ",${DO_STEP}," == *,magma,* ]] || return 0
  if [[ "$REPLACE" != TRUE && -s "$MAGMA_DIR/magma.done" && -s "$MAGMA_PREFIX.genes.out" && -s "$MAGMA_PREFIX.genes.raw" ]]; then
    gwas_post_prune_magma_dir; gwas_post_log "MAGMA exists: $MAGMA_PREFIX.genes.out"; return 0
  fi
  command -v magma >/dev/null 2>&1 || { echo "ERROR: magma not found in PATH" >&2; return 1; }
  gwas_post_prepare_views
  mkdir -p "$MAGMA_DIR"; rm -f "$MAGMA_DIR/magma.done" "$MAGMA_DIR/magma.meta.tsv"
  [[ -z "$MAGMA_N" || "$MAGMA_N" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid MAGMA N: $MAGMA_N" >&2; return 1; }
  if [[ -z "$MAGMA_N" && -s "$GWAS_DIR/$GWAS.magma.N" ]]; then MAGMA_N=$(awk 'NF{print $1;exit}' "$GWAS_DIR/$GWAS.magma.N"); fi
  pval="$GWAS_POST_TMP/$GWAS.pval"
  if [[ -n "$MAGMA_N" ]]; then
    { printf 'SNP\tP\n'; sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2g "$GWAS_POST_MAGMA_ROWS" | awk -F '\t' '$1!=last{print $1 "\t" $2;last=$1}'; } > "$pval"
    narg="N=$MAGMA_N"
  else
    has_usable_n=$(awk -F '\t' '$3~/^[0-9]+([.][0-9]+)?$/&&$3+0>=50{n++}END{print n+0}' "$GWAS_POST_MAGMA_ROWS")
    if (( has_usable_n > 1000 )); then
      { printf 'SNP\tP\tN\n'; sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2g "$GWAS_POST_MAGMA_ROWS" | awk -F '\t' '$1!=last&&$3~/^[0-9]+([.][0-9]+)?$/&&$3+0>=50{print;last=$1}'; } > "$pval"
      narg='ncol=N'
    else
      MAGMA_N="$GWAS_N"
      { printf 'SNP\tP\n'; sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2g "$GWAS_POST_MAGMA_ROWS" | awk -F '\t' '$1!=last{print $1 "\t" $2;last=$1}'; } > "$pval"
      narg="N=$MAGMA_N"
    fi
  fi
  npval=$(( $(wc -l < "$pval") - 1 )); (( npval > 1000 )) || { echo "ERROR: too few MAGMA SNPs: pval=$npval" >&2; return 1; }
  gwas_post_magma_annotation
  nloc="$MAGMA_SNPLOC_N"
  magma --bfile "$MAGMA_REF" synonyms="$SYNONYMS" --pval "$pval" "$narg" --gene-annot "$MAGMA_ANNOT" --out "$MAGMA_PREFIX"
  [[ -s "$MAGMA_PREFIX.genes.out" && -s "$MAGMA_PREFIX.genes.raw" ]] || { echo "ERROR: MAGMA output missing: $MAGMA_PREFIX" >&2; return 1; }
  meta_tmp="$MAGMA_DIR/magma.meta.tsv.tmp.$$"
  printf 'key\tvalue\ngwas\t%s\ngrch\t%s\ngene_loc\t%s\nld_reference\t%s\nsynonyms\t%s\nwindow_kb\t%s\nsnp_loc_n\t%s\npval_n\t%s\nannotation_cache\t%s\nannotation_key\t%s\nsnp_loc_sha256\t%s\n' \
    "$GWAS" "$GRCH" "$GENE_LOC" "$MAGMA_REF" "$SYNONYMS" "$MAGMA_WINDOW" "$nloc" "$npval" "$MAGMA_ANNOT" "$MAGMA_ANNOT_KEY" "$MAGMA_SNPLOC_HASH" > "$meta_tmp"
  mv -f -- "$meta_tmp" "$MAGMA_DIR/magma.meta.tsv"
  gwas_post_prune_magma_dir; date '+%F %T' > "$MAGMA_DIR/magma.done"
  gwas_post_log "MAGMA done: $MAGMA_PREFIX.genes.out"
}

gwas_post_mplot_current() {
  local genes="${MAGMA_PREFIX}.genes.out" signal_input=none
  [[ -z "$ADD_SIGNAL" ]] || signal_input="$ADD_SIGNAL"
  [[ -s "$MH_PNG" && -s "$MH_META" && -s "$MH_FLAG" && "$MH_PNG" -nt "$FINAL" ]] || return 1
  awk -F '\t' -v method="$PLOT_METHOD" -v panel="$ADD_PANEL" -v grch="$GRCH" \
    -v signal="$signal_input" -v match_col="$SIGNAL_MATCH_COL" -v match_value="$SIGNAL_MATCH_VALUE" \
    -v locus_pos="$SIGNAL_LOCUS_POS" -v display_col="$SIGNAL_DISPLAY_COL" -v write_sig="$WRITE_SIG" \
    -v plot_width="$PLOT_WIDTH" -v plot_height="$PLOT_HEIGHT" -v plot_res="$PLOT_RES" '
    $1=="plot_method"&&$2==method{m=1}
    $1=="add_panel"&&$2==panel{p=1}
    $1=="grch"&&$2==grch{g=1}
    $1=="magma_threshold"&&$2+0==2.5e-6{t=1}
    $1=="mplot_style"&&$2==8{s=1}
    $1=="plot_width"&&$2+0==plot_width+0{x=1}
    $1=="plot_height"&&$2==plot_height{y=1}
    $1=="plot_res"&&$2+0==plot_res+0{r=1}
    $1=="write_sig"&&$2==write_sig{w=1}
    $1=="signal_input"&&$2==signal{a=1}
    $1=="signal_match_col"&&$2==match_col{b=1}
    $1=="signal_match_value"&&$2==match_value{c=1}
    $1=="signal_locus_pos"&&$2==locus_pos{d=1}
    $1=="signal_display_col"&&$2==display_col{e=1}
    END{exit !(m&&p&&g&&t&&s&&x&&y&&r&&w&&a&&b&&c&&d&&e)}' "$MH_META" || return 1
  if [[ "$ADD_PANEL" == magma ]]; then
    [[ -s "$genes" && "$MH_PNG" -nt "$genes" ]] || return 1
  fi
  if [[ -n "$ADD_SIGNAL" ]]; then
    [[ -s "$ADD_SIGNAL" && "$MH_PNG" -nt "$ADD_SIGNAL" ]] || return 1
  fi
  if [[ "$WRITE_SIG" == TRUE ]]; then
    [[ -s "$MH_SIG" && "$MH_SIG" -nt "$FINAL" ]] || return 1
    [[ ! -s "$COJO_FILE" || "$MH_SIG" -nt "$COJO_FILE" ]] || return 1
    [[ "$ADD_PANEL" != magma || "$MH_SIG" -nt "$genes" ]] || return 1
  fi
}

gwas_post_append_mplot_flag() {
  local lock_file lock_fd fragment_header aggregate_header
  [[ -s "$MH_FLAG" ]] || { echo "ERROR: mplot flag fragment is missing: $MH_FLAG" >&2; return 1; }
  mkdir -p "$(dirname "$MH_FLAG")" "$(dirname "$MPLOT_FLAG_FILE")"
  command -v flock >/dev/null 2>&1 || { echo "ERROR: flock is required to append $MPLOT_FLAG_FILE" >&2; return 1; }
  lock_file="$(dirname "$MPLOT_FLAG_FILE")/.0flag.lock"
  exec {lock_fd}> "$lock_file"
  flock "$lock_fd"

  fragment_header=$(awk 'NR==1 {sub(/\r$/, ""); print; exit}' "$MH_FLAG")
  if [[ "$fragment_header" != $'GWAS\tFLAG' ]]; then
    echo "ERROR: invalid mplot flag fragment header: $MH_FLAG" >&2
    flock -u "$lock_fd"
    exec {lock_fd}>&-
    return 1
  fi

  if [[ -s "$MPLOT_FLAG_FILE" ]]; then
    aggregate_header=$(awk 'NR==1 {sub(/\r$/, ""); print; exit}' "$MPLOT_FLAG_FILE")
    if [[ "$aggregate_header" != $'GWAS\tFLAG' ]]; then
      echo "ERROR: invalid mplot flag header; refusing to overwrite or append: $MPLOT_FLAG_FILE" >&2
      flock -u "$lock_fd"
      exec {lock_fd}>&-
      return 1
    fi
  else
    printf 'GWAS\tFLAG\n' >> "$MPLOT_FLAG_FILE"
  fi

  # Append only rows whose GWAS is not already present.  The aggregate is never
  # rebuilt, truncated, renamed over, or rewritten, so completed reruns preserve it.
  if ! awk -F '\t' '
    FNR==NR {if (FNR>1 && $1!="") seen[$1]=1; next}
    FNR>1 && $1!="" && $2!="" && !seen[$1]++ {print $1 "\t" $2}
  ' "$MPLOT_FLAG_FILE" "$MH_FLAG" >> "$MPLOT_FLAG_FILE"; then
    echo "ERROR: failed to append mplot flag: $MPLOT_FLAG_FILE" >&2
    flock -u "$lock_fd"
    exec {lock_fd}>&-
    return 1
  fi
  flock -u "$lock_fd"
  exec {lock_fd}>&-
}

gwas_post_mplot() {
  local genes="${MAGMA_PREFIX}.genes.out" meta_tmp panel_input signal_input sig_input cojo_input
  [[ "$DO_STEP" == all || ",${DO_STEP}," == *,mplot,* ]] || return 0
  if [[ "$REPLACE" != TRUE ]] && gwas_post_mplot_current; then
    gwas_post_append_mplot_flag
    gwas_post_log "Manhattan plot exists and is current: $MH_PNG"
    return 0
  fi
  if [[ "$ADD_PANEL" == magma && ! -s "$genes" ]]; then
    echo "ERROR: --add-panel magma requires MAGMA output: $genes" >&2
    echo "ERROR: run the magma,mplot modules together, or create MAGMA output first" >&2
    return 1
  fi
  gwas_post_prepare_views
  [[ -s "$GWAS_POST_MPLOT_INPUT" ]] || { echo "ERROR: Manhattan subset missing" >&2; return 1; }
  mkdir -p "$(dirname "$MH_PNG")" "$(dirname "$MH_META")"
  rm -f -- "$MH_META" "${MH_PNG}.meta.tsv"
  [[ "$WRITE_SIG" != TRUE ]] || rm -f -- "$MH_SIG"
  gwas_post_log "Manhattan plot method=$PLOT_METHOD panel=$ADD_PANEL write-sig=$WRITE_SIG: $GWAS_POST_MPLOT_INPUT -> $MH_PNG"
  Rscript "$MPLOT_R" "$GWAS_POST_MPLOT_INPUT" "$MH_PNG" "$GWAS" "$PLOT_METHOD" "$ADD_PANEL" \
    "$genes" "$PLOT_F" "$CIS_BED" "$MH_PLOT_BED" "$GRCH" "$MH_FLAG" \
    "$ADD_SIGNAL" "$SIGNAL_MATCH_COL" "$SIGNAL_MATCH_VALUE" "$SIGNAL_LOCUS_POS" "$SIGNAL_DISPLAY_COL" \
    "$WRITE_SIG" "$MH_SIG" "$COJO_FILE" "$PLOT_WIDTH" "$PLOT_HEIGHT" "$PLOT_RES"
  [[ -s "$MH_PNG" && -s "$MH_FLAG" ]] || { echo "ERROR: Manhattan plot or flag fragment was not created: $MH_PNG $MH_FLAG" >&2; return 1; }
  [[ "$WRITE_SIG" != TRUE || -s "$MH_SIG" ]] || { echo "ERROR: significant-hit summary was not created: $MH_SIG" >&2; return 1; }
  gwas_post_append_mplot_flag
  panel_input=none
  [[ "$ADD_PANEL" != magma ]] || panel_input="$genes"
  meta_tmp="${MH_META}.tmp.$$"
  signal_input=none
  [[ -z "$ADD_SIGNAL" ]] || signal_input="$ADD_SIGNAL"
  sig_input=none
  cojo_input=none
  if [[ "$WRITE_SIG" == TRUE ]]; then
    sig_input="$MH_SIG"
    [[ ! -s "$COJO_FILE" ]] || cojo_input="$COJO_FILE"
  fi
  printf 'key\tvalue\nplot_method\t%s\nadd_panel\t%s\ngrch\t%s\nmagma_threshold\t2.5e-6\nmplot_style\t8\nplot_width\t%s\nplot_height\t%s\nplot_res\t%s\nwrite_sig\t%s\nsig_output\t%s\ncojo_input\t%s\nsignal_input\t%s\nsignal_match_col\t%s\nsignal_match_value\t%s\nsignal_locus_pos\t%s\nsignal_display_col\t%s\ngwas\t%s\ngwas_input\t%s\nmagma_input\t%s\nflag_fragment\t%s\ncreated\t%s\n' \
    "$PLOT_METHOD" "$ADD_PANEL" "$GRCH" "$PLOT_WIDTH" "$PLOT_HEIGHT" "$PLOT_RES" "$WRITE_SIG" "$sig_input" "$cojo_input" \
    "$signal_input" "$SIGNAL_MATCH_COL" "$SIGNAL_MATCH_VALUE" "$SIGNAL_LOCUS_POS" \
    "$SIGNAL_DISPLAY_COL" "$GWAS" "$FINAL" "$panel_input" "$MH_FLAG" "$(date '+%F %T')" > "$meta_tmp"
  mv -f -- "$meta_tmp" "$MH_META"
}

# Lead extraction now reads the tiny P-threshold view made by the shared scan.
gwas_post_prep_lead_inputs() {
  local refs="$1" suffix=".tmp.$$" assoc ma source="${GWAS_POST_LEAD_VIEW:-$FINAL}"
  rm -f -- "${QC_PREFIX}".*.cojo_skip.log
  while IFS=$'\t' read -r _ _ _ assoc ma _; do rm -f "${assoc}${suffix}" "${ma}${suffix}"; done < "$refs"
  gwas_post_zcat "$source" | awk -v FS='\t' -v OFS='\t' -v refs="$refs" -v suffix="$suffix" -v default_n="$GWAS_N" '
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;return x}
    function get(x){return x in c?$(c[x]):"NA"}
    BEGIN{while((getline line<refs)>0){split(line,a,"\t");chr=a[3];if(chr=="")continue;assoc[chr]=a[4] suffix;ma[chr]=a[5] suffix;print "SNP","CHR","POS","EA","NEA","P">assoc[chr];print "SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","N">ma[chr]}close(refs)}
    NR==1{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}next}
    {s=get("SNP");ch=normchr(get("CHR"));pos=get("POS");ea=get("EA");nea=get("NEA");p=get("P");if(!(ch in assoc)||pos!~/^[0-9]+$/||s=="NA"||!isnum(p))next;print s,ch,pos,ea,nea,p>assoc[ch];b=get("BETA");se=get("SE");if(!isnum(b)||!isnum(se))next;n=get("N");if(!isnum(n))n=default_n;print s,ch,pos,ea,nea,get("EAF"),b,se,p,n>ma[ch]}
  '
  while IFS=$'\t' read -r _ _ _ assoc ma _; do mv -f "${assoc}${suffix}" "$assoc";mv -f "${ma}${suffix}" "$ma";done < "$refs"
}

gwas_post_reference_cache() {
  local source="$1" kind="$2" stat_key key dir cached lock fd tmp
  stat_key=$(stat -c '%n|%s|%Y' "$source")
  key=$(printf '%s|%s' "$kind" "$stat_key" | sha256sum | awk '{print substr($1,1,20)}')
  dir="$LEAD_REF_CACHE/$kind"; cached="$dir/$key.bgz"; lock="$cached.lock"
  mkdir -p "$dir"; exec {fd}> "$lock"; flock "$fd"
  if ! gwas_index_valid "$cached"; then
    tmp="$dir/$key.tmp.$$.bgz"; rm -f -- "$tmp" "$tmp.tbi" "$tmp.csi"
    if [[ "$kind" == pvar ]]; then
      gwas_clean_zcat "$source" | awk -v FS='[ \t]+' -v OFS='\t' '
        /^##/{next} /^#CHROM/{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}print "#CHROM","POS","ID","REF","ALT";next}
        {ch=$(c["CHROM"]);gsub(/^chr/,"",ch);pos=$(c["POS"]);if(ch!=""&&pos~/^[0-9]+$/)print ch,pos,$(c["ID"]),$(c["REF"]),$(c["ALT"])}' |
        { IFS= read -r h; printf '%s\n' "$h"; sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k1,1V -k2,2n; } | bgzip -c > "$tmp"
      tabix -f -s 1 -b 2 -e 2 -S 1 "$tmp"
    else
      gwas_clean_zcat "$source" | awk -v FS='[ \t]+' -v OFS='\t' 'NF>=6{ch=$1;gsub(/^chr/,"",ch);if($4~/^[0-9]+$/)print ch,$2,$3,$4,$5,$6}' |
        sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k1,1V -k4,4n | bgzip -c > "$tmp"
      tabix -f -s 1 -b 4 -e 4 "$tmp"
    fi
    mv -f -- "$tmp" "$cached"; [[ ! -s "$tmp.tbi" ]] || mv -f -- "$tmp.tbi" "$cached.tbi"; [[ ! -s "$tmp.csi" ]] || mv -f -- "$tmp.csi" "$cached.csi"
  fi
  flock -u "$fd"; exec {fd}>&-
  printf '%s\n' "$cached"
}

gwas_post_subset_match_reference() {
  local ref="$1" query="$2" out="$3" kind="$4" source="" ext cached regions
  if [[ "$kind" == pvar ]]; then
    for ext in .pvar .pvar.gz .pvar.bgz .pvar.zst; do [[ -s "${ref}${ext}" ]] && { source="${ref}${ext}"; break; }; done
  else
    for ext in .bim .bim.gz .bim.bgz; do [[ -s "${ref}${ext}" ]] && { source="${ref}${ext}"; break; }; done
  fi
  gwas_post_need_file "$source"
  cached=$(gwas_post_reference_cache "$source" "$kind")
  regions="$GWAS_POST_TMP/reference.$kind.regions.txt"
  awk -v FS='\t' 'NR==1{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);if(h=="CHR")c=i;else if(h=="POS")p=i}next} c&&p&&$p~/^[0-9]+$/{ch=$c;gsub(/^chr/,"",ch);print ch ":" $p "-" $p}' "$query" | sort -u > "$regions"
  if [[ "$kind" == pvar ]]; then printf '#CHROM\tPOS\tID\tREF\tALT\n' > "$out"; else : > "$out"; fi
  if [[ -s "$regions" ]]; then mapfile -t region_args < "$regions"; tabix "$cached" "${region_args[@]}" >> "$out"; fi
  if [[ "$kind" == bim && ! -s "$out" ]]; then printf '0\t.\t0\t0\tN\tN\n' > "$out"; fi
}
