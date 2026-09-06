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

# Skip auxiliary combined panels such as chrXY while enumerating chromosome
# pfiles; lead/MAGMA use the ordinary chr1..chr22/chrX/chrY panels.
ref_clump_pfiles() {
  local ref="$1" f prefix label chr found=FALSE
  if [[ -f "${ref}.pgen" && -f "${ref}.psam" && ( -f "${ref}.pvar" || -f "${ref}.pvar.zst" ) ]]; then
    label=$(basename "$ref"); chr=$(ref_chr "$label") || return 1
    printf '%s %s %s %s\n' "$ref" "$label" "$chr" "${ref}.bad.snp"; return 0
  fi
  if [[ -d "$ref" ]]; then
    while IFS= read -r f; do
      prefix=${f%.pgen}
      [[ -f "${prefix}.psam" && ( -f "${prefix}.pvar" || -f "${prefix}.pvar.zst" ) ]] || continue
      label=$(basename "$prefix"); chr=$(ref_chr "$label") || continue
      printf '%s %s %s %s\n' "$prefix" "$label" "$chr" "${prefix}.bad.snp"; found=TRUE
    done < <(find "$ref" -maxdepth 1 -type f -name 'chr*.pgen' | sort -V)
  else
    while IFS= read -r f; do
      prefix=${f%.pgen}
      [[ -f "${prefix}.psam" && ( -f "${prefix}.pvar" || -f "${prefix}.pvar.zst" ) ]] || continue
      label=$(basename "$prefix"); chr=$(ref_chr "$label") || continue
      printf '%s %s %s %s\n' "$prefix" "$label" "$chr" "${prefix}.bad.snp"; found=TRUE
    done < <(compgen -G "${ref}chr*.pgen" | sort -V)
  fi
  [[ "$found" == TRUE ]] || { echo "ERROR: no chromosome pfile reference found: $ref" >&2; return 1; }
}

gwas_post_bad_snp_cache() {
  local cache="$REF_BAD_DIR/GRCh${GRCH}.bad.snp" lock="${REF_BAD_DIR}/GRCh${GRCH}.lock" fd tmp
  mkdir -p "$REF_BAD_DIR"
  command -v flock >/dev/null 2>&1 || { printf '%s\n' "$cache"; return 0; }
  exec {fd}> "$lock"
  flock "$fd"
  if [[ ! -s "$cache" ]]; then
    tmp="${cache}.tmp.$$"
    {
      while read -r _ _ _ bad; do [[ -s "$bad" ]] && cat -- "$bad"; done < <(ref_clump_pfiles "$REFGEN_CLUMP")
      printf '# cache may be empty\n'
    } | awk 'NF && !seen[$1]++' | sort -u > "$tmp"
    mv -f -- "$tmp" "$cache"
  fi
  flock -u "$fd"; exec {fd}>&-
  printf '%s\n' "$cache"
}

# --small TRUE: filter and standardize in the same raw scan. The old two-pass
# candidate/format path is preserved semantically, including exact-row dedup.
std_small() {
  local src="$1" out="$2" header_out="$3" input_fs bad_snp tmp_tsv tmp audit n_out
  gwas_clean_need_file "$src"; gwas_clean_need_file "$HM3"
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
  bad_snp=$(gwas_post_bad_snp_cache)
  input_fs=$(gwas_clean_detect_fs "$src")
  tmp_tsv="$GWAS_POST_TMP/${GWAS}.small.sorted.tsv"
  tmp="${out}.tmp.$$"
  audit="${QC_PREFIX}.small.audit.tsv"
  gwas_clean_zcat "$src" | awk -v FS="$input_fs" -v OFS='\t' \
    -v hm3_file="$HM3" -v bad_file="$bad_snp" -v pthr="$P_SMALL" \
    -v cis_bed="${CIS_BED:-}" -v cis_name="$GWAS" -v flank="${CIS_FLANK:-0}" \
    -v snp_col="$SNP_col" -v chr_col="$CHR_col" -v pos_col="$POS_col" \
    -v ea_col="$EA_col" -v nea_col="$NEA_col" -v eaf_col="$EAF_col" -v n_col="$N_col" \
    -v beta_col="$BETA_col" -v se_col="$SE_col" -v p_col="$P_col" -v logp_col="$LOG10P_col" \
    -v fill_n="$FILL_N" -v audit="$audit" -v gwas="$GWAS" '
    function get(c,x){x=(c>0?$c:"");gsub(/\r/,"",x);return x}
    function val(c,x){x=get(c);return x==""?"NA":x}
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function validn(x){return x~/^[0-9]+([.][0-9]+)?$/&&x+0>0}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    function incis(chr,pos,i){for(i=1;i<=nr;i++)if(chr==rc[i]&&pos>=rs[i]&&pos<=re[i])return 1;return 0}
    BEGIN{
      logpthr=(pthr+0>0?-log(pthr+0)/log(10):0)
      while((getline x<hm3_file)>0){gsub(/\r/,"",x);split(x,a,/[ \t]+/);if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1}close(hm3_file)
      while((getline x<bad_file)>0){gsub(/\r/,"",x);split(x,a,/[ \t]+/);if(a[1]!=""&&a[1]!~/^#/)bad[a[1]]=1}close(bad_file)
      if(cis_bed!="")while((getline x<cis_bed)>0){gsub(/\r/,"",x);if(x==""||x~/^#/)continue;split(x,a,/[ \t]+/);if(a[4]==cis_name){nr++;rc[nr]=normchr(a[1]);rs[nr]=a[2]-flank;if(rs[nr]<1)rs[nr]=1;re[nr]=a[3]+flank}}close(cis_bed)
    }
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P";next}
    {
      chr=normchr(get(chr_col));pos=get(pos_col);if(chr!~/^[0-9]+$/||chr+0<1||chr+0>25||pos!~/^[0-9]+$/||pos+0<1){badcoord++;next}
      snp=get(snp_col);ea=get(ea_col);nea=get(nea_col);if(snp==""||snp=="NA"||snp=="."){snp=chr ":" pos;if(ea!="")snp=snp ":" ea;if(nea!="")snp=snp ":" nea}
      if(snp in bad){badid++;next}
      p=(p_col>0?get(p_col):"");lp=(logp_col>0?get(logp_col):"");if(p==""&&isnum(lp))p=10^(-lp);if(lp==""&&isnum(p)&&p>0)lp=-log(p)/log(10)
      keep=(snp in hm3)||(isnum(p)&&p+0<=pthr+0)||(nr>0&&incis(chr,pos+0));if(!keep)next
      if(ea=="")ea="NA";if(nea=="")nea="NA";if(p=="")p="NA";if(lp=="")lp="NA"
      n=val(n_col);if(fill_n!=""&&!validn(n)){n=fill_n;nfilled++}else if(validn(n))nkept++
      print snp,chr+0,pos+0,ea,nea,val(eaf_col),n,val(beta_col),val(se_col),p,lp;kept++
    }
    END{print "GWAS\tN_PREFILTER\tN_DROPPED_BAD_COORD\tN_DROPPED_BAD_SNP\tN_EXISTING\tN_FILLED">audit;print gwas "\t" kept+0 "\t" badcoord+0 "\t" badid+0 "\t" nkept+0 "\t" nfilled+0>>audit}
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
  GWAS_POST_MAGMA_ROWS="$GWAS_POST_TMP/${GWAS}.magma.rows.tsv"
  GWAS_POST_LEAD_VIEW="$GWAS_POST_TMP/${GWAS}.lead.tsv"
  raw_plot="$GWAS_POST_TMP/${GWAS}.mplot.tsv"
  plot_tsv="${QC_PREFIX}.mplot.small.gz"
  : > "$GWAS_POST_MAGMA_ROWS"; : > "$GWAS_POST_LEAD_VIEW"; : > "$raw_plot"
  gwas_post_log "single full-GWAS scan for requested module views: magma=$want_magma mplot=$want_mplot lead=$want_lead"
  gwas_post_zcat "$FINAL" | awk -v FS='\t' -v OFS='\t' \
    -v want_magma="$want_magma" -v want_mplot="$want_mplot" -v want_lead="$want_lead" \
    -v magma_out="$GWAS_POST_MAGMA_ROWS" -v plot_out="$raw_plot" -v lead_out="$GWAS_POST_LEAD_VIEW" \
    -v hm3_file="$HM3" -v plot_p='0.001' -v lead_p="$P_LEAD" '
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    BEGIN{if(want_mplot=="TRUE")while((getline x<hm3_file)>0){split(x,a,/[ \t]+/);if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1}close(hm3_file)}
    NR==1{
      for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}
      if(want_magma=="TRUE"&&(!("SNP" in c)||!("P" in c))){print "ERROR: MAGMA view requires SNP/P" > "/dev/stderr";exit 2}
      if(want_mplot=="TRUE")print > plot_out
      if(want_lead=="TRUE")print > lead_out
      next
    }
    {
      snp=("SNP" in c?$(c["SNP"]):"");p=("P" in c?$(c["P"]):"");n=("N" in c?$(c["N"]):"NA")
      if(want_magma=="TRUE"&&snp!=""&&snp!="NA"&&snp!="."&&isnum(p)&&p+0>0&&p+0<=1)print snp,p,n > magma_out
      if(want_mplot=="TRUE"&&((snp in hm3)||(isnum(p)&&p+0<plot_p)))print > plot_out
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
  cache_key="GRCh${GRCH}.window_${window_tag}.ref_${ref_tag}"
  cache_dir="$MAGMA_ANNOT_CACHE/GRCh$GRCH/window_$window_tag/ref_$ref_tag"
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
  local meta="${MH_PNG}.meta.tsv" genes="${MAGMA_PREFIX}.genes.out"
  [[ -s "$MH_PNG" && -s "$meta" && "$MH_PNG" -nt "$FINAL" ]] || return 1
  awk -F '\t' -v method="$PLOT_METHOD" -v panel="$ADD_PANEL" '
    $1=="plot_method"&&$2==method{m=1}
    $1=="add_panel"&&$2==panel{p=1}
    END{exit !(m&&p)}' "$meta" || return 1
  if [[ "$ADD_PANEL" == magma ]]; then
    [[ -s "$genes" && "$MH_PNG" -nt "$genes" ]] || return 1
  fi
}

gwas_post_mplot() {
  local genes="${MAGMA_PREFIX}.genes.out" meta="${MH_PNG}.meta.tsv" meta_tmp panel_input
  [[ "$DO_STEP" == all || ",${DO_STEP}," == *,mplot,* ]] || return 0
  if [[ "$REPLACE" != TRUE ]] && gwas_post_mplot_current; then
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
  mkdir -p "$(dirname "$MH_PNG")"
  rm -f -- "$MH_PNG" "$meta"
  gwas_post_log "Manhattan plot method=$PLOT_METHOD panel=$ADD_PANEL: $GWAS_POST_MPLOT_INPUT -> $MH_PNG"
  Rscript "$MPLOT_R" "$GWAS_POST_MPLOT_INPUT" "$MH_PNG" "$GWAS" "$PLOT_METHOD" "$ADD_PANEL" \
    "$genes" "$PLOT_F" "$CIS_BED" "$MH_PLOT_BED"
  [[ -s "$MH_PNG" ]] || { echo "ERROR: Manhattan plot was not created: $MH_PNG" >&2; return 1; }
  panel_input=none
  [[ "$ADD_PANEL" != magma ]] || panel_input="$genes"
  meta_tmp="${meta}.tmp.$$"
  printf 'key\tvalue\nplot_method\t%s\nadd_panel\t%s\ngwas\t%s\ngwas_input\t%s\nmagma_input\t%s\ncreated\t%s\n' \
    "$PLOT_METHOD" "$ADD_PANEL" "$GWAS" "$FINAL" "$panel_input" "$(date '+%F %T')" > "$meta_tmp"
  mv -f -- "$meta_tmp" "$meta"
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
