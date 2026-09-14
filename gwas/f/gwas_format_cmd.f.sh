#!/usr/bin/env bash
# Sourced by gwas_format.sh: build one self-contained worker command per GWAS.
# Uses entry-point configuration and the helpers in gwas_format.f.sh.
# Keep the worker heredoc literal/escaped as written: expansion happens at plan time.

write_gwas_cmd() {
  local gwas="$1" raw trait_dir gwas_dir source_gwas final cmd awk_snp cis_out qc_prefix clump_dir cojo_dir merged_prefix clump_done cojo_done mh_png mh_meta mh_flag mh_sig mplot_flag_file magma_dir magma_prefix clump_kb
  local pgs_dir pgs_score_file pgs_output pgs_done pgs_meta gwas_pgs_pfile_dir
  local gwas_grch gwas_grch_cache gwas_refGen_clump gwas_refGen_cojo gwas_refGen_id_dir gwas_refGen_keep gwas_gene_loc gwas_mh_plot_bed gwas_hm3_pos detection_input cached_grch rc
  raw=""
  # Only formatting consumes the raw GWAS.  Downstream-only steps use files in
  # the trait's gwas directory, so resolving raw here would rescan the entire
  # project once per trait when no matching raw file exists.
  if has_step format; then
    if [[ -n "$raw_file_arg" ]]; then raw="$raw_file_arg"; else raw=$(raw_file_for_name "$gwas" || true); fi
  fi
  if [[ -n "$dir_clean_arg" ]]; then
    gwas_dir="$dir_clean"
    trait_dir=$(dirname "$gwas_dir")
  else
    trait_dir="$dir_out/$category/$gwas"
    gwas_dir="$trait_dir/gwas"
  fi
  final="$gwas_dir/$gwas.gz"
  if [[ "$liftOver" == "TRUE" ]]; then
    if [[ "$hm3_mode" == FALSE ]]; then source_gwas="$trait_dir/qc/$gwas.source.grch37.gz"; else source_gwas="$gwas_dir/$gwas.hm3.gz"; fi
  else
    source_gwas="$final"
  fi
  # Existing liftOver inputs may still have the old suffix.
  if [[ "$liftOver" == TRUE && "$hm3_mode" == TRUE && ! -s "$source_gwas" ]] &&
     ! has_step format && [[ -s "$gwas_dir/$gwas.small.gz" ]]; then
    source_gwas="$gwas_dir/$gwas.small.gz"
  fi
  awk_snp="$gwas_dir/$gwas.awk.snp"
  cis_out="$gwas_dir/$gwas.cis.gz"
  qc_prefix="$trait_dir/qc/$gwas"
  gwas_grch_cache="${qc_prefix}.grch"
  clump_dir="$gwas_dir/clump"
  cojo_dir="$gwas_dir/cojo"
  merged_prefix="$gwas_dir/$gwas"
  clump_done="$gwas_dir/$gwas.clump.done"
  cojo_done="$gwas_dir/$gwas.cojo.done"
  mh_png="$dir_out/mplot/$gwas.png"
  mh_meta="$dir_out/.project/$category/mplot/$gwas.state"
  mh_flag="$dir_out/.project/$category/mplot/flag/$gwas.tsv"
  mh_sig="$dir_out/mplot/$gwas.sig.txt"
  mplot_flag_file="$dir_out/mplot/0flag.tsv"
  pgs_dir="$trait_dir/pgs"
  pgs_score_file="${merged_prefix}.jma.cojo"
  pgs_output="$pgs_dir/$gwas.pgs.gz"
  pgs_done="$pgs_dir/$gwas.pgs.done"
  pgs_meta="$pgs_dir/$gwas.pgs.meta.tsv"
  if [[ -n "$dir_magma" ]]; then magma_dir="$dir_magma"; else magma_dir="$trait_dir/magma"; fi
  magma_prefix="$magma_dir/$gwas"
  clump_kb=$(( (lead_window + 999) / 1000 ))
  cmd="$dir_cmd/$gwas.cmd"
  mkdir -p "$gwas_dir" "$trait_dir/qc" "$(dirname "$mh_meta")"

  # Keep the standalone PGS fast path intentionally minimal: the two non-empty
  # result/marker files alone define completion.  Skip build detection, pfile
  # validation, gzip tests, metadata reads, and timestamp comparisons.
  if [[ "$replace" != "TRUE" && "$step" == pgs ]] &&
     pgs_output_complete "$pgs_output" "$pgs_done"; then
    prune_pgs_dir "$pgs_dir" "$gwas"
    log "SKIP completed PGS: $gwas"
    rm -f "$cmd"
    return 0
  fi

  # A completed lead result does not need GWAS-build detection or reference
  # validation.  Require independent clump/COJO markers plus their artifacts,
  # except for the explicit no-significant-row case.
  if [[ "$replace" != "TRUE" && "$step" == lead ]] && gzip_ok "$final" &&
     gwas_lead_marker_matches "$clump_done" "$gwas" clump "$p_lead" "$lead_window" "$chrs" &&
     gwas_lead_marker_matches "$cojo_done" "$gwas" cojo "$p_lead" "$lead_window" "$chrs" &&
     [[ -s "$awk_snp" && ! -d "$clump_dir" && ! -d "$cojo_dir" ]] &&
     { { [[ -s "$clump_done" && -s "${merged_prefix}.clumps" && -s "$cojo_done" ]] &&
           [[ -s "${merged_prefix}.jma.cojo" || -s "${merged_prefix}.ldr.cojo" ]]; } ||
       { [[ -s "$clump_done" && -s "$cojo_done" ]] && awk 'NR>1{exit 1}' "$awk_snp"; }; }; then
    log "SKIP completed lead GWAS: $gwas"
    rm -f "$cmd"
    return 0
  fi

  gwas_grch="$grch"
  # With no fixed build (omitted or --grch auto), resolve every GWAS independently
  # from the 39 sentinel rsIDs.  Use RAW before format, SMALL before a standalone
  # liftOver, and the standardized FINAL for downstream-only requests.
  if [[ "$gwas_grch" == auto ]]; then
    if has_step format; then
      detection_input="$raw"
      [[ -s "$detection_input" ]] || { echo "ERROR: missing raw GWAS for GRCh detection: $gwas" >&2; exit 1; }
    elif [[ "$step" == liftover ]]; then
      detection_input="$source_gwas"
      [[ -s "$detection_input" ]] || { echo "ERROR: missing source GWAS for GRCh detection: $gwas" >&2; exit 1; }
    else
      detection_input="$final"
      [[ -s "$detection_input" ]] || { echo "ERROR: missing standardized GWAS for GRCh detection: $gwas" >&2; exit 1; }
    fi
    cached_grch=""
    if [[ "$liftOver" != TRUE && -s "$gwas_grch_cache" && "$gwas_grch_cache" -nt "$detection_input" ]]; then
      cached_grch=$(awk 'NR==1 && ($1==37 || $1==38){print $1}' "$gwas_grch_cache")
    fi
    if [[ -n "$cached_grch" ]]; then
      gwas_grch="$cached_grch"
      log "Reuse cached GRCh$gwas_grch: $gwas"
    elif check_GRCH "$detection_input" >&2; then
      gwas_grch="$CHECK_GRCH_RESULT"
      printf '%s\n' "$gwas_grch" > "${gwas_grch_cache}.tmp.$$"
      mv -f "${gwas_grch_cache}.tmp.$$" "$gwas_grch_cache"
    else
      rc=$?
      echo "ERROR: automatic GRCh detection failed for $gwas; the input must contain usable rsIDs from ${CHECK_GRCH_SNP_LIST:-$dir0/data/ukb/phe/common/snp.lst}, or specify --grch 37/38." >&2
      return "$rc"
    fi
  fi
  gwas_hm3_pos=${hm3_pos//\{grch\}/$gwas_grch}
  gwas_pgs_pfile_dir="${pgs_pfile_dir:-/mnt/i/ukbGen/${gwas_grch}/imp}"
  if has_step pgs; then
    need_pgs_pfiles "$gwas_pgs_pfile_dir"
  fi
  gwas_gene_loc="$gene_loc"
  if [[ -z "$gwas_gene_loc" ]]; then
    [[ "$gwas_grch" == 37 ]] && gwas_gene_loc="$dir0/files/NCBI.37.gene.loc" || gwas_gene_loc="$dir0/files/NCBI.38.gene.loc"
  fi
  gwas_mh_plot_bed=""
  gwas_mh_plot_bed="$mh_plot_bed"
  [[ -n "$gwas_mh_plot_bed" ]] || gwas_mh_plot_bed="$dir0/files/glist.${gwas_grch}.bed"
  if has_step format && [[ "$liftOver" == TRUE && "$gwas_grch" != 37 ]]; then
    echo "ERROR: --liftover TRUE uses the GRCh37-to-GRCh38 chain, but $gwas was detected/configured as GRCh$gwas_grch" >&2
    exit 1
  fi
  gwas_refGen_clump="${refGen_clump:-/mnt/i/refGen/1kg/${gwas_grch}/pfile/}"
  gwas_refGen_cojo="${refGen_cojo:-/mnt/i/refGen/1kg/${gwas_grch}/pfile/${refGen_pop}/}"
  gwas_refGen_id_dir="${refGen_id_dir:-/mnt/i/refGen/1kg/${gwas_grch}/id}"
  gwas_refGen_keep=""
  if has_step lead || wants_magma || { has_step format && [[ "$fill_eaf" == TRUE ]]; }; then
    if [[ "$refGen_pop" != ALL ]]; then
      gwas_refGen_keep="${gwas_refGen_id_dir%/}/${refGen_pop}.id.2col"
      if ! awk 'BEGIN{FS="[ \t]+"}{a=$1;b=$2;gsub(/\r/,"",a);gsub(/\r/,"",b);if(NF!=2||a!=b)bad=1;n++}END{exit bad||n==0?2:0}' "$gwas_refGen_keep"; then
        echo "ERROR: invalid PLINK two-column keep file: $gwas_refGen_keep" >&2
        exit 1
      fi
    fi
    need_refgen_clump "$gwas_refGen_clump"
  fi
  if has_step lead; then
    need_refgen_cojo "$gwas_refGen_cojo"
  fi

  if [[ "$replace" != "TRUE" ]]; then
    case "$step" in
      mplot|thin,mplot|mplot,thin)
        if mplot_output_complete "$mh_png" "$final" "$magma_prefix.genes.out" "$mh_meta" "$gwas_grch" "$mh_flag" "$mplot_flag_file" "$mh_sig" "${merged_prefix}.jma.cojo"; then
          log "SKIP completed Manhattan plot: $gwas panel=$add_panel"
          rm -f "$cmd"
          return 0
        fi
        ;;
      magma)
        if magma_output_complete "$magma_dir" "$magma_prefix"; then
          prune_magma_dir "$magma_dir"
          log "SKIP completed MAGMA: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      cis)
        if cis_output_complete "$cis_out"; then
          log "SKIP completed indexed cis GWAS: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;

    esac
  fi

  cat > "$cmd" <<CMD_TOP
#!/bin/bash
set -euo pipefail
export LC_ALL=C
export PATH=$(q "$dir0/software/bin"):\$PATH

gwas_post_need_file(){
  [[ -s "\$1" ]] || { echo "ERROR: missing or empty file: \$1" >&2; exit 1; }
}
gwas_post_log(){
  printf '[%s] %s\n' "\$(date '+%F %T')" "\$*"
}
gwas_post_zcat(){
  case "\$1" in *.gz|*.bgz) gzip -cd -- "\$1";; *) cat -- "\$1";; esac
}
gwas_post_has_data_rows(){
  [[ -s "\$1" ]] && awk 'NR>1{found=1; exit} END{exit found ? 0 : 1}' "\$1"
}

GWAS=$(q "$gwas")
GWAS_DIR=$(q "$gwas_dir")
GWAS_POST_TMP_ROOT="\$GWAS_DIR/.tmp"
mkdir -p "\$GWAS_POST_TMP_ROOT"
# Isolate temporary files per process and remove them on exit.
for stale_tmp in "\$GWAS_POST_TMP_ROOT"/run.*; do
  [[ -d "\$stale_tmp" ]] || continue
  stale_pid=\${stale_tmp##*.}
  if [[ ! "\$stale_pid" =~ ^[1-9][0-9]*\$ ]] || ! kill -0 "\$stale_pid" 2>/dev/null; then
    rm -rf -- "\$stale_tmp"
  fi
done
GWAS_POST_TMP="\$GWAS_POST_TMP_ROOT/run.\$\$"
mkdir -p "\$GWAS_POST_TMP"
export TMPDIR="\$GWAS_POST_TMP"
gwas_post_cleanup_tmp(){
  rm -rf -- "\$GWAS_POST_TMP"
  rmdir -- "\$GWAS_POST_TMP_ROOT" 2>/dev/null || true
}
trap gwas_post_cleanup_tmp EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
RAW=$(q "$raw")
# SMALL and DO_STEP=small are the unchanged shared 0data.f.sh API.
# The project uses --hm3/HM3_MODE and std_hm3 for the filtering choice.
SMALL=$(q "$source_gwas")
FINAL=$(q "$final")
AWK_SNP=$(q "$awk_snp")
CIS_OUT=$(q "$cis_out")
QC_PREFIX=$(q "$qc_prefix")
CLUMP=$(q "$clump_dir")
COJO=$(q "$cojo_dir")
MERGED=$(q "$merged_prefix")
CLUMP_DONE=$(q "$clump_done")
COJO_DONE=$(q "$cojo_done")
REF_BAD_DIR=$(q "$dir_out/.project/$category/qc/ref_bad")

PHEF=$(q "$phef")
DATA_F=$(q "$data_f")
INDEX_F=$(q "$index_f")
PERF_F=$(q "$perf_f")
PLOT_F=$(q "$plot_f")
MPLOT_R=$(q "$mplot_r")
PLOT_METHOD=self
ADD_PANEL=$(q "$add_panel")
PLOT_WIDTH=$(q "$plot_width")
PLOT_HEIGHT=$(q "$plot_height")
PLOT_RES=$(q "$plot_res")
WRITE_SIG=$(q "$write_sig")
ADD_SIGNAL=$(q "$add_signal")
SIGNAL_MATCH_COL=$(q "$signal_match_col")
SIGNAL_MATCH_VALUE=$(q "$signal_match_value")
SIGNAL_LOCUS_POS=$(q "$signal_locus_pos")
SIGNAL_DISPLAY_COL=$(q "$signal_display_col")
PGS_STEPS=$(q "$step")
LEAD_REF_CACHE=$(q "$dir_out/.project/$category/lead_reference")
MH_PLOT_BED=$(q "$gwas_mh_plot_bed")
HM3=$(q "$hm3_file")
HM3_POS=$(q "$gwas_hm3_pos")
P_HM3=$(q "$p_hm3")
HM3_MODE=$(q "$hm3_mode")
THIN_MODE=$(q "$thin")
THIN_CHR_MAX=$(q "$thin_chr_max")
THIN_R=$(q "$thin_r")
PHE_R=$(q "$phe_r")
THIN_OUT=$(q "${final%.gz}.thin.gz")
MH_PNG=$(q "$mh_png")
MH_META=$(q "$mh_meta")
MH_FLAG=$(q "$mh_flag")
MH_SIG=$(q "$mh_sig")
COJO_FILE=$(q "${merged_prefix}.jma.cojo")
MPLOT_FLAG_FILE=$(q "$mplot_flag_file")
MAGMA_DIR=$(q "$magma_dir")
MAGMA_PREFIX=$(q "$magma_prefix")
MAGMA_ANNOT_CACHE=$(q "$magma_annot_cache")
DELETE_RAW=$(q "$delete_raw")

DO_STEP=$(q "$step")
REPLACE=$(q "$replace")
FILL_EAF=$(q "$fill_eaf")
FILL_N=$(q "$fill_n")
N_TOTAL=$(q "$n_total")
H2_HELPER=$(q "${SCRIPT_PATH%/*}/f/gwas_h2.py")
LIFTOVER_HELPER=$(q "${SCRIPT_PATH%/*}/f/gwas_liftover.py")
H2_SEX=$(q "$h2_sex")
H2_REF_LD=$(q "$h2_ref_ld")
H2_W_LD=$(q "$h2_w_ld")
H2_MERGE_ALLELES=$(q "$h2_merge_alleles")
H2_PYTHON=$(q "$h2_python")
H2_CONDA_ENV=$(q "$h2_conda_env")
H2_REQUESTED=$(q "$step")
DO_LIFTOVER=$(q "$liftOver")
CHAIN=$(q "$chain")
LIFTOVER_BIN=$(q "$liftover_bin")

CIS_BED=$(q "$cis_bed")
CIS_FLANK=$(q "$cis_flank")

P_LEAD=$(q "$p_lead")
LEAD_WINDOW=$(q "$lead_window")
CLUMP_KB=$(q "$clump_kb")
CHRS=$(q "$chrs")
REFGEN_CLUMP=$(q "$gwas_refGen_clump")
REFGEN_POP=$(q "$refGen_pop")
REFGEN_KEEP=$(q "$gwas_refGen_keep")
REFGEN_COJO=$(q "$gwas_refGen_cojo")

GRCH=$(q "$gwas_grch")
MAGMA_REF=$(q "$magma_ref")
GENE_LOC=$(q "$gwas_gene_loc")
SYNONYMS=$(q "$synonyms")
MAGMA_WINDOW=$(q "$magma_window")
MAGMA_N=$(q "$magma_N")
GWAS_N=$(q "$gwas_N")
PGS_DIR=$(q "$pgs_dir")
PGS_SCORE_FILE=$(q "$pgs_score_file")
PGS_OUTPUT=$(q "$pgs_output")
PGS_DONE=$(q "$pgs_done")
PGS_META=$(q "$pgs_meta")
PGS_PFILE_DIR=$(q "${gwas_pgs_pfile_dir%/}")
PGS_THREADS=$(q "$pgs_threads")

# Do not delete trait-wide error files here: another module may be using the
# same trait concurrently.  Each coordinator now owns a module-specific error.

if [[ ",\$DO_STEP," == *,format,* || "\$DO_STEP" == "all" ]]; then
  if [[ ! -s "\$RAW" ]]; then
    echo "ERROR: missing/empty raw GWAS for \$GWAS: \$RAW" >&2
    echo "ERROR: removing incomplete outputs for \$GWAS before re-run." >&2
    find "\$GWAS_DIR" -mindepth 1 -depth ! -path "\$GWAS_DIR/\$GWAS.log" -delete 2>/dev/null || true
    rm -f -- "\${QC_PREFIX}"* "\$MH_PNG" "\$MH_META" "\${MH_PNG}.meta.tsv" "\$MH_SIG"
    exit 1
  fi
fi

source "\$DATA_F"

# Format every raw data row into the standard 11-column schema.  Unlike
# std_hm3(), this intentionally performs no variant filtering.
std_format(){
  local src="\$1" out="\$2" header_out="\$3" tmp input_fs
  gwas_post_need_file "\$src"
  gwas_post_log "Format full GWAS: \$src -> \$out"
  gwas_clean_header_names "\$src" "\$header_out"
  SNP_col=\$(gwas_clean_col "\${SNP_col:-}"); CHR_col=\$(gwas_clean_col "\${CHR_col:-}")
  POS_col=\$(gwas_clean_col "\${POS_col:-}"); EA_col=\$(gwas_clean_col "\${EA_col:-}")
  NEA_col=\$(gwas_clean_col "\${NEA_col:-}"); EAF_col=\$(gwas_clean_col "\${EAF_col:-}")
  N_col=\$(gwas_clean_col "\${N_col:-}"); BETA_col=\$(gwas_clean_col "\${BETA_col:-}")
  SE_col=\$(gwas_clean_col "\${SE_col:-}"); P_col=\$(gwas_clean_col "\${P_col:-}")
  LOG10P_col=\$(gwas_clean_col "\${LOG10P_col:-}")
  [[ "\$SNP_col" -gt 0 && "\$CHR_col" -gt 0 && "\$POS_col" -gt 0 && "\$EA_col" -gt 0 && "\$NEA_col" -gt 0 && "\$BETA_col" -gt 0 && "\$SE_col" -gt 0 && "\$P_col" -gt 0 ]] || {
    echo "ERROR: required GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: \$src" >&2
    exit 1
  }
  input_fs=\$(gwas_clean_detect_fs "\$src"); tmp="\${out}.tmp.\$\$"
  gwas_post_zcat "\$src" | awk -v FS="\$input_fs" -v OFS='\t' \
    -v snp_col="\$SNP_col" -v chr_col="\$CHR_col" -v pos_col="\$POS_col" \
    -v ea_col="\$EA_col" -v nea_col="\$NEA_col" -v eaf_col="\$EAF_col" -v n_col="\$N_col" \
    -v beta_col="\$BETA_col" -v se_col="\$SE_col" -v p_col="\$P_col" -v logp_col="\$LOG10P_col" '
    function get(c, x){x=(c>0 ? \$c : "");gsub(/\r/,"",x);return x}
    function val(c, x){x=get(c); return x=="" ? "NA" : x}
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
    function normchr(x){gsub(/^chr/,"",x); if(x=="X")x="23"; if(x=="Y")x="24"; if(x=="MT"||x=="M")x="25"; return x}
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P"; next}
    {chr=normchr(get(chr_col)); pos=get(pos_col); snp=get(snp_col); ea=get(ea_col); nea=get(nea_col)
     if(snp==""||snp=="NA"||snp=="."){snp=chr":"pos; if(ea!="")snp=snp":"ea; if(nea!="")snp=snp":"nea}
     p=(p_col>0 ? get(p_col) : ""); lp=(logp_col>0 ? get(logp_col) : "")
     if(p=="" && isnum(lp))p=10^(-lp); if(lp=="" && isnum(p) && p>0)lp=-log(p)/log(10)
     if(snp=="")snp="NA"; if(chr=="")chr="NA"; if(pos=="")pos="NA"; if(ea=="")ea="NA"; if(nea=="")nea="NA"
     if(p=="")p="NA"; if(lp=="")lp="NA"
     print snp,chr,pos,ea,nea,val(eaf_col),val(n_col),val(beta_col),val(se_col),p,lp}' | \
    { IFS= read -r format_header; printf '%s\n' "\$format_header"; sort -t \$'\t' -k2,2n -k3,3n; } | \
    gwas_clean_compress > "\$tmp"
  gzip -t "\$tmp"; mv -f "\$tmp" "\$out"
  gwas_post_zcat "\$out" | awk -v g="\$GWAS" 'NR==1{next} END{print g"\t"NR-1}' > "\${QC_PREFIX}.format.nrow.tsv"
}

gwas_post_validate_format_columns(){
  local src="\$1" header_out="\${QC_PREFIX}.header.required.txt"
  gwas_clean_header_names "\$src" "\$header_out"
  SNP_col=\$(gwas_clean_col "\${SNP_col:-}"); CHR_col=\$(gwas_clean_col "\${CHR_col:-}")
  POS_col=\$(gwas_clean_col "\${POS_col:-}"); EA_col=\$(gwas_clean_col "\${EA_col:-}")
  NEA_col=\$(gwas_clean_col "\${NEA_col:-}"); BETA_col=\$(gwas_clean_col "\${BETA_col:-}")
  SE_col=\$(gwas_clean_col "\${SE_col:-}"); P_col=\$(gwas_clean_col "\${P_col:-}")
  [[ "\$SNP_col" -gt 0 && "\$CHR_col" -gt 0 && "\$POS_col" -gt 0 && "\$EA_col" -gt 0 && "\$NEA_col" -gt 0 && "\$BETA_col" -gt 0 && "\$SE_col" -gt 0 && "\$P_col" -gt 0 ]] || {
    echo "ERROR: required raw GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: \$src" >&2
    exit 1
  }
  gwas_post_log "format column map: SNP=\$SNP_col CHR=\$CHR_col POS=\$POS_col EA=\$EA_col NEA=\$NEA_col BETA=\$BETA_col SE=\$SE_col P=\$P_col"
}

gwas_post_fill_missing_fields(){
  local target="\$1" tmp eaf_tmp
  local -a eaf_keep_args=()
  [[ -s "\$target" ]] || { echo "ERROR: formatted GWAS is missing: \$target" >&2; return 1; }

  if [[ -n "\$FILL_N" ]]; then
    tmp="\${target}.fill_n.tmp.\$\$"
    gwas_post_log "fill missing/invalid N with \$FILL_N: \$target"
    gwas_post_zcat "\$target" | awk -v FS='\t' -v OFS='\t' -v fill_n="\$FILL_N" -v audit="\${QC_PREFIX}.fill_n.tsv" '
      function isnum(x){return x ~ /^[0-9]+([.][0-9]+)?\$/ && x+0>0}
      NR==1{for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);c[h]=i}if(!("N" in c)){print "ERROR: standardized GWAS lacks N" > "/dev/stderr";exit 2}print;next}
      {if(isnum(\$(c["N"])))kept++;else{\$(c["N"])=fill_n;filled++}print}
      END{print "STATUS\tN" > audit;print "existing\t" kept+0 >> audit;print "filled\t" filled+0 >> audit}
    ' | gzip -c > "\$tmp"
    gzip -t "\$tmp"; mv -f "\$tmp" "\$target"
  fi

  if [[ "\$FILL_EAF" == TRUE ]]; then
    declare -F match_EAF >/dev/null 2>&1 || gwas_clean_load_phef
    eaf_tmp="\${target}.fill_eaf.tmp.\$\$.gz"
    gwas_post_log "fill missing/invalid EAF from 1KG GRCh\$GRCH: \$target"
    [[ -z "\$REFGEN_KEEP" ]] || eaf_keep_args=(--keep "\$REFGEN_KEEP")
    match_EAF --reference "\$REFGEN_CLUMP" "\${eaf_keep_args[@]}" --output "\$eaf_tmp" \
      --audit "\${QC_PREFIX}.fill_eaf.tsv" "\$target"
    gzip -t "\$eaf_tmp"; mv -f "\$eaf_tmp" "\$target"
  fi
}

gwas_post_prune_magma_dir(){
  [[ -d "\$MAGMA_DIR" ]] || return 0
  find "\$MAGMA_DIR" -mindepth 1 -maxdepth 1 -type f \
    ! -name '*.genes.out' ! -name '*.genes.raw' ! -name '*.log' \
    ! -name 'magma.meta.tsv' ! -name 'magma.done' -delete
}

gwas_post_magma_annotation(){
  local snploc="\${1:-}" cache_key window_tag resource_tag cache_dir annot meta
  local lock_file lock_fd tmp_prefix annot_tmp meta_tmp nloc snploc_hash
  command -v flock >/dev/null 2>&1 || { echo "ERROR: flock not found; required for the shared MAGMA annotation cache" >&2; exit 1; }

  window_tag=\$(printf '%s' "\$MAGMA_WINDOW" | tr -c 'A-Za-z0-9._-' '_')
  resource_tag=\$(printf '%s\n%s\n' "\$MAGMA_REF" "\$GENE_LOC" | sha256sum | awk '{print substr(\$1,1,16)}')
  cache_key="v2.GRCh\${GRCH}.window_\${window_tag}.magma_\${resource_tag}"
  cache_dir="\$MAGMA_ANNOT_CACHE/v2/GRCh\$GRCH/window_\$window_tag/magma_\$resource_tag"
  annot="\$cache_dir/genes.annot"
  meta="\$cache_dir/annotation.meta.tsv"
  mkdir -p "\$(dirname "\$cache_dir")"

  lock_file="\${cache_dir}.lock"
  exec {lock_fd}> "\$lock_file"
  flock "\$lock_fd"
  if [[ ! -s "\$annot" || ! -s "\$meta" ]]; then
    mkdir -p "\$cache_dir"
    if [[ -z "\$snploc" ]]; then
      snploc="\$GWAS_POST_TMP/\$GWAS.snp.loc"
      gwas_post_log "Prepare MAGMA annotation coordinates from \$FINAL"
      gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
        NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}
        {s=\$(c["SNP"]);ch=\$(c["CHR"]);pos=\$(c["POS"]);gsub(/^chr/,"",ch)
         if(s!=""&&s!="NA"&&s!="."&&ch~/^([1-9]|1[0-9]|2[0-3])\$/&&pos~/^[0-9]+\$/)print s,ch,pos}' |
        sort -T "\$GWAS_POST_TMP" -k1,1 -k2,2n -k3,3n -u > "\$snploc"
    fi
    gwas_post_need_file "\$snploc"
    nloc=\$(wc -l < "\$snploc" | tr -d ' ')
    (( nloc > 1000 )) || { echo "ERROR: too few MAGMA annotation SNPs: \$nloc" >&2; return 1; }
    snploc_hash=\$(sha256sum "\$snploc" | awk '{print \$1}')
    tmp_prefix="\$GWAS_POST_TMP/annotation.\$\$"
    rm -f -- "\${tmp_prefix}"*
    gwas_post_log "Build shared MAGMA annotation from GWAS rsID coordinates: \$annot"
    if [[ "\$MAGMA_WINDOW" == "0,0" || "\$MAGMA_WINDOW" == 0 ]]; then
      magma --annotate --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$tmp_prefix"
    else
      magma --annotate window="\$MAGMA_WINDOW" --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$tmp_prefix"
    fi
    gwas_post_need_file "\$tmp_prefix.genes.annot"
    annot_tmp="\${annot}.tmp.\$\$"
    mv -f "\$tmp_prefix.genes.annot" "\$annot_tmp"
    mv -f "\$annot_tmp" "\$annot"
    meta_tmp="\${meta}.tmp.\$\$"
    {
      printf 'key\tvalue\n'
      printf 'grch\t%s\nwindow_kb\t%s\ngene_loc\t%s\nld_reference\t%s\nsource_gwas\t%s\nsnp_loc_n\t%s\nsnp_loc_sha256\t%s\n' \
        "\$GRCH" "\$MAGMA_WINDOW" "\$GENE_LOC" "\$MAGMA_REF" "\$GWAS" "\$nloc" "\$snploc_hash"
    } > "\$meta_tmp"
    mv -f "\$meta_tmp" "\$meta"
  else
    gwas_post_log "Reuse shared MAGMA annotation: \$annot"
  fi
  flock -u "\$lock_fd"
  exec {lock_fd}>&-

  MAGMA_ANNOT="\$annot"
  MAGMA_ANNOT_KEY="\$cache_key"
  MAGMA_SNPLOC_N=\$(awk -F '\t' '\$1=="snp_loc_n"{print \$2;exit}' "\$meta")
  MAGMA_SNPLOC_HASH=\$(awk -F '\t' '\$1=="snp_loc_sha256"{print \$2;exit}' "\$meta")
}

gwas_post_magma(){
  [[ ",\$DO_STEP," == *,magma,* ]] || return 0
  [[ -s "\$FINAL" ]] || { echo "ERROR: missing or empty file: \$FINAL" >&2; exit 1; }

  if [[ "\$GRCH" == "38" ]]; then
    echo " [\$GWAS] GRCh build 38" >&2
  else
    echo "[\$GWAS] GRCh build 37" >&2
  fi
  echo "   gene location : \$GENE_LOC" >&2
  echo "   LD reference  : \$MAGMA_REF" >&2

  if [[ "\$REPLACE" != TRUE && -s "\$MAGMA_DIR/magma.done" && -s "\$MAGMA_PREFIX.genes.out" && -s "\$MAGMA_PREFIX.genes.raw" ]]; then
    gwas_post_prune_magma_dir
    echo "[\$(date '+%F %T')] MAGMA exists: \$MAGMA_PREFIX.genes.out (GRCh\$GRCH)" >&2
    return 0
  fi
  command -v magma >/dev/null 2>&1 || { echo "ERROR: magma not found in PATH" >&2; exit 1; }
  for ext in bed bim fam; do [[ -s "\$MAGMA_REF.\$ext" ]] || { echo "ERROR: missing or empty file: \$MAGMA_REF.\$ext" >&2; exit 1; }; done
  [[ -s "\$GENE_LOC" ]] || { echo "ERROR: missing or empty file: \$GENE_LOC" >&2; exit 1; }
  [[ -s "\$SYNONYMS" ]] || { echo "ERROR: missing or empty file: \$SYNONYMS" >&2; exit 1; }
  mkdir -p "\$MAGMA_DIR"
  # Never leave an old completion marker/meta file behind during a rerun.
  rm -f "\$MAGMA_DIR/magma.done" "\$MAGMA_DIR/magma.meta.tsv"

  if [[ -z "\$MAGMA_N" && -s "\$GWAS_DIR/\$GWAS.magma.N" ]]; then
    MAGMA_N=\$(awk 'NF{print \$1; exit}' "\$GWAS_DIR/\$GWAS.magma.N")
  fi
  if [[ -z "\$MAGMA_N" && ! "\$GWAS_N" =~ ^[1-9][0-9]*\$ ]]; then
    echo "ERROR: invalid default gwas_N: \$GWAS_N" >&2; exit 1
  fi
  if [[ -n "\$MAGMA_N" && ! "\$MAGMA_N" =~ ^[1-9][0-9]*\$ ]]; then
    echo "ERROR: invalid MAGMA sample size for \$GWAS: \$MAGMA_N" >&2; exit 1
  fi

  header=\$(set +o pipefail; gwas_clean_zcat "\$FINAL" | head -1)
  for col in SNP CHR POS P; do
    awk -F '\t' -v c="\$col" '{for(i=1;i<=NF;i++)if(\$i==c)exit 0;exit 1}' <<< "\$header" ||
      { echo "ERROR: \$FINAL lacks required column \$col" >&2; exit 1; }
  done
  snploc="\$GWAS_POST_TMP/\$GWAS.snp.loc"; pval="\$GWAS_POST_TMP/\$GWAS.pval"
  gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
    NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}
    {s=\$(c["SNP"]);ch=\$(c["CHR"]);pos=\$(c["POS"]);gsub(/^chr/,"",ch)
     if(s!=""&&s!="NA"&&s!="."&&ch~/^([1-9]|1[0-9]|2[0-3])\$/&&pos~/^[0-9]+\$/)print s,ch,pos}' |
    sort -k1,1 -k2,2n -k3,3n -u > "\$snploc"

  has_n=FALSE
  awk -F '\t' '{for(i=1;i<=NF;i++)if(\$i=="N")exit 0;exit 1}' <<< "\$header" && has_n=TRUE || true
  if [[ -n "\$MAGMA_N" ]]; then
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"
  elif [[ "\$has_n" == TRUE ]]; then
    { printf 'SNP\tP\tN\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);n=\$(c["N"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1&&n+0>=50)print s,p,n}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg='ncol=N'
  else
    MAGMA_N="\$GWAS_N"
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"
  fi
  nloc=\$(wc -l < "\$snploc" | tr -d ' '); npval=\$((\$(wc -l < "\$pval" | tr -d ' ')-1))
  if (( npval <= 1000 )) && [[ "\$narg" == ncol=N ]]; then
    echo "WARNING: no usable per-SNP N for \$GWAS; falling back to gwas_N=\$GWAS_N" >&2
    MAGMA_N="\$GWAS_N"
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"; npval=\$((\$(wc -l < "\$pval" | tr -d ' ')-1))
  fi
  (( nloc > 1000 && npval > 1000 )) || { echo "ERROR: too few MAGMA SNPs: snploc=\$nloc pval=\$npval" >&2; exit 1; }
  gwas_post_magma_annotation "\$snploc"
  magma --bfile "\$MAGMA_REF" synonyms="\$SYNONYMS" --pval "\$pval" "\$narg" --gene-annot "\$MAGMA_ANNOT" --out "\$MAGMA_PREFIX"
  [[ -s "\$MAGMA_PREFIX.genes.out" ]] || { echo "ERROR: MAGMA output missing: \$MAGMA_PREFIX.genes.out" >&2; exit 1; }
  [[ -s "\$MAGMA_PREFIX.genes.raw" ]] || { echo "ERROR: MAGMA intermediate gene result missing: \$MAGMA_PREFIX.genes.raw" >&2; exit 1; }
  meta_tmp="\$MAGMA_DIR/magma.meta.tsv.tmp.\$\$"
  { printf 'key\tvalue\n'; printf 'gwas\t%s\ngrch\t%s\ngene_loc\t%s\nld_reference\t%s\nsynonyms\t%s\nwindow_kb\t%s\nsnp_loc_n\t%s\npval_n\t%s\nannotation_cache\t%s\nannotation_key\t%s\nsnp_loc_sha256\t%s\n' \
      "\$GWAS" "\$GRCH" "\$GENE_LOC" "\$MAGMA_REF" "\$SYNONYMS" "\$MAGMA_WINDOW" "\$nloc" "\$npval" "\$MAGMA_ANNOT" "\$MAGMA_ANNOT_KEY" "\$MAGMA_SNPLOC_HASH"; } > "\$meta_tmp"
  mv -f "\$meta_tmp" "\$MAGMA_DIR/magma.meta.tsv"
  gwas_post_prune_magma_dir
  date '+%F %T' > "\$MAGMA_DIR/magma.done"
  echo "[\$(date '+%F %T')] MAGMA done: \$MAGMA_PREFIX.genes.out" >&2
}

gwas_post_mplot(){
  [[ ",\$DO_STEP," == *,mplot,* || "\$DO_STEP" == "all" ]] || return 0
  gwas_post_need_file "\$FINAL"
  if [[ "\$REPLACE" != "TRUE" && -s "\$MH_PNG" ]]; then
    gwas_post_log "Manhattan plot exists: \$MH_PNG"
    return 0
  fi
  command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found; required for Manhattan plotting" >&2; exit 1; }
  mkdir -p "\$(dirname "\$MH_PNG")"
  plot_input="\$FINAL"
  if [[ "\$HM3_MODE" == "FALSE" ]]; then
    plot_input="\${QC_PREFIX}.mplot.hm3.gz"
    gwas_post_log "Subset Manhattan input to HM3 or P < 0.001: \$plot_input"
    gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' -v hm3_file="\$HM3" -v pthr='0.001' '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
      BEGIN{while((getline x<hm3_file)>0){split(x,a,/[ \t]+/); if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1} close(hm3_file)}
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i; print; next}
      (("SNP" in c)&&\$(c["SNP"]) in hm3) || (("P" in c)&&isnum(\$(c["P"]))&&\$(c["P"])+0<pthr){print}' | gwas_clean_compress > "\$plot_input"
    gzip -t "\$plot_input"
  fi
  gwas_post_log "Manhattan plot: \$plot_input -> \$MH_PNG"
  Rscript - "\$plot_input" "\$MH_PNG" "\$GWAS" "\$PLOT_F" "\$CIS_BED" "\$MH_PLOT_BED" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]]
output <- args[[2]]
gwas <- args[[3]]
plot_f <- args[[4]]
cis_bed <- if (nzchar(args[[5]])) args[[5]] else NULL
mh_plot_bed <- if (nzchar(args[[6]])) args[[6]] else NULL
source(plot_f)

cis_gene_arg <- NULL
if (!is.null(cis_bed) && file.exists(cis_bed)) {
  cis_rows <- tryCatch(utils::read.table(cis_bed, header=FALSE, comment.char="#", stringsAsFactors=FALSE),
                       error=function(e) NULL)
  if (!is.null(cis_rows) && ncol(cis_rows) >= 4 && any(as.character(cis_rows[[4]]) == gwas))
    cis_gene_arg <- gwas
}

png_args <- list(filename = output, width = 13.333, height = 7.5, units = "in", res = 300)
if (capabilities("cairo")) png_args[["type"]] <- "cairo"
do.call(grDevices::png, png_args)
on.exit(grDevices::dev.off(), add = TRUE)
graphics::par(mar = c(5, 5, 3, 1) + 0.1)
print(mh_plot(input, col=c("gray", "darkgray"), cis_gene=cis_gene_arg, cis_bed=cis_bed, mh_plot_bed=mh_plot_bed,
              cis_color="red", other_top_color="green",
              other_top_max_per_chr=2, locus_size=1e6, main=gwas))
RSCRIPT
  [[ "\$plot_input" == "\$FINAL" ]] || rm -f "\$plot_input"
  [[ -s "\$MH_PNG" ]] || { echo "ERROR: Manhattan plot was not created: \$MH_PNG" >&2; exit 1; }
}

gwas_post_pgs_current(){
  [[ -s "\$PGS_OUTPUT" && -s "\$PGS_DONE" ]]
}

gwas_post_prune_pgs_dir(){
  local f
  [[ -d "\$PGS_DIR" ]] || return 0
  for f in "\$PGS_DIR/\${GWAS}.chr"*; do
    [[ -f "\$f" ]] || continue
    rm -f -- "\$f"
  done
}

gwas_post_prune_empty_pgs_dir(){
  local f
  [[ -d "\$PGS_DIR" ]] || return 0
  for f in "\$PGS_DIR/\${GWAS}.chr"*; do
    [[ -f "\$f" ]] || continue
    case "\$f" in *.log) continue;; esac
    rm -f -- "\$f"
  done
}

gwas_post_pgs(){
  [[ ",\$PGS_STEPS," == *,pgs,* || "\$PGS_STEPS" == "all" ]] || return 0
  if [[ "\$REPLACE" != "TRUE" ]] && gwas_post_pgs_current; then
    if [[ -s "\$PGS_META" ]] && awk -F '\t' '\$1=="matched_variants"&&\$2=="none"{found=1} END{exit !found}' "\$PGS_META"; then
      gwas_post_prune_empty_pgs_dir
    else
      gwas_post_prune_pgs_dir
    fi
    gwas_post_log "PGS exists and is current: \$PGS_OUTPUT"
    return 0
  fi
  gwas_post_need_file "\$PGS_SCORE_FILE"
  gwas_clean_load_phef
  declare -F pgs_plink_calc >/dev/null 2>&1 || { echo "ERROR: pgs_plink_calc is missing from \$PHEF" >&2; exit 1; }
  mkdir -p "\$PGS_DIR"
  rm -f -- "\$PGS_DONE"

  local status_file="\$GWAS_POST_TMP/\$GWAS.pgs.status.tsv"
  local pgs_match_status score_chrs intermediate_status meta_tmp done_tmp done_status=complete
  pgs_plink_calc \
    --input "\$PGS_SCORE_FILE" \
    --pfile-dir "\$PGS_PFILE_DIR" \
    --output "\$PGS_OUTPUT" \
    --label "\$GWAS" \
    --work-dir "\$PGS_DIR" \
    --status-file "\$status_file" \
    --threads "\$PGS_THREADS"
  pgs_match_status=\$(awk -F '\t' '\$1=="matched_variants"{print \$2}' "\$status_file")
  score_chrs=\$(awk -F '\t' '\$1=="chromosomes"{print \$2}' "\$status_file")
  intermediate_status=\$(awk -F '\t' '\$1=="chromosome_intermediates"{print \$2}' "\$status_file")
  [[ "\$pgs_match_status" == scored || "\$pgs_match_status" == none ]] || { echo "ERROR: invalid PGS status: \$pgs_match_status" >&2; exit 1; }
  [[ -n "\$score_chrs" && -n "\$intermediate_status" ]] || { echo "ERROR: incomplete PGS status file: \$status_file" >&2; exit 1; }
  [[ "\$pgs_match_status" != none ]] || done_status=complete_no_matched_variants

  meta_tmp="\${PGS_META}.tmp.\$\$"
  {
    printf 'key\tvalue\n'
    printf 'gwas\t%s\nsource\t%s\neffect_allele\trefA\nscore_weight\tbJ\nscore_stat\tdosage_weighted_sum\nmissing_genotype\tno_mean_imputation\nmatched_variants\t%s\npfile_dir\t%s\nchromosomes\t%s\noutput\t%s\nchromosome_intermediates\t%s\n' \
      "\$GWAS" "\$PGS_SCORE_FILE" "\$pgs_match_status" "\$PGS_PFILE_DIR" "\$score_chrs" "\$PGS_OUTPUT" "\$intermediate_status"
  } > "\$meta_tmp"
  mv -f "\$meta_tmp" "\$PGS_META"
  done_tmp="\${PGS_DONE}.tmp.\$\$"
  printf 'GWAS\tSTATUS\tTIME\n%s\t%s\t%s\n' "\$GWAS" "\$done_status" "\$(date '+%F %T')" > "\$done_tmp"
  mv -f "\$done_tmp" "\$PGS_DONE"
  if [[ "\$pgs_match_status" == none ]]; then
    gwas_post_log "PGS done with no matched variants: \$PGS_OUTPUT (header only; PLINK2 logs retained)"
  else
    gwas_post_log "PGS done: \$PGS_OUTPUT (per-chromosome working files removed from \$PGS_DIR)"
  fi
}

# Load the index and performance overrides after all legacy definitions so
# the optimized implementations replace them deliberately.
gwas_post_magma_annotation_impl=\$(declare -f gwas_post_magma_annotation)
source "\$INDEX_F"
source "\$PERF_F"
# The performance helper's MAGMA override derives annotation IDs from
# REFGEN_CLUMP.  Those IDs are not guaranteed to match MAGMA_REF (for example,
# GRCh38 pvars use CHR:POS:REF:ALT while g1000_eur.bim uses rsIDs).  Keep the
# optimized p-value preparation, but restore the compatible annotation builder.
eval "\$gwas_post_magma_annotation_impl"
unset gwas_post_magma_annotation_impl

if [[ "\$DO_STEP" == "all" ]]; then
  gwas_post_validate_format_columns "\$RAW"
  DO_STEP=small; gwas_clean_run_core
  gwas_post_fill_missing_fields "\$SMALL"
  gwas_post_ensure_index "\$SMALL"
  if [[ "\$DO_LIFTOVER" == TRUE ]]; then
    DO_STEP=liftover; gwas_clean_run_core
  fi
  DO_STEP=all
  gwas_post_ensure_index "\$FINAL"
  gwas_post_thin
  gwas_post_prepare_views
  gwas_post_mplot
  [[ -z "\$CIS_BED" ]] || gwas_clean_make_cis
else
  REQUESTED_STEP="\$DO_STEP"
  if [[ ",\$REQUESTED_STEP," == *,format,* ]]; then
    gwas_post_validate_format_columns "\$RAW"
    DO_STEP=small; gwas_clean_run_core
    gwas_post_fill_missing_fields "\$SMALL"
    gwas_post_ensure_index "\$SMALL"
  fi
  if [[ ",\$REQUESTED_STEP," == *,liftover,* ]]; then
    DO_STEP=liftover; gwas_clean_run_core
  fi
  DO_STEP="\$REQUESTED_STEP"
  gwas_post_ensure_index "\$FINAL"
  gwas_post_thin
  gwas_post_prepare_views
  gwas_post_magma
  gwas_post_mplot
  if [[ ",\$REQUESTED_STEP," == *,cis,* ]]; then gwas_clean_make_cis; fi
  if [[ ",\$REQUESTED_STEP," == *,lead,* ]]; then DO_STEP=lead; else DO_STEP=none; fi
fi

$(declare -f gwas_lead_marker_matches)
$(declare -f gwas_post_run_cojo)

gwas_post_clump_complete(){
  gwas_lead_marker_matches "\$CLUMP_DONE" "\$GWAS" clump "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS" || return 1
  [[ -s "\$CLUMP_DONE" ]] && {
    [[ -s "\${MERGED}.clumps" ]] || gwas_post_lead_awk_has_no_rows ||
      awk -F '\t' 'NR==2&&\$2=="clump"&&\$6=="no_reference_matched_variants"{ok=1}END{exit !ok}' "\$CLUMP_DONE"
  }
}

gwas_post_cojo_complete(){
  gwas_lead_marker_matches "\$COJO_DONE" "\$GWAS" cojo "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS" || return 1
  [[ -s "\$COJO_DONE" ]] && {
    [[ -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]] || gwas_post_lead_awk_has_no_rows ||
      awk -F '\t' 'NR==2&&\$2=="cojo"&&\$6=="no_reference_matched_variants"{ok=1}END{exit !ok}' "\$COJO_DONE"
  }
}

gwas_post_lead_awk_has_no_rows(){
  [[ -s "\$AWK_SNP" ]] || return 1
  awk 'NR>1{found=1} END{exit found ? 1 : 0}' "\$AWK_SNP"
}

gwas_post_lead_done(){
  [[ "\$REPLACE" != "TRUE" ]] || return 1
  [[ -s "\$AWK_SNP" ]] || return 1
  [[ ! -d "\$CLUMP" && ! -d "\$COJO" ]] || return 1
  gwas_post_clump_complete && gwas_post_cojo_complete
}

gwas_post_mark_phase_done(){
  local phase="\$1" status="\${2:-complete}" marker tmp
  case "\$phase" in
    clump) marker="\$CLUMP_DONE" ;;
    cojo) marker="\$COJO_DONE" ;;
    *) echo "ERROR: unknown lead phase: \$phase" >&2; return 2 ;;
  esac
  tmp="\${marker}.tmp.\$\$"
  mkdir -p "\$(dirname "\$marker")"
  {
    printf 'GWAS\tPHASE\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "\$GWAS" "\$phase" "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS" "\$status"
  } > "\$tmp"
  mv -f "\$tmp" "\$marker"
}

# Exclude ambiguous PLINK IDs. Reuse the build/chromosome output until deleted.
gwas_post_filter_clump_multiallelic(){
  local ref="\$1" tag="\$2" assoc="\$3" pvar pvar_real bad_dir bad ref_meta tmp meta_tmp all_ids multi_ids dup_ids filtered n_before n_after n_removed
  if [[ ! -s "\${ref}.pvar" && -s "\${ref}.pvar.zst" ]]; then
    pvar="\${ref}.pvar.zst"
  else
    pvar="\${ref}.pvar"
  fi
  gwas_post_need_file "\$pvar"
  pvar_real=\$(readlink -f -- "\$pvar")
  [[ -n "\$pvar_real" ]] || { echo "ERROR: cannot resolve reference pvar: \$pvar" >&2; exit 1; }
  bad_dir="\$REF_BAD_DIR/GRCh\${GRCH}"
  bad="\$bad_dir/\${tag}.ambiguous.snp"
  ref_meta="\${bad}.reference.tsv"
  mkdir -p "\$bad_dir"

  command -v flock >/dev/null 2>&1 || { echo "ERROR: flock is required for shared lead reference caches" >&2; exit 1; }
  lock_file="\${bad}.lock"
  exec {bad_lock_fd}> "\$lock_file"
  flock "\$bad_lock_fd"
  if [[ ! -e "\$bad" ]]; then
    tmp="\${bad}.tmp.\$\$"
    meta_tmp="\${ref_meta}.tmp.\$\$"
    all_ids="\${tmp}.all"
    multi_ids="\${tmp}.multi"
    dup_ids="\${tmp}.dup"
    : > "\$multi_ids"
    : > "\$dup_ids"
    gwas_clean_zcat "\$pvar" | awk -v FS='\t' -v multi="\$multi_ids" '
      \$1=="#CHROM" {for(i=1;i<=NF;i++){x=\$i; sub(/^#/,"",x); if(x=="ID")id=i; if(x=="ALT")alt=i} next}
      /^##/ {next}
      id>0 && \$id!="" && \$id!="." {print \$id; if(alt>0 && index(\$alt,",")>0) print \$id > multi}
    ' > "\$all_ids"
    sort "\$all_ids" | uniq -d > "\$dup_ids"
    cat "\$multi_ids" "\$dup_ids" 2>/dev/null | sort -u > "\$tmp"
    mv -f "\$tmp" "\$bad"
    printf 'GRCH\tPVAR\n%s\t%s\n' "\$GRCH" "\$pvar_real" > "\$meta_tmp"
    mv -f "\$meta_tmp" "\$ref_meta"
    rm -f "\$all_ids" "\$multi_ids" "\$dup_ids"
  fi
  flock -u "\$bad_lock_fd"
  exec {bad_lock_fd}>&-

  n_before=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  filtered="\${assoc}.filtered.\$\$"
  awk -v FS='\t' -v OFS='\t' '
    FILENAME==ARGV[1] {bad[\$1]=1; next}
    FNR==1 {print; next}
    !(\$1 in bad) {print}
  ' "\$bad" "\$assoc" > "\$filtered"
  mv -f "\$filtered" "\$assoc"
  n_after=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  n_removed=\$((n_before-n_after))
  printf 'GWAS\tCHR\tN_BEFORE\tN_REMOVED_AMBIGUOUS_REF_ID\tN_AFTER\n%s\t%s\t%s\t%s\t%s\n' \
    "\$GWAS" "\$tag" "\$n_before" "\$n_removed" "\$n_after" > "\${QC_PREFIX}.\${tag}.clump_multiallelic.tsv"
  gwas_post_log "clump ambiguous-reference-ID filter \$tag: before=\$n_before removed=\$n_removed after=\$n_after cache=GRCh\${GRCH}/\${tag}.ambiguous.snp"
}

gwas_post_prep_lead_inputs(){
  local refs="\$1" suffix=".tmp.\$\$" tag assoc ma
  rm -f -- "\${QC_PREFIX}".*.cojo_skip.log
  while IFS=\$'\t' read -r _ tag _ assoc ma _; do
    rm -f "\${assoc}\${suffix}" "\${ma}\${suffix}"
  done < "\$refs"
  gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' -v refs="\$refs" \
    -v suffix="\$suffix" -v default_n="\$GWAS_N" -v p_lead="\$P_LEAD" '
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
    function valid(x){return x!="" && x!="NA" && x!="."}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;return x}
    function get(name){return name in c ? \$(c[name]) : "NA"}
    BEGIN{
      while((getline line < refs)>0){split(line,a,"\t");chr=a[3];if(chr=="")continue;
        assoc[chr]=a[4] suffix;ma[chr]=a[5] suffix
        print "SNP","CHR","POS","EA","NEA","P" > assoc[chr]
        print "SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","N" > ma[chr]}
      close(refs)
    }
    NR==1{
      for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);c[h]=i}
      required[1]="SNP";required[2]="CHR";required[3]="POS";required[4]="EA"
      required[5]="NEA";required[6]="BETA";required[7]="SE";required[8]="P"
      for(i=1;i<=8;i++)if(!(required[i] in c)){print "ERROR: required GWAS column is missing: " required[i] > "/dev/stderr";fatal=1}
      if(fatal)exit 2
      next
    }
    {
      snp=get("SNP");chr=normchr(get("CHR"));pos=get("POS");ea=get("EA");nea=get("NEA")
      beta=get("BETA");se=get("SE");p=get("P");eaf=get("EAF");n=get("N")
      # Filter significant variants before the in-memory ID match.
      if(!(chr in assoc)||pos!~/^[0-9]+\$/||pos+0<=0||!valid(snp)||!valid(ea)||!valid(nea)||!isnum(p)||p+0>p_lead)next
      print snp,chr,pos,ea,nea,p > assoc[chr]
      if(!isnum(beta)||!isnum(se))next
      if(!isnum(n))n=default_n
      print snp,chr,pos,ea,nea,eaf,beta,se,p,n > ma[chr]
    }
  '
  while IFS=\$'\t' read -r _ _ _ assoc ma _; do
    mv -f "\${assoc}\${suffix}" "\$assoc";mv -f "\${ma}\${suffix}" "\$ma"
  done < "\$refs"
}

gwas_post_subset_match_reference(){
  local ref="\$1" query="\$2" out="\$3" kind="\$4" ref_file="" ext
  case "\$kind" in
    pvar)
      case "\$ref" in
        *.pvar|*.pvar.gz|*.pvar.bgz|*.pvar.zst) ref_file="\$ref" ;;
        *)
          for ext in .pvar .pvar.gz .pvar.bgz .pvar.zst; do
            [[ -s "\${ref}\${ext}" ]] && { ref_file="\${ref}\${ext}"; break; }
          done
          ;;
      esac
      ;;
    bim)
      case "\$ref" in
        *.bim|*.bim.gz|*.bim.bgz) ref_file="\$ref" ;;
        *)
          for ext in .bim .bim.gz .bim.bgz; do
            [[ -s "\${ref}\${ext}" ]] && { ref_file="\${ref}\${ext}"; break; }
          done
          ;;
      esac
      ;;
    *) echo "ERROR: invalid reference subset kind: \$kind" >&2; return 2 ;;
  esac
  gwas_post_need_file "\$ref_file"

  # Restrict the reference to query coordinates before match_SNP.
  awk -v FS='[ \t]+' -v kind="\$kind" '
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x=="" ? "0" : x}
    NR==FNR{
      if(FNR==1){for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);if(h=="CHR")qchr=i;else if(h=="POS")qpos=i}next}
      if(!qchr||!qpos){fatal=1;next}
      if(\$(qpos)~/^[0-9]+\$/)wanted[normchr(\$(qchr)) SUBSEP \$(qpos)]=1
      next
    }
    kind=="pvar" && /^##/{print;next}
    kind=="pvar" && /^#CHROM/{print;header=1;next}
    kind=="pvar"{
      if((normchr(\$1) SUBSEP \$2) in wanted){print;kept++}
      next
    }
    kind=="bim"{
      if(NF>=4 && (normchr(\$1) SUBSEP \$4) in wanted){print;kept++}
      next
    }
    END{
      if(fatal){print "ERROR: lead query requires CHR and POS columns" > "/dev/stderr";exit 2}
      if(kind=="pvar"&&!header){print "ERROR: invalid .pvar header" > "/dev/stderr";exit 2}
      # Keep an empty BIM subset non-empty so match_SNP can report all query
      # rows as unmatched instead of rejecting the reference argument.
      if(kind=="bim"&&kept==0)print "0\t.\t0\t0\tN\tN"
    }
  ' "\$query" <(gwas_clean_zcat "\$ref_file") > "\$out"
}

# GCTA treats only numeric chromosomes up to --autosome-num as usable.  PLINK
# references commonly encode chrX/chrY as X/Y, which makes GCTA silently load
# just one unusable record from an otherwise complete sex-chromosome BIM.  For
# sex chromosomes, make a small per-trait BED containing only the already
# matched COJO candidates and emit numeric chromosome codes (X=23, Y=24).
gwas_post_prepare_gcta_bfile(){
  local source="\$1" chr="\$2" ma="\$3" out="\$4"
  GCTA_BFILE="\$source"
  GCTA_CHR_ARGS=()
  (( chr > 22 )) || return 0
  command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 is required to prepare chr\$chr for GCTA" >&2; return 1; }
  gwas_post_need_file "\${source}.bed"
  gwas_post_need_file "\${source}.bim"
  gwas_post_need_file "\${source}.fam"
  gwas_post_need_file "\$ma"
  rm -f -- "\${out}.bed" "\${out}.bim" "\${out}.fam" "\${out}.log" "\${out}.err"
  run_tool "\$out" plink2 --bfile "\$source" --extract "\$ma" --make-bed --output-chr 26 --out "\$out"
  gwas_post_need_file "\${out}.bed"
  gwas_post_need_file "\${out}.bim"
  gwas_post_need_file "\${out}.fam"
  GCTA_BFILE="\$out"
  GCTA_CHR_ARGS=(--autosome-num "\$chr")
}

gwas_post_match_lead_inputs(){
  local ref="\$1" cref="\$2" tag="\$3" assoc="\$4" ma="\$5"
  local matched ref_subset n_assoc_before n_assoc_after n_ma_before n_ma_after
  declare -F match_SNP >/dev/null 2>&1 || gwas_clean_load_phef

  n_assoc_before=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  if (( n_assoc_before > 0 )); then
    ref_subset="\$GWAS_POST_TMP/\${tag}.clump.reference.pvar"
    gwas_post_subset_match_reference "\$ref" "\$assoc" "\$ref_subset" pvar
    matched="\${assoc}.matched.\$\$"
    match_SNP --reference "\$ref_subset" --output "\$matched" \
      --audit "\${QC_PREFIX}.\${tag}.clump_match.tsv" \
      --unmatched "\${QC_PREFIX}.\${tag}.clump_unmatched.tsv" "\$assoc"
    awk -v FS='\t' -v OFS='\t' '
      NR==1{next}
      !((\$1) in best) || (\$6+0)<best[\$1]{best[\$1]=\$6+0;line[\$1]=\$1 OFS \$6}
      END{print "SNP","P";for(snp in line)print line[snp]}' "\$matched" > "\${assoc}.tmp.\$\$"
    mv -f "\${assoc}.tmp.\$\$" "\$assoc"; rm -f "\$matched"
    n_assoc_after=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  else
    rm -f "\${QC_PREFIX}.\${tag}.clump_match.tsv" "\${QC_PREFIX}.\${tag}.clump_unmatched.tsv"
    n_assoc_after=0
  fi

  n_ma_before=\$(awk 'NR>1{n++} END{print n+0}' "\$ma")
  if (( n_ma_before > 0 )); then
    ref_subset="\$GWAS_POST_TMP/\${tag}.cojo.reference.bim"
    gwas_post_subset_match_reference "\$cref" "\$ma" "\$ref_subset" bim
    matched="\${ma}.matched.\$\$"
    match_SNP --reference "\$ref_subset" --output "\$matched" \
      --audit "\${QC_PREFIX}.\${tag}.cojo_match.tsv" \
      --unmatched "\${QC_PREFIX}.\${tag}.cojo_unmatched.tsv" "\$ma"
    awk -v FS='\t' -v OFS='\t' '
      NR==1{next}
      !((\$1) in best) || (\$9+0)<best[\$1]{best[\$1]=\$9+0;line[\$1]=\$1 OFS \$4 OFS \$5 OFS \$6 OFS \$7 OFS \$8 OFS \$9 OFS \$10}
      END{print "SNP","A1","A2","freq","b","se","p","N";for(snp in line)print line[snp]}' "\$matched" > "\${ma}.tmp.\$\$"
    mv -f "\${ma}.tmp.\$\$" "\$ma"; rm -f "\$matched"
    n_ma_after=\$(awk 'NR>1{n++} END{print n+0}' "\$ma")
  else
    rm -f "\${QC_PREFIX}.\${tag}.cojo_match.tsv" "\${QC_PREFIX}.\${tag}.cojo_unmatched.tsv"
    n_ma_after=0
  fi

  gwas_post_log "reference-ID match \$tag: clump=\$n_assoc_after/\$n_assoc_before cojo=\$n_ma_after/\$n_ma_before"
}

if run_lead; then
  if gwas_post_lead_done; then
    gwas_post_log "lead SNP discovery exists: \$AWK_SNP"
  else
    lead_t0=\$(date +%s)
    clump_was_done=FALSE
    cojo_was_done=FALSE
    cojo_skipped=FALSE
    clump_has_usable_input=FALSE
    cojo_has_matched_input=FALSE
    if [[ "\$REPLACE" != "TRUE" ]]; then
      gwas_post_clump_complete && clump_was_done=TRUE
      gwas_post_cojo_complete && cojo_was_done=TRUE
    fi
    gwas_post_need_file "\$FINAL"
    labels="\${GWAS_POST_TMP}/labels.\$\$.txt"
    refs="\${GWAS_POST_TMP}/refs.\$\$.tsv"
    : > "\$labels"
    : > "\$refs"

  # 1) Resolve selected reference chromosomes. GWAS coordinates, alleles, and
  # frequencies are never filled from the reference panel.
  while read -r ref lab chr unused; do
    [[ -n "\$ref" ]] || continue
    want_chr "\$chr" "\$CHRS" || continue
    tag=\$(chr_label "\$chr")
    echo "\$tag" >> "\$labels"
    cp="\${CLUMP}/\${tag}"
    jp="\${COJO}/\${tag}"
    assoc="\${cp}.assoc"
    ma="\${jp}.ma"
    cref=\$(cojo_bfile "\$REFGEN_COJO" "\$chr")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "\$ref" "\$tag" "\$chr" "\$assoc" "\$ma" "\$cref" >> "\$refs"
  done < <(ref_clump_pfiles "\$REFGEN_CLUMP")
  gwas_post_check_chromosome_coverage lead "\$refs"
  rm -f -- "\${MERGED}.lead.done"
  if [[ "\$clump_was_done" != "TRUE" ]]; then
    rm -f -- "\$CLUMP_DONE" "\${MERGED}.clumps"
  fi
  if [[ "\$cojo_was_done" != "TRUE" ]]; then
    rm -f -- "\$COJO_DONE" "\${MERGED}.jma.cojo" "\${MERGED}.cma.cojo" "\${MERGED}.ldr.cojo"
  fi
  rm -rf -- "\$CLUMP" "\$COJO"
  mkdir -p "\$CLUMP" "\$COJO"
  gwas_post_log "prepare P<=\$P_LEAD lead inputs from required SNP/CHR/POS/EA/NEA/BETA/SE/P fields: \$FINAL"
  gwas_post_prep_lead_inputs "\$refs"

  # 2) No-LD top SNP per chromosome/window and per-chromosome tool inputs.
  gwas_post_log "awk distance lead SNPs: \${GWAS_POST_LEAD_VIEW:-\$FINAL} -> \$AWK_SNP"
  gwas_post_zcat "\${GWAS_POST_LEAD_VIEW:-\$FINAL}" | awk -v FS='\t' -v OFS='\t' -v pthr="\$P_LEAD" -v win="\$LEAD_WINDOW" -f <(awk_lead) > "\${AWK_SNP}.tmp.\$\$"
  mv -f "\${AWK_SNP}.tmp.\$\$" "\$AWK_SNP"
  awk -v g="\$GWAS" 'NR==1{next} END{print g"\\t"NR-1}' "\$AWK_SNP" > "\${QC_PREFIX}.awk.nrow.tsv"

  # 3) plink2 --clump and gcta --cojo per chromosome.
  while IFS=\$'\t' read -r ref tag chr assoc ma cref; do
    [[ -n "\$ref" ]] || continue
    cp="\${CLUMP}/\${tag}"
    jp="\${COJO}/\${tag}"
    gwas_post_match_lead_inputs "\$ref" "\$cref" "\$tag" "\$assoc" "\$ma"
    # Inputs were already restricted to P<=P_LEAD before reference matching,
    # which bounds memory independently of the full GWAS size.
    cojo_skip_log="\${QC_PREFIX}.\${tag}.cojo_skip.log"
    awk -v FS='\t' -v OFS='\t' '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
      BEGIN{print "STATUS","REASON","SNP","A1","A2","EAF","P"}
      NR>1 && (!isnum(\$4)||\$4+0<0||\$4+0>1){print "SKIP_GCTA","missing_or_invalid_EAF_after_reference_match",\$1,\$2,\$3,\$4,\$7}
    ' "\$ma" > "\${cojo_skip_log}.tmp.\$\$"
    if gwas_post_has_data_rows "\${cojo_skip_log}.tmp.\$\$"; then
      mv -f "\${cojo_skip_log}.tmp.\$\$" "\$cojo_skip_log"
    else
      rm -f "\${cojo_skip_log}.tmp.\$\$" "\$cojo_skip_log"
    fi
    if [[ "\$clump_was_done" == "TRUE" ]]; then
      gwas_post_log "SKIP plink2 chr\$chr: clump phase already complete"
    elif gwas_post_has_data_rows "\$assoc"; then
      gwas_post_filter_clump_multiallelic "\$ref" "\$tag" "\$assoc"
      if gwas_post_has_data_rows "\$assoc"; then
        clump_has_usable_input=TRUE
        pfile_modifier=()
        [[ -s "\${ref}.pvar" || ! -s "\${ref}.pvar.zst" ]] || pfile_modifier=(vzs)
        pfile_keep=()
        [[ -z "\$REFGEN_KEEP" ]] || pfile_keep=(--keep "\$REFGEN_KEEP")
        run_tool "\$cp" plink2 --pfile "\$ref" "\${pfile_modifier[@]}" "\${pfile_keep[@]}" --clump "\$assoc" --clump-snp-field SNP --clump-p-field P --clump-p1 "\$P_LEAD" --clump-kb "\$CLUMP_KB" --out "\$cp"
      else
        gwas_post_log "skip plink2 chr\$chr: all significant IDs are ambiguous in the reference"
      fi
    else
      gwas_post_log "skip plink2 chr\$chr: no SNPs in \$assoc"
    fi

    # 3) gcta --cojo: COJO selection with the chromosome-specific bed/bim/fam reference
    if [[ "\$cojo_was_done" == "TRUE" ]]; then
      gwas_post_log "SKIP gcta chr\$chr: COJO phase already complete"
    elif [[ -s "\${QC_PREFIX}.\${tag}.cojo_skip.log" ]]; then
      cojo_skipped=TRUE
      gwas_post_log "SKIP gcta chr\$chr: significant variant has missing/invalid EAF; see \${QC_PREFIX}.\${tag}.cojo_skip.log"
    elif gwas_post_has_data_rows "\$ma"; then
      gwas_post_prepare_gcta_bfile "\$cref" "\$chr" "\$ma" "\${jp}.gcta_ref"
      # A matched SNP can be monomorphic in the reference population. GCTA
      # 1.94.1 can segfault on its zero genotype variance during COJO selection.
      # This floor excludes fixed alleles while retaining rare 1000G variants.
      gwas_post_run_cojo "\$jp" "\${QC_PREFIX}.\${tag}.cojo_empty" gcta --bfile "\$GCTA_BFILE" "\${GCTA_CHR_ARGS[@]}" \
        --maf 1e-6 --cojo-file "\$ma" --cojo-slct --cojo-p "\$P_LEAD" --out "\$jp"
      if [[ "\$GCTA_COJO_HAS_MATCH" == TRUE ]]; then cojo_has_matched_input=TRUE; fi
    else
      gwas_post_log "skip gcta chr\$chr: no SNPs in \$ma"
    fi
  done < "\$refs"

    if gwas_post_lead_awk_has_no_rows; then
      rm -f "\$labels" "\$refs"
      rm -rf -- "\$CLUMP" "\$COJO"
      gwas_post_mark_phase_done clump no_significant_variants
      gwas_post_mark_phase_done cojo no_significant_variants
    else
      if [[ "\$clump_was_done" != "TRUE" ]]; then
        concat_chr_outputs "\$CLUMP" "\$MERGED" "\$labels" clumps
        if [[ -s "\${MERGED}.clumps" ]]; then
          gwas_post_mark_phase_done clump
        elif [[ "\$clump_has_usable_input" != "TRUE" ]]; then
          gwas_post_log "clump complete with no reference-matched significant variants"
          gwas_post_mark_phase_done clump no_reference_matched_variants
        fi
      fi
      if [[ "\$cojo_was_done" != "TRUE" && "\$cojo_skipped" != "TRUE" ]]; then
        concat_chr_outputs "\$COJO" "\$MERGED" "\$labels" jma.cojo cma.cojo ldr.cojo
        if [[ -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]]; then
          gwas_post_mark_phase_done cojo
        elif [[ "\$cojo_has_matched_input" != "TRUE" ]]; then
          gwas_post_log "COJO complete with no reference-matched significant variants"
          gwas_post_mark_phase_done cojo no_reference_matched_variants
        fi
      fi
      rm -f "\$labels" "\$refs"
      if ! gwas_post_clump_complete; then
        echo "ERROR: clump phase did not create \${MERGED}.clumps; retaining \$CLUMP and \$COJO" >&2
        exit 1
      fi
      if ! gwas_post_cojo_complete; then
        echo "ERROR: COJO phase is incomplete for \$GWAS; GCTA skip details are under \${QC_PREFIX}.*.cojo_skip.log; retaining \$CLUMP and \$COJO" >&2
        exit 1
      fi
      gwas_post_log "delete lead intermediates: \$CLUMP \$COJO \${MERGED}.cma.cojo"
      rm -rf -- "\$CLUMP" "\$COJO"
      rm -f -- "\${MERGED}.cma.cojo"
    fi
    gwas_post_log "lead SNP discovery done in \$((\$(date +%s)-lead_t0)) sec"
  fi
fi

# PGS depends on the genome-wide .jma.cojo created by the lead/COJO phase, so
# it is deliberately invoked after lead even though it is declared after mplot.
gwas_post_pgs

if [[ ",\$H2_REQUESTED," == *,h2,* || "\$H2_REQUESTED" == all ]]; then
  h2_args=()
  [[ -z "\$H2_PYTHON" ]] || h2_args+=(--python "\$H2_PYTHON")
  python3 "\$H2_HELPER" --gwas-file "\$FINAL" --trait "\$GWAS" \
    --output-dir "\$(dirname "\$GWAS_DIR")/h2" --sex "\$H2_SEX" \
    --ref-ld-chr "\$H2_REF_LD" --w-ld-chr "\$H2_W_LD" --merge-alleles "\$H2_MERGE_ALLELES" \
    --conda-env "\$H2_CONDA_ENV" --replace "\$REPLACE" "\${h2_args[@]}"
fi

CMD_TOP

  chmod +x "$cmd"
  echo "$cmd"
}
