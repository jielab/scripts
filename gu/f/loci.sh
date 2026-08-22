#!/usr/bin/env bash
set -euo pipefail



usage() {
  cat <<'EOF'
Usage:
  ./gu.sh [MODE]

Default:
  ./gu.sh
      Run the full GU workflow for all configured traits.

Modes:
  all      Run the staged workflow. For avcf, this runs prep -> ld -> core -> mat -> hap -> phy, then extra QC.
  clean    Remove previous outputs for each configured trait, keeping gu.log when possible.
  prep     Step1: prepare lead SNPs and risk alleles from COJO/GWAS files.
  ld       Step2-Step3: match 1000G variants, calculate LD blocks, and choose loci.
  core     Step4: build core VCFs for selected loci.
  mat      Step5: convert core VCFs into matrix tables.
  hap      Step6: identify candidate inherited haplotypes and write report tables.
  phy      Step7: build PHYLIP files, run phyml/fallback trees, and draw tree plots.
  extra    Summarize existing avcf reports for manuscript QC/sensitivity analyses.
  proxy    Optional UKB haplotype-proxy association using diagnostic variants.

Help:
  ./gu.sh -h
  ./gu.sh --h
  ./gu.sh --help

Main examples:
  method=avcf ./gu.sh
  method=asnp ./gu.sh
  method=avcf traits="bald12" ./gu.sh
  START_STEP=s2 method=avcf traits="bald12" ./gu.sh
  method=avcf ./gu.sh prep
  method=avcf ./gu.sh ld
  method=avcf ./gu.sh core
  method=avcf ./gu.sh mat
  method=avcf ./gu.sh hap
  method=avcf ./gu.sh phy
  method=avcf ./gu.sh extra

Extra QC grid example:
  method=avcf \
    pils_grid="0.1,0.01,0.001,0.0001" \
    freq_grid="0.01,0.025,0.05,0.10" \
    diag_n_grid="1,3,5,10" \
    lineage_prop_grid="0.5,0.7,0.8,0.9" \
    hap_n_grid="5,10,20" \
    ./gu.sh extra

Proxy association example:
  method=avcf \
    pfile_tpl=/mnt/d/data.BIG/ukb/plink2/ukb_imp_chr{CHR} \
    pheno_cov=/mnt/d/analysis/bald/bald_pheno_cov.tsv \
    id_col=eid \
    outcomes="bald bald12 bald13 bald14" \
    binary_outcomes="bald12 bald13 bald14" \
    covars="age array center PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" \
    ./gu.sh proxy

Common variables:
  method=avcf|asnp
  traits="bald bald12 bald13 bald14"
  chrs="1 2 3 X"
  max_cores=16
  job_of_trait=1
  job_in_trait=1
  yri_freq_th=0.05
  diagnostic_delta_th=0.5
  strict_p_ils_th=0.1
  strict_prop_match_risk_th=0.5
  strict_hap_min_n=10

Input/output variables:
  dir0=/mnt/d
  dirscript=/mnt/d/scripts/gu
  dir_ref=/mnt/d/data.BIG/refGen
  dirgwas=/mnt/d/data.BIG/gwas/bald
  dircojo=/mnt/d/data.BIG/gwas/bald/cojo
  dirout_root=/mnt/d/analysis/gu/loci_avcf or /mnt/d/analysis/gu/loci_asnp

Notes:
  - Use environment variables for settings; do not pass settings as --options.
  - Unknown --options are rejected intentionally to keep the interface simple.
  - START_STEP=s2 reuses existing lead/pick.tsv when Step1-Step3 are already complete.
EOF
}

for arg in "$@"; do
  case "$arg" in
    -h|--h|--help)
      usage
      exit 0
      ;;
    --*)
      echo "Unknown option: $arg" >&2
      usage >&2
      exit 1
      ;;
  esac
done
if [[ $# -gt 1 ]]; then
  echo "ERROR: gu.sh accepts at most one MODE, but got: $*" >&2
  usage >&2
  exit 1
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Paths, traits, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir0=${dir0:-/mnt/d}
dirscript=${dirscript:-$dir0/scripts/gu/f}
# Self-contained chrX PAR coordinates; no private helper dependency.
X_PAR1_START_b37=${X_PAR1_START_b37:-60001}; X_PAR1_END_b37=${X_PAR1_END_b37:-2699520}
X_PAR2_START_b37=${X_PAR2_START_b37:-154931044}; X_PAR2_END_b37=${X_PAR2_END_b37:-155260560}
X_PAR1_START_b38=${X_PAR1_START_b38:-10001}; X_PAR1_END_b38=${X_PAR1_END_b38:-2781479}
X_PAR2_START_b38=${X_PAR2_START_b38:-155701383}; X_PAR2_END_b38=${X_PAR2_END_b38:-156030895}
phe_sh=GU_v2_internal_PAR_coordinates

method=${method:-avcf} # avcf | asnp
case "$method" in
	avcf|asnp) ;;
	*) echo "ERROR invalid method=$method; use method=avcf or method=asnp." >&2; exit 1 ;;
esac
dir_ref=${dir_ref:-$dir0/data.BIG/refGen}
[[ -d $dir_ref ]] || { echo "ERROR missing refGen directory: $dir_ref" >&2; exit 1; }
dirarch=${dirarch:-${dir_archaic:-$dir_ref/archaic/GRCH${GRCH:-37}}/$method}
dirmod=${dirmod:-${dir_1kg:-$dir_ref/1kg/GRCH${GRCH:-37}}}
sample_file=${sample_file:-$dirmod/vcf/samples_v3.ALL.panel}

# Shared genomic interval definitions come from 0phe.f.sh.
genome_build=${genome_build:-b37}
export genome_build
for v in X_PAR1_START_${genome_build} X_PAR1_END_${genome_build} X_PAR2_START_${genome_build} X_PAR2_END_${genome_build}; do [[ -n ${!v:-} ]] && export "$v"; done
export PATH="$dir0/software/bin:$PATH"
dirgwas=${dirgwas:-$dir0/data.BIG/gwas/bald}
dircojo=${dircojo:-$dirgwas/cojo}

dirout_root=${dirout_root:-$dir0/analysis/gu/$method}
dirout=
dir_report=
dir_plot=

positive_loci=${positive_loci:-${dir_archaic:-$dir_ref/archaic/GRCH${GRCH:-37}}/gu.loci.bed}
add_positive_loci=${add_positive_loci:-auto}
traits=${traits:-"bald bald12 bald13 bald14"}
chrs=${chrs:-"$(seq 1 22) X"}
ld_r2=${ld_r2:-0.98}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Method switches
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ref_pop=${ref_pop:-ALL} # 📍 ALL | EUR
lead_match=${lead_match:-avcf}   
ld_calc=${ld_calc:-avcf} 
archaic_VCF_match=${archaic_VCF_match:-avcf} 
pick_lineage=${pick_lineage:-avcf} 
hap_filter=${hap_filter:-avcf}  
yri_freq_th=${yri_freq_th:-0.05} # 📍 avcf hap_filter: keep loci with YRI_max_freq <= threshold; >=1 disables this filter
diagnostic_delta_th=${diagnostic_delta_th:-0.5} # 📍 diagnostic archaic allele must be enriched on risk haplotypes by this carry-noncarry frequency delta

strict_p_ils_th=${strict_p_ils_th:-0.1}
strict_min_compared_risk=${strict_min_compared_risk:-2}
strict_min_match_risk=${strict_min_match_risk:-2}
strict_prop_match_risk_th=${strict_prop_match_risk_th:-0.5}
strict_hap_min_n=${strict_hap_min_n:-10}
strict_require_carry_risk=${strict_require_carry_risk:-1}
strict_yri_freq_th=${strict_yri_freq_th:-$yri_freq_th}

# Extra manuscript-validation analyses. These do not rerun the GU haplotype screen;
# they summarize existing report/*.tsv files and optionally test diagnostic haplotype
# proxies in UKB.
extra_outdir=${extra_outdir:-$dirout_root/extra}
pils_grid=${pils_grid:-"0.1,0.01,0.001,0.0001"}
freq_grid=${freq_grid:-"0.01,0.025,0.05,0.10"}
diag_n_grid=${diag_n_grid:-"1,3,5,10"}
lineage_prop_grid=${lineage_prop_grid:-"0.5,0.7,0.8,0.9"}
hap_n_grid=${hap_n_grid:-"5,10,20"}

# Optional UKB haplotype-proxy association. This needs local UKB pfiles and a
# phenotype/covariate table. Run explicitly with: method=avcf ./gu.sh proxy
proxy_outdir=${proxy_outdir:-$extra_outdir/ukb_haplotype_proxy_assoc}
diag_tsv=${diag_tsv:-$extra_outdir/diagnostic_variants_for_haplotype_proxy.tsv}
pfile_tpl=${pfile_tpl:-$dir0/data.BIG/ukb/plink2/ukb_imp_chr{CHR}}
pheno_cov=${pheno_cov:-$dir0/analysis/bald/bald_pheno_cov.tsv}
id_col=${id_col:-eid}
outcomes=${outcomes:-"bald bald12 bald13 bald14"}
binary_outcomes=${binary_outcomes:-"bald12 bald13 bald14"}
covars=${covars:-"age array center PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10"}
plink2=${plink2:-plink2}

phy_input=${phy_input:-avcf} # 📍 only main trees are generated by default

# asnp does not read archaic VCFs; it uses the aSNP haplotype map in refGen/archaic/asnp.
if [[ $method == asnp || ${lead_match:-} == asnp || ${ld_calc:-} == asnp || ${archaic_VCF_match:-} == asnp || ${pick_lineage:-} == asnp || ${hap_filter:-} == asnp || ${phy_input:-} == asnp ]]; then
	asnp_asnp=${asnp_asnp:-${dir_archaic:-$dir_ref/archaic/GRCH${GRCH:-37}}/asnp/aSNPs.haplotypes.v1.tsv} # Yermakovich-style aSNP haplotype map
	asnp_p_th=${asnp_p_th:-1e-9}                                  # default threshold used in the 2026 GBE paper
	asnp_freq_th=${asnp_freq_th:-0.01}                            # archaic allele/haplotype frequency threshold
	asnp_min_asnp=${asnp_min_asnp:-5}                             # require haplotypes with >= this many aSNPs
	asnp_window_kb=${asnp_window_kb:-1000}                        # +/- window around lead SNP for regional scan
	asnp_ld_r2=${asnp_ld_r2:-0.9}                                 # high-LD threshold for network candidate region
	asnp_make_network=${asnp_make_network:-1}                     # 1: try pegas network if 1000G VCF + packages are available
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Parallelism and runtime limits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
max_cores=${max_cores:-16}
job_of_trait=${job_of_trait:-1}
job_in_trait=${job_in_trait:-1}
job_phyml=${job_phyml:-$job_in_trait}
plink_threads=${plink_threads:-1}
internal_threads=${internal_threads:-1}
phyml_cpus=${phyml_cpus:-1}
phyml_retry_cpus=${phyml_retry_cpus:-1}
phyml_boot=${phyml_boot:-0}        # 0 is fast; set phyml_boot=100 for final bootstrap support values
phyml_timeout=${phyml_timeout:-12h}
phyml_full_required=${phyml_full_required:-0}
phyml_run_full=${phyml_run_full:-0}
phyml_nj_fallback=${phyml_nj_fallback:-1}
pvar_scan_timeout=${pvar_scan_timeout:-10m}
plink_timeout=${plink_timeout:-6h}
require_slim_pvar=${require_slim_pvar:-1}
pvar_lookup_debug=${pvar_lookup_debug:-0}
filter_pop=${filter_pop:-YRI}
filter_max_count=${filter_max_count:-1}
start_step=${START_STEP:-${start_step:-s1}}

configured_traits=$traits
read -r -a configured_traits_arr <<< "$configured_traits"
positive_trait=${positive_trait:-${configured_traits_arr[0]}}
traits_arr=()
exec 3>&1
exec 4>&2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Logging and small utilities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
count_lines(){ awk 'NF{n++} END{print n+0}' "$1"; }
data_rows(){
	local file=$1
	[[ -s $file ]] || { echo 0; return; }
	if [[ $file == *.gz ]]; then gzip -dc "$file"; else cat "$file"; fi | awk 'NR>1{n++} END{print n+0}'
}
nrow(){ data_rows "$1"; }
ok(){ [[ -s $1 && $(nrow "$1") -gt 0 ]]; }
chrn(){ [[ $1 == X ]] && echo 23 || echo "$1"; }
chrl(){ [[ $1 == 23 || $1 == X ]] && echo X || echo "$1"; }
clean_msg(){ [[ -f $1 ]] && tail -n 12 "$1" | tr '\t\r\n' '   ' | sed 's/  */ /g' | cut -c1-700; }
current_trait_is_positive(){ [[ ${traits_arr[0]:-} == "$positive_trait" ]]; }
run_timeout_cmd(){
	local limit=$1
	shift
	if [[ -n $limit && $limit != 0 && $limit != none ]] && command -v timeout >/dev/null 2>&1; then
		timeout "$limit" "$@"
	else
		"$@"
	fi
}
pvar_awk(){
	local desc=$1 rc
	shift
	if [[ -n ${pvar_scan_timeout:-} && ${pvar_scan_timeout:-} != 0 && ${pvar_scan_timeout:-} != none ]] && command -v timeout >/dev/null 2>&1; then
		timeout "$pvar_scan_timeout" awk "$@"
	else
		awk "$@"
	fi
	rc=$?
	if (( rc != 0 )); then
		log "ERROR pvar scan failed/timed out rc=$rc timeout=${pvar_scan_timeout:-none}: $desc" >&2
	fi
	return "$rc"
}
check_pvar_format(){
	local pvar=$1 context=${2:-pvar} header nf
	[[ ${require_slim_pvar:-1} == 1 ]] || return 0
	header=$(pvar_awk "$context header $pvar" '/^#CHROM/{print; exit}' "$pvar") || return 1
	if [[ -z $header ]]; then
		log "ERROR malformed pvar header for $context: $pvar"
		return 1
	fi
	nf=$(awk 'BEGIN{FS="\t"} {print NF; exit}' <<< "$header")
	if (( nf > 5 )); then
		log "ERROR fat pvar detected for $context: $pvar has $nf columns. Regenerate a 5-column pvar (#CHROM POS ID REF ALT) or set require_slim_pvar=0."
		return 1
	fi
}

set_trait_context(){
	local t=$1
	traits=$t
	traits_arr=("$t")
	dirout="$dirout_root/$t"
	dir_report="$dirout/report"
	dir_plot="$dirout/plot"
	export GU_REPORT_DIRNAME="$(basename "$dir_report")"
	export GU_PLOT_DIRNAME="$(basename "$dir_plot")"
	mkdir -p "$dirout"
}

open_trait_log(){
	# Reset stdout/stderr to the original descriptors before opening a new per-trait tee.
	exec >&3 2>&4
	exec > >(tee "$dirout/gu.log" >&3) 2>&1
}

export OMP_NUM_THREADS=$internal_threads
export OPENBLAS_NUM_THREADS=$internal_threads
export MKL_NUM_THREADS=$internal_threads
export VECLIB_MAXIMUM_THREADS=$internal_threads
export NUMEXPR_NUM_THREADS=$internal_threads
export R_DATATABLE_NUM_THREADS=$internal_threads

is_asnp(){
	[[ "$method" == "asnp" || " $lead_match $ld_calc $archaic_VCF_match $pick_lineage $hap_filter $phy_input " == *" asnp "* ]]
}

validate_methods(){
	local x v ok=1
	case "$method" in avcf|asnp) ;; *) log "ERROR invalid method=$method; use method=avcf or method=asnp."; ok=0 ;; esac
	for x in lead_match ld_calc archaic_VCF_match pick_lineage hap_filter phy_input; do
		eval "v=\${$x}"
		case "$v" in avcf|asnp) ;; *) log "ERROR invalid $x=$v; use avcf or asnp."; ok=0 ;; esac
	done
	(( ok == 1 )) || exit 1
}

validate_parallelism(){
	local trait_need=$(( job_of_trait * job_in_trait ))
	local phyml_need=$(( job_of_trait * job_phyml * phyml_cpus ))
	if (( trait_need > max_cores )); then
		log "ERROR trait parallelism exceeds max_cores: job_of_trait=$job_of_trait * job_in_trait=$job_in_trait > max_cores=$max_cores"
		exit 1
	fi
	if (( phyml_need > max_cores )); then
		log "ERROR phyml parallelism exceeds max_cores: job_of_trait=$job_of_trait * job_phyml=$job_phyml * phyml_cpus=$phyml_cpus > max_cores=$max_cores"
		exit 1
	fi
	if (( job_in_trait > max_cores )); then
		log "ERROR job_in_trait=$job_in_trait exceeds max_cores=$max_cores"
		exit 1
	fi
}

validate_asnp_inputs(){
	is_asnp || return 0
	[[ -s "$dirscript/loci.R" ]] || { log "ERROR asnp requires $dirscript/loci.R"; exit 1; }
	if [[ ! -s $asnp_asnp ]]; then
		log "ERROR asnp aSNP haplotype map not found: $asnp_asnp"
		log "Put aSNPs.haplotypes.v1.tsv under ${dir_archaic:-$dir_ref/archaic/GRCH${GRCH:-37}}/asnp, or export asnp_asnp before running."
		exit 1
	fi
}

validate_avcf_inputs(){
	is_asnp && return 0
	[[ -d $dirarch ]] || { log "ERROR avcf archaic VCF directory not found: $dirarch"; exit 1; }
	local any_vcf
	any_vcf=$(find "$dirarch" -mindepth 2 -maxdepth 2 -type f -name '*.vcf.gz' -print -quit 2>/dev/null || true)
	[[ -n $any_vcf ]] || { log "ERROR no archaic VCF found under $dirarch"; exit 1; }
}

log_config(){
	validate_methods
	validate_parallelism
	log "Now running GU workflow; traits=$traits; output=$dirout; max_cores=$max_cores"
	log "Methods: method=$method, lead_match=$lead_match, ld_calc=$ld_calc, archaic_VCF_match=$archaic_VCF_match, pick_lineage=$pick_lineage, hap_filter=$hap_filter, yri_freq_th=$yri_freq_th, diagnostic_delta_th=$diagnostic_delta_th, phy_input=$phy_input"
	log "Strict-style report columns: p_ils<$strict_p_ils_th, n_compared_risk>=$strict_min_compared_risk, n_match_risk>=$strict_min_match_risk, prop_match_risk>=$strict_prop_match_risk_th, hap_n>$strict_hap_min_n, require_carry_risk=$strict_require_carry_risk, strict_yri_freq_th=$strict_yri_freq_th"
	if is_asnp; then
		log "asnp: asnp=$asnp_asnp p_th=$asnp_p_th freq_th=$asnp_freq_th min_asnp=$asnp_min_asnp window_kb=$asnp_window_kb ld_r2=$asnp_ld_r2 make_network=$asnp_make_network"
	fi
	log "Inputs: GWAS=$dirgwas; COJO=$dircojo; 1000G=$dirmod; pfile=$dirmod/pfile; vcf=$dirmod/vcf; archaic=$dirarch; sanity_loci=$positive_loci"
	log "Controls: positive_trait=$positive_trait; genome_build=$genome_build; phe_sh=${phe_sh:-NA}"
	validate_avcf_inputs
	validate_asnp_inputs
}

init_output_dirs(){
	mkdir -p "$dirout"/{lead,ld,coreVcf,mat,hap,phy} "$dir_plot" "$dir_report"
}

single_trait_root(){
	[[ ${#traits_arr[@]} -eq 1 && $(basename "$dirout") == "$1" ]]
}
stage_trait_dir(){
	local stage=$1 t=$2
	if single_trait_root "$t"; then
		printf "%s/%s\n" "$dirout" "$stage"
	else
		printf "%s/%s/%s\n" "$dirout" "$stage" "$t"
	fi
}
lead_assoc_file(){
	local t=$1
	if single_trait_root "$t"; then
		printf "%s/lead/lead.assoc\n" "$dirout"
	else
		printf "%s/lead/%s.lead.assoc\n" "$dirout" "$t"
	fi
}
merge_lead_assoc_files(){
	local first=${traits%% *} tmp
	tmp="$dirout/lead/.lead.assoc.$$"
	{ head -n1 "$(lead_assoc_file "$first")"; for t in $traits; do awk 'NR>1' "$(lead_assoc_file "$t")"; done; } > "$tmp" && mv "$tmp" "$dirout/lead/lead.assoc"
}
lead_3col_file(){
	local t=$1
	if single_trait_root "$t"; then
		printf "%s/lead/lead.3col\n" "$dirout"
	else
		printf "%s/lead/%s.lead.3col\n" "$dirout" "$t"
	fi
}
lead_id_sanity_file(){
	local t=$1
	if single_trait_root "$t"; then
		printf "%s/lead/id_sanity.tsv\n" "$dirout"
	else
		printf "%s/lead/%s.id_sanity.tsv\n" "$dirout" "$t"
	fi
}
lead_pick_file(){
	local t=$1 kind=${2:-pick}
	if single_trait_root "$t"; then
		if [[ $kind == single ]]; then
			printf "%s/lead/pick.single.tsv\n" "$dirout"
		else
			printf "%s/lead/pick.tsv\n" "$dirout"
		fi
	else
		if [[ $kind == single ]]; then
			printf "%s/lead/%s.pick.single.tsv\n" "$dirout" "$t"
		else
			printf "%s/lead/%s.pick.tsv\n" "$dirout" "$t"
		fi
	fi
}

id_mode(){ awk 'BEGIN{IGNORECASE=1} $1!="" && $1!="."{n++; if($1~/^rs[0-9]+$/) rs++; else if(toupper($1)~/^(CHR)?([0-9]+|X|Y|MT|M):[0-9]+:[ACGTN]+:[ACGTN,]+$/) cp++; if(n==5000) exit} END{if(n==0) print "missing"; else if(rs/n>.8) print "rsid"; else if(cp/n>.8) print "chrpos"; else print "other"}'; }
trait_rows(){ [[ -s $1 ]] && awk -v t="$2" 'BEGIN{FS="\t"} NR>1 && $1==t{n++} END{print n+0}' "$1" || echo 0; }
count_files(){ find "$1" -name "$2" 2>/dev/null | wc -l; }
valid_hap_locus(){
	local t=$1 id=$2 d
	d=$(stage_trait_dir hap "$t")
	[[ -s $d/$id.region.tsv && -s $d/$id.core.tsv ]] || [[ -s $d/$id.done ]] || return 1
	if [[ $hap_filter == avcf ]]; then
		awk 'NR==1{for(i=1;i<=NF;i++) if($i=="is_diagnostic_archaic") ok=1; exit !ok}' "$d/$id.core.tsv" || return 1
		awk 'NR==1{y=0; for(i=1;i<=NF;i++) if($i=="yri_risk_freq") y=i; if(!y) exit 1; next} y && $y!="" && $y!="NA"{ok=1} END{exit !ok}' "$d/$id.core.tsv"
	fi
}
valid_hap_stage(){
	local tr chr snp bp st en _n _size_bp id
	[[ -s $dir_report/hap_match.tsv && -s $dir_report/region_summary.tsv ]] || return 1
	while IFS=$'\t' read -r tr chr snp bp st en _n _size_bp; do
		id="${chr}.${snp}.${bp}"
		valid_hap_locus "$tr" "$id" || return 1
	done < <(awk 'BEGIN{FS=OFS="\t"} NR>1 && !seen[$1 FS $2 FS $3 FS $4]++{print}' "$dirout/lead/pick.tsv")
	if [[ $hap_filter == avcf ]]; then
		[[ -s $dir_report/inherited_region.tsv ]] || return 1
		[[ -s $dir_report/hap_site_count.tsv && -s $dir_report/core_archaic_match.tsv && -s $dir_report/core_risk.tsv ]] || return 1
		[[ -s $dir_report/yri_filter.tsv ]] || return 1
		awk -v th="$yri_freq_th" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="yri_freq_th") c=i; next} c && sprintf("%.12g",$c+0)==sprintf("%.12g",th+0){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk -v th="$diagnostic_delta_th" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="diagnostic_delta_th") c=i; next} c && sprintf("%.12g",$c+0)==sprintf("%.12g",th+0){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk -v th="$strict_p_ils_th" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="strict_p_ils_th") c=i; next} c && sprintf("%.12g",$c+0)==sprintf("%.12g",th+0){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk -v th="$strict_min_compared_risk" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="strict_min_compared_risk") c=i; next} c && int($c)==int(th){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk -v th="$strict_min_match_risk" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="strict_min_match_risk") c=i; next} c && int($c)==int(th){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk -v th="$strict_prop_match_risk_th" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="strict_prop_match_risk_th") c=i; next} c && sprintf("%.12g",$c+0)==sprintf("%.12g",th+0){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk -v th="$strict_hap_min_n" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="strict_hap_min_n") c=i; next} c && int($c)==int(th){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk -v th="$strict_yri_freq_th" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="strict_yri_freq_th") c=i; next} c && sprintf("%.12g",$c+0)==sprintf("%.12g",th+0){ok=1} END{exit !ok}' "$dir_report/yri_filter.tsv" || return 1
		awk 'NR==1{for(i=1;i<=NF;i++) if($i=="strict_filter_pass") ok=1; exit !ok}' "$dir_report/inherited_region.tsv" || return 1
		awk 'NR==1{for(i=1;i<=NF;i++) if($i=="is_diagnostic_archaic") ok=1; exit !ok}' "$dir_report/core_risk.tsv" || return 1
		awk 'NR==1{y=0; for(i=1;i<=NF;i++) if($i=="yri_risk_freq") y=i; if(!y) exit 1; next} y && $y!="" && $y!="NA"{ok=1} END{exit !ok}' "$dir_report/core_risk.tsv" || return 1
		awk 'NR==1{for(i=1;i<=NF;i++){if($i=="trait") t=1; if($i=="id") id=1; if($i=="lead_snp") snp=1; if($i=="yri_freq") y=1} exit !(t&&id&&snp&&y)}' "$dir_report/inherited_region.tsv" || return 1
	else
		[[ -s $dir_report/candidate_inherited_region.tsv ]] || return 1
		[[ -s $dir_report/filtered_hap_match.tsv && -s $dir_report/inherited_region.tsv ]] || return 1
		[[ -s $dir_report/all.xlsx && -s $dir_report/filtered.xlsx && -s $dir_report/selected.xlsx ]] || return 1
	fi
}
valid_phy_file(){
	local f=$1
	[[ -s $f ]] || return 1
	awk 'NR==1{next} NF && substr($0,11,1)!=" "{bad=1; exit} END{exit bad}' "$f"
}
valid_phy_locus(){
	local t=$1 id=$2 d f
	d=$(stage_trait_dir phy "$t")
	if [[ -e $d/$id.main.phy || -e $d/$id.full.phy ]]; then
		for f in "$d/$id.full.phy" "$d/$id.main.phy"; do
			[[ -s $f && -s ${f%.phy}.meta.tsv ]] || return 1
			valid_phy_file "$f" || return 1
		done
	else
		[[ -s $d/$id.done ]]
	fi
}
valid_phy_stage(){
	local miss invalid n_phy n_png
	n_phy=$(find "$dirout/phy" -name '*.full.phy' 2>/dev/null | wc -l)
	(( n_phy > 0 )) || return 1
	invalid=$(find "$dirout/phy" \( -name '*.full.phy' -o -name '*.main.phy' \) 2>/dev/null | while read -r f; do valid_phy_file "$f" || echo "$f"; done | head -1)
	[[ -z $invalid ]] || return 1
	miss=$(find "$dirout/phy" \( -name '*.full.phy' -o -name '*.main.phy' \) 2>/dev/null | while read -r f; do [[ -s "${f}_phyml_tree.txt" ]] || echo "$f"; done | head -1)
	[[ -z $miss ]] || return 1
	n_png=$(find "$dir_plot" -maxdepth 1 -name 's8_tree_full_*.png' 2>/dev/null | wc -l)
	(( n_png > 0 ))
}
log_trait_count(){
	local label=$1 f=$2 t
	log "summary: $label (all traits; file path $f) $(data_rows "$f")"
	for t in $traits; do log "summary: $label (trait $t; file path $f) $(trait_rows "$f" "$t")"; done
}
cojo_file_for_trait(){
	local t=$1 f
	for f in "$dircojo/$t.jma.cojo" "$dircojo/$t.cma.cojo" "$dircojo/$t.cojo" "$dirgwas/cojo/$t/$t.jma.cojo" "$dirgwas/cojo/$t/$t.cma.cojo"; do
		[[ -s $f ]] && { printf "%s\n" "$f"; return 0; }
	done
	return 1
}
log_process_s1(){
	local t f ncojo nlead npos
	log "process: Step1 input GWAS_dir=$dirgwas COJO_dir=$dircojo"
	for t in $traits; do
		if f=$(cojo_file_for_trait "$t"); then ncojo=$(data_rows "$f"); else f="NA"; ncojo=0; fi
		nlead=$(data_rows "$(lead_assoc_file "$t")")
		log "process: Step1 trait=$t COJO_file=$f COJO_loci=$ncojo lead_assoc_loci=$nlead output=$(lead_assoc_file "$t")"
	done
	nlead=$(data_rows "$dirout/lead/lead.assoc")
	npos=0; [[ -s $positive_loci ]] && npos=$(awk 'NF>=4 && $1 !~ /^#/{n++} END{print n+0}' "$positive_loci")
	log "process: Step1 merged lead loci=$nlead output=$dirout/lead/lead.assoc"
	log "process: Step1 gu.loci.bed input_snp_rows=$npos file=$positive_loci"
	log "FLOW: step=1 name=prepare_lead_loci cojo_loci=$ncojo lead_loci=$nlead positive_control_input=$npos output=$dirout/lead/lead.assoc"
}
log_process_s3(){
	local nlead nfail nok npick nsingle npos
	nlead=$(data_rows "$dirout/lead/lead.assoc")
	nfail=$(data_rows "$dirout/lead/lead_1000G.fail.tsv")
	nok=$(awk 'BEGIN{FS="\t"} NR>1 && $9 ~ /^ok/{k=$1 FS $2 FS $3 FS $4; ok[k]=1} END{print length(ok)+0}' "$dirout/lead/ld_debug.tsv" 2>/dev/null || echo 0)
	npick=$(data_rows "$dirout/lead/pick.tsv")
	nsingle=$(data_rows "$dirout/lead/pick.single.tsv")
	npos=$(data_rows "$dirout/lead/positive_pick.tsv")
	log "process: Step2 1000G matching input_lead_loci=$nlead matched_loci_with_LD=$nok failed_or_unmatched=$nfail debug=$dirout/lead/ld_debug.tsv fail=$dirout/lead/lead_1000G.fail.tsv"
	log "process: Step3 LD block split multi_snp_loci=$npick single_snp_loci=$nsingle pick=$dirout/lead/pick.tsv single=$dirout/lead/pick.single.tsv"
	log "process: Step3 gu.loci.bed loci_written=$npos output=$dirout/lead/positive_pick.tsv"
	log "FLOW: step=2 name=match_1000G_and_compute_LD input_lead_loci=$nlead matched_loci_with_LD=$nok failed_or_unmatched=$nfail"
	log "FLOW: step=3 name=split_LD_blocks multi_snp_loci=$npick single_snp_loci=$nsingle positive_control_loci=$npos"
}
log_process_s4_s5(){
	local npick ncore nmat
	npick=$(data_rows "$dirout/lead/pick.tsv")
	ncore=$([[ -d "$dirout/coreVcf" ]] && find "$dirout/coreVcf" -mindepth 1 -maxdepth 1 -type d | wc -l || echo 0)
	nmat=$([[ -d "$dirout/mat" ]] && find "$dirout/mat" -mindepth 1 -maxdepth 1 -type d | wc -l || echo 0)
	log "process: Step4 core VCF input_loci=$npick core_locus_dirs=$ncore output=$dirout/coreVcf"
	log "process: Step5 matrix input_loci=$npick matrix_locus_dirs=$nmat output=$dirout/mat"
	log "FLOW: step=4 name=build_core_VCF input_loci=$npick core_locus_dirs=$ncore"
	log "FLOW: step=5 name=build_genotype_matrices input_loci=$npick matrix_locus_dirs=$nmat"
}
log_process_s6(){
	local nhap nregion nselected nphy nmap nfate ymethod yth dth ndiag
	nhap=$(data_rows "$dir_report/hap_match.tsv")
	nregion=$(data_rows "$dir_report/region_summary.tsv")
	nselected=$(data_rows "$dir_report/inherited_region.tsv")
	nphy=$(data_rows "$dir_report/phy_region.tsv")
	nmap=$(data_rows "$dir_report/haplotype_sample_map.tsv")
	nfate=$(data_rows "$dir_report/positive_loci_fate.tsv")
	ymethod=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="yri_filter_method") c=i; next} c{print $c; exit}' "$dir_report/yri_filter.tsv" 2>/dev/null || true)
	yth=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="yri_freq_th") c=i; next} c{print $c; exit}' "$dir_report/yri_filter.tsv" 2>/dev/null || true)
	dth=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="diagnostic_delta_th") c=i; next} c{print $c; exit}' "$dir_report/yri_filter.tsv" 2>/dev/null || true)
	ndiag=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="n_yri_filter_sites") c=i; next} c && NR>1{s+=$c} END{print s+0}' "$dir_report/region_summary.tsv" 2>/dev/null || echo 0)
	log "process: Step6 haplotypes hap_rows=$nhap inherited_candidate_loci=$nregion inherited_loci_after_filter=$nselected hap_match=$dir_report/hap_match.tsv inherited=$dir_report/inherited_region.tsv"
	log "process: Step6 YRI diagnostic filter method=${ymethod:-NA} yri_freq_th=${yth:-NA} diagnostic_delta_th=${dth:-NA} diagnostic_sites=$ndiag summary=$dir_report/yri_filter.tsv"
	log "process: Step6 Roman hap labels mapped_to_1KG_copies=$nmap output=$dir_report/haplotype_sample_map.tsv"
	log "process: Step6 gu.loci.bed fate_rows=$nfate output=$dir_report/positive_loci_fate.tsv"
	log "process: Step7 phy input_loci=$nphy input=$dir_report/phy_region.tsv"
	log "FLOW: step=6 name=identify_archaic_haplotypes hap_rows=$nhap matched_haplotype_loci=$nregion selected_after_YRI_filter=$nselected yri_filter_method=${ymethod:-NA} yri_freq_th=${yth:-NA} diagnostic_sites=$ndiag"
}
log_process_s7(){
	local nphy nmain npng
	nphy=$(data_rows "$dir_report/phy_region.tsv")
	nmain=$(find "$dirout/phy" -name '*.main.phy' 2>/dev/null | wc -l)
	npng=$(find "$dir_plot" -maxdepth 1 -name 's8_tree_main_*.png' 2>/dev/null | wc -l)
	log "process: Step7 phy files input_loci=$nphy main_phy_files=$nmain phy_dir=$dirout/phy"
	log "process: Step8 tree plots plotted_loci=$npng plot_dir=$dir_plot"
	log "FLOW: step=7 name=build_phylogeny_input phy_loci=$nphy main_phy_files=$nmain"
	log "FLOW: step=8 name=render_tree_plots plotted_loci=$npng plot_dir=$dir_plot"
}
write_hap_sample_map(){
	local out="$dir_report/haplotype_sample_map.tsv"
	local hap="$dir_report/hap_match.tsv"
	[[ -s $hap ]] || { log "WARN haplotype sample map skipped: missing $hap"; return 0; }
	Rscript "$dirscript/loci.R" hap_sample_map "$dirout" "$sample_file" "$out"
	log "summary: haplotype Roman-label sample map written: $out ($(data_rows "$out") rows)"
}
log_positive_loci_fate(){
	local bed="$positive_loci" out="$dir_report/positive_loci_fate.tsv"
	[[ -s $bed ]] || return 0
	Rscript "$dirscript/loci.R" positive_loci_fate "$dirout" "$bed" "$out"
	log "summary: sanity loci fate table written: $out"
	if [[ -s $out ]]; then
		awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print "SANITY_LOCUS_FATE:",$5,$1,$6}' "$out" | while IFS= read -r line; do log "$line"; done
	fi
}
clean_report_workbooks_only(){
	find "$dir_report" -mindepth 1 -type d -exec rm -rf {} +
}
valid_header(){ [[ -s $1 ]] && awk 'NR==1 && NF>1{ok=1} END{exit !ok}' "$1"; }
valid_chr_ld(){
	local d=$1 nlead=$2
	[[ -s $d/ld.tsv && -s $d/block.tsv ]] || return 1
	valid_header "$d/ld.tsv" && valid_header "$d/block.tsv" || return 1
	if (( nlead == 0 )); then valid_header "$d/ld.tsv" && valid_header "$d/block.tsv"; else [[ $(data_rows "$d/block.tsv") -gt 0 ]]; fi
}
valid_mat_locus(){
	local o=$1 expected_arch=$2 n_arch
	[[ -s $o/kg.tsv && $(count_lines "$o/kg.tsv") -gt 0 ]] || return 1
	[[ -s $o/kg.samples.tsv ]] || return 1
	n_arch=$(find "$o" -maxdepth 1 -type f -name '*.tsv' ! -name 'kg.tsv' ! -name 'kg.samples.tsv' -size +0c 2>/dev/null | wc -l)
	(( n_arch >= expected_arch ))
}
check_pbase(){
	local pbase=$1 assoc=$2 label=$3 flag=$4 m1 m2 e1 e2 ns nv
	require_pfile "$pbase" "$label"
	[[ -f $flag ]] && return
	m1=$(pvar_awk "$label id mode $pbase.pvar" '!/^#/ && $3!=""{print $3; if(++n==5000) exit}' "$pbase.pvar" | id_mode)
	m2=$(awk 'NR>1{print $3}' "$assoc" | id_mode)
	e1=$(pvar_awk "$label id examples $pbase.pvar" '!/^#/ && $3!="." && $3!=""{print $3; if(++n==5) exit}' "$pbase.pvar" | paste -sd, -)
	e2=$(awk 'NR>1{print $3; if(++n==5) exit}' "$assoc" | paste -sd, -)
	ns=$(awk 'NR>1{n++} END{print n+0}' "$pbase.psam")
	nv=NA
	printf "label\tlead_mode\tlead_examples\tG1000_mode\tG1000_examples\tn_samples\tn_variants\tpbase\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$label" "$m2" "$e2" "$m1" "$e1" "$ns" "$nv" "$pbase" > "$flag"
	[[ $m1 == rsid || $m1 == chrpos ]] || log "WARN 1000G pvar ID mode is $m1 for $pbase.pvar; will rely on bp/allele matching and plink2 --set-missing-var-ids"
}
write_fail(){ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" >> "$dirout/lead/lead_1000G.fail.tsv"; }
write_debug(){ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" >> "$dirout/lead/ld_debug.tsv"; }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Command checks and phyml helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
check_cmd(){
	for x in Rscript bcftools bgzip tabix plink2 awk sed sort; do command -v "$x" >/dev/null || { log "ERROR missing command: $x"; exit 1; }; done
	[[ -s "$dirscript/loci.R" ]] || { log "ERROR missing $dirscript/loci.R"; exit 1; }
	grep -q "add_positive_loci" "$dirscript/loci.R" || { log "ERROR $dirscript/loci.R is an old version without add_positive_loci; copy the latest gu.R to $dirscript/ first."; exit 1; }
}
check_phyml(){ command -v phyml >/dev/null || { log "ERROR missing phyml"; exit 1; }; }
run_phyml_cmd(){
	local cpus=$1 f=$2 leave_duplicates=$3 logf=$4 boot=${5:-$phyml_boot} model=${6:-HKY85}
	local cats=${7:-4}
	local args=(-i "$f" -m "$model" -c "$cats" -b "$boot")
	if [[ $model == HKY85 ]]; then
		args+=(-a e -v e)
	fi
	[[ $leave_duplicates == 1 ]] && args+=(--leave_duplicates)
	if [[ $leave_duplicates == 1 ]]; then
		if command -v timeout >/dev/null 2>&1; then
			PHYMLCPUS="$cpus" timeout "$phyml_timeout" phyml "${args[@]}" > "$logf" 2>&1
		else
			PHYMLCPUS="$cpus" phyml "${args[@]}" > "$logf" 2>&1
		fi
	else
		if command -v timeout >/dev/null 2>&1; then
			PHYMLCPUS="$cpus" timeout "$phyml_timeout" phyml "${args[@]}" > "$logf" 2>&1
		else
			PHYMLCPUS="$cpus" phyml "${args[@]}" > "$logf" 2>&1
		fi
	fi
}
run_phyml_nj_fallback(){
	local f=$1 tree=$2 stats=$3 logf=$4
	[[ $phyml_nj_fallback == 1 ]] || return 1
	Rscript - "$f" "$tree" "$stats" > "$logf" 2>&1 <<'RS'
args <- commandArgs(TRUE)
phy <- args[1]
treef <- args[2]
statsf <- args[3]
suppressPackageStartupMessages(library(ape))
x <- read.dna(phy, format = "sequential")
if (nrow(x) < 3L) stop("need at least 3 sequences for fallback tree: ", phy)
d <- dist.dna(x, model = "JC69", pairwise.deletion = TRUE, as.matrix = FALSE)
tr <- tryCatch(nj(d), error = function(e) {
  message("nj failed; retry njs: ", conditionMessage(e))
  njs(d)
})
write.tree(tr, file = treef)
writeLines(c(
  paste("fallback_tree_method", "ape::nj/njs JC69 distance", sep = "\t"),
  paste("source_phy", phy, sep = "\t"),
  paste("n_sequences", nrow(x), sep = "\t"),
  paste("n_sites", ncol(x), sep = "\t")
), statsf)
RS
}
run_phyml_file(){
	local f=$1 leave_duplicates=0 rc=0
	local tree="${f}_phyml_tree.txt" stats="${f}_phyml_stats.txt" logf="${f}.phyml.run.log"
	local boot_stats="${f}_phyml_boot_stats.txt" boot_trees="${f}_phyml_boot_trees.txt" last_log
	[[ $# -ge 2 ]] && leave_duplicates=$2
	[[ -s $f ]] || { log "ERROR missing PHYLIP input: $f"; return 1; }
	if [[ $phy_input == strict ]]; then
		valid_phy_file "$f" || { log "ERROR invalid PHYLIP spacing/label format: $f; remove or regenerate this .phy"; return 1; }
	fi
	if [[ -s $tree ]]; then log "phyml skip existing tree: $f"; return 0; fi
	if [[ -e $tree && ! -s $tree ]]; then log "WARN remove empty phyml tree before rerun: $tree"; rm -f "$tree"; fi
	if [[ -e $stats && ! -s $stats ]]; then log "WARN remove empty phyml stats before rerun: $stats"; rm -f "$stats"; fi
	log "phyml start: $f cpus=$phyml_cpus boot=$phyml_boot timeout=$phyml_timeout"
	last_log="$logf"
	run_phyml_cmd "$phyml_cpus" "$f" "$leave_duplicates" "$logf" "$phyml_boot" "HKY85" || rc=$?
	if [[ $rc -ne 0 && $phyml_retry_cpus -lt $phyml_cpus ]]; then
		log "phyml optimized run unavailable rc=$rc at cpus=$phyml_cpus; retry cpus=$phyml_retry_cpus: $f"
		rm -f "$tree" "$stats" "$boot_stats" "$boot_trees"
		rc=0
		last_log="${logf}.retry${phyml_retry_cpus}"
		run_phyml_cmd "$phyml_retry_cpus" "$f" "$leave_duplicates" "$last_log" "$phyml_boot" "HKY85" || rc=$?
	fi
	if [[ $rc -ne 0 || ! -s $tree ]]; then
		log "WARN phyml HKY85 bootstrap failed or produced empty tree; retry HKY85 boot=0: $f"
		rm -f "$tree" "$stats" "$boot_stats" "$boot_trees"
		rc=0
		last_log="${logf}.boot0"
		run_phyml_cmd "$phyml_retry_cpus" "$f" "$leave_duplicates" "$last_log" 0 "HKY85" || rc=$?
	fi
	if [[ $rc -ne 0 || ! -s $tree ]]; then
		log "WARN phyml HKY85 boot=0 failed or produced empty tree; retry JC69 boot=0: $f"
		rm -f "$tree" "$stats" "$boot_stats" "$boot_trees"
		rc=0
		last_log="${logf}.jc69"
		run_phyml_cmd "$phyml_retry_cpus" "$f" "$leave_duplicates" "$last_log" 0 "JC69" || rc=$?
	fi
	if [[ $rc -ne 0 || ! -s $tree ]]; then
		log "WARN phyml JC69 boot=0 failed or produced empty tree; retry JC69 boot=0 c=1: $f"
		rm -f "$tree" "$stats" "$boot_stats" "$boot_trees"
		rc=0
		last_log="${logf}.jc69.c1"
		run_phyml_cmd "$phyml_retry_cpus" "$f" "$leave_duplicates" "$last_log" 0 "JC69" 1 || rc=$?
	fi
	if [[ $rc -ne 0 || ! -s $tree ]]; then
		log "WARN phyml failed after all model retries; build ape NJ/NJS fallback tree: $f"
		rm -f "$tree" "$stats" "$boot_stats" "$boot_trees"
		rc=0
		last_log="${logf}.nj"
		run_phyml_nj_fallback "$f" "$tree" "$stats" "$last_log" || rc=$?
	fi
	if [[ $rc -ne 0 ]]; then
		log "ERROR phyml failed rc=$rc: $f; see $last_log"
		return 1
	fi
	if [[ ! -s $tree ]]; then
		log "ERROR phyml finished but tree is missing/empty: $tree; see $last_log; $(clean_msg "$last_log")"
		return 1
	fi
	log "phyml done: $f"
	return 0
}

pbase_auto(){
	local c=$1
	if [[ $ref_pop == ALL ]]; then
		[[ -s "$dirmod/pfile/ALL.chr$c.pgen" ]] && echo "$dirmod/pfile/ALL.chr$c" || echo "$dirmod/pfile/chr$c"
	else echo "$dirmod/pfile/${ref_pop}.chr$c"; fi
}
require_pfile(){
	local pfx=$1 context=pfile
	[[ $# -ge 2 ]] && context=$2
	local miss=()
	[[ -s $pfx.pgen ]] || miss+=("$pfx.pgen")
	[[ -s $pfx.pvar ]] || miss+=("$pfx.pvar")
	[[ -s $pfx.psam ]] || miss+=("$pfx.psam")
	if (( ${#miss[@]} > 0 )); then
		log "ERROR missing pfile for $context: prefix=$pfx; missing=${miss[*]}"
		exit 1
	fi
	check_pvar_format "$pfx.pvar" "$context" || exit 1
}
vcf_auto(){
	local c=$1
	if [[ $ref_pop == ALL && -s $dirmod/vcf/chr$c.vcf.gz ]]; then echo "$dirmod/vcf/chr$c.vcf.gz"; else echo "$dirmod/vcf/${ref_pop}.chr$c.vcf.gz"; fi
}
chrX_var(){
	# Usage: chrX_var X_PAR1_START 60001 -> value of X_PAR1_START_b37/b38 if defined.
	local stem=$1 default=$2
	local var="${stem}_${genome_build}"
	printf "%s" "${!var:-$default}"
}
chrX_part(){
	local bp=$1 par=nonPar p1s p1e p2s p2e
	p1s=$(chrX_var X_PAR1_START 60001); p1e=$(chrX_var X_PAR1_END 2699520)
	p2s=$(chrX_var X_PAR2_START 154931044); p2e=$(chrX_var X_PAR2_END 155260560)
	((bp>=p1s && bp<=p1e)) && par=par
	((bp>=p2s && bp<=p2e)) && par=par
	echo "$par"
}
first_pbase(){
	local p
	for p in "$@"; do [[ -s "$p.pgen" && -s "$p.pvar" && -s "$p.psam" ]] && { echo "$p"; return 0; }; done
	# Return the first candidate so require_pfile reports an informative error.
	echo "$1"
}
first_vcf(){
	local f
	for f in "$@"; do [[ -s "$f" ]] && { echo "$f"; return 0; }; done
	echo "$1"
}
pbase_by_bp(){
	local c=$1 bp=0 par
	[[ $# -ge 2 ]] && bp=$2
	if [[ $c == X || $c == 23 ]]; then
		par=$(chrX_part "$bp")
		if [[ $ref_pop == ALL ]]; then
			if [[ $par == par ]]; then
				first_pbase "$dirmod/pfile/ALL.male.chrX.par" "$dirmod/pfile/ALL.chrX.par" "$dirmod/pfile/ALL.chrX" "$dirmod/pfile/chrX"
			else
				first_pbase "$dirmod/pfile/ALL.male.chrX.nonPar" "$dirmod/pfile/ALL.chrX.nonPar" "$dirmod/pfile/ALL.chrX" "$dirmod/pfile/chrX"
			fi
		else
			first_pbase "$dirmod/pfile/${ref_pop}.male.chrX.$par" "$dirmod/pfile/${ref_pop}.chrX.$par" "$dirmod/pfile/${ref_pop}.chrX"
		fi
	else
		pbase_auto "$c"
	fi
}
vcf_by_region(){
	local c=$1 bp=0 par
	[[ $# -ge 2 ]] && bp=$2
	if [[ $c == X || $c == 23 ]]; then
		par=$(chrX_part "$bp")
		if [[ $ref_pop == ALL ]]; then
			if [[ $par == par ]]; then
				first_vcf "$dirmod/vcf/ALL.male.chrX.par.vcf.gz" "$dirmod/vcf/ALL.chrX.par.vcf.gz" "$dirmod/vcf/ALL.chrX.vcf.gz" "$dirmod/vcf/chrX.vcf.gz"
			else
				first_vcf "$dirmod/vcf/ALL.male.chrX.nonPar.vcf.gz" "$dirmod/vcf/ALL.chrX.nonPar.vcf.gz" "$dirmod/vcf/ALL.chrX.vcf.gz" "$dirmod/vcf/chrX.vcf.gz"
			fi
		else
			first_vcf "$dirmod/vcf/${ref_pop}.male.chrX.$par.vcf.gz" "$dirmod/vcf/${ref_pop}.chrX.$par.vcf.gz" "$dirmod/vcf/${ref_pop}.chrX.vcf.gz"
		fi
	else
		vcf_auto "$c"
	fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 1000G ID matching and archaic VCF helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lead_id(){
	local pvar=$1 snp=$2 bp=$3 ea="" oa="" line ch rb id rf al p2 rf2 al2
	[[ $# -ge 4 ]] && ea=$4
	[[ $# -ge 5 ]] && oa=$5
	line=$(pvar_awk "$pvar ID+BP $snp $bp" -v s="$snp" -v p="$bp" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p && $3==s{print $1,$2,$3,$4,$5; exit}' "$pvar") || return 1
	if [[ -z $line && $snp =~ ^[^:]+:([0-9]+)_([ACGT]+)_([ACGT]+)$ ]]; then
		p2=${BASH_REMATCH[1]}; rf2=${BASH_REMATCH[2]}; al2=${BASH_REMATCH[3]}
		line=$(pvar_awk "$pvar chrpos $snp" -v p="$p2" -v r="$rf2" -v a="$al2" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p && toupper($4)==r && toupper($5)==a{print $1,$2,$3,$4,$5; exit}' "$pvar") || return 1
	fi
	if [[ -z $line ]]; then
		line=$(pvar_awk "$pvar ID-only $snp" -v s="$snp" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $3==s{n++; z=$1 FS $2 FS $3 FS $4 FS $5} END{if(n==1) print z}' "$pvar") || return 1
	fi
	if [[ -z $line && -n $ea && -n $oa ]]; then
		line=$(pvar_awk "$pvar allele+BP $snp $bp" -v p="$bp" -v a="$ea" -v b="$oa" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p{r=toupper($4); alt=toupper($5); a=toupper(a); b=toupper(b); if((r==a && alt==b)||(r==b && alt==a)){print $1,$2,$3,$4,$5; exit}}' "$pvar") || return 1
	fi
	if [[ -z $line ]]; then
		line=$(pvar_awk "$pvar unique BP $bp" -v p="$bp" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p{n++; z=$1 FS $2 FS $3 FS $4 FS $5} END{if(n==1) print z}' "$pvar") || return 1
	fi
	[[ -n $line ]] || return 1
	read -r ch rb id rf al <<< "$line"
	[[ $id == "." || -z $id ]] && id="${ch}:${rb}:${rf}:${al}"
	printf "%s\t%s\n" "$rb" "$id"
}

arch_vcf(){
	local a=$1 c=$2
	case "$a" in
		vindija) echo "$dirarch/Vindija/chr${c}_mq25_mapab100.vcf.gz" ;;
		altai) echo "$dirarch/Altai/chr${c}_mq25_mapab100.vcf.gz" ;;
		chagyr) [[ -s $dirarch/Chagyr/chr${c}.noRB.bgz.vcf.gz ]] && echo "$dirarch/Chagyr/chr${c}.noRB.bgz.vcf.gz" || echo "$dirarch/Chagyr/chr${c}.noRB.vcf.gz" ;;
		denisova) echo "$dirarch/Denisova/chr${c}_mq25_mapab100.vcf.gz" ;;
		denisova25) echo "$dirarch/Denisova25/chr${c}.Den25.L35MQ25.B30.map35_100.vcf.gz" ;;
	esac
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step1: prepare COJO/GWAS lead SNPs and risk alleles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s1_prep(){
	check_cmd
	if [[ $lead_match == strict ]]; then
		local need_prep=0 t first
		for t in $traits; do
			[[ -s $(lead_assoc_file "$t") && -s $(lead_3col_file "$t") && -s $(lead_id_sanity_file "$t") ]] || need_prep=1
			[[ -s $dirout/lead/prep.version ]] && grep -qx 'prep_input_v3_id_first_pos_if_missing_id' "$dirout/lead/prep.version" || need_prep=1
			[[ $(data_rows "$(lead_assoc_file "$t")") -gt 0 ]] || need_prep=1
		done
		if (( need_prep == 0 )); then
			log "s1 prep_input skip: existing lead files are complete for traits=$traits"
		else
			log "Now running Step1 prep input: input=$dirgwas; output=$dirout/lead"
			Rscript "$dirscript/loci.R" prep_input_local --dirgwas "$dirgwas" --dirout "$dirout" --dirmod "$dirmod/pfile" --ref_pop "$ref_pop" --traits "${traits_arr[@]}"
			log "s1 prep_input local done: $(wc -l < "$dirout/lead/lead.assoc") lines in lead.assoc"
		fi
		if [[ " $traits " == *" bald "* ]] && ! awk 'NR>1 && $2==3 && $3=="rs35044562"{f=1} END{exit !f}' "$(lead_assoc_file bald)"; then
			awk 'BEGIN{FS=OFS="\t"} NR==1{print; print "bald",3,"rs35044562",45909024,"G","A",1,"G","rs35044562","ID_manual"; next} {print}' "$(lead_assoc_file bald)" > "$dirout/lead/.lead.assoc" && mv "$dirout/lead/.lead.assoc" "$(lead_assoc_file bald)"
			awk 'BEGIN{FS=OFS="\t"} NR==1{print; print 3,"rs35044562",45909024; next} {print}' "$(lead_3col_file bald)" > "$dirout/lead/.lead.3col" && mv "$dirout/lead/.lead.3col" "$(lead_3col_file bald)"
			log "manual add bald rs35044562 at 1000G POS 45909024"
		fi
		merge_lead_assoc_files
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.single.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\treason\tdetail\tpbase\tvid\n" > "$dirout/lead/lead_1000G.fail.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tPAR\tLD_pfile\n" > "$dirout/lead/chrX_PAR.tsv"
		printf "trait\tchr\tsnp\tlead_bp\tpvar_bp\tvid\tmatch\tpbase\tstatus\tn_ld_snp\tmessage\n" > "$dirout/lead/ld_debug.tsv"
	else
		local need_prep=0
		for t in $traits; do
			[[ -s $(lead_assoc_file "$t") && $(data_rows "$(lead_assoc_file "$t")") -gt 0 ]] || need_prep=1
		done
		if [[ -s $dirout/lead/lead.assoc && $need_prep -eq 0 ]]; then
			log "s1 prep_input skip: existing lead files are complete for traits=$traits"
		else
			log "Now running Step1 prep input: input=$dirgwas and $dircojo; output=$dirout/lead"
			Rscript "$dirscript/loci.R" prep_input --dirgwas "$dirgwas" --dircojo "$dircojo" --dirout "$dirout" --traits "${traits_arr[@]}"
			log "s1 done: lead.assoc $(nrow "$dirout/lead/lead.assoc") rows"
		fi
	fi

	# Positive/negative-control loci must be added only once, to positive_trait from
	# the original configured trait list.  In the per-trait workflow traits_arr[0]
	# is the current trait, so using it directly would add the BED to every trait.
	if [[ $add_positive_loci != 0 && -s $positive_loci && ${traits_arr[0]:-} == "$positive_trait" ]]; then
		local first=$positive_trait
		log "s1 add positive-control loci to configured positive_trait only: $first; bed=$positive_loci"
		Rscript "$dirscript/loci.R" add_positive_loci --dirgwas "$dirgwas" --dirmod "$dirmod/pfile" --dirout "$dirout" --bed "$positive_loci" --trait "$first" --ref_pop "$ref_pop"
		merge_lead_assoc_files
	elif [[ $add_positive_loci != 0 && ${traits_arr[0]:-} != "$positive_trait" ]]; then
		log "s1 skip positive-control BED for trait=${traits_arr[0]:-NA}; positive_trait=$positive_trait"
	elif [[ $add_positive_loci != 0 ]]; then
		log "WARN sanity-check loci BED not found or empty: $positive_loci"
	fi
	log_process_s1
}

force_positive_picks(){
	if [[ $add_positive_loci == 0 ]]; then
		return 0
	fi
	if [[ ${traits_arr[0]:-} != "$positive_trait" ]]; then
		return 0
	fi
	if [[ ! -s $positive_loci ]]; then
		log "WARN sanity-check loci BED not found or empty: $positive_loci"
		return 0
	fi
	local first=$positive_trait pickf=$dirout/lead/pick.tsv singlef=$dirout/lead/pick.single.tsv posf=$dirout/lead/positive_pick.tsv tmp
	[[ -s $pickf && -s "$(lead_assoc_file "$first")" ]] || return 0
	tmp=$(mktemp)
	awk -v trait="$first" -v assoc="$(lead_assoc_file "$first")" '
		BEGIN{FS="[[:space:]]+"; OFS="\t"; while((getline < assoc)>0){ if(++arow==1) continue; assoc_chr[$3]=$2; assoc_bp[$3]=$4 }}
		FILENAME==ARGV[1] && NR==FNR { if(NR>1) have[$1 SUBSEP $3]=1; next }
		/^[[:space:]]*($|#)/ { next }
		{
			chr=$1; start=$2+0; end=$3+0; snp=$4
			gsub(/^chr|^CHR/,"",chr)
			if(chr=="X") chr=23
			if(snp=="" || start<=0 || end<=0) next
			if(start>end){ x=start; start=end; end=x }
			bp=(snp in assoc_bp ? assoc_bp[snp] : end)
			out_chr=(snp in assoc_chr ? assoc_chr[snp] : chr)
			key=trait SUBSEP snp
			if(!(key in have)){
				print trait,out_chr,snp,bp,start,end,2,end-start+1
				have[key]=1
				added++
			}
		}
		END{ if(added > 0) printf("%d", added) > "/dev/stderr" }
	' "$pickf" "$positive_loci" > "$tmp" 2>"$tmp.n"
	if [[ -s $tmp ]]; then
		cat "$tmp" >> "$pickf"
		awk 'BEGIN{FS=OFS="\t"} NR==FNR{if(NR>1) forced[$1 FS $3]=1; next} NR==1 || !(($1 FS $3) in forced)' "$pickf" "$singlef" > "$singlef.tmp" && mv "$singlef.tmp" "$singlef"
		awk 'NR==1 || !seen[$1 FS $2 FS $3 FS $4]++' "$pickf" > "$pickf.tmp" && mv "$pickf.tmp" "$pickf"
		log "s3 force positive-control picks from BED: added $(cat "$tmp.n") row(s)"
	fi
	awk -v trait="$first" -v bed="$positive_loci" '
		BEGIN{FS="[[:space:]]+"; OFS="\t"; print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"}
		NR==FNR{ if(NR>1){ key=$1 SUBSEP $3; row[key]=$1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8; seen_snp[key]=1 } next }
		/^[[:space:]]*($|#)/ { next }
		{
			chr=$1; start=$2+0; end=$3+0; snp=$4
			gsub(/^chr|^CHR/,"",chr)
			if(chr=="X") chr=23
			if(snp=="" || start<=0 || end<=0) next
			key=trait SUBSEP snp
			if((key in seen_snp) && !(key in emitted)){ print row[key]; emitted[key]=1 }
		}
	' "$pickf" "$positive_loci" > "$posf"
	log "summary: read $(data_rows "$posf") sanity-check loci from $positive_loci; saved in $posf"
	rm -f "$tmp" "$tmp.n"
}

calculate_ld_avcf(){
	check_cmd
	local t=$1 c=$2 cn pfx tmp assoc nlead snp bp ea oa rb id out lookup line refbp vid _ref _alt n_id
	cn=$(chrn "$c"); tmp=$(stage_trait_dir ld "$t")/chr$c; assoc=$(lead_assoc_file "$t")
	mkdir -p "$tmp"
	awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $2==cn || toupper($2)==c' "$assoc" > "$tmp/lead.assoc"
	nlead=$(data_rows "$tmp/lead.assoc")
	if valid_chr_ld "$tmp" "$nlead"; then log "s2 skip $t chr$c"; return 0; fi
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"
	(( nlead > 0 )) || { printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$tmp/block.tsv"; touch "$tmp/.done"; return 0; }
	log "s2 avcf-LD $t chr$c"
	awk 'BEGIN{FS=OFS="\t"} NR>1{print $3,$4,$5,$6}' "$tmp/lead.assoc" |
	while IFS=$'\t' read -r snp bp ea oa; do
		pfx=$(pbase_by_bp "$c" "$bp")
		require_pfile "$pfx" "$t chr$c $snp"
		lookup=$(pvar_lookup_local "$pfx" "$tmp/lead.assoc" "$tmp")
		line=$(awk -v snp="$snp" -v bp="$bp" 'BEGIN{FS=OFS="\t"} $1==bp && $2==snp{print $1,$2,$3,$4; exit}' "$lookup")
		if [[ -n $line ]]; then
			read -r refbp vid _ref _alt <<< "$line"
		else
			n_id=$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} $2==snp{n++} END{print n+0}' "$lookup")
			if [[ $n_id -eq 1 ]]; then
				read -r refbp vid _ref _alt <<< "$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} $2==snp{print $1,$2,$3,$4; exit}' "$lookup")"
				[[ $refbp != "$bp" ]] && log "WARN $t chr$c $snp: lead_bp=$bp but pvar_bp=$refbp; using pvar_bp"
			else
				line=$(awk -v bp="$bp" -v ea="$ea" -v oa="$oa" 'BEGIN{FS=OFS="\t"} $1==bp{ref=toupper($3); alt=toupper($4); ea=toupper(ea); oa=toupper(oa); if((ref==ea && alt==oa) || (ref==oa && alt==ea)){print $1,$2,$3,$4; exit}}' "$lookup")
				if [[ -n $line ]]; then
					read -r refbp vid _ref _alt <<< "$line"
				else
					line=$(awk -v bp="$bp" 'BEGIN{FS=OFS="\t"} $1==bp{id=$2; r=$3; a=$4; n++} END{if(n==1) print bp,id,r,a}' "$lookup")
					if [[ -n $line ]]; then
						read -r refbp vid _ref _alt <<< "$line"
					else
						log "WARN not in 1000G: $t chr$c $snp $bp"
						continue
					fi
				fi
			fi
		fi
		rb=$refbp; id=$vid; out=$tmp/$snp
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$cn" "$snp" "$rb" "$cn" "$rb" "$id" >> "$tmp/ld.tsv"
		run_timeout_cmd "$plink_timeout" plink2 --pfile "$pfx" --set-missing-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 missing --r2-unphased allow-ambiguous-allele --ld-snp "$id" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --threads "$plink_threads" --memory 4000 --out "$out" > "$out.plink.log" 2>&1 || { log "WARN plink2 failed/timeout: $out.plink.log $(clean_msg "$out.plink.log")"; continue; }
		[[ -s $out.vcor ]] || continue
		awk -v t="$t" -v c="$cn" -v s="$snp" -v b="$rb" -v lead="$id" 'BEGIN{FS="[ \t]+";OFS="\t"} NR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="POS_A") pa=i; if($i=="ID_A") ia=i; if($i=="POS_B") pb=i; if($i=="ID_B") ib=i; if($i~/R2$/) ir=i}; next} pa&&ia&&pb&&ib&&ir{if($ia==lead) print t,c,s,b,c,$pb,$ib,$ir; else if($ib==lead) print t,c,s,b,c,$pa,$ia,$ir}' "$out.vcor" >> "$tmp/ld.tsv"
	done
	awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR==1{next}{k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0; if(!(k in mn)||bp<mn[k]) mn[k]=bp; if(!(k in mx)||bp>mx[k]) mx[k]=bp; if(pos<mn[k]) mn[k]=pos; if(pos>mx[k]) mx[k]=pos; u[k FS bp]=u[k FS pos]=1} END{print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"; for(k in mn){n=0; for(i in u) if(index(i,k FS)==1) n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
	touch "$tmp/.done"
}

pvar_lookup_local(){
	local pbase=$1 lead_assoc=$2 tmp=$3 cache key out nkeys
	cache="$tmp/pvar.$(basename "$pbase").lookup.tsv"
	[[ -s $cache ]] && { printf "%s\n" "$cache"; return; }
	key="$tmp/pvar.$(basename "$pbase").keys.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3,$4,$5,$6}' "$lead_assoc" > "$key"
	nkeys=$(data_rows "$lead_assoc")
	out="$cache.tmp.$$"
	[[ ${pvar_lookup_debug:-0} == 1 ]] && log "pvar lookup start: $(basename "$pbase").pvar keys=$nkeys timeout=${pvar_scan_timeout:-none}" >&2
	if ! pvar_awk "$pbase.pvar lookup" 'BEGIN{FS=OFS="\t"} NR==FNR{bp[$3]=1; id[$2]=1; allele[$3 FS toupper($4) FS toupper($5)]=1; allele[$3 FS toupper($5) FS toupper($4)]=1; next} /^#/{next} (($2 in bp) || ($3 in id) || (($2 FS toupper($4) FS toupper($5)) in allele)){print $2,$3,$4,$5}' "$key" "$pbase.pvar" > "$out"; then
		rm -f "$key" "$out"
		return 1
	fi
	mv "$out" "$cache"
	[[ ${pvar_lookup_debug:-0} == 1 ]] && log "pvar lookup done: $(basename "$pbase").pvar matched=$(wc -l < "$cache")" >&2
	rm -f "$key"
	printf "%s\n" "$cache"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step2: calculate high-LD blocks from 1000G pfiles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calculate_ld_strict(){
	check_cmd
	local t=$1 c=$2 cn assoc tmp nlead snp bp ea oa par pbase flag lookup line refbp vid _ref _alt match n_id outpre plog outvcor msg cmdtxt
	cn=$(chrn "$c"); assoc=$(lead_assoc_file "$t"); tmp=$(stage_trait_dir ld "$t")/chr$c
	mkdir -p "$tmp"
	awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR==1 || $2==cn || toupper($2)==c' "$assoc" > "$tmp/lead.assoc"
	nlead=$(data_rows "$tmp/lead.assoc")
	if valid_chr_ld "$tmp" "$nlead"; then log "LD $t chr$c: skip existing complete result ($nlead lead SNPs)"; return 0; fi
	log "LD $t chr$c: $nlead lead SNPs"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tvid\tpbase\tlead_bp0\tmatch\n" > "$tmp/lead.map.tsv"
	(( nlead > 0 )) || { printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$tmp/block.tsv"; touch "$tmp/.done"; return 0; }
	awk 'BEGIN{FS=OFS="\t"} NR>1{print $3,$4,$5,$6}' "$tmp/lead.assoc" | while IFS=$'\t' read -r snp bp ea oa; do
		[[ -n $snp && -n $bp ]] || continue
		pbase=$(pbase_by_bp "$c" "$bp")
		if [[ $c == X ]]; then par=$(chrX_part "$bp"); echo -e "$t\t23\t$snp\t$bp\t$par\t$pbase" >> "$dirout/lead/chrX_PAR.tsv"; fi
		flag="$tmp/sanity.$(basename "$pbase").tsv"
		check_pbase "$pbase" "$tmp/lead.assoc" "$t chr$c" "$flag"
		lookup=$(pvar_lookup_local "$pbase" "$tmp/lead.assoc" "$tmp")
		line=$(awk -v snp="$snp" -v bp="$bp" 'BEGIN{FS=OFS="\t"} $1==bp && $2==snp{print $1,$2,$3,$4; exit}' "$lookup")
		if [[ -n $line ]]; then read -r refbp vid _ref _alt <<< "$line"; match=ID_BP; else
			n_id=$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} $2==snp{n++} END{print n+0}' "$lookup")
			if [[ $n_id -eq 1 ]]; then
				read -r refbp vid _ref _alt <<< "$(awk -v snp="$snp" 'BEGIN{FS=OFS="\t"} $2==snp{print $1,$2,$3,$4; exit}' "$lookup")"; match=ID_only
				[[ $refbp != "$bp" ]] && log "WARN $t chr$c $snp: lead_bp=$bp but pvar_bp=$refbp; using pvar_bp"
			else
				line=$(awk -v bp="$bp" -v ea="$ea" -v oa="$oa" 'BEGIN{FS=OFS="\t"} $1==bp{ref=toupper($3); alt=toupper($4); ea=toupper(ea); oa=toupper(oa); if((ref==ea && alt==oa) || (ref==oa && alt==ea)){print $1,$2,$3,$4; exit}}' "$lookup")
				if [[ -n $line ]]; then read -r refbp vid _ref _alt <<< "$line"; match=ALLELE_BP; else
					line=$(awk -v bp="$bp" 'BEGIN{FS=OFS="\t"} $1==bp{id=$2; r=$3; a=$4; n++} END{if(n==1) print bp,id,r,a}' "$lookup")
					if [[ -n $line ]]; then read -r refbp vid _ref _alt <<< "$line"; match=UNIQUE_BP; else
						msg="not found by ID+BP, ID-only, allele+BP, or unique BP"
						write_fail "$t" "$cn" "$snp" "$bp" "not_in_1000G_pvar_or_allele_mismatch" "$msg" "$pbase" "NA"
						write_debug "$t" "$cn" "$snp" "$bp" "NA" "NA" "none" "$pbase" "pvar_miss" 0 "$msg"
						continue
					fi
				fi
			fi
		fi
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$cn" "$snp" "$refbp" "$cn" "$refbp" "$vid" >> "$tmp/ld.tsv"
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$t" "$cn" "$snp" "$refbp" "$vid" "$pbase" "$bp" "$match" >> "$tmp/lead.map.tsv"
	done
	awk 'BEGIN{FS=OFS="\t"} NR>1{print $6}' "$tmp/lead.map.tsv" | sort -u | while IFS= read -r pbase; do
		[[ -n $pbase ]] || continue
		outpre="$tmp/$(basename "$pbase").batch"; plog="$outpre.plink2.log"; outvcor="$outpre.vcor"
		awk -v pbase="$pbase" 'BEGIN{FS=OFS="\t"} NR>1 && $6==pbase{print $5}' "$tmp/lead.map.tsv" | sort -u > "$outpre.ld_ids"
		[[ -s $outpre.ld_ids ]] || continue
		log "LD $t chr$c: PLINK start $(basename "$pbase") with $(wc -l < "$outpre.ld_ids") lead IDs"
		cmdtxt="plink2 --pfile $pbase --allow-extra-chr --make-founders --r2-unphased allow-ambiguous-allele --ld-snp-list $outpre.ld_ids --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 $ld_r2 --threads $plink_threads --out $outpre"
		printf "%s\n" "$cmdtxt" > "$outpre.cmd"
		if run_timeout_cmd "$plink_timeout" plink2 --pfile "$pbase" --allow-extra-chr --make-founders --r2-unphased allow-ambiguous-allele --ld-snp-list "$outpre.ld_ids" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --threads "$plink_threads" --out "$outpre" > "$plog" 2>&1; then
			if [[ -s $outvcor ]]; then
				log "LD $t chr$c: PLINK done $(basename "$pbase"); vcor_rows=$(awk 'NR>1{n++} END{print n+0}' "$outvcor")"
				awk -v pbase="$pbase" 'BEGIN{FS="[ \t]+";OFS="\t"} NR==FNR{if(NR>1 && $6==pbase) map[$5]=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4; next} FNR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="ID_A") ia=i; if($i=="POS_B") pb=i; if($i=="ID_B") ib=i; if($i~/R2$/) ir=i}; next} ia&&pb&&ib&&ir&&($ia in map){split(map[$ia],m,SUBSEP); print m[1],m[2],m[3],m[4],m[2],$pb,$ib,$ir}' "$tmp/lead.map.tsv" "$outvcor" >> "$tmp/ld.tsv"
				awk -v pbase="$pbase" -v dbg="$dirout/lead/ld_debug.tsv" 'BEGIN{FS="[ \t]+";OFS="\t"} NR==FNR{if(NR>1 && $6==pbase) map[$5]=$0; next} FNR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="ID_A") ia=i}; next} ia&&($ia in map){n[$ia]++} END{for(id in map){split(map[id],m,FS); print m[1],m[2],m[3],m[7],m[4],m[5],m[8],m[6],"ok_batch",n[id]+0,"plink_ok_batch" >> dbg}}' "$tmp/lead.map.tsv" "$outvcor"
			else
				msg="PLINK finished but $outvcor is missing/empty"
				awk -v pbase="$pbase" -v dbg="$dirout/lead/ld_debug.tsv" -v msg="$msg" 'BEGIN{FS=OFS="\t"} NR>1 && $6==pbase{print $1,$2,$3,$7,$4,$5,$8,$6,"ok_no_vcor",0,msg >> dbg}' "$tmp/lead.map.tsv"
				log "LD $t chr$c: PLINK done $(basename "$pbase") but vcor is empty"
			fi
		else
			msg=$(clean_msg "$plog")
			awk -v pbase="$pbase" -v fail="$dirout/lead/lead_1000G.fail.tsv" -v dbg="$dirout/lead/ld_debug.tsv" -v msg="$msg" 'BEGIN{FS=OFS="\t"} NR>1 && $6==pbase{print $1,$2,$3,$7,"plink2_ld_failed",msg,$6,$5 >> fail; print $1,$2,$3,$7,$4,$5,$8,$6,"plink_failed",0,msg >> dbg}' "$tmp/lead.map.tsv"
			log "FAIL $t chr$c $(basename "$pbase") batch: $msg"
		fi
	done
	awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR==1{next}{k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0; if(!(k in mn)||bp<mn[k]) mn[k]=bp; if(!(k in mx)||bp>mx[k]) mx[k]=bp; if(pos<mn[k]) mn[k]=pos; if(pos>mx[k]) mx[k]=pos; u[k FS bp]=u[k FS pos]=1} END{print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"; for(k in mn){n=0; for(i in u) if(index(i,k FS)==1) n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
	touch "$tmp/.done"
	log "LD $t chr$c: done; ld_rows=$(data_rows "$tmp/ld.tsv") blocks=$(data_rows "$tmp/block.tsv")"
}

trait_ld_complete(){
	local t=$1 d c cn nlead
	d=$(stage_trait_dir ld "$t")
	for c in $chrs; do
		cn=$(chrn "$c")
		nlead=$(awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR>1 && ($2==cn || toupper($2)==c){n++} END{print n+0}' "$(lead_assoc_file "$t")")
		[[ $(data_rows "$d/chr$c/lead.assoc") -eq $nlead ]] || return 1
		valid_chr_ld "$d/chr$c" "$nlead" || return 1
	done
}

merge_trait_ld(){
	local t=$1 d x c f first pickf singlef
	d=$(stage_trait_dir ld "$t")
	pickf=$(lead_pick_file "$t" pick)
	singlef=$(lead_pick_file "$t" single)
	for x in ld block; do
		first=1; : > "$d/$x.tsv"
		for c in $chrs; do
			f=$d/chr$c/$x.tsv; [[ -f $f ]] || continue
			if [[ $first -eq 1 ]]; then cat "$f"; first=0; else awk 'NR>1' "$f"; fi
		done > "$d/$x.tsv"
		awk 'NR==1 || $1!="trait"' "$d/$x.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.$x" && mv "$d/.$x" "$d/$x.tsv"
	done
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$pickf"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$singlef"
	awk 'BEGIN{FS=OFS="\t"} NR>1 && $7>1' "$d/block.tsv" >> "$pickf"
	awk 'BEGIN{FS=OFS="\t"} NR>1 && $7==1' "$d/block.tsv" >> "$singlef"
	log "========== DONE $t: $(awk 'NR>1{n++} END{print n+0}' "$d/block.tsv") blocks; pick=$(data_rows "$pickf"); single=$(data_rows "$singlef") =========="
}

calculate_trait_ld_strict(){
	local t=$1 c d nlead ld_status=0
	log "========== DO $t =========="
	if trait_ld_complete "$t"; then log "LD $t: all chromosome results complete; merge existing results"; merge_trait_ld "$t"; return 0; fi
	for c in $chrs; do
		calculate_ld_strict "$t" "$c" &
		while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do wait -n || ld_status=1; done
	done
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || ld_status=1; done
	(( ld_status == 0 )) || { log "ERROR one or more LD chromosome jobs failed for $t"; return 1; }
	d=$(stage_trait_dir ld "$t")
	for c in $chrs; do
		nlead=$(data_rows "$d/chr$c/lead.assoc")
		valid_chr_ld "$d/chr$c" "$nlead" || { log "ERROR incomplete LD output for $t chr$c; remove $d/chr$c and rerun"; return 1; }
	done
	merge_trait_ld "$t"
}

s2_ld(){
	local t c running=0 trait_status=0
	log "Now running Step2 LD: input=$dirout/lead/lead.assoc; output=$dirout/ld and $dirout/lead/pick.tsv"
	if [[ $ld_calc == strict ]]; then
		for t in $traits; do
			calculate_trait_ld_strict "$t" &
			while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do wait -n || trait_status=1; done
		done
		while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || trait_status=1; done
		(( trait_status == 0 )) || { log "ERROR one or more LD trait jobs failed"; exit 1; }
		if [[ ${#traits_arr[@]} -eq 1 ]]; then
			log_trait_count "candidate high-LD loci for haplotype analysis" "$dirout/lead/pick.tsv"
			log_trait_count "single-SNP lead loci excluded before haplotype analysis" "$dirout/lead/pick.single.tsv"
			log "s2 LD done. Check: $dirout/lead/ld_debug.tsv and $dirout/lead/lead_1000G.fail.tsv"
			return 0
		fi
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.single.tsv"
		for t in $traits; do
			[[ -s "$(lead_pick_file "$t" pick)" ]] && awk 'NR>1' "$(lead_pick_file "$t" pick)" >> "$dirout/lead/pick.tsv"
			[[ -s "$(lead_pick_file "$t" single)" ]] && awk 'NR>1' "$(lead_pick_file "$t" single)" >> "$dirout/lead/pick.single.tsv"
		done
		log_trait_count "candidate high-LD loci for haplotype analysis" "$dirout/lead/pick.tsv"
		log_trait_count "single-SNP lead loci excluded before haplotype analysis" "$dirout/lead/pick.single.tsv"
		log "s2 LD done. Check: $dirout/lead/ld_debug.tsv and $dirout/lead/lead_1000G.fail.tsv"
		return 0
	fi
	for t in $traits; do
		for c in $chrs; do
			limit_jobs "$job_in_trait"
			if [[ $ld_calc == strict ]]; then calculate_ld_strict "$t" "$c" & else calculate_ld_avcf "$t" "$c" & fi
		done
	done
	wait_all
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step3: merge LD blocks and choose multi-SNP loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s3_pick(){
	if [[ -s $dirout/lead/pick.tsv && -s $dirout/lead/pick.single.tsv ]] &&
		valid_header "$dirout/lead/pick.tsv" && valid_header "$dirout/lead/pick.single.tsv"; then
		local complete=1 t
		for t in $traits; do trait_ld_complete "$t" || complete=0; done
		if (( complete == 1 )); then
			force_positive_picks
			log "s3 pick skip: existing pick.tsv and pick.single.tsv are complete"
			log_process_s3
			return 0
		fi
		log "s3 pick refresh: existing pick files found, but LD completeness check needs a rebuild"
	fi
	log "Now running Step3 pick loci: input=$dirout/ld/*/block.tsv; output=$dirout/lead/pick.tsv"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.single.tsv"
	for t in $traits; do
		d=$(stage_trait_dir ld "$t"); mkdir -p "$d"
		for x in ld block; do
			first=1; : > "$d/$x.tsv"
			for c in $chrs; do
				f=$d/chr$c/$x.tsv; [[ -s $f ]] || continue
				[[ $first -eq 1 ]] && { cat "$f"; first=0; } || awk 'NR>1' "$f"
			done > "$d/$x.tsv"
			awk 'NR==1 || $1!="trait"' "$d/$x.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.$x" && mv "$d/.$x" "$d/$x.tsv"
		done
		awk 'BEGIN{FS=OFS="\t"} NR>1 && $7>1' "$d/block.tsv" >> "$dirout/lead/pick.tsv"
		awk 'BEGIN{FS=OFS="\t"} NR>1 && $7==1' "$d/block.tsv" >> "$dirout/lead/pick.single.tsv"
	done
	force_positive_picks
	log "s3 done: pick $(nrow "$dirout/lead/pick.tsv") rows; single $(nrow "$dirout/lead/pick.single.tsv") rows"
	log_process_s3
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step4: build core VCF and genotype matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build_locus_matrix_strict(){
	local tr=$1 chr=$2 snp=$3 bp=$4 st=$5 en=$6 _n=$7 _size_bp=$8
	local cl=$chr lo=$st hi=$en d o kg ad name outname avcf f arc_raw arc_vcf arc_log nvar dot expected_arch
	expected_arch=$(find "$dirarch" -mindepth 1 -maxdepth 1 -type d | wc -l)
	[[ $bp -lt $lo ]] && lo=$bp
	[[ $bp -gt $hi ]] && hi=$bp
	d=$(stage_trait_dir coreVcf "$tr")/${chr}.${snp}.${bp}
	o=$(stage_trait_dir mat "$tr")/${chr}.${snp}.${bp}
	if valid_mat_locus "$o" "$expected_arch"; then log "s4 skip existing mat: $tr ${chr}.${snp}.${bp}"; return 0; fi
	log "s4 mat start: $tr ${chr}.${snp}.${bp} region=$chr:$lo-$hi"
	rm -rf "$d" "$o"; mkdir -p "$d" "$o"
	cl=$(chrl "$chr")
	kg=$(vcf_by_region "$cl" "$bp")
	[[ -s $kg ]] || { log "ERROR missing 1000G VCF for $tr chr$chr $snp: $kg"; return 1; }
	bcftools view -r "$cl:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg" || { log "ERROR bcftools failed: 1000G $tr chr$chr $snp"; return 1; }
	tabix -f -p vcf "$d/kg.vcf.gz" || { log "ERROR tabix failed: $d/kg.vcf.gz"; return 1; }
	bcftools query -l "$d/kg.vcf.gz" > "$o/kg.samples.tsv" || { log "ERROR bcftools query samples failed: $d/kg.vcf.gz"; return 1; }
	for ad in "$dirarch"/*; do
		[[ -d $ad ]] || continue
		name=$(basename "$ad"); outname=$(echo "$name" | tr '[:upper:]' '[:lower:]'); avcf=""
		for f in "$ad"/*chr"$cl"_*.vcf.gz "$ad"/*chr"$cl".*.vcf.gz; do [[ -s $f ]] || continue; avcf=$f; break; done
		[[ -n $avcf ]] || { log "WARN missing archaic VCF: $outname chr$cl"; continue; }
		arc_raw="$d/$outname.raw.vcf.gz"; arc_vcf="$d/$outname.vcf.gz"; arc_log="$d/$outname.project.log"
		bcftools view -r "$cl:$lo-$hi" -Oz -o "$arc_raw" "$avcf" || { log "ERROR bcftools failed: archaic raw $outname $cl:$lo-$hi"; return 1; }
		tabix -f -p vcf "$arc_raw" || { log "ERROR tabix failed: $arc_raw"; return 1; }
		Rscript "$dirscript/loci.R" prep_archaic "$d/kg.vcf.gz" "$arc_raw" "$arc_vcf" "$outname" > "$arc_log" 2>&1 || { log "ERROR archaic projection failed: $outname $cl:$lo-$hi; see $arc_log"; tail -20 "$arc_log"; return 1; }
		nvar=$(bcftools view -H "$arc_vcf" | wc -l)
		dot=$(bcftools view -H "$arc_vcf" | awk 'BEGIN{FS="\t"} $5=="."{n++} END{print n+0}')
		[[ $nvar -gt 0 && $dot -eq 0 ]] || { log "ERROR bad projected archaic VCF: $arc_vcf n=$nvar ALT_dot=$dot"; return 1; }
		rm -f "$arc_raw" "$arc_raw.tbi"
	done
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{for(i=6;i<=NF;i++) if($i=="0" || $i=="1") $i=$i"/"$i; print}' > "$o/kg.tsv" || { log "ERROR bcftools query failed: $d/kg.vcf.gz"; return 1; }
	for f in "$d"/*.vcf.gz; do
		name=$(basename "$f" .vcf.gz)
		[[ $name == kg || $name == *.raw ]] && continue
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$f" | awk 'BEGIN{FS=OFS="\t"} $4!="." && $3~/^[ACGT]$/ && $4~/^[ACGT]$/ {if($5=="0" || $5=="1") $5=$5"/"$5; print}' > "$o/$name.tsv" || { log "ERROR bcftools query failed: $f"; return 1; }
	done
	valid_mat_locus "$o" "$expected_arch" || { log "ERROR incomplete mat output: $o"; return 1; }
	rm -rf "$d"
	log "s4 mat done: $tr ${chr}.${snp}.${bp}; kg_rows=$(count_lines "$o/kg.tsv") archaic_files=$(find "$o" -maxdepth 1 -type f -name '*.tsv' ! -name 'kg.tsv' ! -name 'kg.samples.tsv' | wc -l)"
}

build_trait_matrices_strict(){
	local t=$1 tr chr snp bp st en _n _size_bp status=0 n_trait
	n_trait=$(awk -v t="$t" 'BEGIN{FS="\t"} NR>1 && $1==t{n++} END{print n+0}' "$dirout/lead/pick.tsv")
	log "s4 trait start: $t ($n_trait loci)"
	while IFS=$'\t' read -r tr chr snp bp st en _n _size_bp; do
		build_locus_matrix_strict "$tr" "$chr" "$snp" "$bp" "$st" "$en" "$_n" "$_size_bp" &
		while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do wait -n || status=1; done
	done < <(awk -v t="$t" 'BEGIN{FS=OFS="\t"} NR>1 && $1==t && !seen[$1 FS $2 FS $3 FS $4 FS $5 FS $6]++{print}' "$dirout/lead/pick.tsv")
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
	(( status == 0 )) || { log "ERROR one or more s4 locus jobs failed for $t"; return 1; }
	log "s4 trait done: $t"
}

build_matrices_strict(){
	local t status=0 npick
	npick=$(data_rows "$dirout/lead/pick.tsv")
	(( npick > 0 )) || { log "ERROR no high-LD blocks in $dirout/lead/pick.tsv; inspect $dirout/lead/lead_1000G.fail.tsv and rerun START_STEP=s1"; exit 1; }
	log_trait_count "candidate high-LD loci for haplotype analysis" "$dirout/lead/pick.tsv"
	log "Now running Step4 core/matrix: input=$dirout/lead/pick.tsv; output=$dirout/coreVcf and $dirout/mat"
	for t in $traits; do
		build_trait_matrices_strict "$t" &
		while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do wait -n || status=1; done
	done
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
	(( status == 0 )) || { log "ERROR one or more s4 trait jobs failed"; exit 1; }
	rm -rf "$dirout/coreVcf"
}

build_core_vcf_avcf(){
	check_cmd
	local t=$1 c=$2 cn cl snp bp st en n size lo hi d a f kgvcf
	cn=$(chrn "$c"); cl=$(chrl "$cn")
	log "s4 core $t chr$c"
	awk -v t="$t" -v c="$cn" 'BEGIN{FS=OFS="\t"} NR>1 && $1==t && $2==c{print}' "$dirout/lead/pick.tsv" |
	while IFS=$'\t' read -r t cn snp bp st en n size; do
		lo=$st; hi=$en; [[ $bp -lt $lo ]] && lo=$bp; [[ $bp -gt $hi ]] && hi=$bp
		d=$(stage_trait_dir coreVcf "$t")/${cn}.${snp}.${bp}; mkdir -p "$d"
		kgvcf=$(vcf_by_region "$c" "$bp")
		[[ -s $kgvcf ]] || { log "WARN missing 1000G VCF: $kgvcf"; continue; }
		if [[ ! -s $d/kg.vcf.gz ]]; then bcftools view -r "$cl:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kgvcf"; tabix -f -p vcf "$d/kg.vcf.gz"; fi
		for a in vindija altai chagyr denisova denisova25; do
			f=$(arch_vcf "$a" "$cl"); [[ -s $f ]] || continue
			[[ -s $d/$a.vcf.gz ]] && continue
			if [[ $archaic_VCF_match == strict ]]; then
				raw="$d/$a.raw.vcf.gz"
				bcftools view -r "$cl:$lo-$hi" -Oz -o "$raw" "$f"; tabix -f -p vcf "$raw"
				Rscript "$dirscript/loci.R" prep_archaic "$d/kg.vcf.gz" "$raw" "$d/$a.vcf.gz" "$a"
			else
				bcftools view -r "$cl:$lo-$hi" -Oz -o "$d/$a.vcf.gz" "$f"; tabix -f -p vcf "$d/$a.vcf.gz"
			fi
		done
	done
}

s4_core(){
	local t c
	log "Now running Step4 core VCF: input=$dirout/lead/pick.tsv; output=$dirout/coreVcf"
	if [[ $archaic_VCF_match == strict ]]; then
		build_matrices_strict
		log_process_s4_s5
		return 0
	fi
	for t in $traits; do for c in $chrs; do limit_jobs "$job_in_trait"; build_core_vcf_avcf "$t" "$c" & done; done
	wait_all
	log_process_s4_s5
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step5: convert core VCFs to matrix tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build_matrix_avcf(){
	check_cmd
	local t=$1 c=$2 cn d o a
	cn=$(chrn "$c")
	log "s5 mat $t chr$c"
	for d in "$(stage_trait_dir coreVcf "$t")"/${cn}.*.*; do
		[[ -d $d ]] || continue
		o=$(stage_trait_dir mat "$t")/$(basename "$d"); mkdir -p "$o"
		[[ -s $o/kg.tsv ]] || bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' "$d/kg.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{for(i=6;i<=NF;i++) if($i=="0" || $i=="1") $i=$i"/"$i; print}' > "$o/kg.tsv"
		for a in vindija altai chagyr denisova denisova25; do
			[[ -s $d/$a.vcf.gz && ! -s $o/$a.tsv ]] || continue
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$d/$a.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}{for(i=5;i<=NF;i++) if($i=="0" || $i=="1") $i=$i"/"$i; print}' > "$o/$a.tsv"
		done
	done
}

s5_mat(){
	local t c
	log "Now running Step5 matrix: input=$dirout/coreVcf; output=$dirout/mat"
	if [[ $archaic_VCF_match == strict ]]; then
		log "s5 skip: archaic_VCF_match=strict already wrote mat/<trait>/<locus> in s4"
		return 0
	fi
	for t in $traits; do for c in $chrs; do limit_jobs "$job_in_trait"; build_matrix_avcf "$t" "$c" & done; done
	wait_all
	log_process_s4_s5
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step6: identify inherited haplotypes and summarize reports
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s6_hap(){
	local t tr chr snp bp st en _n _size_bp id status=0 missing=0 matched_loci selected_loci yri_filter_note
	export PICK_LINEAGE_METHOD="$pick_lineage"
	export HAP_FILTER_METHOD="$hap_filter"
	export YRI_FREQ_TH="$yri_freq_th"
	export DIAGNOSTIC_DELTA_TH="$diagnostic_delta_th"
	export STRICT_P_ILS_TH="$strict_p_ils_th"
	export STRICT_MIN_COMPARED_RISK="$strict_min_compared_risk"
	export STRICT_MIN_MATCH_RISK="$strict_min_match_risk"
	export STRICT_PROP_MATCH_RISK_TH="$strict_prop_match_risk_th"
	export STRICT_HAP_MIN_N="$strict_hap_min_n"
	export STRICT_REQUIRE_CARRY_RISK="$strict_require_carry_risk"
	export STRICT_YRI_FREQ_TH="$strict_yri_freq_th"
	export SAMPLE_FILE="$sample_file"
	yri_filter_note="yri_freq_th=$yri_freq_th"
	awk -v th="$yri_freq_th" 'BEGIN{exit !(th+0 >= 1)}' && yri_filter_note="$yri_filter_note; disabled"
	log "Now running Step6 haplotype: input=$dirout/mat; output=$dirout/hap and $dir_report"
	log "process: Step6 YRI filter sample_file=$sample_file yri_freq_th=$yri_freq_th diagnostic_delta_th=$diagnostic_delta_th definition=diagnostic_archaic_allele_fixed_in_high_coverage_archaics_and_absent_or_rare_in_YRI"
	log "process: Step6 strict-style statistics only: p_ils<$strict_p_ils_th n_compared>=$strict_min_compared_risk n_match>=$strict_min_match_risk prop_match>=$strict_prop_match_risk_th hap_n>$strict_hap_min_n carry_risk_required=$strict_require_carry_risk strict_yri_freq_th=$strict_yri_freq_th"
	if valid_hap_stage; then
		log "s6 hap skip: report outputs already complete"
		write_hap_sample_map
		log_positive_loci_fate
		log_process_s6
		if [[ $hap_filter == avcf ]]; then
			matched_loci=$(data_rows "$dir_report/region_summary.tsv")
			selected_loci=$(data_rows "$dir_report/inherited_region.tsv")
			log "step6 finished, $matched_loci loci have matched haplotypes; $selected_loci loci remain after YRI frequency filtering ($yri_filter_note)."
		fi
		log "s5 done: region_summary $(nrow "$dir_report/region_summary.tsv") rows; inherited_region $(nrow "$dir_report/inherited_region.tsv") rows"
		return 0
	fi
	[[ -s "$dirscript/loci.R" ]] || { log "ERROR hap stage requires $dirscript/loci.R"; exit 1; }
	while IFS=$'\t' read -r tr chr snp bp st en _n _size_bp; do
		id="${chr}.${snp}.${bp}"
		[[ -d $(stage_trait_dir mat "$tr")/$id ]] || { log "WARN s6 missing mat for $tr $id; skip hap until mat is available"; continue; }
		if valid_hap_locus "$tr" "$id"; then
			log "s6 hap skip existing locus: $tr $id"
			continue
		fi
		missing=$((missing + 1))
		(
			log "s6 hap start locus: $tr $id"
			Rscript "$dirscript/loci.R" make_hap "$dirout" "$tr" "$id"
			mkdir -p "$(stage_trait_dir hap "$tr")"
			touch "$(stage_trait_dir hap "$tr")/$id.done"
			log "s6 hap done locus: $tr $id"
		) &
		while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do wait -n || status=1; done
	done < <(awk 'BEGIN{FS=OFS="\t"} NR>1 && !seen[$1 FS $2 FS $3 FS $4]++{print}' "$dirout/lead/pick.tsv")
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
	(( status == 0 )) || { log "ERROR one or more s6 hap locus jobs failed"; exit 1; }
	log "s6 hap locus pass complete: missing_or_incomplete=$missing"
	Rscript "$dirscript/loci.R" make_hap "$dirout" merge || { log "ERROR gu.R make_hap merge failed"; exit 1; }
	write_hap_sample_map
	log_positive_loci_fate
	log_process_s6
	if [[ $hap_filter == avcf ]]; then
		matched_loci=$(data_rows "$dir_report/region_summary.tsv")
		selected_loci=$(data_rows "$dir_report/inherited_region.tsv")
		log "step6 finished, $matched_loci loci have matched haplotypes; $selected_loci loci remain after YRI frequency filtering ($yri_filter_note)."
	fi
	if [[ $hap_filter == strict ]]; then
		local candidate_tsv hap_match_tsv filtered_hap_tsv inherited_tsv before_inherited_loci after_inherited_loci before_inherited_hap after_inherited_hap
		candidate_tsv="$dir_report/candidate_inherited_region.tsv"
		hap_match_tsv="$dir_report/hap_match.tsv"
		log_trait_count "candidate inherited loci" "$candidate_tsv"
		log "summary: loci confidence reference (file path $candidate_tsv; columns p_ils,best_lineage,matched_archaics)"
		log "s5 population filter start: pop=$filter_pop max_count=$filter_max_count"
		Rscript "$dirscript/loci.R" filter_hap "$dirout" "$sample_file" "$filter_pop" "$filter_max_count" || { log "ERROR gu.R filter_hap failed"; exit 1; }
		filtered_hap_tsv="$dir_report/filtered_hap_match.tsv"
		inherited_tsv="$dir_report/inherited_region.tsv"
		before_inherited_loci=$(data_rows "$candidate_tsv")
		after_inherited_loci=$(data_rows "$inherited_tsv")
		before_inherited_hap=$(data_rows "$hap_match_tsv")
		after_inherited_hap=$(data_rows "$filtered_hap_tsv")
		log_trait_count "inherited_region.tsv inherited loci" "$inherited_tsv"
		log "FILTER SUMMARY: inherited loci before_filter=$before_inherited_loci after_filter=$after_inherited_loci criteria=${filter_pop}<=$filter_max_count"
		log "FILTER SUMMARY: inherited haplotypes before_filter=$before_inherited_hap after_filter=$after_inherited_hap criteria=${filter_pop}<=$filter_max_count"
		log "FILTER SUMMARY: before-filter locus counts are in $dir_report/all.xlsx sheet inherited_loci_counts"
		log "FILTER SUMMARY: before-filter haplotype counts are in $dir_report/all.xlsx sheet inherited_haplotype_counts"
		log "FILTER SUMMARY: after-filter loci/haplotypes are in $dir_report/filtered.xlsx sheets filtered_loci and filtered_haplotypes"
	fi
	log "s5 done: region_summary $(nrow "$dir_report/region_summary.tsv") rows; inherited_region $(nrow "$dir_report/inherited_region.tsv") rows"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step7: build PHYLIP files, run phyml, and draw trees
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s7_phy(){
	check_cmd; check_phyml
	local f fs=() t id status=0 n miss_tree n_tree_png hapf
	log "Now running Step7 phylogeny: input=$dir_report/phy_region.tsv or inherited_region.tsv; output=$dirout/phy and $dir_plot"
	if [[ $phy_input == strict ]]; then
		if valid_phy_stage; then
			log "s7 phy skip: main PHYLIP, phyml trees, and tree PNGs already complete"
			log_process_s7
			return 0
		fi
		[[ -s "$dirscript/loci.R" ]] || { log "ERROR phy_input=strict requires $dirscript/loci.R"; exit 1; }
		hapf="$dir_report/filtered_hap_match.tsv"; [[ -s $hapf ]] || hapf="$dir_report/hap_match.tsv"
		[[ -s $hapf ]] || { log "ERROR no hap_match input for phy stage"; exit 1; }
		while IFS=$'\t' read -r t id; do
			[[ -n $t && -n $id ]] || continue
			if valid_phy_locus "$t" "$id"; then
				log "s7 make_phy skip existing locus: $t $id"
				continue
			fi
			( log "s7 make_phy start locus: $t $id"; Rscript "$dirscript/loci.R" make_phy "$dirout" "$t" "$id"; mkdir -p "$(stage_trait_dir phy "$t")"; touch "$(stage_trait_dir phy "$t")/$id.done"; log "s7 make_phy done locus: $t $id" ) &
			while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do wait -n || status=1; done
		done < <(awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="trait") it=i; if($i=="id") iid=i}; next} it&&iid{print $it,$iid}' "$hapf" | sort -u)
		while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
		(( status == 0 )) || { log "ERROR one or more s7 make_phy locus jobs failed"; exit 1; }
		local n_main_phy n_full_phy
		n_main_phy=$(count_files "$dirout/phy" '*.main.phy')
		n_full_phy=$(count_files "$dirout/phy" '*.full.phy')
		log "summary: *.main.phy filtered haplotypes (file path $dirout/phy) $n_main_phy"
		log "summary: *.full.phy all phylogeny haplotypes (file path $dirout/phy) $n_full_phy"
		[[ $hap_filter == strict ]] && log "summary: filtering reference (workbook $dir_report/filtered.xlsx; criteria $filter_pop <= $filter_max_count)"
		clean_report_workbooks_only
		if (( n_main_phy + n_full_phy == 0 )); then
			log "WARN no loci available for phylogeny after filtering; no PHYLIP files or tree PNGs will be generated."
			log_process_s7
			log "summary: report workbooks $dir_report/all.xlsx $dir_report/filtered.xlsx $dir_report/selected.xlsx"
			return 0
		fi
		status=0
		local phy_trait_jobs=$(( max_cores / (job_phyml * phyml_cpus) ))
		(( phy_trait_jobs < 1 )) && phy_trait_jobs=1
		for t in $traits; do
			(
				n=$(find "$(stage_trait_dir phy "$t")" \( -name '*.full.phy' -o -name '*.main.phy' \) 2>/dev/null | wc -l)
				log "s6 phyml trait start: $t ($n full/main .phy files)"
				while IFS= read -r f; do
					(
						if [[ -s "${f}_phyml_tree.txt" ]]; then log "s6 skip existing phyml tree: $f"; exit 0; fi
						run_phyml_file "$f" 0 || exit 1
					) &
					while [[ $(jobs -rp | wc -l) -ge $job_phyml ]]; do wait -n || exit 1; done
				done < <(find "$(stage_trait_dir phy "$t")" \( -name '*.full.phy' -o -name '*.main.phy' \) 2>/dev/null | sort)
				while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || exit 1; done
				log "s6 phyml trait done: $t"
			) &
			while [[ $(jobs -rp | wc -l) -ge $phy_trait_jobs ]]; do wait -n || status=1; done
		done
		while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
		(( status == 0 )) || { log "ERROR one or more s6 phyml trait jobs failed or timed out; not marking trait complete"; exit 1; }
		miss_tree=$(find "$dirout/phy" \( -name '*.full.phy' -o -name '*.main.phy' \) | while read -r f; do [[ -s "${f}_phyml_tree.txt" ]] || echo "$f"; done | head -1)
		if [[ -n $miss_tree ]]; then
			log "ERROR some full/main .phy files still lack trees after Step7. First missing: $miss_tree"
			exit 1
		fi
		MAKE_TREE_KEEP_EXISTING=0 MAKE_TREE_TAGS=${make_tree_tags:-main} MAKE_TREE_PPTX=0 Rscript "$dirscript/loci.R" make_tree "$dirout" || { log "ERROR gu.R make_tree failed"; exit 1; }
		n_tree_png=$(find "$dir_plot" -maxdepth 1 \( -name 's8_tree_full_*.png' -o -name 's8_tree_main_*.png' \) | wc -l)
		log "summary: loci phylogeny tree PNG (source *.full.phy and *.main.phy; file path $dir_plot) $n_tree_png"
		[[ $n_tree_png -gt 0 ]] || { log "ERROR no tree PNG generated in $dir_plot"; exit 1; }
		log_process_s7
		clean_report_workbooks_only
		log "summary: report workbooks $dir_report/all.xlsx $dir_report/filtered.xlsx $dir_report/selected.xlsx"
		return 0
	else
		MAKE_PHY_METHOD=avcf Rscript "$dirscript/loci.R" make_phy "$dirout"
		if [[ $phyml_run_full == 1 || $phyml_full_required == 1 ]]; then
			mapfile -t fs < <(find "$dirout/phy" \( -name "*.full.phy" -o -name "*.main.phy" \) | sort)
		else
			mapfile -t fs < <(find "$dirout/phy" -name "*.main.phy" | sort)
		fi
	fi
	if [[ ${#fs[@]} -eq 0 ]]; then
		if [[ -s "$dir_report/phy_region.tsv" && $(data_rows "$dir_report/phy_region.tsv") -eq 0 ]]; then
			log "WARN no loci available for phylogeny; phy_region.tsv has 0 rows. This is allowed when no inherited loci are detected and no sanity BED is provided."
			return 0
		fi
		log "ERROR no .phy files generated; inspect $dir_report/phy_region.tsv and $dir_report/hap_match.tsv"
		exit 1
	fi
	for f in "${fs[@]}"; do
		limit_jobs "$job_phyml"
		(
			if run_phyml_file "$f" 1; then
				exit 0
			fi
			if [[ ${f##*/} == *.full.phy && $phyml_full_required != 1 ]]; then
				log "WARN optional full.phy tree failed; continuing because phyml_full_required=$phyml_full_required: $f"
				exit 0
			fi
			exit 1
		) &
	done
	wait_all || { log "ERROR one or more phyml jobs failed; not marking trait complete"; exit 1; }
	if [[ $phyml_full_required == 1 ]]; then
		miss_tree=$(find "$dirout/phy" \( -name '*.full.phy' -o -name '*.main.phy' \) | while read -r f; do [[ -s "${f}_phyml_tree.txt" ]] || echo "$f"; done | head -1)
	else
		miss_tree=$(find "$dirout/phy" -name '*.main.phy' | while read -r f; do [[ -s "${f}_phyml_tree.txt" ]] || echo "$f"; done | head -1)
	fi
	if [[ -n $miss_tree ]]; then
		log "ERROR missing phyml tree after Step7: $miss_tree"
		exit 1
	fi
	MAKE_TREE_KEEP_EXISTING=0 MAKE_TREE_TAGS=${make_tree_tags:-main} MAKE_TREE_PPTX=0 Rscript "$dirscript/loci.R" make_tree "$dirout" || { log "ERROR gu.R make_tree failed"; exit 1; }
	local n_tree_png
	n_tree_png=$(find "$dir_plot" -maxdepth 1 -name 's8_tree_*.png' 2>/dev/null | wc -l)
	log "summary: loci phylogeny tree PNG (file path $dir_plot) $n_tree_png"
	[[ $n_tree_png -gt 0 ]] || { log "ERROR no tree PNG generated in $dir_plot"; exit 1; }
	log_process_s7
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step8: extra QC/sensitivity tables for manuscript validation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s8_extra(){
	if [[ $method != avcf ]]; then
		log "WARN extra QC is currently designed for avcf outputs; method=$method, skip extra QC"
		return 0
	fi
	if [[ -s "$extra_outdir/locus_screening_summary.tsv" &&
		-s "$extra_outdir/candidate_quality_table.tsv" &&
		-s "$extra_outdir/filter_sensitivity_grid.tsv" &&
		-s "$extra_outdir/empirical_background.tsv" ]]; then
		log "process: extra QC skip; complete sensitivity tables already exist"
		return 0
	fi
	mkdir -p "$extra_outdir"
	log "s8 extra QC start: root=$dirout_root traits=$configured_traits out=$extra_outdir"
	Rscript "$dirscript/loci.R" extra_qc \
		--dirout_root "$dirout_root" \
		--traits "$configured_traits" \
		--outdir "$extra_outdir" \
		--pils_grid "$pils_grid" \
		--freq_grid "$freq_grid" \
		--diag_n_grid "$diag_n_grid" \
		--lineage_prop_grid "$lineage_prop_grid" \
		--hap_n_grid "$hap_n_grid" || { log "ERROR gu.R extra_qc failed"; exit 1; }
	log "summary: extra locus screening table $(nrow "$extra_outdir/locus_screening_summary.tsv") rows; file=$extra_outdir/locus_screening_summary.tsv"
	log "summary: candidate quality table $(nrow "$extra_outdir/candidate_quality_table.tsv") rows; file=$extra_outdir/candidate_quality_table.tsv"
	log "summary: filter sensitivity grid $(nrow "$extra_outdir/filter_sensitivity_grid.tsv") rows; file=$extra_outdir/filter_sensitivity_grid.tsv"
	log "summary: empirical background table $(nrow "$extra_outdir/empirical_background.tsv") rows; file=$extra_outdir/empirical_background.tsv"
	log "summary: diagnostic variants for UKB proxy $(nrow "$extra_outdir/diagnostic_variants_for_haplotype_proxy.tsv") rows; file=$extra_outdir/diagnostic_variants_for_haplotype_proxy.tsv"
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step9: optional UKB haplotype-proxy association
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s9_proxy(){
	if [[ $method != avcf ]]; then
		log "ERROR proxy association is currently designed for avcf outputs; method=$method"
		exit 1
	fi
	[[ -s "$diag_tsv" ]] || { log "WARN missing diagnostic table $diag_tsv; running s8_extra first"; s8_extra; }
	[[ -s "$diag_tsv" ]] || { log "ERROR missing diagnostic table after s8_extra: $diag_tsv"; exit 1; }
	[[ -s "$pheno_cov" ]] || { log "ERROR missing pheno_cov=$pheno_cov"; exit 1; }
	command -v "$plink2" >/dev/null 2>&1 || { log "ERROR missing plink2 executable: $plink2"; exit 1; }
	mkdir -p "$proxy_outdir/scores" "$proxy_outdir/log"
	log "s9 proxy association start: diag=$diag_tsv pheno_cov=$pheno_cov out=$proxy_outdir"
	Rscript "$dirscript/loci.R" proxy_assoc make_scores \
		--diag_tsv "$diag_tsv" \
		--outdir "$proxy_outdir/scores" || { log "ERROR gu.R proxy_assoc make_scores failed"; exit 1; }

	local score_manifest="$proxy_outdir/scores/score_manifest.tsv"
	local sscore_manifest="$proxy_outdir/sscore_files.tsv"
	[[ -s "$score_manifest" ]] || { log "ERROR missing proxy score manifest: $score_manifest"; exit 1; }
	printf "trait\tid\tchr\tsscore\tlead_snp\tlead_bp\n" > "$sscore_manifest"
	while IFS=$'\t' read -r trait id chr score_file lead_snp lead_bp n_score_variants; do
		[[ $trait == trait ]] && continue
		local prefix=${pfile_tpl/\{CHR\}/$chr}
		if [[ ! -s ${prefix}.pgen || ! -s ${prefix}.pvar || ! -s ${prefix}.psam ]]; then
			log "WARN missing pfile for chr$chr: $prefix; skip $trait $id"
			continue
		fi
		local safe_id
		safe_id=$(echo "${trait}_${id}" | tr -c 'A-Za-z0-9_.-' '_')
		local outpref="$proxy_outdir/scores/$safe_id"
		"$plink2" --pfile "$prefix" \
			--set-all-var-ids '@:#:$r:$a' \
			--score "$score_file" 1 2 3 header-read cols=+scoresums \
			--out "$outpref" > "$proxy_outdir/log/$safe_id.plink2.log" 2>&1 || { log "WARN plink2 score failed for $trait $id; see $proxy_outdir/log/$safe_id.plink2.log"; continue; }
		[[ -s "$outpref.sscore" ]] || { log "WARN missing plink2 sscore for $trait $id: $outpref.sscore"; continue; }
		printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$trait" "$id" "$chr" "$outpref.sscore" "$lead_snp" "$lead_bp" >> "$sscore_manifest"
	done < "$score_manifest"
	[[ $(nrow "$sscore_manifest") -gt 0 ]] || { log "ERROR no proxy .sscore files were generated; inspect $proxy_outdir/log"; exit 1; }
	Rscript "$dirscript/loci.R" proxy_assoc fit_models \
		--sscore_manifest "$sscore_manifest" \
		--pheno_cov "$pheno_cov" \
		--id_col "$id_col" \
		--outcomes "$outcomes" \
		--binary_outcomes "$binary_outcomes" \
		--covars "$covars" \
		--outdir "$proxy_outdir" || { log "ERROR gu.R proxy_assoc fit_models failed"; exit 1; }
	log "summary: UKB haplotype-proxy association $(nrow "$proxy_outdir/haplotype_proxy_assoc.tsv") rows; file=$proxy_outdir/haplotype_proxy_assoc.tsv"
}

clean_all(){
	find "$dirout" -mindepth 1 ! -name 'gu.log' -exec rm -rf {} + 2>/dev/null || true
	init_output_dirs
	log "summary: cleaned analysis directory; output=$dirout"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 asnp workflow: reference-based aSNP / archaic haplotype map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://github.com/SillySabertooth/DNA_virus_load_Neand_Intrgoression
run_asnp_current_trait(){
	local mode=${1:-all}
	log_config
	if [[ $mode == prep ]]; then
		s1_prep
		return 0
	fi
	# asnp starts from the same lead-locus/risk-allele preparation, then switches to
	# the 2026 GBE-style reference aSNP / archaic haplotype-map scan.
	s1_prep
	# Generate the same lead-SNP LD summaries when possible; asnp uses these to
	# label whether a candidate aSNP is in high LD with the GWAS lead SNP.
	# If a locus lacks LD output, gu.R still produces asnp reports with LD marked NA.
	s2_ld
	s3_pick
	if [[ $mode == ld ]]; then
		log "summary: asnp LD/prep complete; pick=$dirout/lead/pick.tsv"
		return 0
	fi
	mkdir -p "$dir_report" "$dir_plot" "$dirout/asnp"
	[[ -s "$dirscript/loci.R" ]] || { log "ERROR asnp requires $dirscript/loci.R"; exit 1; }
	if [[ ! -s $asnp_asnp ]]; then
		log "ERROR asnp aSNP haplotype map not found: $asnp_asnp"
		log "Put aSNPs.haplotypes.v1.tsv under ${dir_archaic:-$dir_ref/archaic/GRCH${GRCH:-37}}/asnp, or export asnp_asnp before running."
		exit 1
	fi
	log "Now running asnp workflow: reference aSNP/haplotype-map scan and optional haplotype-network output"
	Rscript "$dirscript/loci.R" asnp \
		--dirgwas "$dirgwas" \
		--dircojo "$dircojo" \
		--dirout "$dirout" \
		--dirmod "$dirmod" \
		--asnp "$asnp_asnp" \
		--traits "${traits_arr[@]}" \
		--p-th "$asnp_p_th" \
		--freq-th "$asnp_freq_th" \
		--min-asnp "$asnp_min_asnp" \
		--window-kb "$asnp_window_kb" \
		--ld-r2 "$asnp_ld_r2" \
		--make-network "$asnp_make_network" || { log "ERROR gu.R asnp failed"; exit 1; }
	log "summary: asnp aSNP hits $(nrow "$dir_report/asnp_aSNP_hits.tsv") rows; file=$dir_report/asnp_aSNP_hits.tsv"
	log "summary: asnp candidate inherited segments $(nrow "$dir_report/asnp_inherited_segments.tsv") rows; file=$dir_report/asnp_inherited_segments.tsv"
	log "summary: asnp region summary $(nrow "$dir_report/asnp_region_summary.tsv") rows; file=$dir_report/asnp_region_summary.tsv"
	log "summary: asnp network inputs/plots are under $dirout/asnp and $dir_plot"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Main workflow and command dispatch
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_current_trait_all(){
	if is_asnp; then
		run_asnp_current_trait all
		return 0
	fi
	log_config
	if [[ $start_step == s1 ]]; then
		s1_prep
		s2_ld
		s3_pick
	elif [[ $start_step == s2 ]]; then
		for f in "$dirout/lead/pick.tsv" "$dirout/lead/lead.assoc"; do [[ -s "$f" ]] || { log "ERROR START_STEP=s2 requires existing file: $f"; exit 1; }; done
		[[ $(data_rows "$dirout/lead/pick.tsv") -gt 0 ]] || { log "ERROR START_STEP=s2 requires non-empty $dirout/lead/pick.tsv; rerun START_STEP=s1"; exit 1; }
		log "Skip s1/s2 LD; reuse existing $dirout/lead/pick.tsv"
	else
		log "ERROR unknown START_STEP=$start_step; use s1 or s2"
		exit 1
	fi
	s4_core
	s5_mat
	s6_hap
	s7_phy
	log_trait_count "candidate high-LD loci for haplotype analysis" "$dirout/lead/pick.tsv"
	log_trait_count "single-SNP lead loci excluded before haplotype analysis" "$dirout/lead/pick.single.tsv"
	log "summary: trait complete $traits; output=$dirout"
	touch "$dirout/.loci_avcf.complete"
}

run_current_trait_mode(){
	local mode=$1
	if is_asnp; then
		case "$mode" in
			all|ld|core|mat|hap|phy) run_asnp_current_trait "$mode" ;;
			prep) log_config; s1_prep ;;
			extra) log "WARN extra QC is currently designed for avcf outputs; method=$method, skip" ;;
			proxy) log "ERROR proxy association is currently designed for avcf outputs; method=$method"; exit 1 ;;
			*) echo "ERROR: unknown mode: $mode" >&2; usage >&2; exit 1 ;;
		esac
		return 0
	fi
	case "$mode" in
		all) run_current_trait_all ;;
		prep) log_config; s1_prep ;;
		ld) log_config; s2_ld; s3_pick ;;
		core) log_config; s4_core ;;
		mat) log_config; s5_mat ;;
		hap) log_config; s6_hap ;;
		phy) log_config; s7_phy ;;
		extra) s8_extra ;;
		proxy) s9_proxy ;;
		*) echo "ERROR: unknown mode: $mode" >&2; usage >&2; exit 1 ;;
	esac
}

run_workflow(){
	local mode=$1
	init_output_dirs
	if [[ $mode == all &&
		-s "$dirout/lead/lead.assoc" &&
		-s "$dirout/lead/pick.tsv" ]] &&
		valid_hap_stage && valid_phy_stage; then
		log "process: loci_avcf skip trait=$traits; lead, haplotype, report, and phylogeny outputs are complete"
		log "ALL DONE: traits=$traits; reused complete output=$dirout"
		return 0
	fi
	run_current_trait_mode "$mode"
	log "ALL DONE: traits=$traits; output=$dirout"
}

run_trait_mode(){
	local t=$1 mode=$2
	set_trait_context "$t"
	open_trait_log
	if [[ -e "$dirout/.loci_avcf.complete" &&
		-s "$dirout/lead/lead.assoc" &&
		-s "$dirout/lead/pick.tsv" &&
		-s "$dir_report/hap_match.tsv" &&
		-s "$dir_report/region_summary.tsv" ]]; then
		log "process: loci_avcf skip trait=$t action=$mode; completion marker and key outputs are valid"
		return 0
	fi
	if [[ -s "$dirout/lead/lead.assoc" &&
		-s "$dirout/lead/pick.tsv" ]] &&
		valid_hap_stage && valid_phy_stage; then
		log "process: loci_avcf skip trait=$t action=$mode; complete outputs already exist"
		touch "$dirout/.loci_avcf.complete"
		return 0
	fi
	case "$mode" in
		all) run_workflow all ;;
		clean) clean_all ;;
		prep|ld|core|mat|hap|phy) run_workflow "$mode" ;;
		*) echo "ERROR: unknown mode: $mode" >&2; usage >&2; exit 1 ;;
	esac
}

run_traits_mode_parallel(){
	local mode=$1 t status=0 status_dir failed
	status_dir="$dirout_root/.trait_status.${mode}.$$"
	mkdir -p "$status_dir"
	for t in "${configured_traits_arr[@]}"; do
		(
			trait_job=$t
			mode_job=$mode
			status_file="$status_dir/${t}.status"
			trap 'rc=$?; printf "%s\t%s\n" "$trait_job" "$rc" > "$status_file"; if (( rc != 0 )); then log "ERROR trait job failed: trait=$trait_job mode=$mode_job rc=$rc"; fi' EXIT
			run_trait_mode "$trait_job" "$mode_job"
		) &
		while [[ $(jobs -rp | wc -l) -ge $job_of_trait ]]; do wait -n || status=1; done
	done
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
	if (( status != 0 )); then
		failed=$(awk -F '\t' '$2 != 0{printf "%s%s(rc=%s)", sep, $1, $2; sep=", "}' "$status_dir"/*.status 2>/dev/null || true)
		[[ -n $failed ]] || failed="unknown"
		rm -f "$status_dir"/*.status 2>/dev/null || true
		rmdir "$status_dir" 2>/dev/null || true
		log "ERROR one or more trait jobs failed in mode=$mode: $failed"
		exit 1
	fi
	rm -f "$status_dir"/*.status 2>/dev/null || true
	rmdir "$status_dir" 2>/dev/null || true
}

run_all_traits_staged(){
	local stage
	for stage in prep ld core mat hap phy; do
		log "STAGE START: $stage for traits=${configured_traits_arr[*]} job_of_trait=$job_of_trait job_in_trait=$job_in_trait"
		run_traits_mode_parallel "$stage"
		log "STAGE DONE: $stage for traits=${configured_traits_arr[*]}"
	done
	log "STAGE START: extra QC for traits=${configured_traits_arr[*]}"
	s8_extra
	log "STAGE DONE: extra QC"
	for stage in "${configured_traits_arr[@]}"; do
		[[ -s "$dirout_root/$stage/lead/lead.assoc" &&
			-s "$dirout_root/$stage/lead/pick.tsv" &&
			-s "$dirout_root/$stage/report/hap_match.tsv" &&
			-s "$dirout_root/$stage/report/region_summary.tsv" ]] &&
			touch "$dirout_root/$stage/.loci_avcf.complete"
	done
}


mode=all
[[ $# -gt 0 ]] && mode=$1
case "$mode" in
	all)
		if is_asnp; then
			run_traits_mode_parallel all
		else
			run_all_traits_staged
		fi
		;;
	clean|prep|ld|core|mat|hap|phy)
		run_traits_mode_parallel "$mode"
		;;
	extra)
		s8_extra
		;;
	proxy)
		s9_proxy
		;;
	*) echo "ERROR: unknown mode: $mode" >&2; usage >&2; exit 1 ;;
esac
