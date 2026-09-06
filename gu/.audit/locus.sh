#!/usr/bin/env bash
set -euo pipefail


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Paths, traits, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dirscript="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$dirscript/comm.inc"
export PATH="$dir0/software/bin:$PATH"
dirgwas=$dir0/data/gwas/bald
dircojo=$dirgwas/cojo

method=simple # 📍 simple  | strict | viome 
dirout=$dir0/analysis/archaic/result/locus; mkdir -p "$dirout"

positive_loci=$dir0/files/gu.loci.bed
add_positive_loci=auto
traits="bald bald12 bald13 bald14"
chrs="$(seq 1 22) X"
ld_r2=0.98


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Method switches
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ref_pop=ALL                 # 📍 ALL | EUR
lead_match=${lead_match:-$method}                 # 📍 simple | strict | viome
ld_calc=${ld_calc:-$method}                       # 📍 simple | strict | viome
archaic_VCF_match=${archaic_VCF_match:-$method}   # 📍 simple | strict | viome
pick_lineage=${pick_lineage:-$method}             # 📍 simple | strict | viome
hap_filter=${hap_filter:-$method}                 # 📍 simple | strict | viome
yri_freq_th=0.05            # 📍 simple hap_filter: keep loci with YRI_max_freq <= threshold; >=1 disables this filter
diagnostic_delta_th=0.5     # 📍 diagnostic archaic allele must be enriched on risk haplotypes by this carry-noncarry frequency delta
phy_input=${phy_input:-$method}                   # 📍 simple | strict | viome

# viome = reference-based aSNP / archaic haplotype map workflow inspired by Rajpara et al. 2026 GBE.
# Override these paths/thresholds here or by environment variables before running bash locus.sh.
if [[ $method == viome || ${lead_match:-} == viome || ${ld_calc:-} == viome || ${archaic_VCF_match:-} == viome || ${pick_lineage:-} == viome || ${hap_filter:-} == viome || ${phy_input:-} == viome ]]; then
	viome_asnp=${viome_asnp:-$dir0/files/aSNPs.haplotypes.v1.tsv}   # Yermakovich-style aSNP haplotype map
	viome_p_th=${viome_p_th:-1e-9}                                  # default threshold used in the 2026 GBE paper
	viome_freq_th=${viome_freq_th:-0.01}                            # archaic allele/haplotype frequency threshold
	viome_min_asnp=${viome_min_asnp:-5}                             # require haplotypes with >= this many aSNPs
	viome_window_kb=${viome_window_kb:-1000}                         # +/- window around lead SNP for regional scan
	viome_ld_r2=${viome_ld_r2:-0.9}                                  # high-LD threshold for network candidate region
	viome_make_network=${viome_make_network:-1}                      # 1: try pegas network if 1000G VCF + packages are available
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Parallelism and runtime limits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
max_cores=16
job_of_trait=1
job_in_trait=12
job_phyml=1
plink_threads=1
internal_threads=1
phyml_cpus=16
phyml_retry_cpus=1
phyml_fallback=1
phyml_boot=100
phyml_timeout=8h
filter_pop=YRI
filter_max_count=1
start_step=s1

read -r -a traits_arr <<< "$traits"
mkdir -p "$dirout"
exec > >(tee "$dirout/locus.log") 2>&1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Logging and small utilities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log(){
	local msg="$*" tag=""
	case "$msg" in
		traits=*|methods:*|START_STEP=*|haplotype\ population*|parallel\ cap:*|paths:*) return 0 ;;
		*"skip existing locus"*|*"s2 skip "*|*"s2 simple-LD "*|*"pick refresh:"*|LD\ *\ chr*|*"s4 core "*|*"s5 mat "*|*"phyml skip existing tree"*|*"skip existing phyml tree"*|*"make_phy skip existing locus"*|*"s4 skip existing mat"*) return 0 ;;
	esac
	case "$msg" in
		ERROR*|FAIL*) tag=$'\U0001F6D1 STOP ' ;;
		WARN*) tag=$'\u26A0 ' ;;
		summary:*|FILTER\ SUMMARY:*|ALL\ DONE|*" done:"*|*" done;"*|s1\ done:*|s2\ LD\ done*|s3\ done:*|s3\ force*|s5\ done:*|*" plotted "*|*"identified "*|*"Added positive"*) tag=$'\u2733 ' ;;
	esac
	echo "[$(date '+%F %T')] ${tag}${msg}"
}
nrow(){ data_rows "$1"; }
ok(){ [[ -s $1 && $(nrow "$1") -gt 0 ]]; }
chrn(){ [[ $1 == X ]] && echo 23 || echo "$1"; }
chrl(){ [[ $1 == 23 || $1 == X ]] && echo X || echo "$1"; }
clean_msg(){ [[ -f $1 ]] && tail -n 12 "$1" | tr '\t\r\n' '   ' | sed 's/  */ /g' | cut -c1-700; }
export OMP_NUM_THREADS=$internal_threads
export OPENBLAS_NUM_THREADS=$internal_threads
export MKL_NUM_THREADS=$internal_threads
export VECLIB_MAXIMUM_THREADS=$internal_threads
export NUMEXPR_NUM_THREADS=$internal_threads
export R_DATATABLE_NUM_THREADS=$internal_threads

is_viome(){
	[[ " $lead_match $ld_calc $archaic_VCF_match $pick_lineage $hap_filter $phy_input " == *" viome "* ]]
}

validate_methods(){
	local x v ok=1
	case "$method" in simple|strict|viome) ;; *) log "ERROR invalid method=$method; use simple | strict | viome"; ok=0 ;; esac
	for x in lead_match ld_calc archaic_VCF_match pick_lineage hap_filter phy_input; do
		eval "v=\${$x}"
		case "$v" in simple|strict|viome) ;; *) log "ERROR invalid $x=$v; use simple | strict | viome"; ok=0 ;; esac
	done
	(( ok == 1 )) || exit 1
}

log_config(){
	validate_methods
	log "Now running locus workflow; traits=$traits; output=$dirout; max_cores=$max_cores"
	log "Methods: method=$method, lead_match=$lead_match, ld_calc=$ld_calc, archaic_VCF_match=$archaic_VCF_match, pick_lineage=$pick_lineage, hap_filter=$hap_filter, yri_freq_th=$yri_freq_th, diagnostic_delta_th=$diagnostic_delta_th, phy_input=$phy_input"
	if is_viome; then
		log "Viome: asnp=$viome_asnp p_th=$viome_p_th freq_th=$viome_freq_th min_asnp=$viome_min_asnp window_kb=$viome_window_kb ld_r2=$viome_ld_r2 make_network=$viome_make_network"
	fi
	log "Inputs: GWAS=$dirgwas; COJO=$dircojo; 1000G=$dirmod; pfile=$dirmod/pfile; vcf=$dirmod/vcf; archaic=$arch0; sanity_loci=$positive_loci"
}

init_output_dirs(){
	mkdir -p "$dirout"/{lead,ld,coreVcf,mat,hap,phy,plot,report}
}

id_mode(){ awk 'BEGIN{IGNORECASE=1} $1!="" && $1!="."{n++; if($1~/^rs[0-9]+$/) rs++; else if(toupper($1)~/^(CHR)?([0-9]+|X|Y|MT|M):[0-9]+:[ACGTN]+:[ACGTN,]+$/) cp++; if(n==5000) exit} END{if(n==0) print "missing"; else if(rs/n>.8) print "rsid"; else if(cp/n>.8) print "chrpos"; else print "other"}'; }
trait_rows(){ [[ -s $1 ]] && awk -v t="$2" 'BEGIN{FS="\t"} NR>1 && $1==t{n++} END{print n+0}' "$1" || echo 0; }
count_files(){ find "$1" -name "$2" 2>/dev/null | wc -l; }
valid_hap_locus(){
	local t=$1 id=$2 d
	d=$dirout/hap/$t
	[[ -s $d/$id.region.tsv && -s $d/$id.core.tsv ]] || [[ -s $d/$id.done ]] || return 1
	if [[ $hap_filter == simple ]]; then
		awk 'NR==1{for(i=1;i<=NF;i++) if($i=="is_diagnostic_archaic") ok=1; exit !ok}' "$d/$id.core.tsv" || return 1
		awk 'NR==1{y=0; for(i=1;i<=NF;i++) if($i=="yri_risk_freq") y=i; if(!y) exit 1; next} y && $y!="" && $y!="NA"{ok=1} END{exit !ok}' "$d/$id.core.tsv"
	fi
}
valid_hap_stage(){
	local tr chr snp bp st en _n _size_bp id
	[[ -s $dirout/report/hap_match.tsv && -s $dirout/report/region_summary.tsv ]] || return 1
	while IFS=$'\t' read -r tr chr snp bp st en _n _size_bp; do
		id="${chr}.${snp}.${bp}"
		valid_hap_locus "$tr" "$id" || return 1
	done < <(awk 'BEGIN{FS=OFS="\t"} NR>1 && !seen[$1 FS $2 FS $3 FS $4]++{print}' "$dirout/lead/pick.tsv")
	if [[ $hap_filter == simple ]]; then
		[[ -s $dirout/report/selected_region.tsv ]] || return 1
		[[ -s $dirout/report/hap_site_count.tsv && -s $dirout/report/core_archaic_match.tsv && -s $dirout/report/core_risk.tsv ]] || return 1
		[[ -s $dirout/report/yri_filter.tsv ]] || return 1
		awk -v th="$yri_freq_th" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="yri_freq_th") c=i; next} c && sprintf("%.12g",$c+0)==sprintf("%.12g",th+0){ok=1} END{exit !ok}' "$dirout/report/yri_filter.tsv" || return 1
		awk -v th="$diagnostic_delta_th" 'BEGIN{FS="\t"; ok=0} NR==1{for(i=1;i<=NF;i++) if($i=="diagnostic_delta_th") c=i; next} c && sprintf("%.12g",$c+0)==sprintf("%.12g",th+0){ok=1} END{exit !ok}' "$dirout/report/yri_filter.tsv" || return 1
		awk 'NR==1{for(i=1;i<=NF;i++) if($i=="is_diagnostic_archaic") ok=1; exit !ok}' "$dirout/report/core_risk.tsv" || return 1
		awk 'NR==1{y=0; for(i=1;i<=NF;i++) if($i=="yri_risk_freq") y=i; if(!y) exit 1; next} y && $y!="" && $y!="NA"{ok=1} END{exit !ok}' "$dirout/report/core_risk.tsv" || return 1
		awk 'NR==1{for(i=1;i<=NF;i++){if($i=="trait") t=1; if($i=="id") id=1; if($i=="lead_snp") snp=1; if($i=="yri_freq") y=1} exit !(t&&id&&snp&&y)}' "$dirout/report/selected_region.tsv" || return 1
	else
		[[ -s $dirout/report/inherited_segments.tsv ]] || return 1
		[[ -s $dirout/report/filtered_hap_match.tsv && -s $dirout/report/filtered_inherited_segments.tsv ]] || return 1
		[[ -s $dirout/report/all.xlsx && -s $dirout/report/filtered.xlsx && -s $dirout/report/selected.xlsx ]] || return 1
	fi
}
valid_phy_file(){
	local f=$1
	[[ -s $f ]] || return 1
	awk 'NR==1{next} NF && substr($0,11,1)!=" "{bad=1; exit} END{exit bad}' "$f"
}
valid_phy_locus(){
	local t=$1 id=$2 f=$dirout/phy/$t/$id.main.phy
	if [[ -e $f ]]; then
		[[ -s $f && -s $dirout/phy/$t/$id.main.meta.tsv ]] && valid_phy_file "$f"
	else
		[[ -s $dirout/phy/$t/$id.done ]]
	fi
}
valid_phy_stage(){
	local miss invalid n_phy n_png
	n_phy=$(find "$dirout/phy" -name '*.main.phy' 2>/dev/null | wc -l)
	(( n_phy > 0 )) || return 1
	invalid=$(find "$dirout/phy" -name '*.main.phy' 2>/dev/null | while read -r f; do valid_phy_file "$f" || echo "$f"; done | head -1)
	[[ -z $invalid ]] || return 1
	miss=$(find "$dirout/phy" -name '*.main.phy' 2>/dev/null | while read -r f; do [[ -s "${f}_phyml_tree.txt" ]] || echo "$f"; done | head -1)
	[[ -z $miss ]] || return 1
	n_png=$(find "$dirout/plot" -maxdepth 1 -name 's8_tree_main_*.png' 2>/dev/null | wc -l)
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
		nlead=$(data_rows "$dirout/lead/$t.lead.assoc")
		log "process: Step1 trait=$t COJO_file=$f COJO_loci=$ncojo lead_assoc_loci=$nlead output=$dirout/lead/$t.lead.assoc"
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
	ncore=$(find "$dirout/coreVcf" -mindepth 2 -maxdepth 3 -type d 2>/dev/null | wc -l)
	nmat=$(find "$dirout/mat" -mindepth 2 -maxdepth 2 -type d 2>/dev/null | wc -l)
	log "process: Step4 core VCF input_loci=$npick core_locus_dirs=$ncore output=$dirout/coreVcf"
	log "process: Step5 matrix input_loci=$npick matrix_locus_dirs=$nmat output=$dirout/mat"
	log "FLOW: step=4 name=build_core_VCF input_loci=$npick core_locus_dirs=$ncore"
	log "FLOW: step=5 name=build_genotype_matrices input_loci=$npick matrix_locus_dirs=$nmat"
}
log_process_s6(){
	local nhap nregion nselected nphy nmap nfate ymethod yth dth ndiag
	nhap=$(data_rows "$dirout/report/hap_match.tsv")
	nregion=$(data_rows "$dirout/report/region_summary.tsv")
	nselected=$(data_rows "$dirout/report/selected_region.tsv")
	nphy=$(data_rows "$dirout/report/phy_region.tsv")
	nmap=$(data_rows "$dirout/report/haplotype_sample_map.tsv")
	nfate=$(data_rows "$dirout/report/positive_loci_fate.tsv")
	ymethod=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="yri_filter_method") c=i; next} c{print $c; exit}' "$dirout/report/yri_filter.tsv" 2>/dev/null || true)
	yth=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="yri_freq_th") c=i; next} c{print $c; exit}' "$dirout/report/yri_filter.tsv" 2>/dev/null || true)
	dth=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="diagnostic_delta_th") c=i; next} c{print $c; exit}' "$dirout/report/yri_filter.tsv" 2>/dev/null || true)
	ndiag=$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="n_yri_filter_sites") c=i; next} c && NR>1{s+=$c} END{print s+0}' "$dirout/report/region_summary.tsv" 2>/dev/null || echo 0)
	log "process: Step6 haplotypes hap_rows=$nhap inherited_candidate_loci=$nregion selected_loci_after_YRI_filter=$nselected hap_match=$dirout/report/hap_match.tsv selected=$dirout/report/selected_region.tsv"
	log "process: Step6 YRI diagnostic filter method=${ymethod:-NA} yri_freq_th=${yth:-NA} diagnostic_delta_th=${dth:-NA} diagnostic_sites=$ndiag summary=$dirout/report/yri_filter.tsv"
	log "process: Step6 Roman hap labels mapped_to_1KG_copies=$nmap output=$dirout/report/haplotype_sample_map.tsv"
	log "process: Step6 gu.loci.bed fate_rows=$nfate output=$dirout/report/positive_loci_fate.tsv"
	log "process: Step7 phy input_loci=$nphy input=$dirout/report/phy_region.tsv"
	log "FLOW: step=6 name=identify_archaic_haplotypes hap_rows=$nhap matched_haplotype_loci=$nregion selected_after_YRI_filter=$nselected yri_filter_method=${ymethod:-NA} yri_freq_th=${yth:-NA} diagnostic_sites=$ndiag"
}
log_process_s7(){
	local nphy nmain npng
	nphy=$(data_rows "$dirout/report/phy_region.tsv")
	nmain=$(find "$dirout/phy" -name '*.main.phy' 2>/dev/null | wc -l)
	npng=$(find "$dirout/plot" -maxdepth 1 -name 's8_tree_main_*.png' 2>/dev/null | wc -l)
	log "process: Step7 phy files input_loci=$nphy main_phy_files=$nmain phy_dir=$dirout/phy"
	log "process: Step8 tree plots plotted_loci=$npng plot_dir=$dirout/plot"
	log "FLOW: step=7 name=build_phylogeny_input phy_loci=$nphy main_phy_files=$nmain"
	log "FLOW: step=8 name=render_tree_plots plotted_loci=$npng plot_dir=$dirout/plot"
}
write_hap_sample_map(){
	local out="$dirout/report/haplotype_sample_map.tsv"
	local hap="$dirout/report/hap_match.tsv"
	[[ -s $hap ]] || { log "WARN haplotype sample map skipped: missing $hap"; return 0; }
	Rscript "$dirscript/locus.R" hap_sample_map "$dirout" "$sample_file" "$out"
	log "summary: haplotype Roman-label sample map written: $out ($(data_rows "$out") rows)"
}
log_positive_loci_fate(){
	local bed="$positive_loci" out="$dirout/report/positive_loci_fate.tsv"
	[[ -s $bed ]] || return 0
	Rscript "$dirscript/locus.R" positive_loci_fate "$dirout" "$bed" "$out"
	log "summary: sanity loci fate table written: $out"
	if [[ -s $out ]]; then
		awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print "SANITY_LOCUS_FATE:",$5,$1,$6}' "$out" | while IFS= read -r line; do log "$line"; done
	fi
}
clean_report_workbooks_only(){
	find "$dirout/report" -mindepth 1 -type d -exec rm -rf {} +
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
	[[ -s $o/kg.tsv && $(data_rows "$o/kg.tsv") -gt 0 ]] || return 1
	[[ -s $o/kg.samples.tsv ]] || return 1
	n_arch=$(find "$o" -maxdepth 1 -type f -name '*.tsv' ! -name 'kg.tsv' ! -name 'kg.samples.tsv' -size +0c 2>/dev/null | wc -l)
	(( n_arch >= expected_arch ))
}
check_pbase(){
	local pbase=$1 assoc=$2 label=$3 flag=$4 m1 m2 e1 e2 ns nv
	require_pfile "$pbase" "$label"
	[[ -f $flag ]] && return
	m1=$(awk '!/^#/ && $3!=""{print $3}' "$pbase.pvar" | id_mode)
	m2=$(awk 'NR>1{print $3}' "$assoc" | id_mode)
	e1=$(awk '!/^#/ && $3!="." && $3!=""{print $3; if(++n==5) exit}' "$pbase.pvar" | paste -sd, -)
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
	for x in Rscript bcftools tabix plink2 awk sed sort; do command -v "$x" >/dev/null || { log "ERROR missing command: $x"; exit 1; }; done
	[[ -s "$dirscript/locus.R" ]] || { log "ERROR missing $dirscript/locus.R"; exit 1; }
	grep -q "add_positive_loci" "$dirscript/locus.R" || { log "ERROR $dirscript/locus.R is an old version without add_positive_loci; copy the latest locus.R to $dirscript/ first."; exit 1; }
}
check_phyml(){ command -v phyml >/dev/null || { log "ERROR missing phyml"; exit 1; }; }
run_phyml_cmd(){
	local cpus=$1 f=$2 leave_duplicates=$3 logf=$4
	if [[ $leave_duplicates == 1 ]]; then
		if command -v timeout >/dev/null 2>&1; then
			PHYMLCPUS="$cpus" timeout "$phyml_timeout" phyml -i "$f" -m HKY85 -c 4 -a e -v e -b "$phyml_boot" --leave_duplicates > "$logf" 2>&1
		else
			PHYMLCPUS="$cpus" phyml -i "$f" -m HKY85 -c 4 -a e -v e -b "$phyml_boot" --leave_duplicates > "$logf" 2>&1
		fi
	else
		if command -v timeout >/dev/null 2>&1; then
			PHYMLCPUS="$cpus" timeout "$phyml_timeout" phyml -i "$f" -m HKY85 -c 4 -a e -v e -b "$phyml_boot" > "$logf" 2>&1
		else
			PHYMLCPUS="$cpus" phyml -i "$f" -m HKY85 -c 4 -a e -v e -b "$phyml_boot" > "$logf" 2>&1
		fi
	fi
}
run_phyml_fallback_cmd(){
	local f=$1 leave_duplicates=$2 logf=$3
	local args=(-i "$f" -m HKY85 -c 1 -a 0 -v 0 -b 0 -o n)
	[[ $leave_duplicates == 1 ]] && args+=(--leave_duplicates)
	if command -v timeout >/dev/null 2>&1; then
		PHYMLCPUS=1 timeout "$phyml_timeout" phyml "${args[@]}" > "$logf" 2>&1
	else
		PHYMLCPUS=1 phyml "${args[@]}" > "$logf" 2>&1
	fi
}
run_phyml_file(){
	local f=$1 leave_duplicates=0 tree="${f}_phyml_tree.txt" stats="${f}_phyml_stats.txt" logf="${f}.phyml.run.log" rc=0
	[[ $# -ge 2 ]] && leave_duplicates=$2
	[[ -s $f ]] || { log "ERROR missing PHYLIP input: $f"; return 1; }
	if [[ $phy_input == strict ]]; then
		valid_phy_file "$f" || { log "ERROR invalid PHYLIP spacing/label format: $f; remove or regenerate this .phy"; return 1; }
	fi
	if [[ -s $tree ]]; then log "phyml skip existing tree: $f"; return 0; fi
	if [[ -e $tree && ! -s $tree ]]; then log "WARN remove empty phyml tree before rerun: $tree"; rm -f "$tree"; fi
	if [[ -e $stats && ! -s $stats ]]; then log "WARN remove empty phyml stats before rerun: $stats"; rm -f "$stats"; fi
	log "phyml start: $f cpus=$phyml_cpus boot=$phyml_boot timeout=$phyml_timeout"
	run_phyml_cmd "$phyml_cpus" "$f" "$leave_duplicates" "$logf" || rc=$?
	if [[ $phy_input == strict && $rc -ne 0 && $phyml_retry_cpus -lt $phyml_cpus ]]; then
		log "WARN phyml failed rc=$rc at cpus=$phyml_cpus; retry cpus=$phyml_retry_cpus: $f"
		rm -f "$tree" "$stats"
		rc=0
		run_phyml_cmd "$phyml_retry_cpus" "$f" "$leave_duplicates" "${logf}.retry${phyml_retry_cpus}" || rc=$?
	fi
	if [[ $phy_input == strict && $rc -ne 0 && $phyml_fallback == 1 ]]; then
		log "WARN phyml optimized run failed rc=$rc; fallback BioNJ/no-optimization/no-bootstrap: $f"
		rm -f "$tree" "$stats"
		rc=0
		run_phyml_fallback_cmd "$f" "$leave_duplicates" "${logf}.fallback" || rc=$?
		if [[ $rc -eq 0 && -s $tree ]]; then
			touch "${f}.phyml.fallback"
			log "WARN phyml fallback tree written without bootstrap support: $tree"
		fi
	fi
	if [[ $rc -ne 0 ]]; then
		log "ERROR phyml failed rc=$rc: $f; see $logf; retry=$(clean_msg "${logf}.retry${phyml_retry_cpus}") fallback=$(clean_msg "${logf}.fallback")"
		return 1
	fi
	if [[ ! -s $tree ]]; then
		log "ERROR phyml finished but tree is missing/empty: $tree; see $logf; $(clean_msg "$logf")"
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
}
vcf_auto(){
	local c=$1
	if [[ $ref_pop == ALL && -s $dirmod/vcf/chr$c.vcf.gz ]]; then echo "$dirmod/vcf/chr$c.vcf.gz"; else echo "$dirmod/vcf/${ref_pop}.chr$c.vcf.gz"; fi
}
chrX_part(){
	local bp=$1 par=nonPar
	((bp>=60001 && bp<=2699520)) && par=par
	((bp>=154931044 && bp<=155260560)) && par=par
	echo "$par"
}
pbase_by_bp(){
	local c=$1 bp=0 par
	[[ $# -ge 2 ]] && bp=$2
	if [[ $c == X || $c == 23 ]]; then
		if [[ $ref_pop == ALL ]]; then
			[[ -s "$dirmod/pfile/ALL.male.chrX.nonPar.pgen" ]] && echo "$dirmod/pfile/ALL.male.chrX.nonPar" || echo "$dirmod/pfile/chrX"
		else par=$(chrX_part "$bp"); echo "$dirmod/pfile/${ref_pop}.male.chrX.$par"; fi
	else
		pbase_auto "$c"
	fi
}
vcf_by_region(){
	local c=$1 lo=0 par
	[[ $# -ge 2 ]] && lo=$2
	if [[ $c == X || $c == 23 ]]; then
		if [[ $ref_pop == ALL ]]; then
			[[ -s "$dirmod/vcf/ALL.male.chrX.nonPar.vcf.gz" ]] && echo "$dirmod/vcf/ALL.male.chrX.nonPar.vcf.gz" || echo "$dirmod/vcf/chrX.vcf.gz"
		else par=$(chrX_part "$lo"); echo "$dirmod/vcf/${ref_pop}.male.chrX.$par.vcf.gz"; fi
	else
		vcf_auto "$c"
	fi
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 1000G ID matching and archaic VCF helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lead_id(){
	local pvar=$1 snp=$2 bp=$3 ea= oa= line ch rb id rf al p2 rf2 al2
	[[ $# -ge 4 ]] && ea=$4
	[[ $# -ge 5 ]] && oa=$5
	line=$(awk -v s="$snp" -v p="$bp" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p && $3==s{print $1,$2,$3,$4,$5; exit}' "$pvar" || true)
	if [[ -z $line && $snp =~ ^[^:]+:([0-9]+)_([ACGT]+)_([ACGT]+)$ ]]; then
		p2=${BASH_REMATCH[1]}; rf2=${BASH_REMATCH[2]}; al2=${BASH_REMATCH[3]}
		line=$(awk -v p="$p2" -v r="$rf2" -v a="$al2" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p && toupper($4)==r && toupper($5)==a{print $1,$2,$3,$4,$5; exit}' "$pvar" || true)
	fi
	[[ -z $line ]] && line=$(awk -v s="$snp" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $3==s{n++; z=$1 FS $2 FS $3 FS $4 FS $5} END{if(n==1) print z}' "$pvar" || true)
	if [[ -z $line && -n $ea && -n $oa ]]; then
		line=$(awk -v p="$bp" -v a="$ea" -v b="$oa" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p{r=toupper($4); alt=toupper($5); a=toupper(a); b=toupper(b); if((r==a && alt==b)||(r==b && alt==a)){print $1,$2,$3,$4,$5; exit}}' "$pvar" || true)
	fi
	[[ -z $line ]] && line=$(awk -v p="$bp" 'BEGIN{FS=OFS="\t"} $1!~/^#/ && $2==p{n++; z=$1 FS $2 FS $3 FS $4 FS $5} END{if(n==1) print z}' "$pvar" || true)
	[[ -n $line ]] || return 1
	read -r ch rb id rf al <<< "$line"
	[[ $id == "." || -z $id ]] && id="${ch}:${rb}:${rf}:${al}"
	printf "%s\t%s\n" "$rb" "$id"
}

arch_vcf(){
	local a=$1 c=$2
	case "$a" in
		vindija) echo "$arch0/Vindija/chr${c}_mq25_mapab100.vcf.gz" ;;
		altai) echo "$arch0/Altai/chr${c}_mq25_mapab100.vcf.gz" ;;
		chagyr) [[ -s $arch0/Chagyr/chr${c}.noRB.bgz.vcf.gz ]] && echo "$arch0/Chagyr/chr${c}.noRB.bgz.vcf.gz" || echo "$arch0/Chagyr/chr${c}.noRB.vcf.gz" ;;
		denisova) echo "$arch0/Denisova/chr${c}_mq25_mapab100.vcf.gz" ;;
		denisova25) echo "$arch0/Denisova25/chr${c}.Den25.L35MQ25.B30.map35_100.vcf.gz" ;;
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
			[[ -s $dirout/lead/$t.lead.assoc && -s $dirout/lead/$t.lead.3col && -s $dirout/lead/$t.id_sanity.tsv ]] || need_prep=1
			[[ -s $dirout/lead/prep.version ]] && grep -qx 'prep_input_v3_id_first_pos_if_missing_id' "$dirout/lead/prep.version" || need_prep=1
			[[ $(data_rows "$dirout/lead/$t.lead.assoc") -gt 0 ]] || need_prep=1
		done
		if (( need_prep == 0 )); then
			log "s1 prep_input skip: existing lead files are complete for traits=$traits"
		else
			log "Now running Step1 prep input: input=$dirgwas; output=$dirout/lead"
			Rscript "$dirscript/locus.R" prep_input_local --dirgwas "$dirgwas" --dirout "$dirout" --dirmod "$dirmod/pfile" --ref_pop "$ref_pop" --traits "${traits_arr[@]}"
			log "s1 prep_input local done: $(wc -l < "$dirout/lead/lead.assoc") lines in lead.assoc"
		fi
		if [[ " $traits " == *" bald "* ]] && ! awk 'NR>1 && $2==3 && $3=="rs35044562"{f=1} END{exit !f}' "$dirout/lead/bald.lead.assoc"; then
			awk 'BEGIN{FS=OFS="\t"} NR==1{print; print "bald",3,"rs35044562",45909024,"G","A",1,"G","rs35044562","ID_manual"; next} {print}' "$dirout/lead/bald.lead.assoc" > "$dirout/lead/.bald.lead.assoc" && mv "$dirout/lead/.bald.lead.assoc" "$dirout/lead/bald.lead.assoc"
			awk 'BEGIN{FS=OFS="\t"} NR==1{print; print 3,"rs35044562",45909024; next} {print}' "$dirout/lead/bald.lead.3col" > "$dirout/lead/.bald.lead.3col" && mv "$dirout/lead/.bald.lead.3col" "$dirout/lead/bald.lead.3col"
			log "manual add bald rs35044562 at 1000G POS 45909024"
		fi
		first=${traits%% *}
		{ head -n1 "$dirout/lead/$first.lead.assoc"; for t in $traits; do awk 'NR>1' "$dirout/lead/$t.lead.assoc"; done; } > "$dirout/lead/lead.assoc"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.single.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\treason\tdetail\tpbase\tvid\n" > "$dirout/lead/lead_1000G.fail.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tPAR\tLD_pfile\n" > "$dirout/lead/chrX_PAR.tsv"
		printf "trait\tchr\tsnp\tlead_bp\tpvar_bp\tvid\tmatch\tpbase\tstatus\tn_ld_snp\tmessage\n" > "$dirout/lead/ld_debug.tsv"
	else
		local need_prep=0
		for t in $traits; do
			[[ -s $dirout/lead/$t.lead.assoc && $(data_rows "$dirout/lead/$t.lead.assoc") -gt 0 ]] || need_prep=1
		done
		if [[ -s $dirout/lead/lead.assoc && $need_prep -eq 0 ]]; then
			log "s1 prep_input skip: existing lead files are complete for traits=$traits"
		else
			log "Now running Step1 prep input: input=$dirgwas and $dircojo; output=$dirout/lead"
			Rscript "$dirscript/locus.R" prep_input --dirgwas "$dirgwas" --dircojo "$dircojo" --dirout "$dirout" --traits "${traits_arr[@]}"
			log "s1 done: lead.assoc $(nrow "$dirout/lead/lead.assoc") rows"
		fi
	fi

	# Positive/negative-control loci are added only to the first trait.
	# The BED parser uses only columns 1-4; optional notes in later columns are ignored.
	if [[ -s $positive_loci && $add_positive_loci != 0 ]]; then
		local first=${traits_arr[0]}
		log "s1 add positive-control loci to first trait only: $first; bed=$positive_loci"
		Rscript "$dirscript/locus.R" add_positive_loci --dirgwas "$dirgwas" --dirmod "$dirmod/pfile" --dirout "$dirout" --bed "$positive_loci" --trait "$first" --ref_pop "$ref_pop"
		{ head -n1 "$dirout/lead/$first.lead.assoc"; for t in $traits; do awk 'NR>1' "$dirout/lead/$t.lead.assoc"; done; } > "$dirout/lead/lead.assoc"
	elif [[ $add_positive_loci != 0 ]]; then
		log "WARN sanity-check loci BED not found or empty: $positive_loci"
	fi
	log_process_s1
}

force_positive_picks(){
	if [[ ! -s $positive_loci || $add_positive_loci == 0 ]]; then
		[[ $add_positive_loci == 0 ]] || log "WARN sanity-check loci BED not found or empty: $positive_loci"
		return 0
	fi
	local first=${traits_arr[0]} pickf=$dirout/lead/pick.tsv singlef=$dirout/lead/pick.single.tsv posf=$dirout/lead/positive_pick.tsv tmp
	[[ -s $pickf && -s "$dirout/lead/$first.lead.assoc" ]] || return 0
	tmp=$(mktemp)
	awk -v trait="$first" -v assoc="$dirout/lead/$first.lead.assoc" '
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

calculate_ld_simple(){
	check_cmd
	local t=$1 c=$2 cn pfx tmp assoc snp bp ea oa got rb id out
	cn=$(chrn "$c"); tmp=$dirout/ld/$t/chr$c; assoc=$dirout/lead/$t.lead.assoc
	mkdir -p "$tmp"
	if ok "$tmp/block.tsv"; then log "s2 skip $t chr$c"; return 0; fi
	printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$tmp/ld.tsv"
	log "s2 simple-LD $t chr$c"
	awk -v c="$cn" 'BEGIN{FS=OFS="\t"} NR>1 && $2==c{print $3,$4,$5,$6}' "$assoc" |
	while IFS=$'\t' read -r snp bp ea oa; do
		pfx=$(pbase_by_bp "$c" "$bp")
		require_pfile "$pfx" "$t chr$c $snp"
		got=$(lead_id "$pfx.pvar" "$snp" "$bp" "$ea" "$oa") || { log "WARN not in 1000G: $t chr$c $snp $bp"; continue; }
		rb=${got%%$'\t'*}; id=${got##*$'\t'}; out=$tmp/$snp
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\n" "$t" "$cn" "$snp" "$rb" "$cn" "$rb" "$id" >> "$tmp/ld.tsv"
		plink2 --pfile "$pfx" --set-missing-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 missing --r2-unphased allow-ambiguous-allele --ld-snp "$id" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --threads "$plink_threads" --memory 4000 --out "$out" > "$out.plink.log" 2>&1 || { log "WARN plink2 failed: $out.plink.log $(clean_msg "$out.plink.log")"; continue; }
		[[ -s $out.vcor ]] || continue
		awk -v t="$t" -v c="$cn" -v s="$snp" -v b="$rb" -v lead="$id" 'BEGIN{FS="[ \t]+";OFS="\t"} NR==1{for(i=1;i<=NF;i++){gsub(/^#/,"",$i); if($i=="POS_A") pa=i; if($i=="ID_A") ia=i; if($i=="POS_B") pb=i; if($i=="ID_B") ib=i; if($i~/R2$/) ir=i}; next} pa&&ia&&pb&&ib&&ir{if($ia==lead) print t,c,s,b,c,$pb,$ib,$ir; else if($ib==lead) print t,c,s,b,c,$pa,$ia,$ir}' "$out.vcor" >> "$tmp/ld.tsv"
	done
	awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld" && mv "$tmp/.ld" "$tmp/ld.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR==1{next}{k=$1 FS $2 FS $3 FS $4; bp=$4+0; pos=$6+0; if(!(k in mn)||bp<mn[k]) mn[k]=bp; if(!(k in mx)||bp>mx[k]) mx[k]=bp; if(pos<mn[k]) mn[k]=pos; if(pos>mx[k]) mx[k]=pos; u[k FS bp]=u[k FS pos]=1} END{print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"; for(k in mn){n=0; for(i in u) if(index(i,k FS)==1) n++; split(k,x,FS); print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1}}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
}

pvar_lookup_local(){
	local pbase=$1 lead_assoc=$2 tmp=$3 cache key
	cache="$tmp/pvar.$(basename "$pbase").lookup.tsv"
	[[ -s $cache ]] && { printf "%s\n" "$cache"; return; }
	key="$tmp/pvar.$(basename "$pbase").keys.tsv"
	awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3,$4,$5,$6}' "$lead_assoc" > "$key"
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{bp[$3]=1; id[$2]=1; allele[$3 FS toupper($4) FS toupper($5)]=1; allele[$3 FS toupper($5) FS toupper($4)]=1; next} /^#/{next} (($2 in bp) || ($3 in id) || (($2 FS toupper($4) FS toupper($5)) in allele)){print $2,$3,$4,$5}' "$key" "$pbase.pvar" > "$cache"
	rm -f "$key"
	printf "%s\n" "$cache"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step2: calculate high-LD blocks from 1000G pfiles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calculate_ld_strict(){
	check_cmd
	local t=$1 c=$2 cn assoc tmp nlead snp bp ea oa par pbase flag lookup line refbp vid _ref _alt match n_id outpre plog outvcor msg cmdtxt
	cn=$(chrn "$c"); assoc=$dirout/lead/$t.lead.assoc; tmp=$dirout/ld/$t/chr$c
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
		if plink2 --pfile "$pbase" --allow-extra-chr --make-founders --r2-unphased allow-ambiguous-allele --ld-snp-list "$outpre.ld_ids" --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 "$ld_r2" --threads "$plink_threads" --out "$outpre" > "$plog" 2>&1; then
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
	local t=$1 d=$dirout/ld/$t c cn nlead
	for c in $chrs; do
		cn=$(chrn "$c")
		nlead=$(awk -v cn="$cn" -v c="$c" 'BEGIN{FS=OFS="\t"} NR>1 && ($2==cn || toupper($2)==c){n++} END{print n+0}' "$dirout/lead/$t.lead.assoc")
		[[ $(data_rows "$d/chr$c/lead.assoc") -eq $nlead ]] || return 1
		valid_chr_ld "$d/chr$c" "$nlead" || return 1
	done
}

merge_trait_ld(){
	local t=$1 d=$dirout/ld/$t x c f first pickf singlef
	pickf=$dirout/lead/$t.pick.tsv
	singlef=$dirout/lead/$t.pick.single.tsv
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
	d=$dirout/ld/$t
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
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.tsv"
		printf "trait\tlead_chr\tlead_snp\tlead_bp\tstart\tend\tn\tsize_bp\n" > "$dirout/lead/pick.single.tsv"
		for t in $traits; do
			[[ -s "$dirout/lead/$t.pick.tsv" ]] && awk 'NR>1' "$dirout/lead/$t.pick.tsv" >> "$dirout/lead/pick.tsv"
			[[ -s "$dirout/lead/$t.pick.single.tsv" ]] && awk 'NR>1' "$dirout/lead/$t.pick.single.tsv" >> "$dirout/lead/pick.single.tsv"
		done
		log_trait_count "candidate high-LD loci for haplotype analysis" "$dirout/lead/pick.tsv"
		log_trait_count "single-SNP lead loci excluded before haplotype analysis" "$dirout/lead/pick.single.tsv"
		log "s2 LD done. Check: $dirout/lead/ld_debug.tsv and $dirout/lead/lead_1000G.fail.tsv"
		return 0
	fi
	for t in $traits; do
		for c in $chrs; do
			limit_jobs "$job_in_trait"
			if [[ $ld_calc == strict ]]; then calculate_ld_strict "$t" "$c" & else calculate_ld_simple "$t" "$c" & fi
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
		d=$dirout/ld/$t; mkdir -p "$d"
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
	expected_arch=$(find "$arch0" -mindepth 1 -maxdepth 1 -type d | wc -l)
	[[ $bp -lt $lo ]] && lo=$bp
	[[ $bp -gt $hi ]] && hi=$bp
	d=$dirout/coreVcf/$tr/${chr}.${snp}.${bp}
	o=$dirout/mat/$tr/${chr}.${snp}.${bp}
	if valid_mat_locus "$o" "$expected_arch"; then log "s4 skip existing mat: $tr ${chr}.${snp}.${bp}"; return 0; fi
	log "s4 mat start: $tr ${chr}.${snp}.${bp} region=$chr:$lo-$hi"
	rm -rf "$d" "$o"; mkdir -p "$d" "$o"
	cl=$(chrl "$chr")
	kg=$(vcf_by_region "$cl" "$lo")
	[[ -s $kg ]] || { log "ERROR missing 1000G VCF for $tr chr$chr $snp: $kg"; return 1; }
	bcftools view -r "$cl:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kg" || { log "ERROR bcftools failed: 1000G $tr chr$chr $snp"; return 1; }
	tabix -f -p vcf "$d/kg.vcf.gz" || { log "ERROR tabix failed: $d/kg.vcf.gz"; return 1; }
	bcftools query -l "$d/kg.vcf.gz" > "$o/kg.samples.tsv" || { log "ERROR bcftools query samples failed: $d/kg.vcf.gz"; return 1; }
	for ad in "$arch0"/*; do
		[[ -d $ad ]] || continue
		name=$(basename "$ad"); outname=$(echo "$name" | tr '[:upper:]' '[:lower:]'); avcf=""
		for f in "$ad"/*chr"$cl"_*.vcf.gz "$ad"/*chr"$cl".*.vcf.gz; do [[ -s $f ]] || continue; avcf=$f; break; done
		[[ -n $avcf ]] || { log "WARN missing archaic VCF: $outname chr$cl"; continue; }
		arc_raw="$d/$outname.raw.vcf.gz"; arc_vcf="$d/$outname.vcf.gz"; arc_log="$d/$outname.project.log"
		bcftools view -r "$cl:$lo-$hi" -Oz -o "$arc_raw" "$avcf" || { log "ERROR bcftools failed: archaic raw $outname $cl:$lo-$hi"; return 1; }
		tabix -f -p vcf "$arc_raw" || { log "ERROR tabix failed: $arc_raw"; return 1; }
		Rscript "$dirscript/locus.R" prep_archaic "$d/kg.vcf.gz" "$arc_raw" "$arc_vcf" "$outname" > "$arc_log" 2>&1 || { log "ERROR archaic projection failed: $outname $cl:$lo-$hi; see $arc_log"; tail -20 "$arc_log"; return 1; }
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
	log "s4 mat done: $tr ${chr}.${snp}.${bp}; kg_rows=$(data_rows "$o/kg.tsv") archaic_files=$(find "$o" -maxdepth 1 -type f -name '*.tsv' ! -name 'kg.tsv' | wc -l)"
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

build_core_vcf_simple(){
	check_cmd
	local t=$1 c=$2 cn cl snp bp st en n size lo hi d a f kgvcf
	cn=$(chrn "$c"); cl=$(chrl "$cn")
	log "s4 core $t chr$c"
	awk -v t="$t" -v c="$cn" 'BEGIN{FS=OFS="\t"} NR>1 && $1==t && $2==c{print}' "$dirout/lead/pick.tsv" |
	while IFS=$'\t' read -r t cn snp bp st en n size; do
		lo=$st; hi=$en; [[ $bp -lt $lo ]] && lo=$bp; [[ $bp -gt $hi ]] && hi=$bp
		d=$dirout/coreVcf/$t/${cn}.${snp}.${bp}; mkdir -p "$d"
		kgvcf=$(vcf_by_region "$c" "$lo")
		[[ -s $kgvcf ]] || { log "WARN missing 1000G VCF: $kgvcf"; continue; }
		if [[ ! -s $d/kg.vcf.gz ]]; then bcftools view -r "$cl:$lo-$hi" -m2 -M2 -v snps -Oz -o "$d/kg.vcf.gz" "$kgvcf"; tabix -f -p vcf "$d/kg.vcf.gz"; fi
		for a in vindija altai chagyr denisova denisova25; do
			f=$(arch_vcf "$a" "$cl"); [[ -s $f ]] || continue
			[[ -s $d/$a.vcf.gz ]] && continue
			if [[ $archaic_VCF_match == strict ]]; then
				raw="$d/$a.raw.vcf.gz"
				bcftools view -r "$cl:$lo-$hi" -Oz -o "$raw" "$f"; tabix -f -p vcf "$raw"
				Rscript "$dirscript/locus.R" prep_archaic "$d/kg.vcf.gz" "$raw" "$d/$a.vcf.gz" "$a"
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
	for t in $traits; do for c in $chrs; do limit_jobs "$job_in_trait"; build_core_vcf_simple "$t" "$c" & done; done
	wait_all
	log_process_s4_s5
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step5: convert core VCFs to matrix tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build_matrix_simple(){
	check_cmd
	local t=$1 c=$2 cn d o a
	cn=$(chrn "$c")
	log "s5 mat $t chr$c"
	for d in "$dirout/coreVcf/$t"/${cn}.*.*; do
		[[ -d $d ]] || continue
		o=$dirout/mat/$t/$(basename "$d"); mkdir -p "$o"
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
	for t in $traits; do for c in $chrs; do limit_jobs "$job_in_trait"; build_matrix_simple "$t" "$c" & done; done
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
	export SAMPLE_FILE="$sample_file"
	yri_filter_note="yri_freq_th=$yri_freq_th"
	awk -v th="$yri_freq_th" 'BEGIN{exit !(th+0 >= 1)}' && yri_filter_note="$yri_filter_note; disabled"
	log "Now running Step6 haplotype: input=$dirout/mat; output=$dirout/hap and $dirout/report"
	log "process: Step6 YRI filter sample_file=$sample_file yri_freq_th=$yri_freq_th diagnostic_delta_th=$diagnostic_delta_th definition=diagnostic_archaic_allele_fixed_in_high_coverage_archaics_and_absent_or_rare_in_YRI"
	if valid_hap_stage; then
		log "s6 hap skip: report outputs already complete"
		write_hap_sample_map
		log_positive_loci_fate
		log_process_s6
		if [[ $hap_filter == simple ]]; then
			matched_loci=$(data_rows "$dirout/report/region_summary.tsv")
			selected_loci=$(data_rows "$dirout/report/selected_region.tsv")
			log "step6 finished, $matched_loci loci have matched haplotypes; $selected_loci loci remain after YRI frequency filtering ($yri_filter_note)."
		fi
		log "s5 done: region_summary $(nrow "$dirout/report/region_summary.tsv") rows; selected_region $(nrow "$dirout/report/selected_region.tsv") rows"
		return 0
	fi
	[[ -s "$dirscript/locus.R" ]] || { log "ERROR hap stage requires $dirscript/locus.R"; exit 1; }
	while IFS=$'\t' read -r tr chr snp bp st en _n _size_bp; do
		id="${chr}.${snp}.${bp}"
		[[ -d $dirout/mat/$tr/$id ]] || { log "WARN s6 missing mat for $tr $id; skip hap until mat is available"; continue; }
		if valid_hap_locus "$tr" "$id"; then
			log "s6 hap skip existing locus: $tr $id"
			continue
		fi
		missing=$((missing + 1))
		(
			log "s6 hap start locus: $tr $id"
			Rscript "$dirscript/locus.R" make_hap "$dirout" "$tr" "$id"
			touch "$dirout/hap/$tr/$id.done"
			log "s6 hap done locus: $tr $id"
		) &
		while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do wait -n || status=1; done
	done < <(awk 'BEGIN{FS=OFS="\t"} NR>1 && !seen[$1 FS $2 FS $3 FS $4]++{print}' "$dirout/lead/pick.tsv")
	while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
	(( status == 0 )) || { log "ERROR one or more s6 hap locus jobs failed"; exit 1; }
	log "s6 hap locus pass complete: missing_or_incomplete=$missing"
	Rscript "$dirscript/locus.R" make_hap "$dirout" merge || { log "ERROR locus.R make_hap merge failed"; exit 1; }
	write_hap_sample_map
	log_positive_loci_fate
	log_process_s6
	if [[ $hap_filter == simple ]]; then
		matched_loci=$(data_rows "$dirout/report/region_summary.tsv")
		selected_loci=$(data_rows "$dirout/report/selected_region.tsv")
		log "step6 finished, $matched_loci loci have matched haplotypes; $selected_loci loci remain after YRI frequency filtering ($yri_filter_note)."
	fi
	if [[ $hap_filter == strict ]]; then
		local inherited_tsv hap_match_tsv filtered_hap_tsv filtered_inherited_tsv before_inherited_loci after_inherited_loci before_inherited_hap after_inherited_hap
		inherited_tsv="$dirout/report/inherited_segments.tsv"
		hap_match_tsv="$dirout/report/hap_match.tsv"
		log_trait_count "inherited_segments.tsv inherited loci" "$inherited_tsv"
		log "summary: loci confidence reference (file path $inherited_tsv; columns p_ils,best_lineage,matched_archaics)"
		log "s5 population filter start: pop=$filter_pop max_count=$filter_max_count"
		Rscript "$dirscript/locus.R" filter_hap "$dirout" "$sample_file" "$filter_pop" "$filter_max_count" || { log "ERROR locus.R filter_hap failed"; exit 1; }
		filtered_hap_tsv="$dirout/report/filtered_hap_match.tsv"
		filtered_inherited_tsv="$dirout/report/filtered_inherited_segments.tsv"
		before_inherited_loci=$(data_rows "$inherited_tsv")
		after_inherited_loci=$(data_rows "$filtered_inherited_tsv")
		before_inherited_hap=$(data_rows "$hap_match_tsv")
		after_inherited_hap=$(data_rows "$filtered_hap_tsv")
		log_trait_count "filtered inherited loci" "$filtered_inherited_tsv"
		log "FILTER SUMMARY: inherited loci before_filter=$before_inherited_loci after_filter=$after_inherited_loci criteria=${filter_pop}<=$filter_max_count"
		log "FILTER SUMMARY: inherited haplotypes before_filter=$before_inherited_hap after_filter=$after_inherited_hap criteria=${filter_pop}<=$filter_max_count"
		log "FILTER SUMMARY: before-filter locus counts are in $dirout/report/all.xlsx sheet inherited_loci_counts"
		log "FILTER SUMMARY: before-filter haplotype counts are in $dirout/report/all.xlsx sheet inherited_haplotype_counts"
		log "FILTER SUMMARY: after-filter loci/haplotypes are in $dirout/report/filtered.xlsx sheets filtered_loci and filtered_haplotypes"
	fi
	log "s5 done: region_summary $(nrow "$dirout/report/region_summary.tsv") rows; selected_region $(nrow "$dirout/report/selected_region.tsv") rows"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step7: build PHYLIP files, run phyml, and draw trees
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s7_phy(){
	check_cmd; check_phyml
	local f fs=() t id status=0 n miss_tree n_tree_png hapf
	log "Now running Step7 phylogeny: input=$dirout/report/phy_region.tsv or selected_region.tsv; output=$dirout/phy and $dirout/plot"
	if [[ $phy_input == strict ]]; then
		if valid_phy_stage; then
			log "s7 phy skip: main PHYLIP, phyml trees, and tree PNGs already complete"
			log_process_s7
			return 0
		fi
		[[ -s "$dirscript/locus.R" ]] || { log "ERROR phy_input=strict requires $dirscript/locus.R"; exit 1; }
		hapf="$dirout/report/filtered_hap_match.tsv"; [[ -s $hapf ]] || hapf="$dirout/report/hap_match.tsv"
		[[ -s $hapf ]] || { log "ERROR no hap_match input for phy stage"; exit 1; }
		while IFS=$'\t' read -r t id; do
			[[ -n $t && -n $id ]] || continue
			if valid_phy_locus "$t" "$id"; then
				log "s7 make_phy skip existing locus: $t $id"
				continue
			fi
			( log "s7 make_phy start locus: $t $id"; Rscript "$dirscript/locus.R" make_phy "$dirout" "$t" "$id"; mkdir -p "$dirout/phy/$t"; touch "$dirout/phy/$t/$id.done"; log "s7 make_phy done locus: $t $id" ) &
			while [[ $(jobs -rp | wc -l) -ge $job_in_trait ]]; do wait -n || status=1; done
		done < <(awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="trait") it=i; if($i=="id") iid=i}; next} it&&iid{print $it,$iid}' "$hapf" | sort -u)
		while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
		(( status == 0 )) || { log "ERROR one or more s7 make_phy locus jobs failed"; exit 1; }
		log "summary: *.main.phy filtered haplotypes (file path $dirout/phy) $(count_files "$dirout/phy" '*.main.phy')"
		[[ $hap_filter == strict ]] && log "summary: filtering reference (workbook $dirout/report/filtered.xlsx; criteria $filter_pop <= $filter_max_count)"
		clean_report_workbooks_only
		status=0
		local phy_trait_jobs=$(( max_cores / (job_phyml * phyml_cpus) ))
		(( phy_trait_jobs < 1 )) && phy_trait_jobs=1
		for t in $traits; do
			(
				n=$(find "$dirout/phy/$t" -name '*.main.phy' 2>/dev/null | wc -l)
				log "s6 phyml trait start: $t ($n main.phy files; tree source=*.main.phy only)"
				while IFS= read -r f; do
					(
						if [[ -s "${f}_phyml_tree.txt" ]]; then log "s6 skip existing phyml tree: $f"; exit 0; fi
						run_phyml_file "$f" 0 || exit 1
					) &
					while [[ $(jobs -rp | wc -l) -ge $job_phyml ]]; do wait -n || exit 1; done
				done < <(find "$dirout/phy/$t" -name '*.main.phy' 2>/dev/null | sort)
				while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || exit 1; done
				log "s6 phyml trait done: $t"
			) &
			while [[ $(jobs -rp | wc -l) -ge $phy_trait_jobs ]]; do wait -n || status=1; done
		done
		while [[ $(jobs -rp | wc -l) -gt 0 ]]; do wait -n || status=1; done
		(( status == 0 )) || log "WARN one or more s6 phyml trait jobs failed or timed out; continuing with completed trees"
		miss_tree=$(find "$dirout/phy" -name '*.main.phy' | while read -r f; do [[ -s "${f}_phyml_tree.txt" ]] || echo "$f"; done | head -1)
		if [[ -n $miss_tree ]]; then
			log "WARN some main.phy files still lack trees; make_tree will continue with completed trees. First missing: $miss_tree"
		fi
		MAKE_TREE_KEEP_EXISTING=1 Rscript "$dirscript/locus.R" make_tree "$dirout" || { log "ERROR locus.R make_tree failed"; exit 1; }
		n_tree_png=$(find "$dirout/plot" -maxdepth 1 -name 's8_tree_main_*.png' | wc -l)
		log "summary: loci phylogeny tree PNG (source *.main.phy; file path $dirout/plot) $n_tree_png"
		[[ $n_tree_png -gt 0 ]] || { log "ERROR no tree PNG generated in $dirout/plot"; exit 1; }
		log_process_s7
		clean_report_workbooks_only
		log "summary: report workbooks $dirout/report/all.xlsx $dirout/report/filtered.xlsx $dirout/report/selected.xlsx"
		return 0
	else
		MAKE_PHY_METHOD=simple Rscript "$dirscript/locus.R" make_phy "$dirout"
		mapfile -t fs < <(find "$dirout/phy" -name "*.phy" | sort)
	fi
	if [[ ${#fs[@]} -eq 0 ]]; then
		if [[ -s "$dirout/report/phy_region.tsv" && $(data_rows "$dirout/report/phy_region.tsv") -eq 0 ]]; then
			log "WARN no loci available for phylogeny; phy_region.tsv has 0 rows. This is allowed when no inherited loci are detected and no sanity BED is provided."
			return 0
		fi
		log "ERROR no .phy files generated; inspect $dirout/report/phy_region.tsv and $dirout/report/hap_match.tsv"
		exit 1
	fi
	for f in "${fs[@]}"; do
		limit_jobs "$job_phyml"
		( run_phyml_file "$f" 1 || exit 1 ) &
	done
	wait_all || log "WARN one or more phyml jobs failed; completed trees remain available under $dirout/phy"
	MAKE_TREE_KEEP_EXISTING=1 Rscript "$dirscript/locus.R" make_tree "$dirout" || { log "ERROR locus.R make_tree failed"; exit 1; }
	local n_tree_png
	n_tree_png=$(find "$dirout/plot" -maxdepth 1 -name 's8_tree_*.png' 2>/dev/null | wc -l)
	log "summary: loci phylogeny tree PNG (file path $dirout/plot) $n_tree_png"
	[[ $n_tree_png -gt 0 ]] || { log "ERROR no tree PNG generated in $dirout/plot"; exit 1; }
	log_process_s7
}

clean_all(){
	find "$dirout" -mindepth 1 ! -name 'locus.log' -exec rm -rf {} + 2>/dev/null || true
	init_output_dirs
	log "summary: cleaned analysis directory; output=$dirout"
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Viome workflow: reference-based aSNP / archaic haplotype map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_viome_current_trait(){
	local mode=${1:-all}
	log_config
	if [[ $mode == prep ]]; then
		s1_prep
		return 0
	fi
	# viome starts from the same lead-locus/risk-allele preparation, then switches to
	# the 2026 GBE-style reference aSNP / archaic haplotype-map scan.
	s1_prep
	# Generate the same lead-SNP LD summaries when possible; viome uses these to
	# label whether a candidate aSNP is in high LD with the GWAS lead SNP.
	# If a locus lacks LD output, locus.R still produces viome reports with LD marked NA.
	s2_ld
	s3_pick
	mkdir -p "$dirout/report" "$dirout/plot" "$dirout/viome"
	[[ -s "$dirscript/locus.R" ]] || { log "ERROR viome requires $dirscript/locus.R"; exit 1; }
	if [[ ! -s $viome_asnp ]]; then
		log "ERROR viome aSNP haplotype map not found: $viome_asnp"
		log "Set viome_asnp=/path/to/aSNPs.haplotypes.v1.tsv in locus.sh or export viome_asnp before running."
		exit 1
	fi
	log "Now running viome workflow: reference aSNP/haplotype-map scan and optional haplotype-network output"
	Rscript "$dirscript/locus.R" viome \
		--dirgwas "$dirgwas" \
		--dircojo "$dircojo" \
		--dirout "$dirout" \
		--dirmod "$dirmod" \
		--asnp "$viome_asnp" \
		--traits "${traits_arr[@]}" \
		--p-th "$viome_p_th" \
		--freq-th "$viome_freq_th" \
		--min-asnp "$viome_min_asnp" \
		--window-kb "$viome_window_kb" \
		--ld-r2 "$viome_ld_r2" \
		--make-network "$viome_make_network" || { log "ERROR locus.R viome failed"; exit 1; }
	log "summary: viome aSNP hits $(nrow "$dirout/report/viome_aSNP_hits.tsv") rows; file=$dirout/report/viome_aSNP_hits.tsv"
	log "summary: viome candidate inherited segments $(nrow "$dirout/report/viome_inherited_segments.tsv") rows; file=$dirout/report/viome_inherited_segments.tsv"
	log "summary: viome region summary $(nrow "$dirout/report/viome_region_summary.tsv") rows; file=$dirout/report/viome_region_summary.tsv"
	log "summary: viome network inputs/plots are under $dirout/viome and $dirout/plot"
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Main workflow and command dispatch
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_current_trait_all(){
	if is_viome; then
		run_viome_current_trait all
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
}

run_current_trait_mode(){
	local mode=$1
	if is_viome; then
		case "$mode" in
			all|ld|core|mat|hap|phy) run_viome_current_trait "$mode" ;;
			prep) log_config; s1_prep ;;
			*) echo "usage: bash locus.sh [all|clean|prep|ld|core|mat|hap|phy]"; exit 1 ;;
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
		*) echo "usage: bash locus.sh [all|clean|prep|ld|core|mat|hap|phy]"; exit 1 ;;
	esac
}

run_workflow(){
	local mode=$1
	init_output_dirs
	run_current_trait_mode "$mode"
	log "ALL DONE: traits=$traits; output=$dirout"
}

mode=all
[[ $# -gt 0 ]] && mode=$1
case "$mode" in
	all) run_workflow all ;;
	clean) clean_all ;;
	prep|ld|core|mat|hap|phy) run_workflow "$mode" ;;
	*) echo "usage: bash locus.sh [all|clean|prep|ld|core|mat|hap|phy]"; exit 1 ;;
esac
