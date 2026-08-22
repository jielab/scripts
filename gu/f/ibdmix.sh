#!/usr/bin/env bash
set -euo pipefail



usage() {
	cat <<'EOF'
Usage:
  ./gu.sh ibdmix [run|report|compare|igv]

Default behavior:
  ./gu.sh ibdmix             Run the genome scan, reports, and IGV export.
  ./gu.sh ibdmix --h         Show this help message only; do not run.

Actions:
  run       Full run. Existing complete expensive chromosome outputs are reused.
  report    Reuse saved scan/summary output and rebuild final calls and reports.
  compare   Reuse final calls and rebuild comparisons and IGV tracks.
  igv       Rebuild IGV files only.
  check     Validate target samples, archaic VCFs, populations, sex, and chrX units without scanning.

Common IBDmix examples:
  dirmod=/mnt/d/data.BIG/refGen/1kg/GRCH37 ./gu.sh ibdmix
  dirmod=/mnt/d/data.BIG/refGen/ukb ./gu.sh ibdmix
  ./gu.sh ibdmix report

Optional loci used in comparison reports:
  Default: /mnt/d/files/ibdmix.loci.bed
  BED-like format: chr start end name
  Example: chrX    66186240    66423003    bald13_rs6525167_core

Common IBDmix variables:
  dir0=/mnt/d
  dir_ref=$dir0/data.BIG/refGen
  dirout=$dir0/analysis/gu/ibdmix
  chrs="1 2 ... 22 X"                 X expands to XPAR, XNONPAR_F, XNONPAR_M
  genome_build=b37                        use PAR/non-PAR coordinates from 0phe.f.sh
  sample_file=FILE                         header must include sample, population, super-population, sex
  target_sample_policy=exact               require all target VCF samples to be represented
  refs="Altai Chagyr Denisova Denisova25 Vindija"
  lod_cut=4
  len_cut=50000
  job_of_chr=2
  job_in_chr=8
  stats_by_group=super_pop
  admixed_african_pops=ACB,ASW
  window_size=1000000

Main outputs:
  final/all_archaic_refs.lod4.len50000.segments.tsv.gz
  final/apparent_neanderthal_AfricanIndividuals.lod4.len50000.segments.tsv.gz
  report/genomewide_segments.raw.tsv.gz
  report/genomewide_segments.reduced.tsv.gz
  report/genomewide_segments.any_ref.reduced.tsv.gz
  report/all_percent.csv
  report/all_length.csv
  report/all.xlsx
  igv/igv_batch_hg19_windows.txt
EOF
}

if [[ $# -gt 0 ]]; then
	case "$1" in
		-h|--h|--help)
			usage
			exit 0
			;;
		--*)
			echo "ERROR unknown option: $1" >&2
			usage >&2
			exit 1
			;;
	esac
	echo "ERROR f/ibdmix.sh is internal; run ./gu.sh ibdmix or ./gu.sh trace." >&2
	exit 1
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Paths, references, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir0=${dir0:-/mnt/d}
dirscript=${dirscript:-$dir0/scripts/gu/f}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Self-contained PAR definitions; GU v2 no longer depends on private 0phe.f.sh.
X_PAR_b37=${X_PAR_b37:-"X:60001-2699520,X:154931044-155260560"}
X_NONPAR_b37=${X_NONPAR_b37:-"X:2699521-154931043"}
X_PAR_b38=${X_PAR_b38:-"X:10001-2781479,X:155701383-156030895"}
X_NONPAR_b38=${X_NONPAR_b38:-"X:2781480-155701382"}
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
count_lines(){ awk 'NF{n++} END{print n+0}' "$1"; }
join_by(){ local sep=$1 first=1 x; shift; for x in "$@"; do (( first )) || printf '%s' "$sep"; printf '%s' "$x"; first=0; done; }
gzip_ok(){ gzip -t "$1" >/dev/null 2>&1; }
n_data_rows(){
	local file=$1
	if [[ $file == *.gz ]]; then gzip -dc "$file"; else cat "$file"; fi | awk 'NR>1{n++} END{print n+0}'
}
data_rows(){ n_data_rows "$1"; }

dir_ref=${dir_ref:-$dir0/data.BIG/refGen}
[[ -d $dir_ref ]] || { log "ERROR missing refGen directory: $dir_ref"; exit 1; }
dirarch=${dirarch:-${dir_archaic:-$dir_ref/archaic/GRCH${GRCH:-37}}/avcf}
dirmod=${dirmod:-${dir_1kg:-$dir_ref/1kg/GRCH${GRCH:-37}}}
sample_file=${sample_file:-$dirmod/vcf/samples_v3.ALL.panel}
dirsoft=${dirsoft:-$dir0/software/IBDmix}
genome_build=${genome_build:-b37}          # b37 or b38; used for PAR/non-PAR coordinates
sample_id_col=${sample_id_col:-}
pop_col=${pop_col:-}
super_pop_col=${super_pop_col:-}
sex_col=${sex_col:-}
target_sample_policy=${target_sample_policy:-exact} # exact: panel must equal target VCF samples; subset: explicit subset

dirout=${dirout:-$dir0/analysis/gu/ibdmix}
dir_report="$dirout/report"
dir_igv="$dirout/igv"

chrs_arg=${chrs:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}
read -r -a chrs <<< "$chrs_arg"
case "$genome_build" in b37|b38) ;; *) echo "ERROR genome_build must be b37 or b38" >&2; exit 1;; esac
case "$target_sample_policy" in exact|subset) ;; *) echo "ERROR target_sample_policy must be exact or subset" >&2; exit 1;; esac
analysis_units=()
for chr0 in "${chrs[@]}"; do
    chr0=${chr0#chr}
    if [[ $chr0 == X || $chr0 == 23 ]]; then
        analysis_units+=(XPAR XNONPAR_F XNONPAR_M)
    else
        analysis_units+=("$chr0")
    fi
done
refs_arg=${refs:-"Altai Chagyr Denisova Denisova25 Vindija"}
read -r -a refs <<< "$refs_arg"
afr_denisova_pops_arg=${afr_denisova_pops:-"ESN GWD LWK MSL YRI"}
read -r -a afr_denisova_pops <<< "$afr_denisova_pops_arg"

job_of_chr=${job_of_chr:-2}
job_in_chr=${job_in_chr:-8}
lod_cut=${lod_cut:-4}
len_cut=${len_cut:-50000}
emit_lod_cut=${emit_lod_cut:-$lod_cut}
minor_allele_count=${minor_allele_count:-1}
archaic_error=${archaic_error:-0.01}
modern_error_max=${modern_error_max:-0.002}
modern_error_proportion=${modern_error_proportion:-2}
stats_by_group=${stats_by_group:-super_pop}
admixed_african_pops=${admixed_african_pops:-ACB,ASW} # comma-separated; separated in report diagnostics
selected_loci=${selected_loci:-$dir0/files/ibdmix.loci.bed}
gu_inherited=${gu_inherited:-}             # optional GU report/inherited_region.tsv for GU-vs-IBDmix comparison
locus_min_cover=${locus_min_cover:-0.5}    # selected-locus match threshold in ibdmix.R
window_size=${window_size:-1000000}        # genome-wide IGV window size in ibdmix.R
igv_genome=${igv_genome:-hg19}             # 1000G b37/UCSC-style chr names load as hg19 in IGV
plot_style=${plot_style:-both} # retained for old command lines; genome-wide PNG plots are disabled

# TRACE starts from pre-computed ARG tree-sequence files and does not use
# archaic genomes as input.
trace_py=${trace_py:-$dirscript/trace.py}
trace_dirarg=${TRACE_DIRARG:-${trace_dirarg:-${dir_1kg:-$dir_ref/1kg/GRCH${GRCH:-37}}/arg}}
trace_dirout=${TRACE_DIROUT:-${trace_dirout:-$dirout/trace}}
trace_sample_map=${TRACE_SAMPLE_MAP:-${trace_sample_map:-}}
trace_include_dir=${TRACE_INCLUDE_DIR:-${trace_include_dir:-}}
trace_genetic_map_dir=${TRACE_GENETIC_MAP_DIR:-${trace_genetic_map_dir:-}}
trace_chrs_arg=${TRACE_CHRS:-${trace_chrs:-$chrs_arg}}
read -r -a trace_chrs <<< "$trace_chrs_arg"
trace_t_archaic=${TRACE_T_ARCHAIC:-${trace_t_archaic:-15000}}
trace_window_size=${TRACE_WINDOW_SIZE:-${trace_window_size:-}}
trace_report_window_size=${TRACE_REPORT_WINDOW_SIZE:-${trace_window_size:-$window_size}}
trace_posterior_threshold=${TRACE_POSTERIOR_THRESHOLD:-${trace_posterior_threshold:-0.9}}
trace_physical_length_threshold=${TRACE_PHYSICAL_LENGTH_THRESHOLD:-${trace_physical_length_threshold:-50000}}
trace_genetic_distance_threshold=${TRACE_GENETIC_DISTANCE_THRESHOLD:-${trace_genetic_distance_threshold:-0.05}}
trace_group_size=${TRACE_GROUP_SIZE:-${trace_group_size:-50}}
trace_job_extract=${TRACE_JOB_EXTRACT:-${trace_job_extract:-1}}
trace_job_infer=${TRACE_JOB_INFER:-${trace_job_infer:-8}}
trace_job_summarize=${TRACE_JOB_SUMMARIZE:-${trace_job_summarize:-8}}
trace_conda_env=${TRACE_CONDA_ENV:-${trace_conda_env:-gu}}

action=${GU_ACTION:-ibdmix_run}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Logging, locks, and runtime checks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir -p "$dirout"/{samples,genotype,raw,summary,combined,final,log,tmp} "$dir_report" "$dir_igv"
exec > >(tee "$dirout/ibdmix.log") 2>&1

summary_header(){ printf "ID\tchrom\tstart\tend\tlength\tslod\tsites\tpositive_lods\tnegative_lods\tsample_set\tsuper_pop\tanc\n"; }
raw_header(){ printf "ID\tchrom\tstart\tend\tslod\tsites\tpositive_lods\tnegative_lods\tmask_and_maf\tin_mask\tmaf_low\tmaf_high\trec_2_0\trec_0_2\n"; }

lock_dir=$dirout/.ibdmix.lock
completed=0
if ! mkdir "$lock_dir" 2>/dev/null; then
	log "ERROR another ibdmix.sh run may be active, or stale lock exists: $lock_dir"
	log "If no ibdmix.sh process is running, remove this lock directory before rerunning."
	exit 1
fi
echo "$$" > "$lock_dir/pid"
cleanup_common(){
	local status=$1
	rm -rf "$lock_dir"
	if (( completed == 0 )); then
		log "process: ibdmix.sh stopped before ALL DONE status=$status; report may be empty because combine/final calls/report did not finish"
		log "process: cleanup tmp directory $dirout/tmp"
		rm -rf "$dirout/tmp"
	fi
}
cleanup_on_exit(){ cleanup_common "$?"; }
cleanup_on_signal(){
	local sig=$1
	log "ERROR received signal $sig"
	cleanup_common 130
	trap - EXIT
	exit 130
}
trap cleanup_on_exit EXIT
trap 'cleanup_on_signal INT' INT
trap 'cleanup_on_signal TERM' TERM

is_trace_action(){
	case "$action" in
		trace_run|trace_infer|trace_segments|trace_report) return 0 ;;
		*) return 1 ;;
	esac
}

if ! is_trace_action; then
	need_cmds=(awk bcftools bedtools find gzip zcat bash Rscript sort md5sum python3 comm)
	for x in "${need_cmds[@]}"; do
		command -v "$x" >/dev/null || { log "ERROR missing command: $x"; exit 1; }
	done
	for x in "$dirsoft/build/src/generate_gt" "$dirsoft/build/src/ibdmix" "$dirsoft/src/summary.sh" "$dirscript/ibdmix.R" "$dirscript/sample_panel.py" "$dirscript/vcf_gt_fix.py" "$sample_file"; do
		[[ -s "$x" ]] || { log "ERROR missing required file: $x"; exit 1; }
	done
	[[ -d "$dirmod/vcf" ]] || { log "ERROR missing modern VCF folder: $dirmod/vcf"; exit 1; }
	if (( job_of_chr < 1 || job_in_chr < 1 )); then
		log "ERROR job_of_chr and job_in_chr must be positive integers"
		exit 1
	fi
else
	for x in awk sort find gzip bash; do
		command -v "$x" >/dev/null || { log "ERROR missing command: $x"; exit 1; }
	done
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 input: discover inputs, validate VCFs, and split 1000G samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ref_dir(){
	case "${1,,}" in
		altai) echo "$dirarch/Altai" ;;
		denisova|denisovan) echo "$dirarch/Denisova" ;;
		*) echo "$dirarch/$1" ;;
	esac
}
ref_label(){
	case "${1,,}" in
		denisovan) echo "Denisova" ;;
		*) echo "$1" ;;
	esac
}
ref_vcf(){
	local ref=$1 chr=$2 d
	d=$(ref_dir "$ref")
	[[ -d $d ]] || return 0
	find "$d" -maxdepth 1 -type f \( -name "*chr${chr}_*.vcf.gz" -o -name "*chr${chr}.*.vcf.gz" -o -name "*chr${chr}.vcf.gz" \) | sort | sed -n '1p'
}
chr_vcf(){
	local chr=$1
	echo "$dirmod/vcf/chr$chr.vcf.gz"
}
pop_super(){
    local pop=$1 x
    if [[ -s $dirout/samples/sample_counts.tsv ]]; then
        x=$(awk -F'	' -v p="$pop" 'NR>1 && $1==p && $4=="all"{print $3; exit}' "$dirout/samples/sample_counts.tsv")
        [[ -n $x ]] && { echo "$x"; return; }
    fi
    case "$pop" in
        CHB|JPT|CHS|CDX|KHV) echo EAS ;;
        CEU|TSI|FIN|GBR|IBS) echo EUR ;;
        YRI|LWK|GWD|MSL|ESN|ASW|ACB) echo AFR ;;
        MXL|PUR|CLM|PEL) echo AMR ;;
        GIH|PJL|BEB|STU|ITU) echo SAS ;;
        *) echo NA ;;
    esac
}
unit_source_chr(){ case "$1" in XPAR|XNONPAR_F|XNONPAR_M) echo X;; *) echo "$1";; esac; }
unit_output_chr(){ case "$1" in XPAR|XNONPAR_F|XNONPAR_M) echo X;; *) echo "$1";; esac; }
unit_sex(){ case "$1" in XNONPAR_F) echo female;; XNONPAR_M) echo male;; *) echo all;; esac; }
unit_region(){
    local unit=$1 var
    case "$unit" in
        XPAR) var=X_PAR_${genome_build}; printf '%s
' "${!var}" ;;
        XNONPAR_F|XNONPAR_M) var=X_NONPAR_${genome_build}; printf '%s
' "${!var}" ;;
        *) printf '
' ;;
    esac
}
vcf_first_contig(){
    python3 - "$1" <<'PYVCF'
import gzip
import sys
path = sys.argv[1]
with open(path, "rb") as raw:
    magic = raw.read(2)
opener = gzip.open if magic == b"\x1f\x8b" else open
with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line and not line.startswith("#"):
            print(line.split("\t", 1)[0].strip())
            break
    else:
        raise SystemExit(f"no VCF records found: {path}")
PYVCF
}
adapt_region_contig(){
    local vcf=$1 region=$2 contig
    [[ -n $region ]] || { echo ""; return; }
    contig=$(vcf_first_contig "$vcf")
    case "$contig" in chrX|X|23) ;; *) echo "ERROR: expected chrX/X/23 contig in $vcf, found $contig" >&2; return 1;; esac
    echo "${region//X:/$contig:}"
}
prepare_samples(){
    local nall require_sex=()
    mkdir -p "$dirout/samples"
    find "$dirout/samples" -maxdepth 1 -type f \( -name '*.txt' -o -name '*.tsv' \) -delete
    printf '%s
' "${analysis_units[@]}" | grep -q '^X' && require_sex=(--require-sex)
    local -a panel_args=(--input "$sample_file" --output "$dirout/samples/panel.normalized.tsv" --list-dir "$dirout/samples")
    [[ -n $sample_id_col ]] && panel_args+=(--sample-col "$sample_id_col")
    [[ -n $pop_col ]] && panel_args+=(--pop-col "$pop_col")
    [[ -n $super_pop_col ]] && panel_args+=(--super-pop-col "$super_pop_col")
    [[ -n $sex_col ]] && panel_args+=(--sex-col "$sex_col")
    python3 "$dirscript/sample_panel.py" "${panel_args[@]}" "${require_sex[@]}"
    nall=$(count_lines "$dirout/samples/ALL.txt")
    (( nall > 0 )) || { log "ERROR no samples written from $sample_file"; exit 1; }
    awk -F'	' 'NR>1 && $4=="all" && $1!="ALL"{print $5}' "$dirout/samples/sample_counts.tsv" > "$dirout/samples/sample_files.list"
    printf 'analysis_unit	source_chr	region	sex
' > "$dirout/samples/analysis_units.tsv"
    local unit
    for unit in "${analysis_units[@]}"; do
        printf '%s	%s	%s	%s
' "$unit" "$(unit_source_chr "$unit")" "$(unit_region "$unit")" "$(unit_sex "$unit")" >> "$dirout/samples/analysis_units.tsv"
    done
    log "process: input sample panel normalized=$dirout/samples/panel.normalized.tsv ALL_samples=$nall populations=$(wc -l < "$dirout/samples/sample_files.list") policy=$target_sample_policy"
    log "process: input X units use genome_build=$genome_build PAR=${X_PAR_b37:-NA}/${X_PAR_b38:-NA} nonPAR=${X_NONPAR_b37:-NA}/${X_NONPAR_b38:-NA}"
    log "process: final Cell-style Denisova subtraction populations=$(join_by , "${afr_denisova_pops[@]}"); genome scan itself includes all target populations for every archaic reference"
}
validate_target_samples(){
    local vcf=$1 label=$2 work=$dirout/tmp/sample_validation
    mkdir -p "$work"
    bcftools query -l "$vcf" | sort -u > "$work/vcf.txt"
    sort -u "$dirout/samples/ALL.txt" > "$work/panel.txt"
    comm -23 "$work/panel.txt" "$work/vcf.txt" > "$work/missing_from_vcf.txt"
    comm -13 "$work/panel.txt" "$work/vcf.txt" > "$work/unrepresented_vcf.txt"
    if [[ -s $work/missing_from_vcf.txt ]]; then
        log "ERROR sample panel contains target IDs absent from $label: $(head -5 "$work/missing_from_vcf.txt" | paste -sd, -)"; return 1
    fi
    if [[ $target_sample_policy == exact && -s $work/unrepresented_vcf.txt ]]; then
        log "ERROR target_sample_policy=exact but $label contains samples absent from sample_file: $(head -5 "$work/unrepresented_vcf.txt" | paste -sd, -)"; return 1
    fi
}
validate_inputs(){
    local unit chr ref vcf miss=0 nmod
    declare -A checked_chr=()
    for unit in "${analysis_units[@]}"; do
        chr=$(unit_source_chr "$unit")
        if [[ -z ${checked_chr[$chr]:-} ]]; then
            checked_chr[$chr]=1
            if [[ -s "$(chr_vcf "$chr")" ]]; then
                nmod=$(bcftools query -l "$(chr_vcf "$chr")" | wc -l)
                log "process: input modern_vcf source_chr=$chr samples=$nmod file=$(chr_vcf "$chr")"
                validate_target_samples "$(chr_vcf "$chr")" "modern target chr$chr VCF" || miss=1
            else
                log "ERROR missing target VCF: $(chr_vcf "$chr")"; miss=1
            fi
            for ref in "${refs[@]}"; do
                vcf=$(ref_vcf "$(ref_label "$ref")" "$chr")
                [[ -s $vcf ]] || { log "ERROR missing archaic VCF ref=$ref source_chr=$chr under $(ref_dir "$ref")"; miss=1; }
            done
        fi
        log "process: input analysis_unit=$unit source_chr=$chr region=$(unit_region "$unit") sex=$(unit_sex "$unit")"
    done
    (( miss == 0 )) || exit 1
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 genotype: prepare chromosome VCFs and generate IBDmix genotype tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prep_modern(){
    local unit=$1 tmp_unit=$2 chr g1k_vcf outvcf mlog region adapted sex sample_opt=() fix_args=()
    chr=$(unit_source_chr "$unit")
    g1k_vcf=$(chr_vcf "$chr")
    outvcf=$tmp_unit/modern.chr$unit.vcf
    mlog=$dirout/log/modern.chr$unit.prep.log
    region=$(unit_region "$unit"); adapted=$(adapt_region_contig "$g1k_vcf" "$region"); sex=$(unit_sex "$unit")
    [[ $sex == female ]] && sample_opt=(-S "$dirout/samples/ALL.female.txt")
    [[ $sex == male ]] && sample_opt=(-S "$dirout/samples/ALL.male.txt")
    if [[ $unit == XNONPAR_M ]]; then fix_args=(--chrom 23 --duplicate-haploid slash)
    elif [[ $unit == XPAR || $unit == XNONPAR_F ]]; then fix_args=(--chrom 23 --require-diploid)
    fi
    echo "[$(date '+%F %T')] bcftools view modern unit=$unit source_chr=$chr region=${adapted:-whole}" > "$mlog"
    local -a view_args=(view --threads 2 -m2 -M2 -v snps)
    [[ -n $adapted ]] && view_args+=(-r "$adapted")
    view_args+=("${sample_opt[@]}" "$g1k_vcf")
    if [[ $unit == XPAR || $unit == XNONPAR_F || $unit == XNONPAR_M ]]; then
        bcftools "${view_args[@]}" 2>> "$mlog" | python3 "$dirscript/vcf_gt_fix.py" "${fix_args[@]}" > "$outvcf" 2>> "$mlog"
    else
        bcftools "${view_args[@]}" 2>> "$mlog" | awk 'BEGIN{OFS="	"} /^#/{print; next} {sub(/^chr/,"",$1); print}' > "$outvcf"
    fi
    [[ -s "$outvcf" ]] || { log "ERROR failed to prepare modern unit=$unit"; return 1; }
}
prep_arch(){
    local unit=$1 ref=$2 avcf=$3 tmp_unit=$4 chr outvcf alog region adapted
    chr=$(unit_source_chr "$unit")
    outvcf=$tmp_unit/$ref.chr$unit.vcf
    alog=$dirout/log/$ref.chr$unit.prep_arch.log
    region=$(unit_region "$unit"); adapted=$(adapt_region_contig "$avcf" "$region")
    echo "[$(date '+%F %T')] bcftools view archaic $ref unit=$unit source_chr=$chr region=${adapted:-whole}" > "$alog"
    local -a view_args=(view --threads 2 -m2 -M2 -v snps)
    [[ -n $adapted ]] && view_args+=(-r "$adapted")
    view_args+=("$avcf")
    if [[ $unit == XPAR || $unit == XNONPAR_F || $unit == XNONPAR_M ]]; then
        bcftools "${view_args[@]}" 2>> "$alog" | python3 "$dirscript/vcf_gt_fix.py" --chrom 23 --duplicate-haploid slash > "$outvcf" 2>> "$alog"
    else
        bcftools "${view_args[@]}" 2>> "$alog" | awk 'BEGIN{OFS="	"} /^#/{print; next} {sub(/^chr/,"",$1); print}' > "$outvcf"
    fi
    [[ -s "$outvcf" ]] || { log "ERROR failed to prepare archaic $ref unit=$unit"; return 1; }
}
check_gt(){
	local gt=$1
	zcat -f "$gt" | awk 'NR==1{nf=NF; next} NF!=nf{print "bad_gt_line",NR,"NF",NF,"expected",nf,"chr",$1,"pos",$2; exit 1}'
}
gt_sample_list(){
	local gt=$1 out=$2
	if [[ ! -s $out || $gt -nt $out ]]; then
		(set +o pipefail; zcat "$gt" | awk 'NR==1{for(i=6;i<=NF;i++) print $i; exit}') | sort -u > "$out"
	fi
}
sample_file_for_gt(){
	local gt=$1 sp=$2 out=$3 stats=$4 gt_samples
	gt_samples=$dirout/samples/gt_headers/$(basename "$gt" .gz).samples.txt
	mkdir -p "$(dirname "$out")" "$(dirname "$gt_samples")"
	gt_sample_list "$gt" "$gt_samples"
	awk 'BEGIN{OFS="\t"} NR==FNR{ok[$1]=1; next} ok[$1]{print; keep++; next} {miss++} END{print keep+0, miss+0 > stats}' stats="$stats" "$gt_samples" "$sp" > "$out"
}

# Generate one archaic-modern genotype table for genotype
generate_gt_one(){
    local unit=$1 ref=$2 tmp_unit=$3 chr avcf mvcf gt glog
    ref=$(ref_label "$ref"); chr=$(unit_source_chr "$unit")
    gt=$dirout/genotype/$ref.chr$unit.gt.txt.gz
    [[ -s $gt ]] && gzip_ok "$gt" && { log "process: genotype generate_gt skip existing ref=$ref unit=$unit gt=$gt"; return 0; }
    rm -f "$gt"
    avcf=$(ref_vcf "$ref" "$chr")
    mvcf=$tmp_unit/modern.chr$unit.vcf
    glog=$dirout/log/$ref.chr$unit.generate_gt.log
    log "process: genotype generate_gt start ref=$ref unit=$unit archaic_vcf=$avcf"
    prep_arch "$unit" "$ref" "$avcf" "$tmp_unit"
    "$dirsoft/build/src/generate_gt" -a "$tmp_unit/$ref.chr$unit.vcf" -m "$mvcf" -o "$tmp_unit/$ref.chr$unit.gt.txt" > "$glog" 2>&1
    gzip -f "$tmp_unit/$ref.chr$unit.gt.txt"
    mv "$tmp_unit/$ref.chr$unit.gt.txt.gz" "$gt"
    check_gt "$gt" >> "$glog" 2>&1 || { log "ERROR malformed gt table ref=$ref unit=$unit gt=$gt"; return 1; }
    rm -f "$tmp_unit/$ref.chr$unit.vcf"
    log "process: genotype generate_gt done ref=$ref unit=$unit gt=$gt"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 genome scan: run IBDmix per population and summarize segment calls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_ibdmix_pop(){
    local unit=$1 ref=$2 pop=$3 outchr gt raw raw_txt raw_tmp sum ilog slog sp sp_run sp_stats super nkeep nmiss
    ref=$(ref_label "$ref"); outchr=$(unit_output_chr "$unit")
    gt=$dirout/genotype/$ref.chr$unit.gt.txt.gz
    raw=$dirout/raw/$ref/$pop.chr$unit.raw.txt.gz
    raw_txt=$dirout/raw/$ref/$pop.chr$unit.raw.txt
    raw_tmp=$dirout/raw/$ref/$pop.chr$unit.raw.txt.gz.tmp.$$
    sum=$dirout/summary/$ref.$pop.chr$unit.segments.tsv.gz
    ilog=$dirout/log/$ref.$pop.chr$unit.ibdmix.log
    slog=$dirout/log/$ref.$pop.chr$unit.summary.log
    sp=$dirout/samples/$pop.txt
    sp_run=$dirout/samples/by_gt/$ref.chr$unit.$pop.txt
    sp_stats=$dirout/samples/by_gt/$ref.chr$unit.$pop.stats
    super=$(pop_super "$pop")
    mkdir -p "$dirout/raw/$ref"
    [[ -s $sp ]] || { log "ERROR missing sample file pop=$pop file=$sp"; return 1; }
    sample_file_for_gt "$gt" "$sp" "$sp_run" "$sp_stats"
    read -r nkeep nmiss < "$sp_stats"
    (( nmiss == 0 )) || log "process: genome scan sex/unit subset ref=$ref pop=$pop unit=$unit keep=$nkeep absent_from_unit=$nmiss"
    if (( nkeep == 0 )); then
        raw_header | gzip -f > "$raw"; summary_header | gzip -f > "$sum"; rm -f "$raw_txt" "$raw_tmp"; return 0
    fi
    if [[ -s $raw ]] && ! gzip_ok "$raw"; then rm -f "$raw"; fi
    if [[ ! -s $raw ]]; then
        echo "[$(date '+%F %T')] ibdmix ref=$ref pop=$pop unit=$unit" > "$ilog"
        rm -f "$raw_txt" "$raw_tmp"
        "$dirsoft/build/src/ibdmix" --genotype <(zcat "$gt") --output "$raw_txt" --sample "$sp_run" \
            --LOD-threshold "$emit_lod_cut" --minor-allele-count-threshold "$minor_allele_count" \
            --archaic-error "$archaic_error" --modern-error-max "$modern_error_max" \
            --modern-error-proportion "$modern_error_proportion" --more-stats >> "$ilog" 2>&1
        [[ -s $raw_txt ]] || { log "ERROR ibdmix wrote empty raw ref=$ref pop=$pop unit=$unit"; return 1; }
        gzip -c "$raw_txt" > "$raw_tmp"; gzip_ok "$raw_tmp" || return 1; mv -f "$raw_tmp" "$raw"; rm -f "$raw_txt"
    fi
    echo "[$(date '+%F %T')] summary ref=$ref pop=$pop unit=$unit len=$len_cut lod=$lod_cut" > "$slog"
    if ! zcat "$raw" | bash "$dirsoft/src/summary.sh" "$len_cut" "$lod_cut" "$pop" - - 2>> "$slog" |
    awk -v ref="$ref" -v pop="$pop" -v super="$super" -v outchr="$outchr" 'BEGIN{OFS="	"} NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "ID","chrom","start","end","length","slod","sites","positive_lods","negative_lods","sample_set","super_pop","anc"; next} $1!="ID"{len=(h["length"] ? $(h["length"]) : $(h["end"])-$(h["start"])); pos=(h["positive_lods"] ? $(h["positive_lods"]) : ""); neg=(h["negative_lods"] ? $(h["negative_lods"]) : ""); print $(h["ID"]),outchr,$(h["start"]),$(h["end"]),len,$(h["slod"]),$(h["sites"]),pos,neg,pop,super,ref}' | gzip -f > "$sum"; then
        log "ERROR summary failed ref=$ref pop=$pop unit=$unit"; return 1
    fi
    log "process: genome scan summary ref=$ref pop=$pop unit=$unit segments=$(n_data_rows "$sum") output=$sum"
}
pop_list_for_ref(){
    # Every archaic reference is scanned against every population in the target panel.
    awk -F'	' 'NR>1 && $4=="all" && $1!="ALL"{print $1}' "$dirout/samples/sample_counts.tsv" | sort -u
}
expected_summary_count(){
    local npop
    npop=$(pop_list_for_ref "$1" | awk 'NF{n++} END{print n+0}')
    echo $(( npop * ${#analysis_units[@]} ))
}
expected_summary_files(){
    local ref=$1 unit pop
    ref=$(ref_label "$ref")
    for unit in "${analysis_units[@]}"; do
        while IFS= read -r pop; do [[ -n $pop ]] && printf '%s/summary/%s.%s.chr%s.segments.tsv.gz
' "$dirout" "$ref" "$pop" "$unit"; done < <(pop_list_for_ref "$ref")
    done
}
genome_scan_complete(){
	local ref expected files f
	[[ -s "$dirout/samples/ALL.txt" ]] || return 1
	if [[ -e "$dirout/.genome_scan.complete" ]]; then
		for ref in "${refs[@]}"; do
			ref=$(ref_label "$ref")
			expected=$(expected_summary_count "$ref")
			mapfile -t files < <(expected_summary_files "$ref")
			(( expected > 0 && ${#files[@]} == expected )) || return 1
            for f in "${files[@]}"; do [[ -s $f ]] && gzip_ok "$f" || return 1; done
		done
		return 0
	fi
	for ref in "${refs[@]}"; do
		ref=$(ref_label "$ref")
		expected=$(expected_summary_count "$ref")
		mapfile -t files < <(expected_summary_files "$ref" | sort)
		(( expected > 0 && ${#files[@]} == expected )) || return 1
		for f in "${files[@]}"; do
			[[ -s $f ]] && gzip_ok "$f" || return 1
			zcat "$f" | awk 'NR==1{ok=($1=="ID" && NF>=10)} END{exit !ok}' || return 1
		done
	done
	touch "$dirout/.genome_scan.complete"
}

# Run genotype and genome scan across all chromosomes
run_ref_unit(){
    local unit=$1 ref=$2 tmp_unit=$3 pop running=0 status=0
    generate_gt_one "$unit" "$ref" "$tmp_unit" || return 1
    while IFS= read -r pop; do
        [[ -n $pop ]] || continue
        run_ibdmix_pop "$unit" "$ref" "$pop" &
        running=$((running + 1))
        if (( running >= job_in_chr )); then wait -n || status=1; running=$((running - 1)); fi
    done < <(pop_list_for_ref "$ref")
    while (( running > 0 )); do wait -n || status=1; running=$((running - 1)); done
    (( status == 0 ))
}
run_unit(){
    local unit=$1 ref tmp_unit status=0 need_gt=0 gt
    log "process: genotype analysis unit start unit=$unit source_chr=$(unit_source_chr "$unit") region=$(unit_region "$unit") sex=$(unit_sex "$unit")"
    tmp_unit=$dirout/tmp/chr$unit
    rm -rf "$tmp_unit"; mkdir -p "$tmp_unit"
    for ref in "${refs[@]}"; do
        ref=$(ref_label "$ref"); gt=$dirout/genotype/$ref.chr$unit.gt.txt.gz
        [[ -s $gt ]] && gzip_ok "$gt" || need_gt=1
    done
    if (( need_gt == 1 )); then prep_modern "$unit" "$tmp_unit" || { rm -rf "$tmp_unit"; return 1; }
    else log "process: genotype modern prep skip unit=$unit because all genotype files already exist"; fi
    for ref in "${refs[@]}"; do run_ref_unit "$unit" "$ref" "$tmp_unit" || status=1; done
    rm -rf "$tmp_unit"
    (( status == 0 )) || { log "ERROR analysis unit failed unit=$unit"; return 1; }
    log "process: genome scan analysis unit done unit=$unit"
}
run_genome_scan(){
    local unit running=0 status=0
    if genome_scan_complete; then log "process: genome scan skip; all expected summary files are complete and readable"; return 0; fi
    rm -rf "$dirout/tmp"; mkdir -p "$dirout/tmp"
    prepare_samples; validate_inputs
    for unit in "${analysis_units[@]}"; do
        run_unit "$unit" &
        running=$((running + 1))
        if (( running >= job_of_chr )); then wait -n || status=1; running=$((running - 1)); fi
    done
    while (( running > 0 )); do wait -n || status=1; running=$((running - 1)); done
    rmdir "$dirout/tmp" 2>/dev/null || true
    (( status == 0 )) || { log "ERROR one or more analysis units failed"; exit 1; }
    touch "$dirout/.genome_scan.complete"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 combine: combine per-population segment summaries by archaic reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
combined_ref_file(){
	local ref=$1
	ref=$(ref_label "$ref")
	printf "%s/combined/%s.lod%s.len%s.segments.tsv.gz\n" "$dirout" "$ref" "$lod_cut" "$len_cut"
}

all_refs_file(){
	printf "%s/final/all_archaic_refs.lod%s.len%s.segments.tsv.gz\n" "$dirout" "$lod_cut" "$len_cut"
}

final_cell_style_file(){
	printf "%s/final/Altai_minus_AFR_Denisova.lod%s.len%s.cell_style.tsv\n" "$dirout" "$lod_cut" "$len_cut"
}

final_neanderthal_file(){
	printf "%s/final/apparent_neanderthal_AfricanIndividuals.lod%s.len%s.segments.tsv.gz\n" "$dirout" "$lod_cut" "$len_cut"
}

combine_ref(){
	local ref=$1 out files n expected
	ref=$(ref_label "$ref")
	out=$(combined_ref_file "$ref")
	mapfile -t files < <(expected_summary_files "$ref" | sort)
	(( ${#files[@]} > 0 )) || { log "ERROR no summary files for ref=$ref"; exit 1; }
	expected=$(expected_summary_count "$ref")
	if (( ${#files[@]} != expected )); then
		log "ERROR incomplete summary files for ref=$ref observed=${#files[@]} expected=$expected"
		log "ERROR do not summarize from partial output; rerun GU_ACTION=ibdmix_run or delete output and run ./gu.sh ibdmix"
		exit 1
	fi
	if [[ ${skip_completed:-0} == 1 &&
		-e "$dirout/.combined.complete" &&
		-s $out ]]; then
		log "process: combine skip ref=$ref; complete output=$out"
		return 0
	fi
	zcat "${files[@]}" | awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} $1!="ID"{print}' | gzip -f > "$out"
	n=$(n_data_rows "$out")
	log "process: combine combine ref=$ref summary_files=${#files[@]} rows=$n output=$out"
}

combine_scan_outputs(){
	local ref
	[[ -s "$dirout/samples/ALL.txt" ]] || prepare_samples
	for ref in "${refs[@]}"; do combine_ref "$ref"; done
	touch "$dirout/.combined.complete"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 final calls: subtract African Denisova-like segments from Altai calls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build_final_calls(){
    local alt den final_bed final_std all_refs final_rows all_rows den_afr_rows afr_csv
    alt=$(combined_ref_file Altai); den=$(combined_ref_file Denisova)
    final_bed=$(final_cell_style_file); final_std=$(final_neanderthal_file); all_refs=$(all_refs_file)
    if [[ ${skip_completed:-0} == 1 && -e "$dirout/.final_calls.complete" && -s $all_refs && -s $final_bed && -s $final_std ]]; then
        log "process: final calls skip; complete outputs already exist"; return 0
    fi
    [[ -s $alt && -s $den ]] || { log "ERROR final subtraction needs $alt and $den"; exit 1; }
    zcat "$dirout"/combined/*.lod${lod_cut}.len${len_cut}.segments.tsv.gz | awk 'BEGIN{FS=OFS="	"} NR==1{print; next} $1!="ID"{print}' | gzip -f > "$all_refs"
    all_rows=$(n_data_rows "$all_refs")
    afr_csv=$(join_by , "${afr_denisova_pops[@]}")
    den_afr_rows=$(zcat "$den" | awk -F'	' -v pops="$afr_csv" 'BEGIN{n=split(pops,a,","); for(i=1;i<=n;i++) ok[a[i]]=1} NR>1 && ok[$10]{c++} END{print c+0}')
    printf "chrom	start	end	length	slod	sites	positive_lods	negative_lods	sample_set	super_pop	anc	ID
" > "$final_bed"
    bedtools subtract \
        -a <(zcat "$alt" | awk 'BEGIN{FS=OFS="	"} NR>1{for(i=2;i<=NF;i++) printf "%s	",$i; print $1}' | sort -k1,1 -k2,2n -k3,3n) \
        -b <(zcat "$den" | awk -v pops="$afr_csv" 'BEGIN{FS=OFS="	"; n=split(pops,a,","); for(i=1;i<=n;i++) ok[a[i]]=1} NR>1 && ok[$10]{print $2,$3,$4}' | sort -k1,1 -k2,2n -k3,3n) |
    awk -v len="$len_cut" 'BEGIN{FS=OFS="	"} {$4=$3-$2; if($4>=len) print}' >> "$final_bed"
    awk 'BEGIN{FS=OFS="	"} NR==1{printf "%s",$NF; for(i=1;i<NF;i++) printf OFS "%s",$i; print ""; next} {printf "%s",$NF; for(i=1;i<NF;i++) printf OFS "%s",$i; print ""}' "$final_bed" | gzip -f > "$final_std"
    final_rows=$(n_data_rows "$final_std")
    log "process: final Cell-style filter Altai_rows=$(n_data_rows "$alt") Denisova_all_target_rows=$(n_data_rows "$den") Denisova_AFR_mask_rows=$den_afr_rows AFR_pops=$afr_csv all_ref_rows=$all_rows final_rows=$final_rows"
    touch "$dirout/.final_calls.complete"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 report: build summary statistics and plots for comparison with IBDmix 2020 Cell
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_selected_loci_bed(){
	local out="$dir_report/selected_loci.bed"
	mkdir -p "$dir_report"
	rm -f "$out"
	if [[ -n ${gu_inherited:-} && -s $gu_inherited ]]; then
		awk 'BEGIN{FS=OFS="\t"}
			NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
			{
				chr=(h["lead_chr"] ? $(h["lead_chr"]) : ""); if(chr=="23") chr="X"; sub(/^chr/,"",chr)
				st=(h["core_start"] ? $(h["core_start"])-1 : (h["lead_bp"] ? $(h["lead_bp"])-1 : -1))
				en=(h["core_end"] ? $(h["core_end"]) : (h["lead_bp"] ? $(h["lead_bp"]) : st+1))
				if(chr=="" || st<0 || en<=st) next
				name=(h["trait"]?$(h["trait"]):"GU") "|" (h["id"]?$(h["id"]):"NA") "|" (h["lead_snp"]?$(h["lead_snp"]):"NA") "|" (h["best_lineage"]?$(h["best_lineage"]):"NA")
				gsub(/[^A-Za-z0-9._|:-]/,"_",name)
				print "chr" chr, st, en, name
			}' "$gu_inherited" | sort -k1,1V -k2,2n -k3,3n > "$out"
	elif [[ -n ${selected_loci:-} && -s $selected_loci ]]; then
		awk 'BEGIN{FS=OFS="\t"} $1 !~ /^#/ && NF>=3{chr=$1; sub(/^chr/,"",chr); name=(NF>=4?$4:"selected_locus"); print "chr"chr,$2,$3,name}' "$selected_loci" | sort -k1,1V -k2,2n -k3,3n > "$out"
	fi
	[[ -s $out ]] && printf "%s\n" "$out"
}

build_ibdmix_reports(){
	local report_in selected_bed selected_args=()
	if [[ ${skip_completed:-0} == 1 &&
		-e "$dirout/.report.complete" &&
		-s "$dir_report/all.xlsx" &&
		-s "$dir_report/genomewide_segments.raw.tsv.gz" &&
		-s "$dir_report/genomewide_segments.reduced.tsv.gz" &&
		-s "$dir_report/genomewide_segments.any_ref.reduced.tsv.gz" ]]; then
		log "process: report skip; complete workbook and segment tables already exist"
		return 0
	fi
	report_in=$(all_refs_file)
	mkdir -p "$dir_report"
	rm -f "$dir_report/input.check.tsv"
	if [[ ! -s $report_in ]] || (( $(n_data_rows "$report_in") == 0 )); then
		log "WARN report report skipped because all-reference calls have 0 rows"
		return 0
	fi
	selected_bed=$(write_selected_loci_bed || true)
	if [[ -n ${selected_bed:-} && -s $selected_bed ]]; then
		selected_args=(--selected_loci "$selected_bed")
	fi
	log "process: report report start sample_file=$sample_file stats_by_group=$stats_by_group selected_loci=${selected_bed:-NA} gu_inherited=${gu_inherited:-NA} report_input=$report_in output_dir=$dir_report"
	Rscript "$dirscript/ibdmix.R" \
		--in "$report_in" \
		--outdir "$dir_report" \
		--sample_keep "$dirout/samples/ALL.txt" \
		--sample_file "$dirout/samples/panel.normalized.tsv" \
        --sample_id_col sample \
		--stats_by_group "$stats_by_group" \
        --genome_build "$genome_build" \
		--window_size "$window_size" \
		--locus_min_cover "$locus_min_cover" \
		--plot_style "$plot_style" \
		--admixed_african_pops "$admixed_african_pops" \
		"${selected_args[@]}" || { log "ERROR ibdmix.R report failed"; exit 1; }
	log "process: report report complete workbook=$dir_report/all.xlsx all_percent=$dir_report/all_percent.csv all_length=$dir_report/all_length.csv"
	log "process: report diagnostic plots all_percent=$dir_report/all_percent.png afr_sensitivity=$dir_report/all_percent_afr_sensitivity.png by_pop=$dir_report/all_percent_by_pop.png denominator=$dir_report/all_percent_denom_sensitivity.png overlap=$dir_report/overlap_reduction.png"
	log "process: report diagnostic tables diagnostics=$dir_report/ibdmix_diagnostics.tsv overlap=$dir_report/overlap_reduction_diagnostics.csv cross_ref=$dir_report/cross_reference_overlap_diagnostics.csv"
	log "process: report reduced segments=$dir_report/genomewide_segments.reduced.tsv.gz any_ref_segments=$dir_report/genomewide_segments.any_ref.reduced.tsv.gz raw_segments=$dir_report/genomewide_segments.raw.tsv.gz"
	if [[ -s "$dir_report/ibdmix_diagnostics.tsv" ]]; then
		awk 'BEGIN{FS="\t"} NR>1{print "process: report DIAG " $1 "=" $2}' "$dir_report/ibdmix_diagnostics.tsv"
	fi
	touch "$dirout/.report.complete"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 IGV: generate IGV input files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
index_track(){
	local f=$1 preset=${2:-bed}
	[[ -s $f ]] || return 0
	if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
		bgzip -f "$f"
		tabix -f -p "$preset" "$f.gz" || log "WARN tabix failed for $f.gz"
		printf "%s.gz\n" "$f"
	else
		gzip -f "$f"
		log "WARN bgzip/tabix not found; wrote $f.gz without tabix index. IGV can load it, but genomic navigation may be slower."
		printf "%s.gz\n" "$f"
	fi
}

make_segment_bed(){
	local in_tsv=$1 out_bed=$2 label=$3
	[[ -s $in_tsv ]] || { log "WARN IGV segment BED skipped: missing $in_tsv"; return 0; }
	zcat -f "$in_tsv" |
	awk -v label="$label" 'BEGIN{FS=OFS="\t"}
		NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
		{
			chr=$(h["chrom"]); sub(/^chr/, "", chr); chr="chr" chr
			start=$(h["start"]); end=$(h["end"])
			id=$(h["ID"])
			anc=(h["anc"] ? $(h["anc"]) : "NA")
			pop=(h["sample_set"] ? $(h["sample_set"]) : "NA")
			super=(h["super_pop"] ? $(h["super_pop"]) : "NA")
			len=(h["seg_len"] ? $(h["seg_len"]) : (h["length"] ? $(h["length"]) : end-start))
			slod=(h["slod"] ? $(h["slod"]) : 0)
			sites=(h["sites"] ? $(h["sites"]) : "NA")
			poslod=(h["positive_lods"] ? $(h["positive_lods"]) : "NA")
			neglod=(h["negative_lods"] ? $(h["negative_lods"]) : "NA")
			nraw=(h["n_raw_segments"] ? $(h["n_raw_segments"]) : "NA")
			bp_rm=(h["bp_removed_by_reduce"] ? $(h["bp_removed_by_reduce"]) : "NA")
			score=int(slod*10); if(score<0) score=0; if(score>1000) score=1000
			rgb="40,90,210"
			if(anc=="Altai") rgb="210,40,40"
			else if(anc=="Chagyr") rgb="230,120,20"
			else if(anc=="Vindija") rgb="30,70,220"
			else if(anc=="Denisova") rgb="130,40,170"
			else if(anc=="Denisova25") rgb="0,145,145"
			else if(anc=="ANY_REF") rgb="80,80,80"
			name=label "|sample=" id "|pop=" pop "|super_pop=" super "|ref=" anc "|len=" len "|slod=" slod "|sites=" sites "|pos_lods=" poslod "|neg_lods=" neglod
			if(nraw!="NA") name=name "|n_raw=" nraw "|bp_removed=" bp_rm
			gsub(/[^A-Za-z0-9._|:=+;,\/-]/, "_", name)
			# BED9: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb.
			# In IGV, mouse-over/click shows the name field, which contains sample ID and segment statistics.
			print chr, start, end, name, score, ".", start, end, rgb
		}' |
	sort -k1,1V -k2,2n -k3,3n > "$out_bed"
	index_track "$out_bed" bed >/dev/null
}

make_any_ref_percent_bedgraph(){
	local win_csv=$1 out_bedgraph=$2
	[[ -s $win_csv ]] || { log "WARN IGV percent bedGraph skipped: missing $win_csv"; return 0; }
	awk -F',' -v OFS='\t' '
		NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
		$(h["anc"])=="ANY_REF"{
			chr=$(h["chrom"]); sub(/^chr/, "", chr); chr="chr" chr
			print chr, $(h["win_start"])-1, $(h["win_end"]), $(h["burden_percent"])
		}
	' "$win_csv" | sort -k1,1V -k2,2n -k3,3n > "$out_bedgraph"
	if [[ ! -s $out_bedgraph ]]; then
		rm -f "$out_bedgraph"
		log "WARN IGV percent bedGraph skipped: no ANY_REF rows in $win_csv"
		return 0
	fi
	index_track "$out_bedgraph" bed >/dev/null
}

make_group_any_ref_percent_bedgraphs(){
	local win_csv=$1 outdir=$2 group_cols_arg=$3
	[[ -s $win_csv ]] || { log "WARN IGV group percent bedGraph skipped: missing $win_csv"; return 0; }
	[[ -n ${group_cols_arg:-} ]] || return 0
	awk -F',' -v OFS='\t' -v out="$outdir" -v groups="$group_cols_arg" '
		BEGIN{ng=split(groups,g,",")}
		NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
		$(h["anc"])=="ANY_REF"{
			tag=""
			for(i=1;i<=ng;i++){
				if(!(g[i] in h)) next
				v=$(h[g[i]])
				gsub(/[^A-Za-z0-9._-]/, "_", v)
				tag=tag (tag=="" ? v : "." v)
			}
			if(tag=="") next
			chr=$(h["chrom"]); sub(/^chr/, "", chr); chr="chr" chr
			file=out "/percent." tag ".bedGraph"
			print chr, $(h["win_start"])-1, $(h["win_end"]), $(h["burden_percent"]) >> file
		}
	' "$win_csv"
	local f sorted
	for f in "$outdir"/percent.*.bedGraph; do
		[[ -e $f ]] || continue
		sorted="$f.sorted"
		sort -k1,1V -k2,2n -k3,3n "$f" > "$sorted"
		mv "$sorted" "$f"
		index_track "$f" bed >/dev/null
	done
}

copy_selected_loci_for_igv(){
	local igv=$1 src="$dir_report/selected_loci.bed"
	[[ -s $src ]] || return 0
	awk 'BEGIN{FS=OFS="\t"}
		$1 !~ /^#/ && NF>=3{
			chr=$1; sub(/^chr/,"",chr); chr="chr" chr
			name=(NF>=4 ? $4 : "selected_locus")
			gsub(/[^A-Za-z0-9._|:=+;,\/-]/,"_",name)
			print chr,$2,$3,name,1000,".",$2,$3,"0,0,0"
		}' "$src" | sort -k1,1V -k2,2n -k3,3n > "$igv/selected_loci.bed"
	index_track "$igv/selected_loci.bed" bed >/dev/null
}

wsl_to_windows_path(){
	local p=$1 drive rest
	if [[ "$p" =~ ^/mnt/([A-Za-z])/(.*)$ ]]; then
		drive=${BASH_REMATCH[1]}
		rest=${BASH_REMATCH[2]//\//\\}
		printf "%s:\\%s" "${drive^^}" "$rest"
	else
		printf "%s" "$p"
	fi
}

write_igv_batches(){
	local igv=$1 f abs win_abs locus type desc
	locus=${igv_locus:-chr3:45859651-45909024}

	: > "$igv/igv_tracks.txt"
	: > "$igv/igv_tracks_windows.txt"
	printf "file\ttype\tdescription\n" > "$igv/igv_manifest.tsv"

	for f in "$igv"/*.bed.gz "$igv"/*.bedGraph.gz; do
		[[ -e $f ]] || continue
		case "$f" in *.tbi) continue ;; esac
		abs=$(readlink -f "$f" 2>/dev/null || printf "%s" "$f")
		win_abs=$(wsl_to_windows_path "$abs")
		printf "%s\n" "$abs" >> "$igv/igv_tracks.txt"
		printf "%s\n" "$win_abs" >> "$igv/igv_tracks_windows.txt"
		case "$f" in
			*.bedGraph.gz) type="bedGraph"; desc="ANY_REF 1 Mb carrier percent track" ;;
			*) type="BED9"; desc="segment feature track; sample ID and statistics are in the feature name" ;;
		esac
		printf "%s\t%s\t%s\n" "$(basename "$f")" "$type" "$desc" >> "$igv/igv_manifest.tsv"
	done

	{
		printf "new\n"
		printf "genome %s\n" "$igv_genome"
		printf "maxPanelHeight 5000\n"
		for f in "$igv"/*.bed.gz "$igv"/*.bedGraph.gz; do
			[[ -e $f ]] || continue
			case "$f" in *.tbi) continue ;; esac
			abs=$(readlink -f "$f" 2>/dev/null || printf "%s" "$f")
			printf "load %s\n" "$abs"
		done
		printf "goto %s\n" "$locus"
		# Expand feature tracks so overlapping sample-level segments are shown as separate rows.
		printf "expand\n"
	} > "$igv/igv_batch_${igv_genome}.txt"

	{
		printf "new\n"
		printf "genome %s\n" "$igv_genome"
		printf "maxPanelHeight 5000\n"
		for f in "$igv"/*.bed.gz "$igv"/*.bedGraph.gz; do
			[[ -e $f ]] || continue
			case "$f" in *.tbi) continue ;; esac
			abs=$(readlink -f "$f" 2>/dev/null || printf "%s" "$f")
			win_abs=$(wsl_to_windows_path "$abs")
			printf "load %s\n" "$win_abs"
		done
		printf "goto %s\n" "$locus"
		# Expand feature tracks so overlapping sample-level segments are shown as separate rows.
		printf "expand\n"
	} > "$igv/igv_batch_${igv_genome}_windows.txt"
}

build_igv_files(){
	local igv="$dir_igv"
	if [[ ${skip_completed:-0} == 1 &&
		-e "$dirout/.igv.complete" &&
		-s "$igv/igv_batch_${igv_genome}_windows.txt" &&
		-s "$igv/segments.raw.bed.gz" &&
		-s "$igv/segments.raw.bed.gz.tbi" ]]; then
		log "process: IGV skip; complete tracks already exist"
		return 0
	fi
	mkdir -p "$igv"
	rm -f "$igv"/* 2>/dev/null || true

	# Segment tracks.  Load these in IGV and set display mode to Expanded, or use the generated batch script.
	make_segment_bed "$dir_report/genomewide_segments.raw.tsv.gz" "$igv/segments.raw.bed" "raw"
	make_segment_bed "$dir_report/genomewide_segments.reduced.tsv.gz" "$igv/segments.reduced.bed" "reduced"
	make_segment_bed "$dir_report/genomewide_segments.any_ref.reduced.tsv.gz" "$igv/segments.any_ref.reduced.bed" "ANY_REF_reduced"

	# Compact percent tracks.  Only ANY_REF is exported: overall plus optional stats_by_group values such as AFR/AMR/EAS/EUR/SAS.
	make_any_ref_percent_bedgraph "$dir_report/all_gw.csv" "$igv/percent.bedGraph"
	make_group_any_ref_percent_bedgraphs "$dir_report/all_gw_by_group.csv" "$igv" "$stats_by_group"

	copy_selected_loci_for_igv "$igv"
	write_igv_batches "$igv"

	log "process: IGV IGV tracks written to $igv"
	log "process: IGV output names use no ibdmix_ prefix; percent tracks are ANY_REF only"
	log "process: IGV recommended Windows IGV batch=$igv/igv_batch_${igv_genome}_windows.txt locus=${igv_locus:-chr3:45859651-45909024} display=expanded"
	touch "$dirout/.igv.complete"
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 TRACE setup: discover ARG files and prepare TRACE sample map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trace_activate_env(){
	# Preferred for reproducibility: TRACE_CONDA_ENV=trace ./gu.sh trace ./ibdmix.sh
	if [[ -n ${trace_conda_env:-} ]]; then
		# gu.sh normally activates this environment before dispatching here.  Do not
		# activate it a second time: conda's deactivate hooks are not nounset-safe.
		if [[ ${CONDA_DEFAULT_ENV:-} == "$trace_conda_env" ]]; then
			log "process: TRACE using active conda environment: $trace_conda_env"
			return 0
		fi
		if command -v conda >/dev/null 2>&1; then
			eval "$(conda shell.bash hook)" || true
			set +u
			local activate_status=0
			conda activate "$trace_conda_env" || activate_status=$?
			set -u
			(( activate_status == 0 )) || { log "ERROR could not conda activate TRACE_CONDA_ENV=$trace_conda_env"; exit 1; }
			log "process: TRACE conda environment activated: $trace_conda_env"
			return 0
		fi
		for conda_sh in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh" "/opt/conda/etc/profile.d/conda.sh"; do
			if [[ -s $conda_sh ]]; then
				source "$conda_sh"
				set +u
				local activate_status=0
				conda activate "$trace_conda_env" || activate_status=$?
				set -u
				(( activate_status == 0 )) || { log "ERROR could not conda activate TRACE_CONDA_ENV=$trace_conda_env"; exit 1; }
				log "process: TRACE conda environment activated: $trace_conda_env"
				return 0
			fi
		done
		log "ERROR TRACE_CONDA_ENV=$trace_conda_env was set, but conda was not found"
		exit 1
	fi

}

trace_need_cmds(){
	local x miss=0
	trace_activate_env
	for x in awk sort find gzip python3 trace-extract trace-infer trace-summarize; do
		if ! command -v "$x" >/dev/null 2>&1; then
			log "ERROR TRACE missing command: $x"
			miss=1
		fi
	done
	[[ -s "$sample_file" ]] || { log "ERROR TRACE missing sample_file=$sample_file"; miss=1; }
	[[ -s "$trace_py" ]] || { log "ERROR TRACE missing helper script: $trace_py"; miss=1; }
	[[ -d "$trace_dirarg" ]] || { log "ERROR TRACE missing ARG directory: $trace_dirarg"; miss=1; }
	(( miss == 0 )) || exit 1
}

trace_chrom_label(){ echo "chr$1"; }

trace_tree_files_for_chr(){
	local chr=$1
	find "$trace_dirarg" -type f \( -name '*.tsz' -o -name '*.trees' \) 2>/dev/null |
	awk -v c="$chr" '
		BEGIN{cc=tolower(c)}
		{
			x=tolower($0)
			if (x ~ "/chr" cc "/" || x ~ "(^|[^[:alnum:]])chr" cc "([^[:alnum:]]|$)") print
		}' | sort
}

trace_include_file_for_chr(){
	local chr=$1 f
	local include_root=$trace_include_dir
	# ARG coordinate chunks carry an automatically generated accessibility BED
	# so the unused span outside each chunk cannot become an inferred segment.
	if [[ -z $include_root && $chr == *_part* ]]; then include_root=$trace_dirarg; fi
	[[ -n $include_root && -d $include_root ]] || return 1
	for f in "$include_root/chr$chr.bed" "$include_root/chr$chr.include.bed" "$include_root/$chr.bed"; do
		[[ -s $f ]] && { echo "$f"; return 0; }
	done
	return 1
}

trace_genetic_map_for_chr(){
	local chr=$1 f
	[[ -n $trace_genetic_map_dir && -d $trace_genetic_map_dir ]] || return 1
	for f in "$trace_genetic_map_dir/genetic_map_hg38_chr$chr.txt" "$trace_genetic_map_dir/genetic_map_hg19_chr$chr.txt" "$trace_genetic_map_dir/chr$chr.txt" "$trace_genetic_map_dir/chr$chr.map"; do
		[[ -s $f ]] && { echo "$f"; return 0; }
	done
	return 1
}

trace_write_manifest(){
	local chr f i n total=0 mode0="" mode
	mkdir -p "$trace_dirout"/{manifest,samples,extract,datafiles,infer,summary,final,report,log,tmp}
	: > "$trace_dirout/manifest/tree_files.tsv"
	printf "chr\tchrom\tposterior_index\ttree_file\n" > "$trace_dirout/manifest/tree_files.tsv"
	for chr in "${trace_chrs[@]}"; do
		i=0
		while IFS= read -r f; do
			[[ -n $f ]] || continue
			i=$((i + 1))
			printf "%s\tchr%s\t%s\t%s\n" "$chr" "$chr" "$i" "$f" >> "$trace_dirout/manifest/tree_files.tsv"
		done < <(trace_tree_files_for_chr "$chr")
		n=$i
		if (( n == 0 )); then
			log "ERROR TRACE no .trees/.tsz found for chr$chr under $trace_dirarg"
			exit 1
		fi
		mode=$([[ $n -gt 1 ]] && echo posterior || echo single)
		[[ -z $mode0 ]] && mode0=$mode
		if [[ $mode != "$mode0" ]]; then
			log "WARN TRACE mixed single/posterior chromosome inputs detected; trace-infer will use --data-files for all chromosomes"
		fi
		if (( n > 1 )) && [[ -z $trace_window_size ]]; then
			log "ERROR TRACE chr$chr has $n posterior tree files. Set TRACE_WINDOW_SIZE, e.g. TRACE_WINDOW_SIZE=10000 ./gu.sh trace ./ibdmix.sh"
			exit 1
		fi
		total=$((total + n))
		log "process: TRACE setup ARG input chr=$chr tree_files=$n mode=$mode"
	done
	log "process: TRACE setup ARG manifest written files=$total output=$trace_dirout/manifest/tree_files.tsv"
}

trace_first_tree_file(){
	awk 'BEGIN{FS="\t"} NR==2{print $4; exit}' "$trace_dirout/manifest/tree_files.tsv"
}

trace_prepare_sample_map(){
	local first mapout
	mapout=$trace_dirout/samples/trace_sample_map.tsv
	mkdir -p "$trace_dirout/samples"
	if [[ -n $trace_sample_map ]]; then
		[[ -s $trace_sample_map ]] || { log "ERROR TRACE_SAMPLE_MAP does not exist or is empty: $trace_sample_map"; exit 1; }
		python3 "$trace_py" normalize-sample-map --sample-map "$trace_sample_map" --sample-file "$sample_file" --out "$mapout"
	else
		first=$(trace_first_tree_file)
		python3 "$trace_py" auto-sample-map --tree-file "$first" --sample-file "$sample_file" --out "$mapout"
	fi
	awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} $1 ~ /^[0-9]+$/ {print}' "$mapout" > "$trace_dirout/samples/trace_sample_map.clean.tsv"
	mv "$trace_dirout/samples/trace_sample_map.clean.tsv" "$mapout"
	awk 'BEGIN{FS=OFS="\t"} NR>1{print $1}' "$mapout" > "$trace_dirout/samples/tree_nodes.txt"
	awk 'BEGIN{FS=OFS="\t"} NR>1{n++; if($4=="" || $4=="NA") miss++} END{print "process: TRACE setup sample map tree_nodes=" n " missing_1kg_sample_annotation=" miss+0}' "$mapout" | while read -r x; do log "$x"; done
	log "process: TRACE setup sample map output=$mapout"
}

trace_split_groups(){
	local list="$trace_dirout/samples/tree_nodes.txt" gdir="$trace_dirout/samples/groups"
	rm -rf "$gdir"; mkdir -p "$gdir"
	awk -v n="$trace_group_size" -v out="$gdir" 'NF{g=int((NR-1)/n)+1; print $1 >> sprintf("%s/group%04d.nodes", out, g)}' "$list"
	find "$gdir" -type f -name 'group*.nodes' | sort > "$trace_dirout/samples/group_files.list"
	while IFS= read -r f; do paste -sd, "$f" > "${f%.nodes}.comma"; done < "$trace_dirout/samples/group_files.list"
	log "process: TRACE setup split tree nodes into groups group_size=$trace_group_size groups=$(count_lines "$trace_dirout/samples/group_files.list")"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Trace input: trace-extract ARG observations by chromosome and node group
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trace_run_extract_one(){
	local chr=$1 group_file=$2 gname=$3 status=0 f idx out prefix opt inc chrom datafile
	chrom=$(trace_chrom_label "$chr")
	datafile=$trace_dirout/datafiles/$gname.$chrom.data.txt
	: > "$datafile"
	inc=""
	if f=$(trace_include_file_for_chr "$chr"); then inc="--include-regions $f --chrom $chrom"; fi
	while IFS=$'\t' read -r _chr _chrom idx f; do
		[[ $_chr == "$chr" ]] || continue
		prefix=$trace_dirout/extract/$chrom/$gname.$chrom.p$idx
		out=$prefix.npz
		mkdir -p "$(dirname "$prefix")"
		if [[ ! -s $out ]]; then
			opt=""
			[[ -n $trace_window_size ]] && opt="--window-size $trace_window_size"
			log "process: Trace input trace-extract start chr=$chr group=$gname posterior=$idx tree=$f"
			if ! trace-extract --tree-file "$f" -t "$trace_t_archaic" --individuals "$(cat "${group_file%.nodes}.comma")" $opt $inc -o "$prefix" > "$trace_dirout/log/extract.$gname.$chrom.p$idx.log" 2>&1; then
				log "ERROR Trace input trace-extract failed chr=$chr group=$gname posterior=$idx log=$trace_dirout/log/extract.$gname.$chrom.p$idx.log"
				status=1
				continue
			fi
		fi
		[[ -s $out ]] || { log "ERROR Trace input missing trace-extract npz=$out"; status=1; continue; }
		echo "$out" >> "$datafile"
	done < "$trace_dirout/manifest/tree_files.tsv"
	(( status == 0 )) || return 1
	log "process: Trace input trace-extract done chr=$chr group=$gname npz_files=$(count_lines "$datafile") datafile=$datafile"
}

trace_extract(){
	local chr gf gname running=0 status=0
	trace_write_manifest
	trace_prepare_sample_map
	trace_split_groups
	while IFS= read -r gf; do
		gname=$(basename "$gf" .nodes)
		for chr in "${trace_chrs[@]}"; do
			trace_run_extract_one "$chr" "$gf" "$gname" &
			running=$((running + 1))
			if (( running >= trace_job_extract )); then wait -n || status=1; running=$((running - 1)); fi
		done
	done < "$trace_dirout/samples/group_files.list"
	while (( running > 0 )); do wait -n || status=1; running=$((running - 1)); done
	(( status == 0 )) || { log "ERROR one or more Trace input trace-extract jobs failed"; exit 1; }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Trace genotype: trace-infer archaic posterior probabilities per haplotype
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trace_datafile_string_for_group(){
	local gname=$1 chr arr=()
	for chr in "${trace_chrs[@]}"; do arr+=("$trace_dirout/datafiles/$gname.$(trace_chrom_label "$chr").data.txt"); done
	join_comma "${arr[@]}"
}

trace_chrom_string(){
	local chr arr=()
	for chr in "${trace_chrs[@]}"; do arr+=("$(trace_chrom_label "$chr")"); done
	join_comma "${arr[@]}"
}

trace_gmap_string(){
	local chr f arr=()
	for chr in "${trace_chrs[@]}"; do
		if f=$(trace_genetic_map_for_chr "$chr"); then arr+=("$f"); else return 1; fi
	done
	join_comma "${arr[@]}"
}

trace_run_infer_one(){
	local node=$1 gname=$2 outprefix=$trace_dirout/infer/hap$node gmaps_arg="" gf chroms dfs
	chroms=$(trace_chrom_string)
	dfs=$(trace_datafile_string_for_group "$gname")
	if gf=$(trace_gmap_string); then gmaps_arg="--genetic-maps $gf"; fi
	if [[ -s "$outprefix.$(trace_chrom_label "${trace_chrs[-1]}").xss.npz" ]]; then
		log "process: Trace genotype trace-infer skip existing node=$node"
		return 0
	fi
	log "process: Trace genotype trace-infer start node=$node group=$gname"
	if ! trace-infer -i "$node" --data-files "$dfs" --chroms "$chroms" $gmaps_arg -o "$outprefix" > "$trace_dirout/log/infer.hap$node.log" 2>&1; then
		log "ERROR Trace genotype trace-infer failed node=$node log=$trace_dirout/log/infer.hap$node.log"
		return 1
	fi
	log "process: Trace genotype trace-infer done node=$node"
}

trace_infer(){
	local gf gname node running=0 status=0
	[[ -s "$trace_dirout/samples/group_files.list" ]] || { log "ERROR missing TRACE groups; run ./gu.sh trace"; exit 1; }
	while IFS= read -r gf; do
		gname=$(basename "$gf" .nodes)
		while IFS= read -r node; do
			[[ -n $node ]] || continue
			trace_run_infer_one "$node" "$gname" &
			running=$((running + 1))
			if (( running >= trace_job_infer )); then wait -n || status=1; running=$((running - 1)); fi
		done < "$gf"
	done < "$trace_dirout/samples/group_files.list"
	while (( running > 0 )); do wait -n || status=1; running=$((running - 1)); done
	(( status == 0 )) || { log "ERROR one or more Trace genotype trace-infer jobs failed"; exit 1; }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Trace genome scan: trace-summarize archaic inheritance segments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trace_xss_string_for_node(){
	local node=$1 chr arr=()
	for chr in "${trace_chrs[@]}"; do arr+=("$trace_dirout/infer/hap$node.$(trace_chrom_label "$chr").xss.npz"); done
	join_comma "${arr[@]}"
}

trace_run_summarize_one(){
	local node=$1 outprefix=$trace_dirout/summary/hap$node files chroms
	files=$(trace_xss_string_for_node "$node")
	chroms=$(trace_chrom_string)
	if [[ -s "$outprefix.summary.txt" ]]; then
		log "process: Trace genome scan trace-summarize skip existing node=$node"
		return 0
	fi
	log "process: Trace genome scan trace-summarize start node=$node"
	if ! trace-summarize -f "$files" -c "$chroms" \
		--posterior-threshold "$trace_posterior_threshold" \
		--physical-length-threshold "$trace_physical_length_threshold" \
		--genetic-distance-threshold "$trace_genetic_distance_threshold" \
		-o "$outprefix" > "$trace_dirout/log/summarize.hap$node.log" 2>&1; then
		log "ERROR Trace genome scan trace-summarize failed node=$node log=$trace_dirout/log/summarize.hap$node.log"
		return 1
	fi
	log "process: Trace genome scan trace-summarize done node=$node segments=$(data_rows "$outprefix.summary.txt")"
}

trace_segments(){
	local node running=0 status=0
	[[ -s "$trace_dirout/samples/tree_nodes.txt" ]] || { log "ERROR missing TRACE tree_nodes.txt; run ./gu.sh trace"; exit 1; }
	while IFS= read -r node; do
		[[ -n $node ]] || continue
		trace_run_summarize_one "$node" &
		running=$((running + 1))
		if (( running >= trace_job_summarize )); then wait -n || status=1; running=$((running - 1)); fi
	done < "$trace_dirout/samples/tree_nodes.txt"
	while (( running > 0 )); do wait -n || status=1; running=$((running - 1)); done
	(( status == 0 )) || { log "ERROR one or more Trace genome scan trace-summarize jobs failed"; exit 1; }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Trace combine: combine TRACE segments and build secondary-analysis reports
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trace_combine_outputs(){
	python3 "$trace_py" combine --root "$trace_dirout" --sample-file "$sample_file" --window-size "$trace_report_window_size" --stats-by-group "$stats_by_group" --genome-build GRCh37
	log "process: Trace combine combine output haplotype_segments=$trace_dirout/final/trace_haplotype_segments.tsv.gz rows=$(gzip -dc "$trace_dirout/final/trace_haplotype_segments.tsv.gz" | awk 'NR>1{n++} END{print n+0}')"
	log "process: Trace combine combine output person_segments=$trace_dirout/final/trace_person_segments.tsv.gz"
	log "process: Trace combine report=$trace_dirout/report/trace_summary.tsv all_percent=$trace_dirout/report/trace_all_percent.csv"
}

run_trace_all(){
	trace_need_cmds
	local window_label include_label gmap_label sample_label
	window_label=$trace_window_size
	include_label=$trace_include_dir
	gmap_label=$trace_genetic_map_dir
	sample_label=$trace_sample_map
	[[ -n $window_label ]] || window_label=none
	[[ -n $include_label ]] || include_label=none
	[[ -n $gmap_label ]] || gmap_label=none
	[[ -n $sample_label ]] || sample_label=auto_from_tree_metadata
	log "TRACE secondary run output=$trace_dirout"
	log "TRACE ARG directory=$trace_dirarg"
	log "TRACE sample_file=$sample_file"
	log "TRACE chromosomes=${trace_chrs[*]}"
	log "TRACE t_archaic=$trace_t_archaic window_size=$window_label"
	log "TRACE thresholds posterior=$trace_posterior_threshold physical_bp=$trace_physical_length_threshold genetic_cM=$trace_genetic_distance_threshold"
	log "TRACE parallel job_extract=$trace_job_extract job_infer=$trace_job_infer job_summarize=$trace_job_summarize group_size=$trace_group_size"
	log "TRACE optional include_dir=$include_label genetic_map_dir=$gmap_label sample_map=$sample_label"
	trace_extract
	trace_infer
	trace_segments
	trace_combine_outputs
}

run_trace_dispatch(){
	local trace_action=$1
	trace_need_cmds
	case "$trace_action" in
		trace_run) run_trace_all ;;
		trace_infer) log "process: TRACE reuse extracted observations"; trace_infer; trace_segments; trace_combine_outputs ;;
		trace_segments) log "process: TRACE reuse inference output"; trace_segments; trace_combine_outputs ;;
		trace_report) trace_combine_outputs ;;
		*) log "ERROR unknown TRACE action=$trace_action"; exit 1 ;;
	esac
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Main workflow and command dispatch
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_summary_and_igv(){
	combine_scan_outputs
	build_final_calls
	build_ibdmix_reports
	build_igv_files
}

if ! is_trace_action; then
	log "IBDmix run output=$dirout"

	run_manifest=$dirout/run.config.tsv
	current_manifest=$dirout/tmp.run.config.$$
	{
    printf 'dirmod\t%s\n' "$dirmod"
    printf 'sample_file\t%s:%s:%s\n' "$sample_file" "$(stat -c %s "$sample_file")" "$(stat -c %Y "$sample_file")"
    printf 'analysis_units\t%s\n' "${analysis_units[*]}"
    printf 'genome_build\t%s\n' "$genome_build"
    printf 'refs\t%s\n' "${refs[*]}"
    printf 'lod_cut\t%s\nlen_cut\t%s\n' "$lod_cut" "$len_cut"
    } > "$current_manifest"
	if [[ -s $run_manifest ]] && ! cmp -s "$current_manifest" "$run_manifest"; then
    log "WARN run configuration changed; removing stale genotype/scan/report outputs to prevent cross-target reuse"
    rm -rf "$dirout/genotype" "$dirout/raw" "$dirout/summary" "$dirout/combined" "$dirout/final" "$dirout/report" "$dirout/igv"
    mkdir -p "$dirout"/{genotype,raw,summary,combined,final,report,igv}
    rm -f "$dirout"/.genome_scan.complete "$dirout"/.combined.complete "$dirout"/.final_calls.complete "$dirout"/.report.complete
	fi
	mv -f "$current_manifest" "$run_manifest"
	log "modern target VCF folder=$dirmod/vcf"
	log "sample_file=$sample_file"
	log "requested chromosomes=${chrs[*]} analysis_units=${analysis_units[*]} genome_build=$genome_build"
	log "IBDmix options LOD=$emit_lod_cut minor_allele_count=$minor_allele_count archaic_error=$archaic_error modern_error_max=$modern_error_max modern_error_proportion=$modern_error_proportion"
	log "summary/filter thresholds LOD=$lod_cut length=$len_cut window_size=$window_size locus_min_cover=$locus_min_cover; final=Altai minus Denisova calls from African populations; genomewide_pngs=disabled_use_IGV"
	log "diagnostic settings stats_by_group=$stats_by_group admixed_african_pops=$admixed_african_pops; all_percent includes ANY_REF cross-reference-deduplicated burden"
	log "refs=${refs[*]} job_of_chr=$job_of_chr job_in_chr=$job_in_chr action=$action"
else
	log "TRACE run output=$trace_dirout chromosomes=${trace_chrs[*]} action=$action"
fi

case "$action" in
	ibdmix_run)
		skip_completed=1
		run_genome_scan
		run_summary_and_igv
		;;
	ibdmix_report)
		log "process: reuse saved genome-scan outputs under $dirout"
		run_summary_and_igv
		;;
	ibdmix_compare)
		log "process: reuse final calls under $dirout/final"
		build_ibdmix_reports
		build_igv_files
		;;
	ibdmix_check)
        prepare_samples
        validate_inputs
        log "INPUT CHECK PASSED: target_samples=$(count_lines "$dirout/samples/ALL.txt") populations=$(pop_list_for_ref Altai | wc -l) analysis_units=${analysis_units[*]}"
        missing_stages=()
        [[ -e "$dirout/.genome_scan.complete" ]] || missing_stages+=(genome_scan)
        [[ -e "$dirout/.combined.complete" ]] || missing_stages+=(combined)
        [[ -e "$dirout/.final_calls.complete" ]] || missing_stages+=(final_calls)
        [[ -e "$dirout/.report.complete" ]] || missing_stages+=(report)
        if ((${#missing_stages[@]})); then
            log "PIPELINE INCOMPLETE: missing_stages=${missing_stages[*]}; run ./gu.sh ibdmix to resume"
        else
            log "PIPELINE COMPLETE: genome_scan combined final_calls report"
        fi
        ;;
	ibdmix_igv)
		log "process: generate IGV inputs only; reuse report outputs under $dir_report"
		build_igv_files
		;;
	trace_run|trace_infer|trace_segments|trace_report)
		run_trace_dispatch "$action"
		;;
	*)
		log "ERROR unknown action=$action"
		exit 1
		;;
esac

if [[ $action == ibdmix_check ]]; then
    log "CHECK DONE: output=$dirout log=$dirout/ibdmix.log"
else
    log "ALL DONE: output=$dirout log=$dirout/ibdmix.log igv=$dir_igv"
fi
completed=1
