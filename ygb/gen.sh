#!/bin/bash -l
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd)/output_dirs.sh"
script_output_directory_args "${BASH_SOURCE[0]}" "$@" || { return 2 2>/dev/null || exit 2; }
set -- "${SCRIPT_OUTPUT_DIRECTORY_ARGS[@]}"

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
./gen.sh  # 默认 pgs -> snp -> apoe；必须补齐各自输入参数
./gen.sh MODULE  # all|pgs|snp|apoe|liftOver|make-pfile|lift-bed
./gen.sh pgs --input FILE --label NAME  # 必填，无通用文件默认值

Modules:
    make-pfile  Convert VCF/BGEN/BED to PGEN
    lift-bed    Convert compressed BED interval coordinates
    liftOver  Convert genotype coordinates
    pgs       Calculate polygenic scores
    snp       Extract selected variants
    apoe      Derive APOE genotypes

Conversion options (make-pfile / lift-bed):
    --dir-input DIR --dir-pfile DIR --input-format auto|vcf|bgen|bed
    --grch 37|38 --out-files pgen,bed --split-sex TRUE|FALSE
    --chr CSV --input-prefix PREFIX --sample-file FILE --replace TRUE|FALSE
    make-pfile defaults to ${DIR_GEN}/1kg/37/vcf and sibling pfile/.
    lift-bed requires --dir-in DIR --dir-out DIR --chain-file FILE.
    --bgzip-threads N  BED compression threads [4]
    --plan-only --job-hours N --queue NAME

Main options:
    --dir-in DIR --dir-out DIR  Input/output
    --dir-gen DIR --data NAME  Genotype dataset
    --dir-runtime DIR  liftOver logs/tmp root [OUTPUT_DIR]; keeps runtime files out of genotype storage
    --chain-file FILE  Liftover chain
    --data NAME  liftOver: process NAME.pvar with matching .pgen/.psam (e.g. ancient)
    --chr CHR  Process only chrCHR.pvar (e.g. 20); default: all chromosomes
    --ref-fasta FILE  Indexed GRCh38 FASTA [${DIR_GEN}/1hgp/38/GRCH38.fasta.gz]
    liftOver preserves row count/order and all fields except POS; failures get POS=-1.
    Original .pgen/.psam are reused through relative symlinks.
    --input FILE --label NAME  PGS input and label
    --export-allele FILE  SNP export counted alleles (variant ID and allele)
    --dir-pgs-pfile DIR  PGS genotypes [${DIR_GEN}/ygb/37/imp]
    --submit local|hpc|auto  Auto-detect HPC/local using shared 0f/0os.sh [auto]; optional override
    --write-cmd / --dry-run  Write the .cmd only; --run-cmd FALSE is equivalent
    --queue NAME --slots N --memory-gb N --walltime HHH:MM  LSF resources
    --threads N  Analysis threads [4]; slots round up by the shared CPU/RAM policy
    --dir-jobs DIR  Command/log directory [OUTPUT_DIR/logs]
    --wait TRUE|FALSE  Wait for verified LSF completion [FALSE]
    Final outputs: snp -> DIR_GEN/ygb/37/imp/snp.raw.gz; pgs -> DIR_GEN/ygb/37/pgs/LABEL.pgs.txt
    APOE output: DIR_GEN/ygb/37/hap/apoe.hap (also in its work directory).
    Logs/intermediates: ANALYSIS_ROOT/ygb/gen/{snp,pgs/LABEL,apoe}; --dir-out overrides the work directory.
    SNP list defaults to $DIR_YGB/common/snp.lst; column 2 + counted allele names the exported columns.
    --background  Legacy alias for --submit hpc
    YGB_PYTHON  Python override [SHARED_PYTHON from 0f/0os.sh]
    -h, --help  Show options and examples below

Permanent output directories (legacy options remain accepted):
  --dir-out-pfile DIR
                                  make-pfile: pfiles [input VCF sibling pfile/].
  --dir-out-lift_bed DIR
                                  lift-bed: destination directory (required).
  --dir-out-liftover DIR
                                  liftover: destination directory (required).
  --dir-out-pgs DIR
                                  Published PGS tables [${DIR_GEN}/ygb/37/pgs].
  --dir-out-snp DIR
                                  Published SNP tables [${DIR_GEN}/ygb/37/imp].

Examples:
    ./gen.sh make-pfile --input-format bgen --dir-input "${DIR_GEN}/ygb/37/hap/raw" --dir-pfile "${DIR_GEN}/ygb/37/hap" --grch 37
    ./gen.sh lift-bed --dir-in "${DIR_GEN}/archaic/37/bed" --dir-out "${DIR_GEN}/archaic/38/bed" --chain-file "${DIR_FILES}/liftOver/hg19ToHg38.over.chain.gz"
    source ../0f/0os.sh
    ./gen.sh pgs --input "$DIR_GWAS/bmi.ref" --label bmi --write-cmd
    ./gen.sh pgs --input "$DIR_GWAS/bmi.ref" --label bmi
    ./gen.sh snp --dir-gen "$DIR_GEN/ygb/37" --data imp --snp-file "$DIR_FILES/Wang_MY.bed" --snp-file-type bed --dir-out "$ANALYSIS_ROOT/ygb" --label Wang
PUBLIC_ENTRY_HELP
    exit 0
fi
# These modules retain their shared project workflow and resource options.
run_conversion_module() (
# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
submit=${submit:-${SUBMIT:-auto}}
source "$_paths_directory/0f/0os.sh"
unset _paths_directory

source "$SCRIPTS_ROOT/0f/project_entry.sh"
project_workflow_init "$@" || exit $?
set -- "${WORKFLOW_ARGS[@]}"

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "${SCRIPTS_ROOT}/0f/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
source "${SCRIPTS_ROOT}/0f/console.sh"
set -euo pipefail

dir0=${dir0:-${DIR_DAT}}
dirfunc=${dirfunc:-$SCRIPTS_ROOT/0f}


usage() { bash "${BASH_SOURCE[0]}" --help; }

# 🚩 Shared helpers
die() {
	echo "ERROR: $*" >&2
	exit 1
}

is_true() {
	case "${1^^}" in
		TRUE | T | YES | Y | 1) return 0 ;;
		FALSE | F | NO | N | 0) return 1 ;;
		*) die "expected TRUE or FALSE, got: $1" ;;
	esac
}

split_csv() {
	local input=$1
	local -n output_ref=$2
	IFS=',' read -r -a output_ref <<<"$input"
	((${#output_ref[@]} > 0)) || die "empty comma-separated option"
}

pfile_exists() {
	local prefix=$1
	[[ -s ${prefix}.pgen && -s ${prefix}.psam && -s ${prefix}.pvar ]]
}

pfile_input_args() {
	local prefix=$1
	local -n output_ref=$2
	[[ -s ${prefix}.pvar ]] || die "missing plain .pvar for: $prefix"
	output_ref=(--pfile "$prefix")
}

bed_exists() {
	local prefix=$1
	[[ -s ${prefix}.bed && -s ${prefix}.bim && -s ${prefix}.fam ]]
}

run_make_pfile() {
	local output=$1
	shift
	if pfile_exists "$output" && ! is_true "$replace"; then
		echo "SKIP pgen: $output already has .pvar (use --replace TRUE)"
		return
	fi
	"$plink2_bin" "$@" --make-pgen psam-cols=+fid --out "$output"
	pfile_exists "$output" || die "incomplete pfile output: $output"
}

run_make_pfile_replacing() {
	# Dynamic scoping makes run_make_pfile see this local override. This is
	# used only to repair an existing pfile with obsolete chromosome content.
	local replace=TRUE
	run_make_pfile "$@"
}

write_sex_update() {
	local output=$1 tmp=${1}.tmp.$$
	awk 'BEGIN{OFS="\t"}
        NR==1 {next}
        {
            s=tolower($4)
            if (s=="male" || s=="m") s=1
            else if (s=="female" || s=="f") s=2
            if (s!=1 && s!=2) {
                print "ERROR: invalid sex for sample " $1 ": " $4 > "/dev/stderr"
                bad=1; next
            }
            print $1,$1,s
        }
        END {exit bad}
    ' "$sample_file" >"$tmp"
	mv -f -- "$tmp" "$output"
}

normalize_chromosomes() {
	local input=$1 raw chr
	local -A seen=()
	local values=()
	split_csv "$input" values
	selected_chrs=()
	for raw in "${values[@]}"; do
		chr=${raw#chr}
		chr=${chr#CHR}
		chr=${chr^^}
		if [[ $chr == AUTOSOME || $chr == AUTOSOMES ]]; then
			for chr in {1..22}; do seen[$chr]=1; done
			continue
		fi
		[[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == XY || $chr == Y || $chr == MT ]] ||
			die "invalid chromosome: $raw"
		seen[$chr]=1
	done
	mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

discover_vcf_chromosomes() {
	local file base chr
	local -A seen=()
	shopt -s nullglob
	for file in "${input_prefix}"*.vcf.gz; do
		base=${file#"$input_prefix"}
		chr=${base%.vcf.gz}
		chr=${chr^^}
		[[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == Y || $chr == MT ]] && seen[$chr]=1
	done
	shopt -u nullglob
	((${#seen[@]} > 0)) || die "no chromosome VCFs found with prefix: $input_prefix"
	mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

discover_bed_chromosomes() {
	local file chr prefix
	local -A seen=()
	shopt -s nullglob
	for file in "${input_prefix}"*.bed; do
		chr=${file#"$input_prefix"}
		chr=${chr%.bed}
		chr=${chr^^}
		[[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == XY || $chr == Y || $chr == MT ]] || continue
		prefix=${file%.bed}
		[[ -s ${prefix}.bim && -s ${prefix}.fam ]] ||
			die "incomplete BED/BIM/FAM input fileset: $prefix"
		seen[$chr]=1
	done
	shopt -u nullglob
	((${#seen[@]} > 0)) || die "no complete BED/BIM/FAM inputs found with prefix: $input_prefix"
	mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

bgen_file_for_chr() {
	local chr=$1 exact=${input_prefix}${1}.bgen
	local candidates=()
	if [[ -s $exact ]]; then
		printf '%s\n' "$exact"
		return
	fi
	shopt -s nullglob
	candidates=("${input_prefix}${chr}"_*.bgen)
	shopt -u nullglob
	((${#candidates[@]} == 1)) ||
		die "expected exactly one BGEN for chr$chr under $input_dir; found ${#candidates[@]}"
	printf '%s\n' "${candidates[0]}"
}

bgen_sample_file() {
	local bgen=$1 same_prefix=${1%.bgen}.sample
	local candidates=()
	if [[ -s $same_prefix ]]; then
		printf '%s\n' "$same_prefix"
		return
	fi
	shopt -s nullglob
	candidates=("${bgen%.bgen}"_s*.sample)
	shopt -u nullglob
	((${#candidates[@]} == 1)) ||
		die "expected exactly one SAMPLE file for $bgen; found ${#candidates[@]}"
	printf '%s\n' "${candidates[0]}"
}

discover_bgen_chromosomes() {
	local chr bgen sample exact
	local candidates=()
	local -A seen=()
	for chr in {1..22} X XY Y MT; do
		exact=${input_prefix}${chr}.bgen
		if [[ -s $exact ]]; then
			seen[$chr]=1
			continue
		fi
		shopt -s nullglob
		candidates=("${input_prefix}${chr}"_*.bgen)
		shopt -u nullglob
		((${#candidates[@]} <= 1)) ||
			die "multiple BGEN inputs match chr$chr and prefix $input_prefix"
		((${#candidates[@]} == 1)) && seen[$chr]=1
	done
	((${#seen[@]} > 0)) || die "no BGEN inputs found with prefix: $input_prefix"
	mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
	for chr in "${selected_chrs[@]}"; do
		bgen=$(bgen_file_for_chr "$chr")
		sample=$(bgen_sample_file "$bgen")
		[[ -s $sample ]] || die "missing SAMPLE file for chr$chr"
	done
}

detect_input_format() {
	local files=()
	case $input_format in
		vcf | bgen | bed) return ;;
		auto) ;;
		*) die "--input-format must be auto, vcf, bgen, or bed" ;;
	esac
	shopt -s nullglob
	files=("${input_prefix}"*.vcf.gz)
	if ((${#files[@]})); then
		input_format=vcf
		shopt -u nullglob
		return
	fi
	files=("${input_prefix}"*.bgen)
	if ((${#files[@]})); then
		input_format=bgen
		shopt -u nullglob
		return
	fi
	files=("${input_prefix}"*.bed)
	if ((${#files[@]})); then
		input_format=bed
		shopt -u nullglob
		return
	fi
	shopt -u nullglob
	die "cannot detect VCF, BGEN, or BED input under $input_dir"
}

discover_input_chromosomes() {
	case $input_format in
		vcf) discover_vcf_chromosomes ;;
		bgen) discover_bgen_chromosomes ;;
		bed) discover_bed_chromosomes ;;
	esac
}

discover_pfile_chromosomes() {
	local file base chr
	local -A seen=()
	shopt -s nullglob
	for file in "$pfile_dir"/chr*.pgen; do
		base=${file##*/}
		chr=${base#chr}
		chr=${chr%.pgen}
		chr=${chr^^}
		[[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == Y || $chr == MT ]] || continue
		pfile_exists "$pfile_dir/chr$chr" && seen[$chr]=1
	done
	shopt -u nullglob
	((${#seen[@]} > 0)) || die "no complete chr*.pgen/.pvar/.psam filesets under $pfile_dir"
	mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

detect_genome_build() {
	local vcf=$1 line detected=""
	while IFS= read -r line; do
		case $line in
			*assembly=b37* | *GRCh37* | *hs37d5* | *hg19* | *length=155270560*)
				detected=37
				break
				;;
			*assembly=b38* | *GRCh38* | *hg38* | *length=156040895*)
				detected=38
				break
				;;
			\#CHROM*) break ;;
		esac
	done < <(phe_zcat "$vcf")
	[[ $detected == 37 || $detected == 38 ]] ||
		die "cannot infer GRCh build from VCF header: $vcf (use --grch 37 or 38)"
	printf '%s\n' "$detected"
}

set_par_coordinates() {
	local build=$1
	# X PAR variables are supplied by the sourced phenotype.sh library.
	# shellcheck disable=SC2154
	if [[ $build == 37 ]]; then
		x_par_head_end=$X_PAR1_END_b37
		x_par_tail_start=$X_PAR2_START_b37
		y_par1_start=10001
		y_par1_end=2649520
		y_par2_start=59034050
		y_par2_end=59363566
		y_length=59373566
	else
		x_par_head_end=$X_PAR1_END_b38
		x_par_tail_start=$X_PAR2_START_b38
		y_par1_start=10001
		y_par1_end=2781479
		y_par2_start=56887903
		y_par2_end=57217415
		y_length=57227415
	fi
}

pvar_file() {
	local prefix=$1
	[[ -s ${prefix}.pvar ]] || die "missing plain .pvar for: $prefix"
	printf '%s\n' "${prefix}.pvar"
}

read_pvar() {
	cat -- "$1"
}

validate_pvar_chromosomes() {
	local prefix=$1 expected=$2 pvar
	pvar=$(pvar_file "$prefix")
	if ! pvar_chromosomes_match "$pvar" "$expected"; then
		die "$prefix failed the post-generation pure-$expected validation"
	fi
}

pvar_chromosomes_match() {
	local pvar=$1 expected=$2 build=${chrX_build:-${grch:-}}
	if [[ $expected == X ]]; then
		[[ $build == 37 || $build == 38 ]] || die "X validation requires --grch 37 or 38"
		set_par_coordinates "$build"
	fi
	read_pvar "$pvar" | awk -v expected="$expected" -v head="${x_par_head_end:-0}" -v tail="${x_par_tail_start:-0}" '
        /^#/ {next}
        {
            count++
            if (expected == "X" && ($1 != "X" || $2 <= head || $2 >= tail)) bad=1
            if (expected == "Y" && $1 != "Y") bad=1
            if (expected == "PAR" && $1 != "XY" && $1 != "PAR1" && $1 != "PAR2") bad=1
        }
        END {if (!count) exit 2; exit bad}
    '
}

pvar_contains_par() {
	local pvar=$1
	read_pvar "$pvar" | awk '
        /^#/ {next}
        $1 == "PAR1" || $1 == "PAR2" {found=1}
        END {exit !found}
    '
}

validate_y_nonpar() {
	local prefix=$1 pvar
	pvar=$(pvar_file "$prefix")
	if ! read_pvar "$pvar" | awk \
		-v p1s="$y_par1_start" -v p1e="$y_par1_end" \
		-v p2s="$y_par2_start" -v p2e="$y_par2_end" '
        /^#/ {next}
        {
            count++
            if ($1 != "Y" || ($2 >= p1s && $2 <= p1e) ||
                ($2 >= p2s && $2 <= p2e)) bad=1
        }
        END {if (!count) exit 2; exit bad}
    '; then
		die "$prefix contains Y-PAR variants; use --replace TRUE after checking old outputs"
	fi
}

make_chrX_male_stats() {
	local prefix=$pfile_dir/chrX.male pvar summary_tmp
	local variant_count min_pos max_pos sample_count nonmale_count hethap_count
	local stats_complete=0 stats_build=unknown
	local pfile_args=()
	validate_pvar_chromosomes "$prefix" X
	if [[ -n ${chrX_build:-} ]]; then
		stats_build=GRCh$chrX_build
	elif [[ -n $grch ]]; then
		stats_build=GRCh$grch
	fi

	if [[ -s ${prefix}.afreq && -s ${prefix}.gcount && -s ${prefix}.smiss &&
		-s ${prefix}.basic-stats.txt ]]; then
		stats_complete=1
	fi
	if ((stats_complete)) && ! is_true "$replace"; then
		echo "SKIP chrX.male statistics: reports already exist (use --replace TRUE)"
		return
	fi

	echo "STATS chrX.male (--freq, --geno-counts, --missing sample-only)"
	pfile_input_args "$prefix" pfile_args
	"$plink2_bin" "${pfile_args[@]}" --freq --out "$prefix"
	"$plink2_bin" "${pfile_args[@]}" --geno-counts --out "$prefix"
	"$plink2_bin" "${pfile_args[@]}" --missing sample-only \
		'scols=maybefid,maybesid,nmiss,nmisshh,hethap,nobs,fmiss,fmisshh' \
		--out "$prefix"
	[[ -s ${prefix}.afreq && -s ${prefix}.gcount && -s ${prefix}.smiss ]] ||
		die "incomplete chrX.male statistics under: $pfile_dir"

	pvar=$(pvar_file "$prefix")
	read -r variant_count min_pos max_pos < <(
		read_pvar "$pvar" | awk '
            /^#/ {next}
            NR == 1 || !count {min=$2}
            {count++; max=$2}
            END {print count+0, min+0, max+0}
        '
	)
	read -r sample_count nonmale_count < <(
		awk '
            /^#/ {
                for (i=1; i<=NF; i++) if ($i == "SEX") sex_col=i
                next
            }
            NF {count++; if (!sex_col || $sex_col != 1) nonmale++}
            END {print count+0, nonmale+0}
        ' "${prefix}.psam"
	)
	hethap_count=$(awk '
        NR == 1 {
            for (i=1; i<=NF; i++) if ($i == "HETHAP_CT") col=i
            next
        }
        col {sum += $col}
        END {if (!col) print "NA"; else print sum+0}
    ' "${prefix}.smiss")
	[[ $variant_count -gt 0 && $sample_count -gt 0 && $nonmale_count -eq 0 ]] ||
		die "chrX.male validation failed: variants=$variant_count samples=$sample_count nonmale=$nonmale_count"
	[[ $hethap_count == 0 ]] ||
		die "chrX.male validation failed: heterozygous-haploid count=$hethap_count"

	summary_tmp=${prefix}.basic-stats.txt.tmp.$$
	printf '%s\n' \
		"pfile=${prefix}" \
		"build=${stats_build}" \
		"chromosome=X" \
		"par_regions_present=FALSE" \
		"male_samples=${sample_count}" \
		"nonmale_samples=${nonmale_count}" \
		"variants=${variant_count}" \
		"min_position=${min_pos}" \
		"max_position=${max_pos}" \
		"plink2_encoding=male_chrX_haploid" \
		"heterozygous_haploid_calls=${hethap_count}" \
		"frequency_report=${prefix}.afreq" \
		"genotype_count_report=${prefix}.gcount" \
		"sample_missingness_report=${prefix}.smiss" \
		>"$summary_tmp"
	mv -f -- "$summary_tmp" "${prefix}.basic-stats.txt"
}

validate_psam_sex() {
	local prefix=$1 expected_code=$2 expected_label=$3
	if ! awk -v expected="$expected_code" '
        /^#/ {
            for (i=1; i<=NF; i++) if ($i == "SEX") sex_col=i
            next
        }
        NF {count++; if (!sex_col || $sex_col != expected) bad++}
        END {if (!count || bad) exit 1}
    ' "${prefix}.psam"; then
		die "$prefix failed the $expected_label-only sample validation"
	fi
}

make_chrX_sex_pfiles() {
	local input=$pfile_dir/chrX output sex keep_flag sex_code
	local pfile_args=()
	validate_pvar_chromosomes "$input" X
	pfile_input_args "$input" pfile_args
	for sex in male female; do
		output=$pfile_dir/chrX.$sex
		if [[ $sex == male ]]; then
			keep_flag=--keep-males
			sex_code=1
		else
			keep_flag=--keep-females
			sex_code=2
		fi
		if pfile_exists "$output" && ! is_true "$replace"; then
			echo "SKIP pgen: $output already has .pvar (use --replace TRUE)"
		else
			echo "MAKE pgen chrX.$sex ($sex-only pure non-PAR X)"
			if [[ $sex == male ]]; then
				run_make_pfile "$output" "${pfile_args[@]}" "$keep_flag" \
					--set-invalid-haploid-missing
			else
				run_make_pfile "$output" "${pfile_args[@]}" "$keep_flag"
			fi
		fi
		validate_pvar_chromosomes "$output" X
		validate_psam_sex "$output" "$sex_code" "$sex"
	done
	make_chrX_male_stats
}

make_chrX_pfiles() {
	local vcf=${input_prefix}X.vcf.gz x_output=$pfile_dir/chrX
	local xy_output=$pfile_dir/chrXY no_par_marker=$pfile_dir/chrXY.no-PAR.txt
	local build_source x_pvar xy_pvar need_import=0 x_valid=0 xy_valid=0
	[[ -s $vcf ]] || die "missing VCF: $vcf"

	if [[ -n $grch ]]; then
		chrX_build=$grch
		build_source="--grch"
	else
		chrX_build=$(detect_genome_build "$vcf")
		build_source="VCF header"
	fi
	set_par_coordinates "$chrX_build"

	if pfile_exists "$x_output"; then
		x_pvar=$(pvar_file "$x_output")
		if pvar_chromosomes_match "$x_pvar" X; then
			x_valid=1
		else
			echo "REBUILD pgen: existing chrX contains PAR or non-X variants"
		fi
	fi
	if pfile_exists "$xy_output"; then
		xy_pvar=$(pvar_file "$xy_output")
		if pvar_chromosomes_match "$xy_pvar" PAR; then
			xy_valid=1
		else
			echo "REBUILD pgen: existing chrXY contains non-PAR variants"
		fi
	fi

	((x_valid)) || need_import=1
	if ((! xy_valid)) && [[ ! -s $no_par_marker ]]; then need_import=1; fi
	is_true "$replace" && need_import=1
	if ((need_import)); then
		(
			local tmp_dir raw raw_pvar raw_has_par=0
			local raw_args=()
			tmp_dir=$(mktemp -d "$pfile_dir/.chrX.make-pfile.XXXXXX")
			trap 'rm -rf -- "$tmp_dir"' EXIT
			raw=$tmp_dir/chrX.split
			echo "Using GRCh$chrX_build for chrX PAR splitting ($build_source): $vcf"
			run_make_pfile "$raw" \
				--vcf "$vcf" --double-id --update-sex "$sex_update" \
				--split-par "$x_par_head_end" "$x_par_tail_start" \
				--set-missing-var-ids '@:#:$r:$a' \
				--new-id-max-allele-len 100 missing
			raw_pvar=$(pvar_file "$raw")
			if pvar_contains_par "$raw_pvar"; then raw_has_par=1; fi
			pfile_input_args "$raw" raw_args

			if ((! x_valid)) || is_true "$replace"; then
				echo "MAKE pgen chrX (non-PAR only)"
				if pfile_exists "$x_output"; then
					run_make_pfile_replacing "$x_output" "${raw_args[@]}" --chr X
				else
					run_make_pfile "$x_output" "${raw_args[@]}" --chr X
				fi
			else
				echo "SKIP pgen: existing chrX is already pure non-PAR X"
			fi

			if ((raw_has_par)); then
				rm -f -- "$no_par_marker"
				if ((! xy_valid)) || is_true "$replace"; then
					echo "MAKE pgen chrXY (PAR1 and PAR2)"
					if pfile_exists "$xy_output"; then
						run_make_pfile_replacing "$xy_output" "${raw_args[@]}" --chr PAR1 PAR2
					else
						run_make_pfile "$xy_output" "${raw_args[@]}" --chr PAR1 PAR2
					fi
				else
					echo "SKIP pgen: existing chrXY is already PAR1/PAR2 only"
				fi
			else
				printf 'source=%s\nbuild=GRCh%s\npar_variants=0\n' \
					"$vcf" "$chrX_build" >"$no_par_marker"
				echo "SKIP pgen chrXY: source chrX VCF has no PAR1/PAR2 variants"
			fi
		)
	else
		echo "SKIP pgen: existing chrX/chrXY outputs are already valid"
	fi
	validate_pvar_chromosomes "$x_output" X
	if pfile_exists "$xy_output"; then
		[[ ! -s $no_par_marker ]] ||
			die "chrXY exists even though the source VCF has no PAR variants"
		validate_pvar_chromosomes "$xy_output" PAR
	else
		[[ -s $no_par_marker ]] || die "chrXY was not generated and no no-PAR marker was written"
	fi
	if is_true "$split_sex"; then make_chrX_sex_pfiles; fi
}

make_chrY_pfile() {
	local vcf=${input_prefix}Y.vcf.gz output=$pfile_dir/chrY
	local chrY_build build_source
	[[ -s $vcf ]] || die "missing VCF: $vcf"
	if [[ -n $grch ]]; then
		chrY_build=$grch
		build_source="--grch"
	else
		chrY_build=$(detect_genome_build "$vcf")
		build_source="VCF header"
	fi
	set_par_coordinates "$chrY_build"
	if pfile_exists "$output" && ! is_true "$replace"; then
		echo "SKIP pgen: $output already has .pvar (use --replace TRUE)"
	else
		(
			local tmp_dir nonpar_bed
			tmp_dir=$(mktemp -d "$pfile_dir/.chrY.make-pfile.XXXXXX")
			trap 'rm -rf -- "$tmp_dir"' EXIT
			nonpar_bed=$tmp_dir/chrY.nonpar.bed1
			printf 'Y\t1\t%d\nY\t%d\t%d\nY\t%d\t%d\n' \
				"$((y_par1_start - 1))" \
				"$((y_par1_end + 1))" "$((y_par2_start - 1))" \
				"$((y_par2_end + 1))" "$y_length" >"$nonpar_bed"
			echo "Using GRCh$chrY_build Y-PAR boundaries ($build_source): $vcf"
			echo "MAKE pgen chrY (non-PAR only)"
			run_make_pfile "$output" \
				--vcf "$vcf" --double-id --update-sex "$sex_update" --chr Y \
				--extract bed1 "$nonpar_bed" --set-invalid-haploid-missing \
				--set-missing-var-ids '@:#:$r:$a' \
				--new-id-max-allele-len 100 missing
		)
	fi
	validate_pvar_chromosomes "$output" Y
	validate_y_nonpar "$output"
}

make_vcf_pgen_one() {
	local chr=$1 vcf=${input_prefix}${1}.vcf.gz output=$pfile_dir/chr${1}
	local args=(--vcf "$vcf" --double-id --update-sex "$sex_update"
		--set-missing-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing)
	if [[ $chr == X ]]; then
		make_chrX_pfiles
		return
	fi
	if [[ $chr == Y ]]; then
		make_chrY_pfile
		return
	fi
	if pfile_exists "$output" && ! is_true "$replace"; then
		echo "SKIP pgen: $output already has .pvar (use --replace TRUE)"
		return
	fi
	[[ -s $vcf ]] || die "missing VCF: $vcf"
	echo "MAKE pgen chr$chr"
	run_make_pfile "$output" "${args[@]}"
}

make_bgen_pgen_one() {
	local chr=$1 output=$pfile_dir/chr${1} bgen sample
	if [[ $chr == X ]]; then
		[[ $grch == 37 || $grch == 38 ]] || die "BGEN X import requires --grch 37 or 38 and pre-split non-PAR X"
	fi
	bgen=$(bgen_file_for_chr "$chr")
	sample=$(bgen_sample_file "$bgen")
	if pfile_exists "$output" && ! is_true "$replace"; then
		echo "SKIP pgen: $output already has .pvar (use --replace TRUE)"
	else
		echo "MAKE pgen chr$chr from BGEN"
		run_make_pfile "$output" --bgen "$bgen" ref-first --sample "$sample" \
			--oxford-single-chr "$chr"
	fi
}

make_bed_pgen_one() {
	local chr=$1 input=${input_prefix}${1} output=$pfile_dir/chr${1}
	if [[ $chr == X ]]; then
		[[ $grch == 37 || $grch == 38 ]] || die "BED X import requires --grch 37 or 38 and pre-split non-PAR X"
		set_par_coordinates "$grch"
		awk -v head="$x_par_head_end" -v tail="$x_par_tail_start" '$4 <= head || $4 >= tail {bad=1} END {exit bad}' "${input}.bim" || die "BED chrX contains PAR coordinates; split PAR before import"
	fi
	[[ -s ${input}.bed && -s ${input}.bim && -s ${input}.fam ]] ||
		die "incomplete BED/BIM/FAM input fileset: $input"
	if pfile_exists "$output" && ! is_true "$replace"; then
		echo "SKIP pgen: $output already has .pvar (use --replace TRUE)"
	else
		echo "MAKE pgen chr$chr from BED/BIM/FAM"
		run_make_pfile "$output" --bfile "$input" --indiv-sort none
	fi
}

make_pgen_one() {
	local chr=$1 output=$pfile_dir/chr${1}
	case $input_format in
		vcf)
			make_vcf_pgen_one "$chr"
			return
			;;
		bgen) make_bgen_pgen_one "$chr" ;;
		bed) make_bed_pgen_one "$chr" ;;
	esac
	case $chr in
		X)
			validate_pvar_chromosomes "$output" X
			if is_true "$split_sex"; then make_chrX_sex_pfiles; fi
			;;
		XY) validate_pvar_chromosomes "$output" PAR ;;
		Y) validate_pvar_chromosomes "$output" Y ;;
	esac
}

make_bed_prefix() (
	local label=$1 input=$2 tmp
	local pfile_args=()
	local temporary
	temporary=$(mktemp -d "$pfile_dir/.${label}.make-bed.XXXXXX")
	trap 'rm -rf -- "$temporary"' EXIT
	trap 'exit 130' INT
	trap 'exit 143' TERM
	trap 'exit 129' HUP
	tmp=$temporary/data
	pfile_exists "$input" || die "missing pfile for $label: $input (request --out-files pgen,bed to create it)"
	if bed_exists "$input" && ! is_true "$replace"; then
		echo "SKIP bed: $input already exists (use --replace TRUE)"
		return
	fi
	rm -f -- "$tmp".bed "$tmp".bim "$tmp".fam "$tmp".log "$tmp".nosex
	echo "MAKE bed  $label (multiallelic variants excluded)"
	pfile_input_args "$input" pfile_args
	"$plink2_bin" "${pfile_args[@]}" \
		--max-alleles 2 --output-chr 26 --make-bed --out "$tmp"
	[[ -s $tmp.bed && -s $tmp.bim && -s $tmp.fam ]] || die "incomplete bed output for $label"
	mv -f -- "$tmp.bed" "$input.bed"
	mv -f -- "$tmp.bim" "$input.bim"
	mv -f -- "$tmp.fam" "$input.fam"
	rm -f -- "$tmp.log" "$tmp.nosex"
)

make_bed_one() {
	local chr=$1
	if [[ $chr == X ]]; then
		validate_pvar_chromosomes "$pfile_dir/chrX" X
		make_bed_prefix chrX "$pfile_dir/chrX"
		if pfile_exists "$pfile_dir/chrXY"; then
			validate_pvar_chromosomes "$pfile_dir/chrXY" PAR
			make_bed_prefix chrXY "$pfile_dir/chrXY"
		else
			echo "SKIP bed chrXY: no chrXY pfile (source chrX VCF had no PAR variants)"
		fi
		if is_true "$split_sex"; then
			local sex
			for sex in male female; do
				validate_pvar_chromosomes "$pfile_dir/chrX.$sex" X
				make_bed_prefix "chrX.$sex" "$pfile_dir/chrX.$sex"
			done
		fi
	else
		make_bed_prefix "chr$chr" "$pfile_dir/chr$chr"
	fi
}


# 🚩 Genotype conversions
run_make_module() {
	local value
	local -A requested=()
	local outputs=()
	split_csv "$out_files_csv" outputs
	for value in "${outputs[@]}"; do
		value=${value,,}
		[[ $value == pfile ]] && value=pgen
		[[ $value == pgen || $value == bed ]] || die "invalid --out-files value: $value"
		requested[$value]=1
	done

	if [[ -n $chrs_csv ]]; then
		normalize_chromosomes "$chrs_csv"
	elif [[ -n ${requested[pgen]:-} ]]; then
		detect_input_format
		discover_input_chromosomes
	else
		discover_pfile_chromosomes
	fi
	if is_true "$split_sex"; then
		local found_x=0
		for value in "${selected_chrs[@]}"; do [[ $value == X ]] && found_x=1; done
		((found_x)) || die "--split-sex TRUE requires chrX to be selected"
		if [[ -z ${requested[pgen]:-} ]] &&
			{ ! pfile_exists "$pfile_dir/chrX.male" ||
				! pfile_exists "$pfile_dir/chrX.female"; }; then
			echo "sex-split chrX pfiles are missing; adding pgen generation before requested bed output"
			requested[pgen]=1
		fi
	fi
	sex_update=""
	if [[ -n ${requested[pgen]:-} ]]; then
		detect_input_format
		if [[ $input_format == vcf ]]; then
			sample_file=${sample_file:-$(dirname -- "$input_dir")/samples.txt}
			[[ -s $sample_file ]] || die "missing sample metadata: $sample_file"
			IFS=$' \t' read -r -a header < <(head -n 1 "$sample_file" | tr -d '\r')
			expected=(sample fatherID motherID sex pop super_pop)
			[[ "${header[*]}" == "${expected[*]}" ]] ||
				die "unexpected header in $sample_file; expected: ${expected[*]}"
			sex_update=$pfile_dir/ALL.update-sex.txt
			write_sex_update "$sex_update"
		fi
	fi
	echo "make-pfile input-format=$input_format input-prefix=$input_prefix dir-input=$input_dir dir-pfile=$pfile_dir chr=$(
		IFS=,
		echo "${selected_chrs[*]}"
	)"
	if [[ -n ${requested[pgen]:-} ]]; then
		for value in "${selected_chrs[@]}"; do make_pgen_one "$value"; done
	fi
	if [[ -n ${requested[bed]:-} ]]; then
		for value in "${selected_chrs[@]}"; do make_bed_one "$value"; done
	fi
	echo "DONE: make-pfile under $pfile_dir"
}

run_lift_bed_module() {
	local input rel output output_dir unmapped_output unmapped_dir
	local tmp_root chain_for_liftover input_bed mapped_bed unmapped_bed
	local output_tmp unmapped_tmp chrom_style_file chrom_style
	local -a bed_files=()
	local command_name

	for command_name in awk bgzip find gzip grep liftOver mktemp sort; do
		command -v "$command_name" >/dev/null 2>&1 || die "missing command: $command_name"
	done
	[[ $bgzip_threads =~ ^[1-9][0-9]*$ ]] ||
		die "--bgzip-threads must be a positive integer"
	[[ -d $dir_in ]] || die "input directory does not exist: $dir_in"
	[[ -s $chain_file ]] || die "chain file does not exist or is empty: $chain_file"

	dir_in=${dir_in%/}
	dir_out=${dir_out%/}
	[[ $dir_out != "$dir_in" && $dir_out != "$dir_in"/* ]] ||
		die "--dir-out must not be --dir-in or a directory below it"

	mapfile -d '' bed_files < <(
		find "$dir_in" -type f -name '*.bed.gz' -print0 | sort -z
	)
	((${#bed_files[@]} > 0)) || die "no *.bed.gz files found below: $dir_in"

	mkdir -p "$dir_out"
	tmp_root=$(mktemp -d "$dir_out/.lift-bed.tmp.XXXXXX")
	trap '[[ -n ${tmp_root:-} && -d $tmp_root ]] && rm -rf -- "$tmp_root"' EXIT
	chain_for_liftover=$chain_file
	if [[ $chain_file == *.gz ]]; then
		chain_for_liftover=$tmp_root/map.over.chain
		gzip -cd -- "$chain_file" >"$chain_for_liftover"
		[[ -s $chain_for_liftover ]] || die "failed to decompress chain file: $chain_file"
	fi

	input_bed=$tmp_root/input.bed
	mapped_bed=$tmp_root/mapped.bed
	unmapped_bed=$tmp_root/unmapped.bed
	chrom_style_file=$tmp_root/chrom.style
	echo "lift-bed files=${#bed_files[@]} dir-in=$dir_in dir-out=$dir_out chain=$chain_file"

	for input in "${bed_files[@]}"; do
		rel=${input#"$dir_in"/}
		output=$dir_out/$rel
		output_dir=$(dirname -- "$output")
		unmapped_output=$dir_out/.unmapped/${rel%.bed.gz}.unmapped.bed.gz
		unmapped_dir=$(dirname -- "$unmapped_output")

		if [[ -s $output ]] && gzip -t -- "$output" 2>/dev/null && ! is_true "$replace"; then
			echo "SKIP lift-bed: $rel already exists (use --replace TRUE)"
			continue
		fi

		mkdir -p "$output_dir" "$unmapped_dir"
		: >"$chrom_style_file"
		# UCSC chains use chr1-style names, while these archaic BED masks use
		# 1-style names. Record the input convention, add chr for liftOver,
		# then restore the original convention in both result files.
		gzip -cd -- "$input" |
			awk -v style_file="$chrom_style_file" 'BEGIN {OFS="\t"}
                /^#/ || /^track([[:space:]]|$)/ || /^browser([[:space:]]|$)/ {print; next}
                NF == 0 {next}
                NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 > $3 {
                    print "invalid BED row at input line " NR > "/dev/stderr"
                    bad=1; next
                }
                {
                    if (!seen++) print ($1 ~ /^chr/ ? "chr" : "bare") > style_file
                    if ($1 !~ /^chr/) $1="chr" $1
                    print
                }
                END {exit bad}
            ' >"$input_bed"
		[[ -s $chrom_style_file ]] || die "no BED records found in: $input"
		chrom_style=$(<"$chrom_style_file")

		: >"$mapped_bed"
		: >"$unmapped_bed"
		echo "LIFT: $rel"
		liftOver "$input_bed" "$chain_for_liftover" "$mapped_bed" "$unmapped_bed"

		output_tmp=$output.tmp.$$
		if [[ $chrom_style == bare ]]; then
			awk 'BEGIN {OFS="\t"}
                /^#/ || /^track([[:space:]]|$)/ || /^browser([[:space:]]|$)/ {print; next}
                {$1=substr($1, 4); print}
            ' "$mapped_bed" | bgzip -@ "$bgzip_threads" -c >"$output_tmp"
		else
			bgzip -@ "$bgzip_threads" -c "$mapped_bed" >"$output_tmp"
		fi
		gzip -t -- "$output_tmp" || die "invalid compressed output for: $rel"
		mv -f -- "$output_tmp" "$output"

		if grep -qv '^#' "$unmapped_bed"; then
			unmapped_tmp=$unmapped_output.tmp.$$
			if [[ $chrom_style == bare ]]; then
				awk 'BEGIN {OFS="\t"}
                    /^#/ {print; next}
                    {$1=substr($1, 4); print}
                ' "$unmapped_bed" | bgzip -@ "$bgzip_threads" -c >"$unmapped_tmp"
			else
				bgzip -@ "$bgzip_threads" -c "$unmapped_bed" >"$unmapped_tmp"
			fi
			gzip -t -- "$unmapped_tmp" || die "invalid unmapped output for: $rel"
			mv -f -- "$unmapped_tmp" "$unmapped_output"
		else
			rm -f -- "$unmapped_output"
		fi
	done

	rm -rf -- "$tmp_root"
	tmp_root=""
	trap - EXIT
	echo "DONE: lift-bed under $dir_out"
}

[[ $# -gt 0 ]] || {
	usage >&2
	exit 2
}
case $1 in
	-h | --help)
		usage
		exit 0
		;;
	make-pfile | lift-bed)
		module=$1
		shift
		;;
	*) die "unknown module '$1'; expected make-pfile or lift-bed" ;;
esac

grch=""
input_name_prefix=chr
input_format=auto
out_files_csv=pgen
chrs_csv=""
root=""
sample_file=""
dir_input=""
dir_vcf=""
dir_pfile=""
dir_in=""
dir_out=""
chain_file=""
input_dir=""
pfile_dir=""
replace=FALSE
split_sex=FALSE
bgzip_threads=4

while (($#)); do
	case $1 in
		--grch)
			grch=${2:?missing value for --grch}
			shift 2
			;;
		--input-prefix | --vcf-prefix)
			input_name_prefix=${2:?missing value for --input-prefix}
			shift 2
			;;
		--input-format)
			input_format=${2:?missing value for --input-format}
			input_format=${input_format,,}
			shift 2
			;;
		--out-files)
			out_files_csv=${2:?missing value for --out-files}
			shift 2
			;;
		--chr)
			chrs_csv=${2:?missing value for --chr}
			shift 2
			;;
		--split-sex)
			split_sex=${2:?missing value for --split-sex}
			is_true "$split_sex" || true
			shift 2
			;;
		--chrX-male | --chrx-male)
			die "--chrX-male was replaced by --split-sex"
			;;
		--root)
			root=${2:?missing value for --root}
			shift 2
			;;
		--bgzip-threads)
			bgzip_threads=${2:?missing value for --bgzip-threads}
			shift 2
			;;
		--sample-file)
			sample_file=${2:?missing value for --sample-file}
			shift 2
			;;
		--dir-input)
			dir_input=${2:?missing value for --dir-input}
			shift 2
			;;
		--dir-vcf)
			dir_vcf=${2:?missing value for --dir-vcf}
			shift 2
			;;
		--dir-pfile | --pfile-dir)
			dir_pfile=${2:?missing value for --dir-pfile}
			shift 2
			;;
		--dir-in)
			dir_in=${2:?missing value for --dir-in}
			shift 2
			;;
		--dir-out)
			dir_out=${2:?missing value for --dir-out}
			shift 2
			;;
		--chain-file)
			chain_file=${2:?missing value for --chain-file}
			shift 2
			;;
		--replace)
			replace=${2:?missing value for --replace}
			is_true "$replace" || true
			shift 2
			;;
		-h | --help)
			usage
			exit 0
			;;
		*) die "unknown option: $1" ;;
	esac
done

if [[ $module == make-pfile ]]; then
	[[ -z $grch || $grch == 37 || $grch == 38 ]] ||
		die "make-pfile --grch must be 37 or 38 when specified"
	[[ -z $root ]] || die "make-pfile does not accept --root; use --dir-input and --dir-pfile"
	if [[ -n $dir_vcf ]]; then
		[[ -z $dir_input ]] || die "use only one of --dir-input and --dir-vcf"
		[[ $input_format == auto || $input_format == vcf ]] ||
			die "--dir-vcf conflicts with --input-format $input_format"
		dir_input=$dir_vcf
		input_format=vcf
	fi
	if [[ -z $dir_input && -z $dir_pfile ]]; then
		default_build=${grch:-37}
		dir_input=${DIR_GEN}/1kg/$default_build/vcf
		dir_pfile=${DIR_GEN}/1kg/$default_build/pfile
		input_format=vcf
	elif [[ -z $dir_input ]]; then
		dir_input=$(dirname -- "${dir_pfile%/}")/vcf
		[[ $input_format != auto ]] || input_format=vcf
	elif [[ -z $dir_pfile ]]; then
		dir_pfile=$(dirname -- "${dir_input%/}")/pfile
		if [[ $(basename -- "$(dirname -- "$(dirname -- "${dir_input%/}")")") == ygb ]]; then
			dir_pfile=$(dirname -- "${dir_input%/}")/imp
		fi
	fi
	input_dir=${dir_input%/}
	pfile_dir=${dir_pfile%/}
elif [[ $module == lift-bed ]]; then
	[[ -z $grch ]] || die "lift-bed does not accept --grch; the chain file determines direction"
	[[ -z $root && -z $dir_input && -z $dir_vcf && -z $dir_pfile && -z $sample_file ]] ||
		die "lift-bed accepts --dir-in, --dir-out, --chain-file, --bgzip-threads, and --replace"
	[[ -n $dir_in ]] || die "lift-bed requires --dir-in DIR"
	[[ -n $dir_out ]] || die "lift-bed requires --dir-out DIR"
	[[ -n $chain_file ]] || die "lift-bed requires --chain-file FILE"
	run_lift_bed_module
	exit 0
fi
phe_sh=${phe_sh:-$dirfunc/phenotype.sh}
plink2_bin=${PLINK2_BIN:-plink2}

command -v "$plink2_bin" >/dev/null 2>&1 || die "missing PLINK2 executable: $plink2_bin"
gen_conversion_plink2_tool=$plink2_bin
gen_conversion_plink2() {
    # PLINK otherwise sizes its workspace from the entire node. Stay within
    # this job's budget while preserving all conversion/filter arguments.
    local memory_mb
    memory_mb=$(python3 -c 'import os; print(max(1, int(float(os.environ.get("WORKFLOW_MEMORY_GIB", "32")) * 1024 * .8)))')
    "$gen_conversion_plink2_tool" --threads "${OMP_NUM_THREADS:-4}" --memory "$memory_mb" "$@"
}
plink2_bin=gen_conversion_plink2
command -v awk >/dev/null 2>&1 || die "missing command: awk"
[[ -s $phe_sh ]] || die "missing chrX definitions: $phe_sh"
# shellcheck source=/dev/null
source "$phe_sh"
mkdir -p "$pfile_dir"

if [[ $input_name_prefix == */* ]]; then input_prefix=$input_name_prefix; else input_prefix=$input_dir/$input_name_prefix; fi

selected_chrs=()


# 🚩 Module dispatch
case $module in
	make-pfile) run_make_module ;;
esac
)

case ${1:-} in
    make-pfile|lift-bed) run_conversion_module "$@"; exit $? ;;
esac

# Plan every module as a .cmd before local execution or LSF submission.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "$SCRIPT_DIR/../0f/phenotype_jobs.sh"
ygb_job_init "$@" || exit $?
set -- "${YGB_ARGS[@]}"
source "$shared_dir/cli_dirs.sh"
script_directory_args "$SCRIPT_DIR/gen.sh" "$@"
set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
worker_args=("$@")
set +H
set -eo pipefail

usage() {
	bash "${BASH_SOURCE[0]}" --help
}

init_gen_env() {
	dir0=${DIR0:-${DIR_BIG}}
	local phe_f="$shared_dir/phenotype.sh"

	if [[ -f "$phe_f" ]]; then
		PHE_F=$phe_f
	elif [[ -f "$dir0/scripts/0f/phenotype.sh" ]]; then
		PHE_F=$dir0/scripts/0f/phenotype.sh
	else
		echo "ERROR: cannot find phenotype.sh beside gen.sh or under $dir0/scripts/0f" >&2
		exit 1
	fi
	source "$PHE_F"
}

to_wsl_path() {
	local path=$1
	if [[ "$path" =~ ^[A-Za-z]:[\\/].* ]]; then
		command -v wslpath >/dev/null 2>&1 || {
			echo "ERROR: Windows path supplied but wslpath is unavailable: $path" >&2
			exit 1
		}
		wslpath -u "$path"
	else
		printf '%s\n' "$path"
	fi
}

pvar_stream() {
	cat -- "$1"
}

liftover_output_is_complete() {
	# The log is written last and acts as the completion marker.  Validate the
	# marker, dependencies, links, and PVAR rows before resuming past it.
	local input_pvar=$1 input_pgen=$2 input_psam=$3 output_pvar=$4 output_log=$5
	local output_dir output_name output_stem
	output_dir=$(dirname "$output_pvar")
	output_name=$(basename "$output_pvar")
	output_stem=${output_name%.pvar}

	[[ -s "$output_pvar" && -s "$output_log" ]] || return 1
	[[ "$output_pvar" -nt "$input_pvar" ]] || return 1
	[[ "$output_pvar" -nt "$chain_file" ]] || return 1
	[[ "$output_pvar" -nt "$ref_fasta" && "$output_pvar" -nt "$ref_fasta.fai" ]] || return 1
	[[ $ref_fasta != *.gz && $ref_fasta != *.bgz && $ref_fasta != *.bgzf ]] || [[ -s $ref_fasta.gzi && "$output_pvar" -nt "$ref_fasta.gzi" ]] || return 1
	[[ -L "$output_dir/$output_stem.pgen" && "$output_dir/$output_stem.pgen" -ef "$input_pgen" ]] || return 1
	[[ -L "$output_dir/$output_stem.psam" && "$output_dir/$output_stem.psam" -ef "$input_psam" ]] || return 1

	awk -F '\t' -v expected_input="$input_pvar" -v expected_output="$output_pvar" '
        $1=="input" {
            input_count++
            input_ok=($2==expected_input)
        }
        $1=="output" {
            output_count++
            output_ok=($2==expected_output)
        }
        $1=="total_variants" && $2 ~ /^[0-9]+$/ {
            total_count++
            total=$2+0
        }
        $1=="unique_forward_same_chr_length" && $2 ~ /^[0-9]+$/ {
            candidate_count++
            candidate=$2+0
        }
        $1=="verified_ref" && $2 ~ /^[0-9]+$/ {
            ok_count++
            ok=$2+0
        }
        $1=="failed_pos_minus_1" && $2 ~ /^[0-9]+$/ {
            fail_count++
            fail=$2+0
        }
        END {
            valid=(input_count==1 && input_ok && output_count==1 && output_ok &&
                   total_count==1 && candidate_count==1 && ok_count==1 && fail_count==1 &&
                   total>0 && candidate<=total && ok<=candidate && fail==total-ok)
            exit !valid
        }
    ' "$output_log" || return 1

	awk -v expected="$(awk -F '\t' '$1=="total_variants" {print $2}' "$output_log")" '
        /^##/ {next}
        /^#CHROM[\t ]/ {header=1; next}
        !/^#/ {if (NF<5) bad=1; n++}
        END {exit !(header && !bad && n==expected)}
    ' "$output_pvar"
}

run_liftover() {
	if [[ -n "${bim_file:-}" ]]; then
		python3 "$SCRIPT_DIR/f/gen.py" lift-bim "$(to_wsl_path "$bim_file")" \
			"$(to_wsl_path "${dir_out:?ERROR: --dir-out required}")" \
			"$(to_wsl_path "${chain_file:?ERROR: --chain-file required}")" \
			"$(to_wsl_path "${ref_fasta:-${DIR_GEN}/1hgp/38/GRCH38.fasta.gz}")"
		return
	fi
	# Keep every PVAR row in its original position.  Only POS changes; rows that
	# cannot be mapped and verified against the target reference receive POS=-1.
	set -o pipefail
	: "${dir_in:?ERROR: liftOver requires --dir-in}"
	: "${dir_out:?ERROR: liftOver requires --dir-out}"
	: "${chain_file:?ERROR: liftOver requires --chain-file}"

	dir_in=$(to_wsl_path "$dir_in")
	dir_out=$(to_wsl_path "$dir_out")
	chain_file=$(to_wsl_path "$chain_file")
	ref_fasta=$(to_wsl_path "${ref_fasta:-${DIR_GEN}/1hgp/38/GRCH38.fasta.gz}")

	[[ -d "$dir_in" ]] || {
		echo "ERROR: input directory not found: $dir_in" >&2
		exit 1
	}
	[[ -f "$chain_file" ]] || {
		echo "ERROR: chain file not found: $chain_file" >&2
		exit 1
	}
	[[ -f "$ref_fasta" ]] || {
		echo "ERROR: GRCh38 FASTA not found: $ref_fasta" >&2
		exit 1
	}
	[[ -f "$ref_fasta.fai" ]] || {
		echo "ERROR: FASTA index not found: $ref_fasta.fai" >&2
		exit 1
	}
	[[ $ref_fasta != *.gz && $ref_fasta != *.bgz && $ref_fasta != *.bgzf ]] || [[ -s $ref_fasta.gzi ]] || {
		echo "ERROR: BGZF index not found: $ref_fasta.gzi" >&2
		exit 1
	}

	local lift_over=${LIFTOVER_BIN:-}
	if [[ -z "$lift_over" ]]; then
		if command -v liftOver >/dev/null 2>&1; then
			lift_over=$(command -v liftOver)
		elif [[ -x ${DIR_SOFTWARE}/bin/liftOver ]]; then
			lift_over=${DIR_SOFTWARE}/bin/liftOver
		else
			echo "ERROR: liftOver not found (set LIFTOVER_BIN if needed)" >&2
			exit 1
		fi
	fi
	if ! command -v bedtools >/dev/null 2>&1 && [[ -x $DIR_WORK/.conda/envs/gu/bin/bedtools ]]; then
		export PATH="$DIR_WORK/.conda/envs/gu/bin:$PATH"
	fi
	command -v bedtools >/dev/null 2>&1 || {
		echo "ERROR: bedtools not found" >&2
		exit 1
	}
	command -v realpath >/dev/null 2>&1 || {
		echo "ERROR: realpath not found" >&2
		exit 1
	}

	mkdir -p "$dir_out"
	dir_in=$(cd "$dir_in" && pwd -P)
	dir_out=$(cd "$dir_out" && pwd -P)
	chain_file=$(cd "$(dirname "$chain_file")" && pwd -P)/$(basename "$chain_file")
	ref_fasta=$(cd "$(dirname "$ref_fasta")" && pwd -P)/$(basename "$ref_fasta")
	[[ "$dir_in" != "$dir_out" ]] || {
		echo "ERROR: --dir-in and --dir-out must differ" >&2
		exit 1
	}

	local -a pvars
	if [[ -n ${snp_data:-} ]]; then
		[[ $snp_data =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo 'ERROR: --data must be a path-free pfile basename' >&2; exit 2; }
		[[ -z ${lift_chr:-} ]] || { echo 'ERROR: choose --data or --chr, not both' >&2; exit 2; }
		pvars=("$dir_in/$snp_data.pvar")
	elif [[ -n ${lift_chr:-} ]]; then
		pvars=("$dir_in/chr$lift_chr.pvar")
	else
		mapfile -t pvars < <(find "$dir_in" -maxdepth 1 -type f -name 'chr*.pvar' -printf '%p\n' | sort -V)
	fi
	((${#pvars[@]} > 0)) || { echo "ERROR: no chr*.pvar files found; use --data NAME for a named pfile" >&2; exit 1; }
	for pvar in "${pvars[@]}"; do
		[[ -f $pvar ]] || { echo "ERROR: missing $pvar" >&2; exit 1; }
	done

	local runtime_root
	runtime_root=$(to_wsl_path "${dir_runtime:-$dir_out}")
	mkdir -p "$runtime_root/logs" "$runtime_root/tmp"
	runtime_root=$(cd "$runtime_root" && pwd -P)
	export TMPDIR="$runtime_root/tmp"
	local lift_work_dir=
	cleanup_liftover_work() {
		if [[ -n "$lift_work_dir" && -d "$lift_work_dir" ]]; then
			find "$lift_work_dir" -maxdepth 1 -type f -delete
			rmdir "$lift_work_dir" 2>/dev/null || true
		fi
	}
	trap cleanup_liftover_work EXIT
	trap 'exit 130' INT TERM

	local pvar name stem pgen psam out_pvar out_log
	local source_bed mapped_bed unmapped_bed sorted_bed candidate_bed verified_map out_tmp out_log_tmp
	local n_total n_candidate n_ok n_fail rel target
	for pvar in "${pvars[@]}"; do
		name=$(basename "$pvar")
		stem=${name%.pvar}
		pgen="$dir_in/$stem.pgen"
		psam="$dir_in/$stem.psam"
		[[ -e "$pgen" ]] || {
			echo "ERROR: missing companion file: $pgen" >&2
			exit 1
		}
		[[ -e "$psam" ]] || {
			echo "ERROR: missing companion file: $psam" >&2
			exit 1
		}

		out_pvar="$dir_out/$name"
		out_log="$runtime_root/logs/$name.liftOver.log"
		# Completed legacy outputs remain reusable without rewriting their receipt.
		if [[ ! -e $out_log && -e $dir_out/$name.liftOver.log ]]; then
			out_log="$dir_out/$name.liftOver.log"
		fi
		if liftover_output_is_complete "$pvar" "$pgen" "$psam" "$out_pvar" "$out_log"; then
			echo "SKIP $name (complete)"
			continue
		fi
		if [[ -e "$out_pvar" || -e "$out_log" ]]; then
			echo "REBUILD $name (existing result is incomplete or stale)"
		fi

		lift_work_dir="$runtime_root/tmp/.liftOver.${stem}.$$"
		mkdir "$lift_work_dir"
		source_bed="$lift_work_dir/source.bed"
		mapped_bed="$lift_work_dir/mapped.bed"
		unmapped_bed="$lift_work_dir/unmapped.bed"
		sorted_bed="$lift_work_dir/mapped.sorted.bed"
		candidate_bed="$lift_work_dir/candidate.bed"
		verified_map="$lift_work_dir/verified.tsv"
		out_tmp="$lift_work_dir/$name"
		out_log_tmp="$lift_work_dir/$name.liftOver.log"

		echo "LIFT $name"
		# The BED name carries row number, expected chromosome, and original REF.
		# It survives liftOver and avoids loading millions of variants into memory.
		n_total=$(pvar_stream "$pvar" | awk -v bed="$source_bed" '
            BEGIN { OFS="\t" }
            /^#/ { next }
            {
                row++
                chr=$1
                sub(/^chr/, "", chr)
                if (chr==23) chr="X"
                else if (chr==24) chr="Y"
                else if (chr==25) chr="XY"
                else if (chr==26 || chr=="M") chr="MT"
                target_chr="chr" chr
                ref=toupper($4)
                if ($2 ~ /^[0-9]+$/ && $2>0 && ref ~ /^[ACGTN]+$/) {
                    start=$2-1
                    print target_chr,start,start+length(ref),row "|" target_chr "|" ref,0,"+" > bed
                }
            }
            END { print row+0 }
        ')
		((n_total > 0)) || {
			echo "ERROR: no variant rows in $pvar" >&2
			exit 1
		}
		: >"$mapped_bed"
		: >"$unmapped_bed"
		: >"$sorted_bed"
		if [[ -s "$source_bed" ]]; then
			"$lift_over" "$source_bed" "$chain_file" "$mapped_bed" "$unmapped_bed" >/dev/null
			LC_ALL=C sort -k4,4n -k2,2n "$mapped_bed" >"$sorted_bed"
		fi

		# Accept exactly one forward-strand mapping to the expected chromosome,
		# with the same interval length.  REF is verified in the next pass.
		: >"$candidate_bed"
		awk -v out="$candidate_bed" '
            BEGIN { OFS="\t" }
            function load_record(    a,n) {
                n=split($4,a,"|")
                id=a[1]+0
                expected_chr=a[2]
                ref=a[3]
                tchr=$1; tstart=$2; tend=$3
                ok=(n==3 && $6=="+" && $1==expected_chr && $2>=0 && $3>$2 && $3-$2==length(ref))
            }
            function emit_group() {
                if (count==1 && ok) {
                    print tchr,tstart,tend,id "|" (tstart+1) "|" ref > out
                }
            }
            NR==1 { load_record(); count=1; next }
            { split($4,current,"|") }
            current[1]+0==id { count++; next }
            { emit_group(); load_record(); count=1 }
            END { if (NR) emit_group() }
        ' "$sorted_bed"

		: >"$verified_map"
		if [[ -s "$candidate_bed" ]]; then
			bedtools getfasta -fi "$ref_fasta" -bed "$candidate_bed" -nameOnly -tab 2>/dev/null |
				awk -F '\t' 'BEGIN{OFS="\t"} {
                    n=split($1,a,"|")
                    if (n==3 && toupper($2)==toupper(a[3])) print a[1],a[2]
                }' >"$verified_map"
		fi
		n_candidate=$(wc -l <"$candidate_bed")
		n_ok=$(wc -l <"$verified_map")
		n_fail=$((n_total - n_ok))

		# verified_map is sorted by row number.  Consume it as a stream while
		# rereading the PVAR, so memory use stays bounded for imputed chromosomes.
		pvar_stream "$pvar" | awk -v map="$verified_map" '
            BEGIN {
                OFS="\t"
                if ((getline line < map)>0) { split(line,a,"\t"); next_id=a[1]+0; next_pos=a[2]+0 }
                else next_id=0
            }
            /^#/ { print; next }
            {
                row++
                if (row==next_id) {
                    $2=next_pos
                    if ((getline line < map)>0) { split(line,a,"\t"); next_id=a[1]+0; next_pos=a[2]+0 }
                    else next_id=0
                } else {
                    $2=-1
                }
                print
            }
            END {
                close(map)
                if (row==0 || next_id!=0) exit 42
            }
        ' >"$out_tmp"
		[[ $(cat -- "$out_tmp" | awk '!/^#/{n++} END{print n+0}') -eq "$n_total" ]] || {
			echo "ERROR: output row-count check failed for $name" >&2
			exit 1
		}
		mv -f "$out_tmp" "$out_pvar"

		local companion_name
		for target in "$pgen" "$psam"; do
			companion_name=$(basename "$target")
			if [[ -e "$dir_out/$companion_name" && ! -L "$dir_out/$companion_name" ]]; then
				echo "ERROR: refusing to replace non-symlink: $dir_out/$companion_name" >&2
				exit 1
			fi
			rel=$(realpath --relative-to="$dir_out" "$target")
			ln -sfn -- "$rel" "$dir_out/$companion_name"
		done

		{
			printf 'input\t%s\n' "$pvar"
			printf 'output\t%s\n' "$out_pvar"
			printf 'total_variants\t%s\n' "$n_total"
			printf 'unique_forward_same_chr_length\t%s\n' "$n_candidate"
			printf 'verified_ref\t%s\n' "$n_ok"
			printf 'failed_pos_minus_1\t%s\n' "$n_fail"
			printf 'pgen_link\t%s\n' "$(readlink "$dir_out/$stem.pgen")"
			printf 'psam_link\t%s\n' "$(readlink "$dir_out/$stem.psam")"
		} >"$out_log_tmp"
		mv -f "$out_log_tmp" "$out_log"
		echo "OK $name: total=$n_total lifted=$n_ok POS=-1=$n_fail"

		cleanup_liftover_work
		lift_work_dir=
	done
	trap - EXIT INT TERM
	echo "OK: lifted PVARs and relative pfile links written to: $dir_out"
}

publish_gen_result() {
    local source=$1 target=$2 action=${3:-move}
    mkdir -p "$(dirname "$target")"
    if [[ $(realpath -m "$source") != "$(realpath -m "$target")" ]]; then
        # Stage beside the target so publication is atomic even across filesystems.
        local staged
        staged=$(mktemp "${target}.tmp.XXXXXX")
        if ! cp -- "$source" "$staged"; then rm -f -- "$staged"; return 1; fi
        mv -f -- "$staged" "$target"
        [[ $action == keep ]] || rm -f -- "$source"
    fi
    echo "OK: final data written to: $target"
}

run_pgs() {
	init_gen_env
	: "${pgs_input:?ERROR: pgs requires --input FILE}"
	: "${label:?ERROR: pgs requires --label LABEL}"
	pgs_input=$(to_wsl_path "$pgs_input")
	pgs_pfile_dir=$(to_wsl_path "${pgs_pfile_dir:-${DIR_GEN}/ygb/37/imp}")
	dir_out=$(to_wsl_path "${dir_out:-$ANALYSIS_ROOT/ygb/gen/pgs/$label}")
	pgs_threads=${pgs_threads:-$YGB_THREADS}
	[[ -s "$pgs_input" ]] || {
		echo "ERROR: PGS score input is missing/empty: $pgs_input" >&2
		exit 1
	}
	[[ "$pgs_threads" =~ ^[1-9][0-9]*$ ]] || {
		echo "ERROR: --pgs-threads must be a positive integer" >&2
		exit 1
	}
	declare -F pgs_plink_calc >/dev/null 2>&1 || {
		echo "ERROR: pgs_plink_calc is missing from $PHE_F" >&2
		exit 1
	}

	mkdir -p "$dir_out/jobs" "$dir_out/tmp" "$dir_out/qc" "$dir_out/logs"
	local cmd_file="$dir_out/jobs/$label.pgs.worker.sh"
	{
		printf '%s\n' '#!/bin/bash -l' 'set +H' 'set -euo pipefail'
		printf 'score_file=%q\n' "$pgs_input"
		printf 'pfile_dir=%q\n' "$pgs_pfile_dir"
		printf 'out_dir=%q\n' "$dir_out"
		printf 'label=%q\n' "$label"
		printf 'threads=%q\n' "$pgs_threads"
		printf 'phe_f=%q\n' "$PHE_F"
		cat <<'EOF2'

source "$phe_f"
declare -F pgs_plink_calc >/dev/null 2>&1 || { echo "ERROR: pgs_plink_calc is missing from $phe_f" >&2; exit 1; }
collect_pgs_logs() {
    local status=$? log
    trap - EXIT
    for log in "$out_dir/tmp/$label".chr*.log; do
        [[ ! -f $log ]] || mv -f -- "$log" "$out_dir/logs/"
    done
    exit "$status"
}
trap collect_pgs_logs EXIT
pgs_plink_calc \
    --input "$score_file" \
    --pfile-dir "$pfile_dir" \
    --output "$out_dir/$label.pgs.txt" \
    --label "$label" \
    --work-dir "$out_dir/tmp" \
    --status-file "$out_dir/qc/$label.pgs.status.tsv" \
    --threads "$threads"
echo "OK: PGS written to: $out_dir/$label.pgs.txt"
EOF2
	} >"$cmd_file"
	chmod +x "$cmd_file"

	bash "$cmd_file"
	publish_gen_result "$dir_out/$label.pgs.txt" "${YGB_PGS_OUTPUT_DIR:-$dir_gen/pgs}/$label.pgs.txt"
}

run_snp() {


	# Extract SNP data
	: "${dir_gen:?ERROR: snp requires --dir-gen}"
	: "${snp_data:?ERROR: snp requires --data imp}"
	: "${snp_file:?ERROR: snp requires --snp-file}"
	: "${snp_file_type:?ERROR: snp requires --snp-file-type bed|snp}"
	: "${dir_out:?ERROR: snp requires --dir-out}"
	: "${label:?ERROR: snp requires --label}"
	[[ "$snp_data" == imp ]] || {
		echo "ERROR: SNP calculation requires --data imp" >&2
		exit 1
	}
	[[ "$snp_file_type" == bed || "$snp_file_type" == snp ]] || {
		echo "ERROR: --snp-file-type must be bed or snp" >&2
		exit 1
	}
	if [[ "$snp_file_type" == bed ]]; then
		flank=${flank:-10000}
		[[ "$flank" =~ ^[0-9]+$ ]] || {
			echo "ERROR: --flank must be a nonnegative integer" >&2
			exit 1
		}
	elif [[ "${flank_set:-N}" == Y ]]; then
		echo "ERROR: --flank is only valid with --snp-file-type bed" >&2
		exit 1
	fi
	[[ -d "$dir_gen/$snp_data" ]] || {
		echo "ERROR: missing data directory: $dir_gen/$snp_data" >&2
		exit 1
	}
	[[ -f "$snp_file" ]] || {
		echo "ERROR: missing SNP file: $snp_file" >&2
		exit 1
	}

	dir_gen=$(cd "$dir_gen" && pwd)
	snp_file=$(cd "$(dirname "$snp_file")" && pwd)/$(basename "$snp_file")
	mkdir -p "$dir_out"
	dir_out=$(cd "$dir_out" && pwd)
	local alias_args=""
	if [[ "$snp_file_type" == snp ]]; then
		printf -v alias_args " --alias-file %q" "$snp_file"
	fi
	local export_allele_args=""
	if [[ -n ${export_allele_file:-} ]]; then
		[[ -s "$export_allele_file" ]] || { echo "Missing export allele file: $export_allele_file" >&2; exit 1; }
		export_allele_file=$(realpath -- "$export_allele_file")
		printf -v export_allele_args " --export-allele %q" "$export_allele_file"
	fi
	work_dir="$dir_out/tmp/${label}.snp-work"
	mkdir -p "$work_dir"

	if [[ "$snp_file_type" == bed ]]; then
		range_file="$work_dir/$label.flank${flank}.bed"
		awk -v flank="$flank" 'BEGIN{OFS="\t"} /^[[:space:]]*#/ || NF<3 {next} {s=$2-flank; if(s<0)s=0; print $1,s,$3+flank,(NF>=4?$4:".")}' "$snp_file" >"$range_file"
		[[ -s "$range_file" ]] || {
			echo "ERROR: no BED intervals found in $snp_file" >&2
			exit 1
		}
		extract_args="--extract bed0 \"$range_file\""
		mapfile -t chromosomes < <(awk '{c=$1; sub(/^chr/,"",c); print c}' "$range_file" | sort -uV)
	else
		extract_file="$work_dir/$label.extract.snp"
		awk 'BEGIN{FS="[[:space:]]+"}
             /^[[:space:]]*#/ || /^[[:space:]]*$/ {next}
             !seen && toupper($1) ~ /^(SNP|RSID|VARIANT|ID)$/ {seen=1; next}
             {seen=1; print $1}' "$snp_file" >"$extract_file"
		[[ -s "$extract_file" ]] || {
			echo "ERROR: no variant IDs found in $snp_file" >&2
			exit 1
		}
		if [[ $(sort "$extract_file" | uniq -d | wc -l) -gt 0 ]]; then
			echo "ERROR: duplicate variant IDs found in first column of $snp_file" >&2
			sort "$extract_file" | uniq -d | head -20 >&2
			exit 1
		fi
		extract_args="--extract \"$extract_file\""
			mapfile -t chromosomes < <(find -L "$dir_gen/imp" -maxdepth 1 -type f -name 'chr*.pgen' -printf '%f\n' | sed -nE 's/^chr([0-9]+|X|Y|XY|MT)\.pgen$/\1/p' | sort -uV)
	fi
	((${#chromosomes[@]} > 0)) || {
		echo "ERROR: no chromosome data files found" >&2
		exit 1
	}

	for chr in "${chromosomes[@]}"; do
			[[ -f "$dir_gen/$snp_data/chr$chr.pgen" && -f "$dir_gen/$snp_data/chr$chr.pvar" && -f "$dir_gen/$snp_data/chr$chr.psam" ]] || {
				echo "ERROR: incomplete $snp_data pfiles for chr$chr" >&2
				exit 1
			}
	done

	cmd_file="$work_dir/$label.worker.sh"
	{
		cat <<EOF2
#!/bin/bash -l
set +H
set -e

scratch=\$(mktemp -d "$work_dir/run.XXXXXX")
trap 'rm -rf -- "\$scratch"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP
cd "\$scratch"
rm -f merge-list.txt
for chr in ${chromosomes[*]}; do
        input_args="--pfile \"$dir_gen/$snp_data/chr\$chr\""
    set +e
    eval plink2 --threads "$YGB_THREADS" --memory "$YGB_PLINK_MEMORY_MB" \$input_args $extract_args --make-pgen --out "chr\$chr"
    plink_status=\$?
    set -e
    if (( plink_status != 0 )); then
        if grep -Eq 'No variants remaining|No variants loaded|0 variants remaining' "chr\$chr.log"; then
            echo "WARN: no variants extracted from chr\$chr" >&2
            rm -f "chr\$chr.pgen" "chr\$chr.pvar" "chr\$chr.psam"
            continue
        fi
        echo "ERROR: PLINK extraction failed for chr\$chr (exit \$plink_status)" >&2
        exit "\$plink_status"
    fi
    if awk '!/^#/{found=1} END{exit !found}' "chr\$chr.pvar"; then
        echo "chr\$chr" >> merge-list.txt
    else
        rm -f "chr\$chr.pgen" "chr\$chr.pvar" "chr\$chr.psam"
        echo "WARN: no variants extracted from chr\$chr" >&2
    fi
done
[[ -s merge-list.txt ]] || { echo "ERROR: no variants were extracted" >&2; exit 1; }
if [[ \$(wc -l < merge-list.txt) -eq 1 ]]; then
    merged=\$(head -1 merge-list.txt)
else
    merged=merged
    plink2 --threads "$YGB_THREADS" --memory "$YGB_PLINK_MEMORY_MB" --pmerge-list merge-list.txt --make-pgen --out "\$merged"
    while read -r part; do rm -f "\$part.pgen" "\$part.pvar" "\$part.psam"; done < merge-list.txt
fi
python3 "$SCRIPT_DIR/f/gen.py" export --output "$dir_out/$label.raw.gz" $alias_args -- plink2 --threads "$YGB_THREADS" --memory "$YGB_PLINK_MEMORY_MB" --pfile "\$merged" --maj-ref --export A $export_allele_args
EOF2
	} >"$cmd_file"
	chmod +x "$cmd_file"

	bash "$cmd_file"
	publish_gen_result "$dir_out/$label.raw.gz" "${YGB_SNP_OUTPUT_DIR:-$dir_gen/imp}/$label.raw.gz"
}

run_apoe() {


	# Extract APOE haplotype
	# e1: C-T; e2: T-T; e3: T-C; e4: C-C
	dir0=${DIR_WORK}
	dirout=$(realpath -m -- "${dir_out:-$ANALYSIS_ROOT/ygb/gen/apoe}")
	dir_gen=${dir_gen:-$DIR_GEN/ygb/37}
	mkdir -p "$dirout/tmp/apoe" "$dirout/jobs" "$dirout/logs"

	cat >"$dirout/jobs/apoe.worker.sh" <<EOF2
#!/bin/bash -l
set +H
set -e

workdir=\$(mktemp -d "$dirout/tmp/apoe/run.XXXXXX")
trap 'rm -rf -- "\$workdir"' EXIT
cd "\$workdir"
printf '%s\n' rs429358 rs7412 > apoe.snps
plink2 --threads "$YGB_THREADS" --memory "$YGB_PLINK_MEMORY_MB" --pfile "$dir_gen/hap/chr19" --extract apoe.snps --export vcf id-paste=iid bgz --out apoe
python3 "$SCRIPT_DIR/f/gen.py" apoe apoe.vcf.gz apoe.hap
awk 'NR>1{print \$2}' apoe.hap | sort | uniq -c > apoe.hap.count.txt
cat apoe.hap.count.txt
mv -f apoe.hap apoe.hap.count.txt "$dirout/"
mv -f apoe.log "$dirout/logs/"
EOF2
	chmod +x "$dirout/jobs/apoe.worker.sh"

	bash "$dirout/jobs/apoe.worker.sh"
	publish_gen_result "$dirout/apoe.hap" "$dir_gen/hap/apoe.hap" keep
}

pos_args=()
flank=
flank_set=N
while [[ $# -gt 0 ]]; do
	case "$1" in
		-h | --help)
			usage
			exit 0
			;;
		--dir-gen)
			dir_gen=${2:?ERROR: --dir-gen needs a value}
			shift 2
			;;
		--data)
			snp_data=${2:?ERROR: --data needs a value}
			shift 2
			;;
		--snp-file)
			snp_file=${2:?ERROR: --snp-file needs a value}
			shift 2
			;;
		--export-allele)
			export_allele_file=${2:?ERROR: --export-allele needs a file}
			shift 2
			;;
		--snp-file-type)
			snp_file_type=${2:?ERROR: --snp-file-type needs a value}
			shift 2
			;;
		--flank)
			flank=${2:?ERROR: --flank needs a value}
			flank_set=Y
			shift 2
			;;
		--chr)
			lift_chr=${2:?ERROR: --chr needs a value}
			lift_chr=${lift_chr#chr}
			[[ $lift_chr =~ ^([1-9]|1[0-9]|2[0-6]|X|X\.male|Y|XY|MT|M)$ ]] || { echo "ERROR: invalid --chr" >&2; exit 2; }
			shift 2
			;;
		--bim-file)
			bim_file=${2:?ERROR: --bim-file needs a value}
			shift 2
			;;
		--dir-runtime)
			dir_runtime=${2:?ERROR: --dir-runtime needs a value}
			shift 2
			;;
		--dir-in)
			dir_in=${2:?ERROR: --dir-in needs a value}
			shift 2
			;;
		--input | --pgs-input)
			pgs_input=${2:?ERROR: $1 needs a value}
			shift 2
			;;
		--pgs-pfile-dir)
			pgs_pfile_dir=${2:?ERROR: --pgs-pfile-dir needs a value}
			shift 2
			;;
		--pgs-threads)
			pgs_threads=${2:?ERROR: --pgs-threads needs a value}
			shift 2
			;;
		--dir-out)
			dir_out=${2:?ERROR: --dir-out needs a value}
			shift 2
			;;
		--chain-file)
			chain_file=${2:?ERROR: --chain-file needs a value}
			shift 2
			;;
		--ref-fasta)
			ref_fasta=${2:?ERROR: --ref-fasta needs a value}
			shift 2
			;;
		--label)
			label=${2:?ERROR: --label needs a value}
			shift 2
			;;
		--*)
			echo "Unknown option: $1" >&2
			usage >&2
			exit 1
			;;
		*)
			pos_args+=("$1")
			shift
			;;
	esac
done
if ((${#pos_args[@]} > 1)); then
	echo "ERROR: too many positional arguments: ${pos_args[*]}" >&2
	usage >&2
	exit 1
fi

task=${pos_args[0]:-all}
task=$(echo "$task" | tr '[:upper:]' '[:lower:]')

case "$task" in pgs|liftover|snp|apoe|all) ;; *) echo "ERROR: unknown task: $task" >&2; exit 2 ;; esac
if [[ -n ${label:-} && ! $label =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo 'ERROR: --label must be a path-free name' >&2; exit 2
fi
if [[ -n ${pgs_threads:-} ]]; then
    [[ $pgs_threads =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: --pgs-threads must be positive' >&2; exit 2; }
    (( YGB_THREADS >= pgs_threads )) || YGB_THREADS=$pgs_threads
fi
dir_gen=${dir_gen:-$DIR_GEN/ygb/37}
case "$task" in
    snp)
        snp_data=${snp_data:-imp}
        label=${label:-snp}
        dir_out=${dir_out:-$ANALYSIS_ROOT/ygb/gen/snp}
        ;;
    pgs)
        dir_out=${dir_out:-$ANALYSIS_ROOT/ygb/gen/pgs/${label:?ERROR: pgs requires --label}}
        ;;
esac
# Default SNP list: first column is the variant ID, second is its output alias.
case "$task" in
    snp|all)
        snp_file=${snp_file:-$DIR_YGB/common/snp.lst}
        snp_file_type=${snp_file_type:-snp}
        ;;
esac
# Validate required arguments before dispatch; data scans and tools run in the worker.
case "$task" in
    pgs|all) : "${pgs_input:?ERROR: pgs requires --input}" "${label:?ERROR: pgs requires --label}" ;;
esac
case "$task" in
    snp|all) : "${snp_data:?ERROR: snp requires --data}" "${snp_file:?ERROR: snp requires --snp-file}"         "${snp_file_type:?ERROR: snp requires --snp-file-type}" "${dir_out:?ERROR: snp requires --dir-out}" "${label:?ERROR: snp requires --label}" ;;
    liftover) : "${dir_out:?ERROR: liftOver requires --dir-out}" "${chain_file:?ERROR: liftOver requires --chain-file}"
        [[ -n ${bim_file:-} ]] || : "${dir_in:?ERROR: liftOver requires --dir-in}" ;;
esac
dir_gen=${dir_gen:-$DIR_GEN/ygb/37}
# SNP and PGS calculations must use imputed pgen inputs, never typed inputs.
case "$task" in
    snp|all)
        [[ "$snp_data" == imp ]] || { echo 'ERROR: SNP calculation requires --data imp' >&2; exit 2; }
        ;;
esac
case "$task" in
    pgs|all)
        pgs_pfile_dir=$(to_wsl_path "${pgs_pfile_dir:-$dir_gen/imp}")
        pgs_resolved=$(realpath -m -- "$pgs_pfile_dir")
        [[ $(basename -- "$pgs_resolved") == imp && "$pgs_resolved" != */typ/* ]] || {
            echo 'ERROR: PGS calculation requires an imp/ pgen directory' >&2; exit 2;
        }
        ;;
esac
job_root=${dir_out:-$ANALYSIS_ROOT/ygb/gen/$task}
if [[ $task == pgs && -z ${dir_out:-} ]]; then job_root=$job_root/$label; fi
ygb_dispatch "gen.$task${label:+.$label}" "$job_root" "$SCRIPT_DIR/gen.sh" "${worker_args[@]}"
ygb_worker_runtime
if [[ $task == pgs || $task == all ]]; then
    PGS_PLINK_MEMORY_MB=${PGS_PLINK_MEMORY_MB:-8192}
    [[ $PGS_PLINK_MEMORY_MB =~ ^[1-9][0-9]*$ ]] && (( PGS_PLINK_MEMORY_MB <= YGB_MEMORY_GIB * 1024 )) || {
        echo 'ERROR: PGS_PLINK_MEMORY_MB exceeds the job RAM budget or is invalid' >&2; exit 2;
    }
    export PGS_PLINK_MEMORY_MB
fi

case "$task" in
	pgs) run_pgs ;;
	liftover) run_liftover ;;
	snp) run_snp ;;
	apoe) run_apoe ;;
	all)
		run_pgs
		run_snp
		run_apoe
		;;
	*)
		echo "ERROR: unknown task: $task" >&2
		usage >&2
		exit 1
		;;
esac
