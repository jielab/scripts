#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# CPU planning flags for the analysis command below:
#   --submit hpc --queue 33 --cores 40 --memory-gb 160 --hours 24 --plan-only
# New HPC plans auto-round slots. Keep tool threads within the allocated total.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
source "$_paths_directory/0f/0os.sh"
unset _paths_directory


# Resolve paths and dispatch before any scientific work or input discovery.
GWAS_ENTRY_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$GWAS_ENTRY_ROOT/f/runtime.sh"
gwas_runtime_init "gwas.sh" "$@" || exit $?
set -- "${GWAS_ARGS[@]}"
source "$GWAS_ENTRY_ROOT/../0f/console.sh"
set -euo pipefail
set +H


# 🚩 Command-line help
# Independent GWAS modules. Statistical/data conversion helpers live in f/;
# parameter parsing, genotype selection and analysis commands are defined here.
usage() {
	cat <<HELP
Modules:
    sparseGRM  Prepare the pruned typed panel (add --sparse-grm FALSE for REGENIE only)
    acount     Prepare shared imputed allele counts and variant lists
    regenie    REGENIE Step 1 + Step 2; HPC submission is the default

Defaults from 0f/0os.sh:
    phenotype: $GWAS_PHENOTYPE_FILE
    continuous covariates: age,sex,PC1,PC2,PC3,PC4
    typed/imputed: $DIR_GEN/ukb/37/typ and $DIR_GEN/ukb/37/imp
    output: $GWAS_ANALYSIS_ROOT/self/common
    HPC queue: 33; per phenotype: 6 workers x 4 threads, 160 GiB RAM (40 reserved slots)

Main options:
    --grch 37|38                  Build; also sets REGENIE PAR boundaries
    --pheno FILE --pheno-name CSV Phenotype table (FID IID plus named columns)
    --type auto|qt|ordinal|bt|t2e  ordinal = linear trend on numeric ordered scores,
                                 via --qt --force-qt; NOT ordinal logistic regression
    --event-col NAME              Event indicator for t2e (0=censored, 1=event)
    --sex all|male|female         Genotype-sex restriction in BOTH REGENIE steps
    --covariates CSV              Continuous covariates; default listed above
    --categorical-covariates CSV  Dummy-coded categorical covariates; default none
    --covariate-file FILE         Shared covariate table; default phenotype table
    --continuous-covariate-file FILE  Explicit file: all non-ID columns by default
    --categorical-covariate-file FILE Explicit file: all non-ID columns by default
    --chr CSV                    REGENIE default 1-22,X; subsets require explicit --chr
    --threads N --jobs N          Threads/process and concurrent Step 2 chromosomes [4 / 6]
    --step2-threads N             Step 2 threads; defaults to --threads
    --pheno-jobs N               HPC: concurrent phenotype jobs [all requested]
                                 Names: regenie.[type].[pheno]
    --dir-output DIR              Use separate output roots for different models
    --replace TRUE               Recompute changed settings in an existing result folder
    --run-cmd FALSE               Write internal commands without running them
    --plan-only                  Write HPC submission plan; do not submit
    -h, --help                    Show copyable commands for the current path table

Examples (run from this copied scripts tree; default: submit to HPC):
    cd "$GWAS_ENTRY_ROOT"

    # Only if typ.prune.bed/bim/fam have not already been generated:
    bash gwas.sh sparseGRM --grch 37 --chr 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 --sparse-grm FALSE --wait TRUE

    # Three independent HPC jobs; covariates age,sex,PC1-PC4.
    bash gwas.sh regenie --grch 37 --pheno-name bmi,height,ldl --type qt --dir-output "$GWAS_ANALYSIS_ROOT/ukb_regenie"

    # Ordered bald scores: male-only linear trend across scores 1-4.
    bash gwas.sh regenie --grch 37 --pheno-name bald --type ordinal --sex male --dir-output "$GWAS_ANALYSIS_ROOT/ukb_regenie"

    # Binary contrasts: 0/1, logistic GWAS with approximate Firth correction.
    bash gwas.sh regenie --grch 37 --pheno-name bald12,bald13,bald14 --type bt --sex male --dir-output "$GWAS_ANALYSIS_ROOT/ukb_regenie"

    # Time to event: positive duration + cvd_cad.Yt2e (0=censored, 1=event).
    bash gwas.sh regenie --grch 37 --pheno-name cvd_cad.t2e --type t2e --event-col cvd_cad.Yt2e --dir-output "$GWAS_ANALYSIS_ROOT/ukb_regenie"

    # Add continuous tdi and categorical center; separate output for this model.
    bash gwas.sh regenie --grch 37 --pheno-name bmi,height,ldl --type qt \\
        --covariates age,sex,PC1,PC2,PC3,PC4,tdi --categorical-covariates center \\
        --dir-output "$GWAS_ANALYSIS_ROOT/ukb_center_tdi"

    # Three separate CIMT files; use ALL continuous and categorical columns in them.
    # Current headers: phe1,phe2; age,PC1-PC10; sex,batch.
    bash gwas.sh regenie --grch 37 --pheno "$DIR_UKB/cimt.phe" --pheno-name phe1,phe2 --type qt \\
        --continuous-covariate-file "$DIR_UKB/cimt.cont.cov" \\
        --categorical-covariate-file "$DIR_UKB/cimt.cat.cov" \\
        --dir-output "$GWAS_ANALYSIS_ROOT/cimt"

chrX contract:
    Reads the same imputed chrX.pgen/pvar/psam prefix; never substitutes a male-only copy.
    Uses PLINK sex (1=male, 2=female), build-specific PAR boundaries, and default
    REGENIE dosage compensation (male non-PAR dosage on the 0-2 scale).
    No --skip-dosage-comp, genotype recoding, or automatic PAR removal is applied.
    The X metadata audit checks selected samples, sex concordance and chromosome labels.
    Valid PGEN dosage/hardcall encoding is an input requirement; this audit does not
    inspect or repair genotypes. Existing results with different settings require
    a separate --dir-output or --replace TRUE.

HELP
	gwas_runtime_usage
}

die() {
	echo "ERROR: $*" >&2
	exit 2
}
bool_word() {
	case "${1^^}" in TRUE) echo TRUE ;; FALSE) echo FALSE ;; *) die 'Use TRUE or FALSE' ;; esac
}
tool_path() {
	if command -v "$1" >/dev/null 2>&1; then
		command -v "$1"
	elif [[ -x $GWAS_SOFTWARE_ROOT/bin/$1 ]]; then
		echo "$GWAS_SOFTWARE_ROOT/bin/$1"
	else echo "$1"; fi
}


# 🚩 Defaults and argument parsing
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
helper="$here/f/gwas.py"
module=${1:-}
# Normalize saved commands to the public module names before configuring jobs.
case "$module" in
	prep_gwas) module=sparseGRM ;;
	account | prep_regenie | filter_variants) module=acount ;;
	run_regenie) module=regenie ;;
esac
case "$module" in
	-h | --help | '')
		usage
		exit 0
		;;
	sparseGRM | acount | run_plink2 | regenie | run_saige) shift ;;
	*) die "Unknown module: $module (see --help)" ;;
esac
grch=37 threads=4 memory=24000 min_mac=100 jobs=6
step2_threads=''
global_maf=1e-4
[[ "$module" != regenie ]] || min_mac=5
typed_dir=''
imputed_dir='' phenotype_file=$GWAS_PHENOTYPE_FILE covariate_file=''
continuous_covariate_file='' categorical_covariate_file=''
covariates_set=FALSE categorical_covariates_set=FALSE
phenotypes=height
type=auto event_col='' covariates=age,sex,PC1,PC2,PC3,PC4 categorical_covariates=''
sex=all
output_dir=$GWAS_ANALYSIS_ROOT/self/common
chromosomes=''
extract='' keep='' sparse_grm=TRUE run=TRUE replace=FALSE
plink2=$(tool_path plink2)
plink=$(tool_path plink)
regenie=$(tool_path regenie)
rscript=$(tool_path Rscript)
saige_dir=$GWAS_SOFTWARE_ROOT/SAIGE/extdata
saige_rscript=$rscript
for saige_env in "${SAIGE_ENV_DIR:-$DIR_WORK/.conda/envs/saige}" "$HOME/anaconda3/envs/saige" "$HOME/miniforge3/envs/saige"; do
	if [[ -f "$saige_env/bin/createSparseGRM.R" ]]; then
		saige_dir="$saige_env/bin"
		export SAIGE_ENV_DIR="$saige_env"
		saige_rscript="$here/f/gwas.f.sh"
		break
	fi
done


# 🚩 Command-line arguments
while (($#)); do
	case "$1" in -h | --help)
		usage
		exit 0
		;;
	esac
	(($# >= 2)) || die "$1 requires a value"
	# run/replace are read through indirect expansion in save_jobs.
	# shellcheck disable=SC2034
	case "$1" in
		--grch) grch=$2 ;; --threads) threads=$2 ;; --memory) memory=$2 ;;
		--step2-threads) step2_threads=$2 ;;
		--sex) sex=$2 ;;
		--typ-dir) typed_dir=$2 ;; --imputed-dir) imputed_dir=$2 ;;
		--jobs) jobs=$2 ;;
		--pheno | --phenotype-file) phenotype_file=$2 ;;
		--pheno-name | --phenotypes | --traits) phenotypes=$2 ;;
		--covariate-file) covariate_file=$2 ;; --type) type=$2 ;; --event-col) event_col=$2 ;;
		--continuous-covariate-file) continuous_covariate_file=$2 ;;
		--categorical-covariate-file) categorical_covariate_file=$2 ;;
		--covariates) covariates=$2; covariates_set=TRUE ;;
		--categorical-covariates) categorical_covariates=$2; categorical_covariates_set=TRUE ;;
		--output-dir) output_dir=$2 ;; --chr) chromosomes=$2 ;;
		--extract) extract=$2 ;; --keep) keep=$2 ;; --min-mac) min_mac=$2 ;;
		--global-maf) global_maf=$2 ;;
		--plink2) plink2=$2 ;; --plink) plink=$2 ;; --regenie) regenie=$2 ;; --rscript) rscript=$2 ;;
		--saige-rscript) saige_rscript=$2 ;; --saige-dir) saige_dir=$2 ;;
		--sparse-grm) sparse_grm=$(bool_word "$2") ;;
		--run-cmd) run=$(bool_word "$2") ;;
		--replace) replace=$(bool_word "$2") ;;
		*) die "Unknown option: $1" ;;
	esac
	shift 2
done
[[ "$grch" == 37 || "$grch" == 38 ]] || die '--grch must be 37 or 38'
case "$type" in auto | qt | ordinal | bt | t2e) ;; *) die 'Invalid --type' ;; esac
case "$sex" in all|male|female) ;; *) die '--sex must be all, male or female' ;; esac
[[ $module == regenie || $sex == all ]] || die '--sex is supported by the regenie module'
step2_threads=${step2_threads:-$threads}
for number in threads step2_threads memory min_mac jobs; do
	[[ ${!number} =~ ^[0-9]+$ ]] || die "Invalid --${number//_/-}"
	printf -v "$number" '%d' "$((10#${!number}))"
done
((threads > 0 && step2_threads > 0 && jobs > 0 && memory >= 640 && min_mac > 0)) || die 'Invalid threads/step2-threads/jobs/memory/min-mac'
python3 - "$global_maf" <<'PY'
import sys
try:
    assert 0 < float(sys.argv[1]) <= 0.5
except (ValueError, AssertionError):
    sys.exit('ERROR: --global-maf must be > 0 and <= 0.5')
PY
typed_dir=${typed_dir:-$DIR_GEN/ukb/$grch/typ}
imputed_dir=${imputed_dir:-$DIR_GEN/ukb/$grch/imp}
selection_dir="$imputed_dir/.${module}"
input_dir="$imputed_dir"
if [[ "$module" == sparseGRM ]]; then
	input_dir="$typed_dir"
	selection_dir="$typed_dir/.prep"
elif [[ "$module" == regenie || "$module" == run_* ]]; then
	selection_dir="$output_dir/.${module}"
fi
mkdir -p "$selection_dir"
selection_log="$selection_dir/chromosomes.log"
[[ $module != regenie ]] || selection_log="$selection_dir/chromosomes.${LSB_JOBID:-$$}.log"
printf 'CHR\tSTATUS\tREASON\n' >"$selection_log"
# Genome-wide REGENIE must not silently omit chromosomes still being transferred.
if [[ $module == regenie && ( -z "$chromosomes" || "$chromosomes" == all ) ]]; then
	chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X
fi
if [[ -z "$chromosomes" || "$chromosomes" == all ]]; then
	selected=()
	for chr in {1..22} X; do
		if [[ -e "$input_dir/chr$chr.pgen" || -e "$input_dir/chr$chr.pvar" || -e "$input_dir/chr$chr.psam" ]]; then
			selected+=("$chr")
		else
			printf '%s\tNOT_PRESENT\tno genotype input at %s/chr%s\n' "$chr" "$input_dir" "$chr" >>"$selection_log"
		fi
	done
	((${#selected[@]})) || die "No chromosome inputs found in $input_dir; see $selection_log"
	chromosomes=$(
		IFS=,
		echo "${selected[*]}"
	)
fi
read -r -a chrs <<<"${chromosomes//,/ }"
((${#chrs[@]})) || die 'Empty chromosome list'
declare -A seen_chr=()
for chr in "${chrs[@]}"; do
	[[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2]|X|Y)$ ]] || die "Invalid chromosome: $chr"
	[[ ! ${seen_chr[$chr]+yes} ]] || die "Duplicate chromosome: $chr"
	seen_chr[$chr]=1
	printf '%s\tSELECTED\trequested or discovered input\n' "$chr" >>"$selection_log"
	for ext in pgen pvar psam; do
		[[ -e "$input_dir/chr$chr.$ext" ]] || printf '%s\tINPUT_MISSING\t%s/chr%s.%s\n' "$chr" "$input_dir" "$chr" "$ext" >>"$selection_log"
	done
done
mapfile -t chrs < <(printf '%s\n' "${chrs[@]}" | sort -V)
chromosomes=$(
	IFS=,
	echo "${chrs[*]}"
)
covariate_file=${covariate_file:-$phenotype_file}
if [[ -n "$continuous_covariate_file" && "$covariates_set" == FALSE ]]; then covariates='*'; fi
if [[ -n "$categorical_covariate_file" && "$categorical_covariates_set" == FALSE ]]; then categorical_covariates='*'; fi
[[ -n "$phenotypes" ]] || die 'Empty phenotype list'
[[ "$phenotypes" != ,* && "$phenotypes" != *, && "$phenotypes" != *,,* ]] || die 'Empty trait name in --pheno-name'


# Validate the complete request before counts/model fitting. Outer HPC planning
# still does no data reads; this executes on the allocated worker only.
if [[ $run == TRUE ]]; then
	required_inputs=()
	for chr in "${chrs[@]}"; do
		for ext in pgen pvar psam; do required_inputs+=("$input_dir/chr$chr.$ext"); done
	done
	if [[ $module == regenie ]]; then
		required_inputs+=("$typed_dir/typ.prune.bed" "$typed_dir/typ.prune.bim" "$typed_dir/typ.prune.fam" "$phenotype_file"
			"${continuous_covariate_file:-$covariate_file}")
		[[ -z $categorical_covariates ]] || required_inputs+=("${categorical_covariate_file:-$covariate_file}")
		[[ -z $keep ]] || required_inputs+=("$keep")
		[[ -z $extract ]] || required_inputs+=("$extract")
	fi
	missing_inputs=0
	declare -A checked_inputs=()
	for required_input in "${required_inputs[@]}"; do
		[[ ! ${checked_inputs[$required_input]+yes} ]] || continue
		checked_inputs[$required_input]=1
		if [[ ! -r $required_input || ! -s $required_input ]]; then
			printf 'INPUT_MISSING_OR_EMPTY: %s\n' "$required_input" >&2
			missing_inputs=$((missing_inputs + 1))
		fi
		if [[ -e $required_input.filepart ]]; then
			printf 'INPUT_TRANSFER_IN_PROGRESS: %s.filepart\n' "$required_input" >&2
			missing_inputs=$((missing_inputs + 1))
		fi
	done
	((missing_inputs == 0)) || die "Input preflight failed ($missing_inputs problems). Finish syncing the matching .pgen/.pvar/.psam files before resubmitting. For an intentional subset, specify --chr. No GWAS tools were started."
fi


# 🚩 Command-plan utilities. The Python helper only handles JSON/checkpoints,
# logged execution, and PGEN file views; it does not define GWAS commands.
jobs_file=$(mktemp)
trap 'rm -f -- "$jobs_file"' EXIT
cmd=() inputs=() outputs=() pfile_args=() geno_inputs=()
add_job() {
	python3 "$helper" plan-node "$1" "${2:-}" "${3:-}" \
		--inputs "${inputs[@]}" --outputs "${outputs[@]}" -- "${cmd[@]}" >>"$jobs_file"
}
save_jobs() {
	local key config=()
	for key in module grch threads memory jobs typed_dir imputed_dir phenotype_file covariate_file \
		continuous_covariate_file categorical_covariate_file \
		phenotypes type event_col covariates categorical_covariates sex output_dir chromosomes \
		extract keep min_mac global_maf plink2 plink regenie rscript saige_rscript saige_dir sparse_grm run replace; do
		config+=("$key=${!key}")
	done
	[[ "$module" != regenie ]] || config+=("step2_threads=$step2_threads")
	python3 "$helper" finalize-plan "$1" "$jobs_file" "${config[@]}"
}
genotype() {
	pfile_args=(--pfile "$1")
	geno_inputs=("$1.pgen" "$1.pvar" "$1.psam")
}


# 🚩 Shared typed-genotype pruning and sparse GRM
sparseGRM() {
	local root cache chr src dst merged prefix actual
	local beds=() merge_inputs=() base=("$plink2" --threads "$threads" --memory "$memory")
	root=$(realpath -m -- "$typed_dir")
	cache="$root/.prep"
	mkdir -p "$cache"
	: >"$jobs_file"
	printf '6\t25000000\t35000000\tMHC\n' >"$cache/mhc.range"
	for chr in "${chrs[@]}"; do
		src="$root/chr$chr"
		dst="$cache/chr$chr.prune"
		genotype "$src"
		# k=0 explicitly preserves the fixed HWE cutoff for the shared GRM panel.
		cmd=("${base[@]}" "${pfile_args[@]}" --maf .01 --mac 100 --max-alleles 2
			--geno .1 --hwe 1e-15 0 --indep-pairwise 1000kb 0.1 --out "$dst")
		[[ "$chr" != 6 ]] || cmd+=(--exclude range "$cache/mhc.range")
		inputs=("${geno_inputs[@]}")
		outputs=("$dst.prune.in")
		add_job "chr$chr.prune"
		cmd=("${base[@]}" "${pfile_args[@]}" --extract "$dst.prune.in" --make-bed --out "$dst")
		inputs=("${geno_inputs[@]}" "$dst.prune.in")
		outputs=("$dst.bed" "$dst.bim" "$dst.fam")
		add_job "chr$chr.bed"
		beds+=("$dst")
		merge_inputs+=("${outputs[@]}")
	done
	merged="$root/typ.prune"
	printf '%s\n' "${beds[@]:1}" >"$cache/merge.list"
	cmd=("$plink" --bfile "${beds[0]}" --make-bed --out "$merged" --threads "$threads" --memory "$memory")
	((${#beds[@]} <= 1)) || cmd+=(--merge-list "$cache/merge.list")
	inputs=("${merge_inputs[@]}")
	outputs=("$merged.bed" "$merged.bim" "$merged.fam")
	add_job merge
	if [[ "$sparse_grm" == TRUE ]]; then
		prefix="$root/typ"
		actual="${prefix}_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
		cmd=("$saige_rscript" "$saige_dir/createSparseGRM.R" "--plinkFile=$merged" "--nThreads=$threads"
			--numRandomMarkerforSparseKin=2000 --relatednessCutoff=0.125 "--outputPrefix=$prefix")
		cmd=(python3 "$helper" create-grm "$actual" "$prefix" "${cmd[@]}")
		inputs=("$merged.bed" "$merged.bim" "$merged.fam")
		outputs=("$prefix.sparseGRM.mtx" "$prefix.sparseGRM.mtx.sampleIDs.txt")
		add_job sparseGRM
	fi
	save_jobs "$cache"
}

# Shared imputed counts are independent of phenotypes and association --keep.
# Reuse existing counts, including those produced before checkpoint tracking.
# Only an explicit acount --replace TRUE may rebuild shared counts.
acount() (
	local replace=$replace
	[[ "$module" != regenie ]] || replace=FALSE
	local module=acount cache="$imputed_dir/.regenie-filter" chr prefix
	local -a cmd inputs outputs pfile_args geno_inputs
	mkdir -p "$cache"
	# Lock before inspecting/writing the shared plan and counts. Phenotype jobs
	# wait only for preparation, then run their models independently. The subshell
	# closes the descriptor before REGENIE Step 1 starts (also on failure).
	local acount_lock
	exec {acount_lock}>"$cache/prepare.lock"
	flock -x "$acount_lock"
	: >"$jobs_file"
	for chr in "${chrs[@]}"; do
		prefix="$imputed_dir/chr$chr"
		genotype "$prefix"
		if [[ -s "$prefix.acount" && "$replace" != TRUE ]]; then
			printf 'REUSE %s\n' "$prefix.acount"
		else
			cmd=("$plink2" --threads "$threads" --memory "$memory" "${pfile_args[@]}"
				--nonfounders --freq counts --out "$prefix")
			inputs=("${geno_inputs[@]}")
			outputs=("$prefix.acount")
			add_job "chr$chr.acount"
		fi
		# Accept extraction lists created by the documented standalone AWK loop.
		# Its summary records the MAF threshold; changed counts/rules require rebuilding.
		if [[ "$replace" != TRUE && -s "$prefix.acount" && -s "$prefix.extract" && -s "$prefix.extract.counts.tsv" &&
			! "$prefix.acount" -nt "$prefix.extract" && ! "$here/f/gwas.extract.awk" -nt "$prefix.extract" ]] &&
			LC_ALL=C awk -v maf="$global_maf" '
				NR == 1 { header = ($1 == "TOTAL" && $2 == "RETAINED" && $3 == "MIN_MAF" && NF == 3) }
				NR == 2 { valid = (NF == 3 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $2 > 0 && $2 <= $1 && $3+0 == maf+0) }
				END { exit !(NR == 2 && header && valid) }
			' "$prefix.extract.counts.tsv"; then
			printf 'SKIP existing %s (MAF >= %s)\n' "$prefix.extract" "$global_maf"
			continue
		fi
		cmd=(bash "$here/f/gwas.f.sh" extract "$prefix.acount" "$prefix.extract" "$global_maf")
		inputs=("$prefix.acount" "$here/f/gwas.f.sh" "$here/f/gwas.extract.awk")
		outputs=("$prefix.extract" "$prefix.extract.counts.tsv")
		add_job "chr$chr.extract"
	done
	save_jobs "$cache"
)


# 🚩 Shared phenotype setup and output collection
setup_trait() {
	skip_trait=FALSE
	trait=$1
	trait_type=$type
	[[ -n "$trait" && "$trait" != . && "$trait" != .. && "$trait" != *[/\\]* && "$trait" != *$'\n'* && "$trait" != *$'\r'* ]] || die 'Invalid trait name'
	if [[ "$trait_type" == auto ]]; then
		case "$trait" in
			bmi | height | ldl | *.adu) trait_type=qt ;; bald) trait_type=ordinal ;;
			bald12 | bald13 | bald14 | cvd_cad.Yt2e) trait_type=bt ;; cvd_cad.t2e) trait_type=t2e ;;
			*) die "Specify --type for $trait" ;;
		esac
	fi
	root="$(realpath -m -- "$output_dir")/$trait/gwas"
	work="$root/$method"
	mkdir -p "$work"
	event=${event_col:-${trait%.t2e}.Yt2e}
	phe="$work/analysis.phe"
	covars="$phe.covars"
	pruned="$typed_dir/typ.prune"
	grm="$typed_dir/typ.sparseGRM.mtx"
	keep_args=()
	extract_args=()
	raw=()
	: >"$jobs_file"
	[[ -z "$keep" ]] || keep_args=(--keep "$keep")
	[[ -z "$extract" ]] || extract_args=(--extract "$extract")
	if [[ "$trait_type" == t2e && "$method" != regenie ]]; then
		printf '%s\n' "$method: time-to-event is unsupported; use regenie --type t2e" | tee "$work/UNSUPPORTED.txt"
		skip_trait=TRUE
		return 0
	fi
	cmd=("$rscript" "$here/f/gwas.inputs.R" "$phenotype_file" "$covariate_file" "$trait" "$trait_type" "$event" "$covariates" "$categorical_covariates" "$phe")
	inputs=("$phenotype_file" "$covariate_file" "$here/f/gwas.inputs.R")
	if [[ -n "$continuous_covariate_file" || -n "$categorical_covariate_file" ]]; then
		cmd+=("$continuous_covariate_file" "$categorical_covariate_file")
		# Only the covariate files actually used should be required and fingerprinted.
		inputs=("$phenotype_file" "${continuous_covariate_file:-$covariate_file}" "$here/f/gwas.inputs.R")
		if [[ -n "$categorical_covariates" ]]; then
			inputs+=("${categorical_covariate_file:-$covariate_file}")
		fi
	fi
	outputs=("$phe" "$covars")
	add_job phenotype
	common_inputs=("$phe" "$covars")
	[[ -z "$keep" ]] || common_inputs+=("$keep")
}
association_inputs() {
	genotype "$prefix"
	inputs=("${geno_inputs[@]}" "${common_inputs[@]}")
	[[ -z "$extract" ]] || inputs+=("$extract")
}
finish_trait() {
	local final="$root/$trait.$method.gz"
	# REGENIE merge writes P to 3 decimal places, using scientific notation below 0.001.
	cmd=(python3 "$here/f/gwas.py" merge "$method" "$final" "${raw[@]}")
	inputs=("${raw[@]}" "$here/f/gwas.py")
	outputs=("$final")
	add_job merge-results
	cmd=(python3 "$helper" build-metadata "$final" "$grch")
	inputs=("$final")
	outputs=("$final.grch")
	add_job build-metadata
	save_jobs "$work"
}


# 🚩 PLINK2 linear / logistic association
run_plink2() {
	method=plink2
	for trait in "${trait_list[@]}"; do
		setup_trait "$trait"
		[[ "$skip_trait" != TRUE ]] || continue
		for chr in "${chrs[@]}"; do
			prefix="$imputed_dir/chr$chr"
			dest="$work/chr$chr"
			association_inputs
			glm=(hide-covar allow-no-covars omit-ref no-x-sex 'cols=+omitted,+a1freq,+nobs')
			result="$dest.$trait.glm.linear"
			if [[ "$trait_type" == bt ]]; then
				glm+=(firth-fallback)
				result="$dest.$trait.glm.logistic.hybrid"
			fi
			cmd=("$plink2" --threads "$threads" --memory "$memory" "${pfile_args[@]}" --glm "${glm[@]}"
				--pheno "$phe" --pheno-name "$trait" --no-psam-pheno --mac "$min_mac" --max-alleles 2 --out "$dest" "${keep_args[@]}" "${extract_args[@]}")
			[[ "$trait_type" != bt ]] || cmd+=(--1)
			[[ "$chr" != X ]] || cmd+=(--xchr-model 2)
			outputs=("$result")
			add_job "chr$chr.association" "$covars" "$method"
			raw+=("$result")
		done
		finish_trait
	done
}


# 🚩 REGENIE model fitting and per-chromosome association
regenie() {
	method=regenie
	acount
	for trait in "${trait_list[@]}"; do
		setup_trait "$trait"
		[[ "$skip_trait" != TRUE ]] || continue
		# Validate sex metadata before fitting any model; never rewrite genotypes.
		cmd=(python3 "$here/f/gwas.audit.py" --phe "$phe" --fam "$pruned.fam"
			--build "$grch" --sex "$sex" --output "$work/genotype.qc.json")
		inputs=("${common_inputs[@]}" "$pruned.fam" "$here/f/gwas.audit.py")
		if [[ ${seen_chr[X]+yes} ]]; then
			cmd+=(--psam "$imputed_dir/chrX.psam" --pvar "$imputed_dir/chrX.pvar")
			inputs+=("$imputed_dir/chrX.psam" "$imputed_dir/chrX.pvar")
		fi
		[[ -z "$keep" ]] || cmd+=(--keep "$keep")
		outputs=("$work/genotype.qc.json")
		add_job genotype-audit
		flags=(--qt)
		phflags=(--phenoCol "$trait")
		case "$trait_type" in
			bt) flags=(--bt) ;;
			ordinal) flags=(--qt --force-qt) ;;
			t2e)
				flags=(--t2e)
				phflags=(--phenoColList "$trait" --eventColList "$event")
				;;
		esac
		common=("$regenie" --threads "$threads" --bsize 1000 --phenoFile "$phe" "${phflags[@]}" "${flags[@]}" --strict --par-region "b$grch" "${keep_args[@]}")
		[[ $sex == all ]] || common+=(--sex-specific "$sex")
		step1="$work/step1"
		cmd=("${common[@]}" --step 1 --bed "$pruned" --lowmem --lowmem-prefix "$work/tmp" --out "$step1")
		inputs=("${common_inputs[@]}" "$pruned.bed" "$pruned.bim" "$pruned.fam")
		outputs=("${step1}_pred.list" "${step1}_1.loco")
		add_job null-model "$covars" "$method"
		common=("$regenie" --threads "$step2_threads" --bsize 1000 --phenoFile "$phe" "${phflags[@]}" "${flags[@]}" --strict --par-region "b$grch" "${keep_args[@]}")
		[[ $sex == all ]] || common+=(--sex-specific "$sex")
		for chr in "${chrs[@]}"; do
			prefix="$imputed_dir/chr$chr"
			dest="$work/chr$chr"
			genotype "$prefix"
			shared_extract="$prefix.extract"
			if [[ -n "$extract" ]]; then
				shared_extract="$work/chr$chr.extract"
				cmd=(bash "$here/f/gwas.f.sh" extract --intersect "$extract" "$prefix.extract" "$shared_extract")
				inputs=("$extract" "$prefix.extract" "$here/f/gwas.f.sh")
				outputs=("$shared_extract")
				add_job "chr$chr.extract"
			fi
			association_inputs
			inputs+=("${step1}_pred.list" "${step1}_1.loco" "$shared_extract")
			cmd=("${common[@]}" --step 2 --pgen "$prefix" --pred "${step1}_pred.list" --minMAC "$min_mac" --out "$dest" --extract "$shared_extract")
			[[ "$trait_type" != bt ]] || cmd+=(--firth --approx --pThresh 0.01)

			result="${dest}_$trait.regenie"
			outputs=("$result")
			add_job "chr$chr.association" "$covars" "$method"
			raw+=("$result")
		done
		finish_trait
	done
}


# 🚩 SAIGE sparse null model and dosage-based association
run_saige() (
	mkdir -p "$output_dir"
	GU_SAIGE_CACHE_DIR=$(mktemp -d "$output_dir/.saige-dosage.XXXXXX")
	export GU_SAIGE_CACHE_DIR
	trap 'rm -rf -- "$GU_SAIGE_CACHE_DIR"' EXIT
	trap 'exit 130' INT
	trap 'exit 143' TERM
	trap 'exit 129' HUP
	method=saige
	[[ -z "$keep" ]] || die 'SAIGE --keep: provide a subset --pheno file instead'
	for trait in "${trait_list[@]}"; do
		setup_trait "$trait"
		[[ "$skip_trait" != TRUE ]] || continue
		step1="$work/step1"
		saige_type=quantitative
		[[ "$trait_type" != bt ]] || saige_type=binary
		cmd=("$saige_rscript" "$saige_dir/step1_fitNULLGLMM.R" "--plinkFile=$pruned"
			"--phenoFile=$phe" "--phenoCol=$trait" --sampleIDColinphenoFile=IID "--traitType=$saige_type"
			--invNormalize=FALSE --useSparseGRMtoFitNULL=TRUE --skipVarianceRatioEstimation=FALSE
			--IsOverwriteVarianceRatioFile=TRUE --LOCO=FALSE "--nThreads=$threads" "--outputPrefix=$step1"
			"--sparseGRMFile=$grm" "--sparseGRMSampleIDFile=$grm.sampleIDs.txt")
		inputs=("${common_inputs[@]}" "$pruned.bed" "$pruned.bim" "$pruned.fam" "$grm" "$grm.sampleIDs.txt")
		outputs=("$step1.rda" "$step1.varianceRatio.txt")
		add_job null-model "$covars" "$method"
		for chr in "${chrs[@]}"; do
			prefix="$imputed_dir/chr$chr"
			dest="$work/chr$chr"
			vcf="$dest.dosage.vcf.gz"
			genotype "$prefix"
			# Force DS output, retaining dosage for SAIGE instead of hard-call BED.
			cmd=("$plink2" --threads "$threads" --memory "$memory" "${pfile_args[@]}" --mac "$min_mac"
				--max-alleles 2 --export vcf bgz vcf-dosage=DS-force id-paste=iid --out "$dest.dosage" "${extract_args[@]}")
			inputs=("${geno_inputs[@]}")
			[[ -z "$extract" ]] || inputs+=("$extract")
			outputs=("$vcf")
			add_job "chr$chr.export"
			cmd=("$(tool_path tabix)" -f -C -p vcf "$vcf")
			inputs=("$vcf")
			outputs=("$vcf.csi")
			add_job "chr$chr.index"
			result="$dest.saige.txt"
			association_inputs
			inputs+=("$vcf" "$vcf.csi" "$step1.rda" "$step1.varianceRatio.txt")
			cmd=("$saige_rscript" "$saige_dir/step2_SPAtests.R" "--vcfFile=$vcf" "--vcfFileIndex=$vcf.csi" --vcfField=DS "--chrom=$chr"
				"--GMMATmodelFile=$step1.rda" "--varianceRatioFile=$step1.varianceRatio.txt"
				"--SAIGEOutputFile=$result" --LOCO=FALSE "--minMAC=$min_mac" --is_output_moreDetails=TRUE --is_overwrite_output=TRUE)
			[[ "$trait_type" != bt ]] || cmd+=(--is_Firth_beta=TRUE --pCutoffforFirth=0.01)
			outputs=("$result")
			add_job "chr$chr.association"
			raw+=("$result")
		done
		finish_trait
	done
)


# 🚩 Dispatch one independent module, never a sequence of GWAS methods.
IFS=, read -r -a trait_list <<<"$phenotypes"
if [[ "$module" != sparseGRM && "$module" != acount && "$run" == TRUE ]]; then
	# Check every requested trait before starting any expensive model fit.
	python3 - "$phenotype_file" "$phenotypes" <<'PY'
import gzip
import sys
path, requested = sys.argv[1:]
opener = gzip.open if path.endswith('.gz') else open
with opener(path, 'rt') as handle:
    columns = handle.readline().split()
missing = [name for name in requested.split(',') if name not in columns]
if missing:
    sys.exit('ERROR: Missing phenotype columns in ' + path + ': ' + ', '.join(missing)
             + '. Generate the requested phenotypes before running GWAS.')
PY
fi
"$module"
