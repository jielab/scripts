#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    compare  Effect-size and allele-frequency comparisons
             Mixed or unknown builds match by SNP ID and alleles
    ldsc     Heritability and genetic correlations
    shiny    Manhattan tracks and LD viewer

Main options:
    --submit auto|local|hpc       Execution mode [auto]
    --dir-project DIR / --gwas-files CSV  Project or explicit GWAS files
    --category NAME --anchor NAME  Project category and reference trait
    --grch 37|38 --require-grch 37|38  Build alignment/requirement
    --dir-output DIR  Result directory
    --p-threshold N --significant first|either|both  Selection
    --run-munge TRUE|FALSE  Prepare LDSC caches [TRUE]
    --ref-ld-chr PREFIX --w-ld-chr PREFIX  LD references
    --port N --host IP  Shiny address [3841 / 127.0.0.1]
    -h, --help  Show options and examples below

Examples:
    source ../0f/0os.sh
    ./compare.sh compare --dir-project "$DIR_GWAS/4grid" --anchor height.EUR --p-threshold 1e-6
    ./compare.sh ldsc --project-dir "$DIR_GWAS/4grid" --category common --run-munge FALSE --run-h2 TRUE --run-rg TRUE --run TRUE --output-dir "$ANALYSIS_ROOT/gwas/4grid/ldsc"
PUBLIC_ENTRY_HELP
    exit 0
fi
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
gwas_runtime_init "compare.sh" "$@" || exit $?
set -- "${GWAS_ARGS[@]}"
source "$GWAS_ENTRY_ROOT/../0f/console.sh"
set -euo pipefail
# Defaults, CLI and module commands live here. Project discovery/build alignment
# analysis functions live in f/ and the Shiny application lives in shiny/.


# 🚩 Command-line help
usage() {
	bash "${BASH_SOURCE[0]}" --help
	gwas_runtime_usage
}


# 🚩 Defaults and module selection
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
helper="$here/f/compare.py"
# shellcheck source=f/compare.f.sh
source "$here/f/compare.f.sh"

module=${1:-}
case "$module" in
	compare | ldsc | shiny) shift ;;
	-h | --help | '')
		usage
		exit 0
		;;
	*)
		echo "ERROR: unknown module: $module (expected compare, ldsc or shiny; see --help)" >&2
		exit 1
		;;
esac

project_dir='' gwas_files='' category=common anchor=''
grch='' require_grch='' check_only=FALSE
output_dir="$GWAS_ANALYSIS_ROOT/$module"

# Comparison settings; these values are passed explicitly to the R helper.
labels='' compare_beta=TRUE compare_eaf=TRUE
p_threshold=5e-8 significant=first

# LDSC settings; no fallback sample size is assumed.
ldsc_software_dir=$GWAS_SOFTWARE_ROOT/ldsc
conda=$(type -P conda || true)
conda=${conda:-$HOME/anaconda3/bin/conda}
conda_env=ldsc ldsc_python=''
merge_alleles=$GWAS_DATA_ROOT/LD/ldsc/hm3/w_hm3.snplist
ref_ld_chr=$GWAS_DATA_ROOT/LD/ldsc/1000G/1000G_Phase3_ldscores/LDscore.
w_ld_chr=$GWAS_DATA_ROOT/LD/ldsc/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
sample_size='' missing_n=error
run_munge=TRUE run_h2=TRUE run_rg=TRUE run=TRUE

# Interactive block viewer. All settings remain visible at the entry point.
block_dir=$GWAS_DATA_ROOT/LD/block
ld_dir=$GWAS_DATA_ROOT/LD/csx ld_python="${GWAS_LD_PYTHON:-python3}" ld_max_snps=300
port=3841 host=127.0.0.1 launch_browser=FALSE max_points=12000 races=''
chain_dir=$GWAS_FILES_ROOT/liftOver liftover_bin=$GWAS_SOFTWARE_ROOT/bin/liftOver
reference_dir=$GWAS_DATA_ROOT/gen/1hgp gene_dir=$GWAS_FILES_ROOT


# 🚩 Command-line arguments
while (($#)); do
	case "$1" in
		-h | --help)
			usage
			exit 0
			;;
		# Preserve the --option=value form supported by the former Python entry point.
		--*=*)
			option=${1%%=*}
			value=${1#*=}
			shift
			set -- "$option" "$value" "$@"
			;;
	esac
	compare_need_arg_value "$1" "${2-}"
	case "$1" in
		--project-dir | --dir-gwas) project_dir=$2 ;;
		--gwas-files) gwas_files=$2 ;;
		--category) category=$2 ;;
		--anchor) anchor=$2 ;;
		--grch) grch=$2 ;;
		--require-grch) require_grch=$2 ;;
		--output-dir | --dir-out) output_dir=$2 ;;
		--check-only) check_only=$2 ;;
		*)
			case "$module:$1" in
				compare:--labels | shiny:--labels) labels=$2 ;;
				compare:--compare-beta) compare_beta=$2 ;;
				compare:--compare-EAF | compare:--compare-eaf) compare_eaf=$2 ;;
				compare:--p-threshold | shiny:--p-threshold) p_threshold=$2 ;;
				compare:--significant) significant=$2 ;;
				ldsc:--ldsc-software-dir) ldsc_software_dir=$2 ;;
				ldsc:--conda) conda=$2 ;;
				ldsc:--conda-env) conda_env=$2 ;;
				ldsc:--python) ldsc_python=$2 ;;
				ldsc:--merge-alleles) merge_alleles=$2 ;;
				ldsc:--ref-ld-chr) ref_ld_chr=$2 ;;
				ldsc:--w-ld-chr) w_ld_chr=$2 ;;
				ldsc:--N) sample_size=$2 ;;
				ldsc:--missing-n) missing_n=$2 ;;
				ldsc:--run-munge) run_munge=$2 ;;
				ldsc:--run-h2) run_h2=$2 ;;
				ldsc:--run-rg) run_rg=$2 ;;
				ldsc:--run) run=$2 ;;
				shiny:--races) races=$2 ;;
				shiny:--block-dir) block_dir=$2 ;;
				shiny:--ld-dir) ld_dir=$2 ;;
				shiny:--ld-python) ld_python=$2 ;;
				shiny:--ld-max-snps) ld_max_snps=$2 ;;
				shiny:--port) port=$2 ;;
				shiny:--host) host=$2 ;;
				shiny:--launch-browser) launch_browser=$2 ;;
				shiny:--max-points) max_points=$2 ;;
				shiny:--chain-dir) chain_dir=$2 ;;
				shiny:--liftover-bin) liftover_bin=$2 ;;
				shiny:--reference-dir) reference_dir=$2 ;;
				shiny:--gene-dir) gene_dir=$2 ;;
				*) compare_die "Unknown $module option: $1 (see --help)" ;;
			esac
			;;
	esac
	shift 2
done
gwas_compare_validate_options


# 🚩 Build the module command from the settings above
cmd=(python3 "$helper" "$module" --category "$category"
	--output-dir "$output_dir" --check-only "$check_only")
if [[ -n "$project_dir" ]]; then
	cmd+=(--project-dir "$project_dir")
else
	cmd+=(--gwas-files "$gwas_files")
fi
[[ -z "$anchor" ]] || cmd+=(--anchor "$anchor")
[[ -z "$grch" ]] || cmd+=(--grch "$grch")
[[ -z "$require_grch" ]] || cmd+=(--require-grch "$require_grch")

case "$module" in
	compare)
		cmd+=(--compare-beta "$compare_beta" --compare-EAF "$compare_eaf"
			--p-threshold "$p_threshold" --significant "$significant")
		[[ -z "$labels" ]] || cmd+=(--labels "$labels")
		;;
	ldsc)
		cmd+=(--ldsc-software-dir "$ldsc_software_dir" --conda "$conda" --conda-env "$conda_env"
			--merge-alleles "$merge_alleles" --ref-ld-chr "$ref_ld_chr" --w-ld-chr "$w_ld_chr"
			--missing-n "$missing_n" --run-munge "$run_munge" --run-h2 "$run_h2"
			--run-rg "$run_rg" --run "$run")
		[[ -z "$ldsc_python" ]] || cmd+=(--python "$ldsc_python")
		[[ -z "$sample_size" ]] || cmd+=(--N "$sample_size")
		;;
	shiny)
		cmd+=(--ld-dir "$ld_dir" --ld-python "$ld_python" --ld-max-snps "$ld_max_snps"
			--block-dir "$block_dir" --port "$port" --host "$host"
			--launch-browser "$launch_browser" --max-points "$max_points"
			--p-threshold "$p_threshold" --chain-dir "$chain_dir" --liftover-bin "$liftover_bin"
			--reference-dir "$reference_dir" --gene-dir "$gene_dir")
		[[ -z "$labels" ]] || cmd+=(--labels "$labels")
		[[ -z "$races" ]] || cmd+=(--races "$races")
		;;
esac
exec "${cmd[@]}"
