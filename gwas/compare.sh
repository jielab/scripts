#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# A foreground Shiny server must keep its listening URL and live diagnostics.
source /mnt/d/scripts/0f/console.sh
set -euo pipefail
# Defaults, CLI and module commands live here. Project discovery/build alignment
# analysis functions live in f/ and the Shiny application lives in shiny/.


# 🚩 Command-line help
usage() {
	cat <<'HELP'
Modules:
    compare  Effect-size and allele-frequency comparisons
    ldsc     Heritability and genetic correlations
    shiny    Manhattan tracks and LD viewer

Main options:
    --dir-project DIR / --gwas-files CSV  Project or explicit GWAS files
    --category NAME --anchor NAME  Project category and reference trait
    --grch 37|38 --require-grch 37|38  Build alignment/requirement
    --dir-output DIR  Result directory
    --p-threshold N --significant first|either|both  Selection
    --run-munge TRUE|FALSE  Prepare LDSC caches [TRUE]
    --ref-ld-chr PREFIX --w-ld-chr PREFIX  LD references
    --port N --host IP  Shiny address [3841 / 127.0.0.1]
    -h, --help  Help; viewer recipe and module options: 0.README.md

Example:
    cd /mnt/d/scripts/gwas
    gwas_root=/mnt/f/gwas/4grid/common

    # 1. All 15 GWAS: height.EUR versus the other 14, using only SNPs with
    # P <= 5e-8 in height.EUR. Inputs must share a build for CHR/POS matching.
    # For mixed builds, add --grch 38 to align cached comparison copies.
    # Output: 14 *.BETA.png and 14 *.EAF.png, plus harmonized data and QC tables.
    ./compare.sh compare \
      --project-dir "${gwas_root%/*}" --category "${gwas_root##*/}" --anchor height.EUR \
      --compare-beta TRUE --compare-EAF TRUE \
      --p-threshold 5e-8 --significant first \
      --output-dir /mnt/d/analysis/gwas/4grid/compare

    # 2. All 15 GWAS: h2 and all 105 pairwise genetic correlations, with rg.png.
    # Reuse <trait>/gwas/<trait>.sumstats.gz beside each <trait>.gz; do not munge.
    # Every adjacent sumstats file must already exist and contain SNP A1 A2 Z N.
    # The configured LD-score reference is shared by all inputs; it is not
    # selected by ancestry. Override --ref-ld-chr/--w-ld-chr when appropriate.
    ./compare.sh ldsc \
      --project-dir "${gwas_root%/*}" --category "${gwas_root##*/}" \
      --run-munge FALSE --run-h2 TRUE --run-rg TRUE --run TRUE \
      --output-dir /mnt/d/analysis/gwas/4grid/ldsc

    # 3. R Shiny: the same six height tracks (EUR, AFR, EAS, HIS, SAS, ALL).
    # Generate .thin.gz first using format.sh thin (see its examples).
    # Hover over a plot and use the camera button to download a PNG.
    trait=height
    gwas_files=()
    for race in EUR AFR EAS HIS SAS ALL; do
      gwas_files+=("$gwas_root/$trait.$race/gwas/$trait.$race.thin.gz")
    done
    gwas_csv=$(IFS=,; echo "${gwas_files[*]}")
    ./compare.sh shiny \
      --gwas-files "$gwas_csv" \
      --block-dir /mnt/f/LD/block --ld-dir /mnt/f/LD/csx --port 3841 \
      --output-dir "/mnt/d/analysis/gwas/4grid/shiny/$trait"

    ./compare.sh compare --dir-project /mnt/f/gwas/4grid --anchor height.EUR --p-threshold 1e-6
    ./compare.sh compare --dir-project /mnt/f/gwas/4grid
HELP
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
output_dir="/mnt/d/analysis/gwas/$module"

# Comparison settings; these values are passed explicitly to the R helper.
labels='' compare_beta=TRUE compare_eaf=TRUE
p_threshold=5e-8 significant=first

# LDSC settings; no fallback sample size is assumed.
ldsc_software_dir=/mnt/d/software/ldsc
conda=$(type -P conda || true)
conda=${conda:-$HOME/anaconda3/bin/conda}
conda_env=ldsc ldsc_python=''
merge_alleles=/mnt/f/LD/ldsc/hm3/w_hm3.snplist
ref_ld_chr=/mnt/f/LD/ldsc/1000G/1000G_Phase3_ldscores/LDscore.
w_ld_chr=/mnt/f/LD/ldsc/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
sample_size='' missing_n=error
run_munge=TRUE run_h2=TRUE run_rg=TRUE run=TRUE

# Interactive block viewer. All settings remain visible at the entry point.
block_dir=/mnt/f/LD/block
ld_dir=/mnt/f/LD/csx ld_python="$HOME/anaconda3/bin/python" ld_max_snps=300
port=3841 host=127.0.0.1 launch_browser=FALSE max_points=12000 races=''
chain_dir=/mnt/d/files/liftOver liftover_bin=/mnt/d/software/bin/liftOver
reference_dir=/mnt/f/gen/1hgp gene_dir=/mnt/d/files


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
