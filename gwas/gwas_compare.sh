#!/usr/bin/env bash
# A foreground Shiny server must keep its listening URL and live diagnostics.
if [[ ${1:-} != shiny ]]; then source /mnt/d/scripts/0f/console.sh; fi
set -euo pipefail
# Defaults, CLI and module commands live here. Project discovery/build alignment
# analysis functions live in f/ and the Shiny application lives in shiny/.


# 🚩 Command-line help
usage(){ cat <<'HELP'
Usage: gwas_compare.sh compare|ldsc|shiny [options]

Shared options:
  --project-dir, --dir-gwas DIR  Discover <DIR>/<category>/<trait>/gwas/<trait>.gz
  --gwas-files A.gz,B.gz        Explicit inputs; use this OR a project directory
  --category common            Project category (default: common)
  --anchor NAME                Put this project trait first
  --grch 37|38                 Compare build (default: input build) / initial shiny build (default: 37)
  --require-grch 37|38         Require finalized inputs in this build
  --output-dir, --dir-out DIR  Default: /mnt/d/analysis/gwas/<module>
  --check-only TRUE|FALSE      Write input manifest without starting a module (FALSE)

compare options (first GWAS versus each follower):
  --labels A,B                 Optional labels matching input order
  --mplot TRUE|FALSE           Manhattan comparison (TRUE)
  --compare-beta TRUE|FALSE    Effect-size comparison (TRUE)
  --compare-EAF TRUE|FALSE     Allele-frequency comparison (TRUE)
  --p-threshold 5e-8           Significant-variant threshold
  --significant first|either|both  Significance selector (first)

ldsc options:
  --ldsc-software-dir DIR      LDSC checkout (default: /mnt/d/software/ldsc)
  --conda PATH --conda-env NAME Conda executable and environment (ldsc)
  --python PATH               Explicit LDSC interpreter; bypass conda
  --merge-alleles FILE         HM3 allele list
  --ref-ld-chr PREFIX          Reference LD scores; preserve trailing / or dot
  --w-ld-chr PREFIX            Regression weights
  --N NUMBER                  Explicit fallback when input has no N column
  --missing-n error|skip       Missing sample-size policy (error)
  --run-munge TRUE|FALSE       Prepare missing/stale sumstats caches (TRUE)
  --run-h2 TRUE|FALSE          Estimate heritability (TRUE)
  --run-rg TRUE|FALSE          Estimate genetic correlations (TRUE)
  --run TRUE|FALSE             Execute LDSC; FALSE writes commands only (TRUE)

shiny options (interactive Shiny Manhattan tracks and LD block boundaries):
  --labels CSV                Optional track labels in input order
  --races CSV                 Optional ancestry labels; default: filename suffix
  --block-dir DIR             [race].[37|38].bed files (default: /mnt/e/refLD/block)
  --ld-dir DIR                PRS-CSx HDF5 root (default: /mnt/e/refLD/csx)
  --ld-python PATH            Python with h5py/numpy/pandas (default: ~/anaconda3/bin/python)
  --ld-max-snps 300           Shared display SNP limit for the five LD panels (2–1000)
  --port 3841 --host 127.0.0.1 Local Shiny server address
  --launch-browser FALSE      Open a browser automatically (TRUE|FALSE)
  --max-points 12000           Maximum displayed points per track; zoom for detail
  --p-threshold 5e-8           Significance guide line
  --chain-dir DIR             UCSC chains (default: /mnt/d/files/liftOver)
  --liftover-bin PATH          Default: /mnt/d/software/bin/liftOver
  --reference-dir DIR          Indexed FASTA (default: /mnt/e/refGen/fasta)
  --gene-dir DIR               glist.37.bed / glist.38.bed (default: /mnt/d/files)
  shiny prefers adjacent .thin.gz files when available; compare/ldsc use full GWAS.
  All inputs must share one source build. The GRCh selector converts cached
  viewing copies together; source GWAS files are never rewritten.

Examples:
cd /mnt/d/scripts/gwas

./gwas_compare.sh compare \
  --project-dir /mnt/d/data/gwas/main --category common --anchor bald0 --require-grch 38 \
  --mplot TRUE --compare-beta TRUE --compare-EAF TRUE \
  --output-dir /mnt/d/analysis/gwas/compare

./gwas_compare.sh ldsc --dir-gwas /mnt/d/data/gwas/main --category common --anchor bald0 --require-grch 38 --missing-n skip --dir-out /mnt/d/analysis/gwas/main/ldsc

# 4grid: six formatted height GWAS, with EUR first for pairwise comparisons.
# List the files explicitly so later LDL/T2DM outputs are not included.
gwas_root=/mnt/d/data/gwas/4grid/common
trait=height
gwas_files=()
for ancestry in EUR AFR EAS HIS SAS ALL; do
  gwas_files+=("$gwas_root/$trait.$ancestry/gwas/$trait.$ancestry.gz")
done
gwas_csv=$(IFS=,; echo "${gwas_files[*]}")

./gwas_compare.sh compare \
  --gwas-files "$gwas_csv" --labels EUR,AFR,EAS,HIS,SAS,ALL \
  --require-grch 37 \
  --mplot TRUE --compare-beta TRUE --compare-EAF TRUE \
  --p-threshold 5e-8 --significant first \
  --output-dir "/mnt/d/analysis/gwas/4grid/compare/$trait"

# LDSC uses the configured LD-score reference for all six inputs; it does not
# select ancestry-specific references from the filenames.
./gwas_compare.sh ldsc \
  --gwas-files "$gwas_csv" --require-grch 37 \
  --missing-n error --run-munge TRUE --run-h2 TRUE --run-rg TRUE --run TRUE \
  --output-dir "/mnt/d/analysis/gwas/4grid/ldsc/$trait"

# R Shiny: generate .thin.gz first using gwas_format.sh thin (see its examples).
trait=height
gwas_root=/mnt/d/data/gwas/4grid/common
gwas_files=()
for race in EUR AFR EAS HIS SAS ALL; do
  gwas_files+=("$gwas_root/$trait.$race/gwas/$trait.$race.thin.gz")
done
gwas_csv=$(IFS=,; echo "${gwas_files[*]}")
./gwas_compare.sh shiny \
  --gwas-files "$gwas_csv" \
  --block-dir /mnt/e/refLD/block --ld-dir /mnt/e/refLD/csx --port 3841 \
  --output-dir "/mnt/d/analysis/gwas/4grid/shiny/$trait"
HELP
}


# 🚩 Defaults and module selection
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
helper="$here/f/gwas_compare_project.py"
# shellcheck source=f/gwas_compare.f.sh
source "$here/f/gwas_compare.f.sh"

module=${1:-}
case "$module" in
  compare|ldsc|shiny) shift ;;
  -h|--help|'') usage; exit 0 ;;
  *) echo "ERROR: unknown module: $module (expected compare, ldsc or shiny; see --help)" >&2; exit 1 ;;
esac

project_dir='' gwas_files='' category=common anchor=''
grch='' require_grch='' check_only=FALSE
output_dir="/mnt/d/analysis/gwas/$module"

# Comparison settings; these values are passed explicitly to the R helper.
labels='' mplot=TRUE compare_beta=TRUE compare_eaf=TRUE
p_threshold=5e-8 significant=first

# LDSC settings; no fallback sample size is assumed.
ldsc_software_dir=/mnt/d/software/ldsc
conda=$(type -P conda || true)
conda=${conda:-$HOME/anaconda3/bin/conda}
conda_env=ldsc ldsc_python=''
merge_alleles=/mnt/e/refLD/ldsc/hm3/w_hm3.snplist
ref_ld_chr=/mnt/e/refLD/ldsc/1000G/1000G_Phase3_ldscores/LDscore.
w_ld_chr=/mnt/e/refLD/ldsc/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
sample_size='' missing_n=error
run_munge=TRUE run_h2=TRUE run_rg=TRUE run=TRUE

# Interactive block viewer. All settings remain visible at the entry point.
block_dir=/mnt/e/refLD/block
ld_dir=/mnt/e/refLD/csx ld_python="$HOME/anaconda3/bin/python" ld_max_snps=300
port=3841 host=127.0.0.1 launch_browser=FALSE max_points=12000 races=''
chain_dir=/mnt/d/files/liftOver liftover_bin=/mnt/d/software/bin/liftOver
reference_dir=/mnt/e/refGen/fasta gene_dir=/mnt/d/files


# 🚩 Command-line arguments
while (( $# )); do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    # Preserve the --option=value form supported by the former Python entry point.
    --*=*) option=${1%%=*}; value=${1#*=}; shift; set -- "$option" "$value" "$@" ;;
  esac
  compare_need_arg_value "$1" "${2-}"
  case "$1" in
    --project-dir|--dir-gwas) project_dir=$2 ;;
    --gwas-files) gwas_files=$2 ;;
    --category) category=$2 ;;
    --anchor) anchor=$2 ;;
    --grch) grch=$2 ;;
    --require-grch) require_grch=$2 ;;
    --output-dir|--dir-out) output_dir=$2 ;;
    --check-only) check_only=$2 ;;
    *)
      case "$module:$1" in
        compare:--labels|shiny:--labels) labels=$2 ;;
        compare:--mplot) mplot=$2 ;;
        compare:--compare-beta) compare_beta=$2 ;;
        compare:--compare-EAF|compare:--compare-eaf) compare_eaf=$2 ;;
        compare:--p-threshold|shiny:--p-threshold) p_threshold=$2 ;;
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
    cmd+=(--mplot "$mplot" --compare-beta "$compare_beta" --compare-EAF "$compare_eaf"
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
