# Shared FLT environment. Execute through flt.sh or flt.paper.sh.
entry=$(cd -- "$(dirname -- "$0")" && pwd -P)/$(basename -- "$0")
project_dir=$(dirname -- "$entry")
export FLT_CODE="$project_dir/f"
export FLT_SHARED="${FLT_SHARED:-$(dirname -- "$project_dir")/0f}"
if [[ -n ${FLT_PYTHON:-} ]]; then export PYTHON_BIN="${PYTHON_BIN:-$FLT_PYTHON}"; fi
source "$FLT_SHARED/0os.sh"
workflow_args
HELPERS=$FLT_SHARED
source "$project_dir/../diet/f/runtime_env.sh" "$@"
export FLT_RSCRIPT="${FLT_RSCRIPT:-$DIET_RSCRIPT}" FLT_PYTHON="${FLT_PYTHON:-$DIET_PYTHON}"
export FLT_R_LIBRARY="${FLT_R_LIBRARY:-$DIR_SOFTWARE/flt/R-library}"
export FLT_BWMR="${FLT_BWMR:-$DIR_SOFTWARE/flt/sources/BWMR_updated.R}"
export FLT_RESOURCE_DIR="${FLT_RESOURCE_DIR:-$DIR_SOFTWARE/flt/resources}"
# User-authorized policy: unavailable optional reference resources are reported and skipped.
export FLT_SKIP_MISSING_RESOURCES="${FLT_SKIP_MISSING_RESOURCES:-TRUE}"
export FLT_PROJECT="${FLT_PROJECT:-${FLT_OUT:-$ANALYSIS_ROOT/flt}}" FLT_YGB="${FLT_YGB:-$DIR_YGB}" FLT_GWAS="${FLT_GWAS:-$DIR_GWAS}"
export FLT_DICTIONARY="${FLT_DICTIONARY:-$DIR_FILES/diet/flt.xlsx}"
export FLT_MRI_DICTIONARY="${FLT_MRI_DICTIONARY:-$DIR_FILES/diet/MRI_analysis_dictionary.csv}"
export FLT_MRI_TOOLS="${FLT_MRI_TOOLS:-${DOSE_TOOLS_DIR:-$DIR_SOFTWARE/diet/.tools}}"
default_atlas="$ANALYSIS_ROOT/diet/dose/.tools/mri_atlases"
if [[ ! -d $default_atlas ]]; then default_atlas="${DOSE_RUNTIME_DIR:-$HOME/.local/share/dose}/mri_atlases"; fi
export FLT_MRI_ATLAS="${FLT_MRI_ATLAS:-${DOSE_MRI_ATLAS_DIR:-$default_atlas}}"
if [[ -z ${FLT_LDREF:-} ]]; then
    FLT_LDREF=$DIR_LD
    if [[ ! -d $FLT_LDREF && -d $DIR_DATA/ldref ]]; then FLT_LDREF=$DIR_DATA/ldref; fi
fi
export FLT_LDREF
export FLT_MR_LD_PREFIX="${FLT_MR_LD_PREFIX:-${MR_LD_PREFIX:-$DIR_GEN/1kg/37/openGWAS/EUR}}"
export FLT_PLINK="${FLT_PLINK:-${PLINK_BIN:-$(command -v plink || echo "$DIR_SOFTWARE/bin/plink")}}"
export FLT_SMR="${FLT_SMR:-$(command -v smr || echo "$DIR_SOFTWARE/bin/smr")}"
if [[ ! -x "$FLT_SMR" && -x "$FLT_RESOURCE_DIR/smr/smr" ]]; then
    FLT_SMR="$FLT_RESOURCE_DIR/smr/smr"
fi
if [[ ! -x "$FLT_SMR" && $FLT_SMR == "$DIR_SOFTWARE/bin/smr" && -x "$DIR_SOFTWARE/linux/smr/smr-1.3.1-linux-x86_64/smr" ]]; then
    FLT_SMR="$DIR_SOFTWARE/linux/smr/smr-1.3.1-linux-x86_64/smr"
fi
if [[ ! -x "$FLT_PLINK" && $FLT_PLINK == "$DIR_SOFTWARE/bin/plink" && -x "$DIR_SOFTWARE/linux/plink/plink" ]]; then
    FLT_PLINK="$DIR_SOFTWARE/linux/plink/plink"
fi
export FLT_EQTLGEN="${FLT_EQTLGEN:-$FLT_GWAS/eqtl/eqtlgen_v1/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense}"
export FLT_GTEX="${FLT_GTEX:-$FLT_GWAS/eqtl/gtex_v8_whole_blood/Whole_Blood/Whole_Blood}"
export MR_LD_PREFIX="$FLT_MR_LD_PREFIX"
export R_DATATABLE_NUM_THREADS="$WORKFLOW_CORES" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export FLT_EXTERNAL_ROOT="${FLT_EXTERNAL_ROOT:-$DIR_FILES/flt}"
export FLT_EXTERNAL_PLACO="${FLT_EXTERNAL_PLACO:-$FLT_EXTERNAL_ROOT/03placo}"
export FLT_PHEWAS_ONLINE="${FLT_PHEWAS_ONLINE:-false}"
export FLT_PHEWAS_INPUT="${FLT_PHEWAS_INPUT:-$FLT_EXTERNAL_PLACO/phewas_cache}"
export FLT_FUMA_DIR="${FLT_FUMA_DIR:-$FLT_EXTERNAL_PLACO/FUMA_job787285}"
export FLT_ATLAS_DIR="${FLT_ATLAS_DIR:-$FLT_EXTERNAL_PLACO/GWAS_Atlas}"
export FLT_METASCAPE_RESULTS="${FLT_METASCAPE_RESULTS:-$FLT_EXTERNAL_PLACO/metascape/Enrichment_GO/_FINAL_GO.csv}"
if ! command -v "$FLT_PYTHON" >/dev/null 2>&1; then
    export FLT_PYTHON; FLT_PYTHON=$(command -v python3)
fi
exec "$FLT_PYTHON" "$FLT_CODE/workflow.py" "$(basename -- "$entry")" "$@"
