#!/usr/bin/env bash
set -eo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd); F=$ROOT/f
[[ -s "$ROOT/gu.env" ]] && source "$ROOT/gu.env"

usage(){ cat <<'HELP'
GU v7 — archaic introgression workflow

Usage:
  ./gu.sh loci_avcf [run]
      GWAS-locus analysis using modern and archaic VCFs.
  ./gu.sh loci_asnp [run]
      GWAS-locus screening/QC using a published archaic SNP/haplotype map.

  ./gu.sh ibdmix [run|report|compare|igv|check]
  ./gu.sh arg [prepare|build|check]
  ./gu.sh trace [run|extract|infer|segments|report]
  ./gu.sh as3 [run|check|prep]

  ./gu.sh normalize
      Normalize outputs and precompute browser-ready BED/link files for R Shiny.
  ./gu.sh shiny
      Open the local GU results browser.

  ./gu.sh ukb inspect-hap|make-panel|batches|hap-vcf|hap-arg-vcf|inspect-typed

Notes:
  ACTION defaults to run for loci_avcf, loci_asnp, ibdmix, trace, and as3;
  ACTION defaults to build for arg.
  ./install.sh installs/repairs the 'gu' Conda environment, software, and references.
  ./install.sh --check checks software, environment, and references without changing them.
  Set GU_ENV_NAME only when using a differently named Conda environment.
  TRACE can reuse an ARG *backend/scaffold*, but target haplotypes must be represented
  in the final tree sequence. A finished 1KG-only ARG cannot directly call UKB samples.
HELP
}

# Help is dependency-free and must remain available before Conda is installed.
case "${1:-help}" in help|-h|--help) usage; exit 0;; esac

# Every method uses the single GU environment. Activation is local to this process;
# the caller's shell environment is unchanged when gu.sh exits.
ENV_NAME=${GU_ENV_NAME:-gu}
if [[ ${CONDA_DEFAULT_ENV:-} != "$ENV_NAME" ]]; then
  env_activated=0
  # First try the Conda already configured in this shell. If it belongs to a
  # different installation and cannot see GU, continue through the known bases.
  if command -v conda >/dev/null 2>&1; then
    if [[ $(type -t conda) != function ]]; then
      current_conda_base=$(conda info --base 2>/dev/null || true)
      [[ -n $current_conda_base && -s $current_conda_base/etc/profile.d/conda.sh ]] && source "$current_conda_base/etc/profile.d/conda.sh"
    fi
    if conda activate "$ENV_NAME" 2>/dev/null; then env_activated=1; fi
  fi
  if (( env_activated == 0 )); then
    for conda_base in "${GU_CONDA_BASE:-}" "$HOME/anaconda3" "$HOME/miniconda3" /opt/conda; do
      [[ -n $conda_base && -s $conda_base/etc/profile.d/conda.sh ]] || continue
      source "$conda_base/etc/profile.d/conda.sh"
      if conda activate "$ENV_NAME" 2>/dev/null; then env_activated=1; break; fi
    done
  fi
  if (( env_activated == 0 )); then
    echo "ERROR: Conda environment '$ENV_NAME' is not installed or not discoverable." >&2
    echo "Run: ./install.sh" >&2
    echo "Or, for an existing differently named environment: GU_ENV_NAME=<name> ./gu.sh ${1:-help}" >&2
    exit 1
  fi
fi
# Do not let a host R_LIBS path shadow Conda packages built for this R version.
unset R_LIBS R_LIBS_USER
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
unset BCFTOOLS_PLUGINS
set -u

METHOD=${1:-help}; ACTION=${2:-}
[[ $# -gt 0 ]] && shift; [[ $# -gt 0 ]] && shift || true
GU_DATA_ROOT=${GU_DATA_ROOT:-/mnt/d}; GU_REF_ROOT=${GU_REF_ROOT:-$GU_DATA_ROOT/data.BIG/refGen}; GU_ANALYSIS_ROOT=${GU_ANALYSIS_ROOT:-$GU_DATA_ROOT/analysis/gu}; GU_SOFT=${GU_SOFT:-$GU_DATA_ROOT/software/gu}
GU_BUILD=${GU_BUILD:-37}; GU_TARGET=${GU_TARGET:-1kg}; GU_TARGET_ROOT=${GU_TARGET_ROOT:-$GU_REF_ROOT/1kg/GRCH$GU_BUILD}; GU_TARGET_VCF_DIR=${GU_TARGET_VCF_DIR:-$GU_TARGET_ROOT/vcf}; GU_ARCHAIC_ROOT=${GU_ARCHAIC_ROOT:-$GU_REF_ROOT/archaic/GRCH$GU_BUILD}; GU_SAMPLE_PANEL=${GU_SAMPLE_PANEL:-$GU_TARGET_VCF_DIR/samples_v3.ALL.panel}
IBDMIX_RUNTIME=${IBDMIX_RUNTIME:-$GU_SOFT/ibdmix}; [[ -d $IBDMIX_RUNTIME ]] || IBDMIX_RUNTIME=$GU_SOFT/IBDmix
AS3_RUNTIME=${AS3_RUNTIME:-$GU_SOFT/as3}; [[ -d $AS3_RUNTIME ]] || AS3_RUNTIME=$GU_SOFT/ArchaicSeeker3.0
export GU_DATA_ROOT GU_REF_ROOT GU_ANALYSIS_ROOT GU_SOFT GU_BUILD GU_TARGET GU_TARGET_ROOT GU_TARGET_VCF_DIR GU_ARCHAIC_ROOT GU_SAMPLE_PANEL
export dir0="$GU_DATA_ROOT" dir_ref="$GU_REF_ROOT" dir_archaic="$GU_ARCHAIC_ROOT" dir_1kg="$GU_TARGET_ROOT" GRCH="$GU_BUILD"

refresh_rshiny(){
  local rshiny=${GU_RSHINY_DIR:-$GU_ANALYSIS_ROOT/Rshiny}
  mkdir -p "$rshiny"
  python3 "$F/normalize_results.py" --analysis-root "$GU_ANALYSIS_ROOT" --output-dir "$rshiny" \
    --build "GRCh$GU_BUILD" --reciprocal-overlap "${GU_CATALOG_RECIP_OVERLAP:-0.5}"
}

run_loci(){
  local m=$1 mode=${ACTION:-run}; [[ $mode == run ]] && mode=all
  local archaic_input=$GU_ARCHAIC_ROOT
  [[ $m == asnp ]] && archaic_input=$GU_ARCHAIC_ROOT/asnp
  env method="$m" dir0="$GU_DATA_ROOT" dir_ref="$GU_REF_ROOT" dir_archaic="$GU_ARCHAIC_ROOT" dirarch="$archaic_input" dirmod="$GU_TARGET_ROOT" sample_file="$GU_SAMPLE_PANEL" dirscript="$F" genome_build="b$GU_BUILD" dirout_root="$GU_ANALYSIS_ROOT/loci_${m}" bash "$F/loci.sh" "$mode"
  [[ $mode == all ]] && refresh_rshiny
}

run_ibdmix(){
  local a=${ACTION:-run} ga; case "$a" in run)ga=ibdmix_run;;report)ga=ibdmix_report;;compare)ga=ibdmix_compare;;igv)ga=ibdmix_igv;;check)ga=ibdmix_check;;*) echo "bad IBDmix action $a" >&2; exit 2;;esac
  env GU_ACTION="$ga" dir0="$GU_DATA_ROOT" dir_ref="$GU_REF_ROOT" dir_archaic="$GU_ARCHAIC_ROOT" dirarch="$GU_ARCHAIC_ROOT" dirmod="$GU_TARGET_ROOT" sample_file="$GU_SAMPLE_PANEL" dirscript="$F" dirsoft="$IBDMIX_RUNTIME" GRCH="$GU_BUILD" genome_build="b$GU_BUILD" dirout="${IBDMIX_OUT:-$GU_ANALYSIS_ROOT/ibdmix}" refs="${IBDMIX_REFS:-Altai Chagyr Vindija Denisova}" lod_cut="${IBDMIX_LOD:-4}" len_cut="${IBDMIX_MIN_BP:-50000}" job_of_chr="${IBDMIX_JOB_OF_CHR:-2}" job_in_chr="${IBDMIX_JOB_IN_CHR:-8}" chrs="${GU_CHRS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X}" bash "$F/ibdmix.sh"
  [[ $a != check ]] && refresh_rshiny
}

run_as3(){
  local a=${ACTION:-run}
  export AS3_CHRS=${AS3_CHRS:-${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"}}
  export AS3_TARGET_VCF_DIR=${AS3_TARGET_VCF_DIR:-$GU_TARGET_VCF_DIR}
  export AS3_TARGET_PANEL_FILE=${AS3_TARGET_PANEL_FILE:-$GU_SAMPLE_PANEL}
  export AS3_AFR_REFERENCE_PANEL_FILE=${AS3_AFR_REFERENCE_PANEL_FILE:-$GU_SAMPLE_PANEL}
  export AS3_DATA_OUT=${AS3_DATA_OUT:-$GU_ARCHAIC_ROOT/as3/preprocessed}
  export AS3_OUT=${AS3_OUT:-$GU_ANALYSIS_ROOT/as3}
  if [[ $a == prep ]]; then bash "$F/as3_prep.sh" --dir-as3 "$GU_ARCHAIC_ROOT" --out "$AS3_DATA_OUT" --chr "$AS3_CHRS"; return; fi
  # Preparation is explicit because AS3 requires build-matched FASTA/masks/archaic refs.
  [[ -s $AS3_DATA_OUT/manifest.tsv ]] || { echo "AS3 prepared manifest absent; running prep first."; bash "$F/as3_prep.sh" --dir-as3 "$GU_ARCHAIC_ROOT" --out "$AS3_DATA_OUT" --chr "$AS3_CHRS"; }
  env GU_ACTION="as3_${a}" AS3_DATA_IN="$AS3_DATA_OUT" AS3_RUNTIME="$AS3_RUNTIME" bash "$F/as3.sh"
}

case "$METHOD" in
  help|-h|--help) usage ;;
  loci_avcf) run_loci avcf ;;
  loci_asnp) run_loci asnp ;;
  ibdmix) run_ibdmix ;;
  arg) exec bash "$F/arg.sh" "${ACTION:-build}" "$@" ;;
  trace)
    if [[ ${ACTION:-run} == run && ${GU_ARG_BACKEND:-tsinfer} == tsinfer ]]; then
      if ! GU_CHRS="${TRACE_CHRS:-${GU_CHRS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X}}" bash "$F/arg.sh" check >/dev/null 2>&1; then
        echo "TRACE ARG missing; building full-chromosome tsinfer exploratory ARGs first." >&2
        GU_CHRS="${TRACE_CHRS:-${GU_CHRS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X}}" bash "$F/arg.sh" build
      fi
    fi
    exec bash "$F/trace.sh" "${ACTION:-run}" "$@" ;;
  as3) run_as3 ;;
  normalize)
    refresh_rshiny ;;
  shiny)
    export GU_SQLITE=${GU_SQLITE:-${GU_RSHINY_DIR:-$GU_ANALYSIS_ROOT/Rshiny}/gu.sqlite}
    exec Rscript -e "shiny::runApp('$ROOT/shiny', host=Sys.getenv('GU_SHINY_HOST','127.0.0.1'), port=as.integer(Sys.getenv('GU_SHINY_PORT','3838')), launch.browser=interactive())" ;;
  ukb) exec bash "$F/ukb.sh" "$ACTION" "$@" ;;
  *) echo "Unknown method: $METHOD" >&2; usage; exit 2 ;;
esac
