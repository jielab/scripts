#!/usr/bin/env bash
# Conda activate/deactivate hooks in common geospatial packages are not nounset-safe.
set -eo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ENV_NAME=${GU_ENV_NAME:-gu}
SOFT=${GU_SOFT:-/mnt/d/software/gu}
STATE=$ROOT/.install.state
YML=$ROOT/environment.yml
IBDMIX_DIR=${IBDMIX_RUNTIME:-$SOFT/ibdmix}
AS3_DIR=${AS3_RUNTIME:-$SOFT/as3}
TRACE_DIR=${TRACE_RUNTIME:-$SOFT/trace}
DATA_ROOT=${GU_DATA_ROOT:-/mnt/d}
REF_ROOT=${GU_REF_ROOT:-$DATA_ROOT/data.BIG/refGen}
TARGET_ROOT=${GU_TARGET_ROOT:-$REF_ROOT/1kg/GRCH37}
TARGET_VCF_DIR=${GU_TARGET_VCF_DIR:-$TARGET_ROOT/vcf}
SAMPLE_PANEL=${GU_SAMPLE_PANEL:-$TARGET_VCF_DIR/samples_v3.ALL.panel}
ARCHAIC_ROOT=${GU_ARCHAIC_ROOT:-$REF_ROOT/archaic/GRCH37}

if ! command -v conda >/dev/null 2>&1; then
  for conda_base in "${GU_CONDA_BASE:-}" "$HOME/anaconda3" "$HOME/miniconda3" /opt/conda; do
    [[ -n $conda_base && -s $conda_base/etc/profile.d/conda.sh ]] || continue
    source "$conda_base/etc/profile.d/conda.sh"
    break
  done
fi
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda is unavailable; install Miniconda/Anaconda or set GU_CONDA_BASE" >&2; exit 1; }

die(){ echo "ERROR: $*" >&2; exit 1; }
have_env(){ conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$ENV_NAME"; }
activate(){ source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$ENV_NAME"; }
yml_hash(){ sha256sum "$YML" | awk '{print $1}'; }
as3_runtime_ok(){
  local exe
  for exe in ArchaicSeeker3.1-mamba ArchaicSeeker3.0-mamba ArchaicSeeker3-mamba; do
    [[ -s "$AS3_DIR/$exe" ]] && return 0
  done
  return 1
}

check(){
  local bad=0 cmd plugin_help
  [[ ${CONDA_DEFAULT_ENV:-} == "$ENV_NAME" ]] || activate
  # Host-level R libraries can contain packages compiled for a different R ABI.
  # GU must use the packages installed in its Conda environment.
  unset R_LIBS R_LIBS_USER
  export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
  unset BCFTOOLS_PLUGINS
  for cmd in python python3 Rscript bcftools bgzip tabix bedtools samtools plink2 bgenix phyml cmake; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "MISSING command: $cmd" >&2; bad=1; }
  done
  plugin_help=$(bcftools plugin -l 2>/dev/null || true)
  grep -Fx fill-tags <<< "$plugin_help" >/dev/null || { echo "MISSING bcftools plugin: fill-tags" >&2; bad=1; }
  python - <<'PY' || bad=1
import importlib
for name in ('numpy','pandas','pysam','cyvcf2','pybedtools','zarr','tskit','tszip',
             'tsinfer','tsdate','bio2zarr','torch','mamba_ssm','causal_conv1d'):
    importlib.import_module(name)
import importlib.metadata, torch
assert importlib.metadata.version('torch') == torch.__version__.split('+')[0], (
    'torch package metadata does not match the imported runtime')
PY
  for cmd in trace-extract trace-infer trace-summarize; do command -v "$cmd" >/dev/null 2>&1 || { echo "MISSING command: $cmd" >&2; bad=1; }; done
  for cmd in data.table ape pegas writexl shiny bslib plotly DT DBI RSQLite; do
    Rscript -e "quit(status=if(requireNamespace('$cmd', quietly=TRUE)) 0 else 1)" >/dev/null 2>&1 || { echo "MISSING or unloadable R package: $cmd" >&2; bad=1; }
  done
  for cmd in ibdmix generate_gt gt_lods; do
    [[ -x "$IBDMIX_DIR/build/src/$cmd" ]] || { echo "MISSING IBDmix executable: $IBDMIX_DIR/build/src/$cmd" >&2; bad=1; }
  done
  as3_runtime_ok || { echo "MISSING AS3 runtime: $AS3_DIR" >&2; bad=1; }
  (( bad == 0 ))
}

find_archaic_vcf(){
  local ref=$1 chr=$2
  find "$ARCHAIC_ROOT/$ref" -maxdepth 1 -type f \( -name "*chr${chr}_*.vcf.gz" -o -name "*chr${chr}.*.vcf.gz" -o -name "*chr${chr}.vcf.gz" \) -print -quit 2>/dev/null
}

check_references(){
  local bad=0 chr ref vcf ext
  [[ -s $SAMPLE_PANEL ]] || { echo "MISSING 1KG sample panel: $SAMPLE_PANEL" >&2; bad=1; }
  for chr in {1..22} X; do
    vcf=$TARGET_VCF_DIR/chr${chr}.vcf.gz
    [[ -s $vcf && -s $vcf.tbi ]] || { echo "MISSING indexed 1KG GRCh37 VCF: $vcf" >&2; bad=1; }
    for ext in pgen pvar psam; do
      [[ -s $TARGET_ROOT/pfile/chr${chr}.$ext ]] || { echo "MISSING 1KG PLINK file: $TARGET_ROOT/pfile/chr${chr}.$ext" >&2; bad=1; }
    done
  done
  for ref in Altai Chagyr Vindija Denisova; do
    for chr in {1..22} X; do
      vcf=$(find_archaic_vcf "$ref" "$chr")
      [[ -n $vcf && -s $vcf && -s $vcf.tbi ]] || { echo "MISSING indexed archaic GRCh37 VCF: $ref chr$chr under $ARCHAIC_ROOT/$ref" >&2; bad=1; }
    done
  done
  (( bad == 0 ))
}

prepare_references(){
  mkdir -p "$TARGET_VCF_DIR"
  if [[ ! -s $SAMPLE_PANEL ]]; then
    wget -O "$SAMPLE_PANEL" \
      https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
  fi
  env MODERN_ROOT="$TARGET_ROOT" ARCHAIC_ROOT="$ARCHAIC_ROOT" sample_file="$SAMPLE_PANEL" GRCH=37 \
    bash "$ROOT/f/gen.clean.sh"
}

command -v conda >/dev/null 2>&1 || die "conda is not installed or not on PATH"
mkdir -p "$SOFT"

if [[ ${1:-} == --check ]]; then
  check && check_references && echo OK
  exit
fi
[[ $# == 0 ]] || die "unknown option: $1"

# A matching successful installation makes subsequent calls a fast health check.
if [[ -s $STATE ]] && grep -Fxq "$(yml_hash)" "$STATE" && have_env && check && check_references; then
  echo OK
  exit 0
fi

if have_env; then
  conda env update -n "$ENV_NAME" -f "$YML" --prune
else
  conda env create -n "$ENV_NAME" -f "$YML"
fi
activate

# Install PyTorch separately. A failed/aborted pip transaction can leave torch's
# metadata present while files are missing; --ignore-installed repairs that state
# without asking pip to uninstall the broken copy first.
if ! python -c 'import torch; assert torch.__version__.split("+")[0] == "2.8.0"' >/dev/null 2>&1; then
  python -m pip install --ignore-installed "torch==2.8.0"
fi
python -m pip install "torchmetrics"

# TRACE has no PyTorch dependency. Installing it separately keeps its source build
# independent of the large PyTorch wheel transaction and reuses the local checkout.
if [[ ! -d $TRACE_DIR/.git ]]; then
  git clone https://github.com/YulinZhang9806/trace.git "$TRACE_DIR"
fi
if ! command -v trace-extract >/dev/null 2>&1 || ! command -v trace-infer >/dev/null 2>&1 || ! command -v trace-summarize >/dev/null 2>&1; then
  python -m pip install --no-build-isolation "$TRACE_DIR"
fi

if [[ ! -d $IBDMIX_DIR/.git ]]; then
  git clone --depth 1 https://github.com/PrincetonUniversity/IBDmix.git "$IBDMIX_DIR"
fi
if [[ -x "$IBDMIX_DIR/build/src/ibdmix" && -x "$IBDMIX_DIR/build/src/generate_gt" && -x "$IBDMIX_DIR/build/src/gt_lods" ]]; then
  echo "IBDmix: reuse existing build under $IBDMIX_DIR/build"
else
  # CMake caches absolute paths. --fresh safely regenerates a cache when an
  # existing checkout/build was moved from another directory.
  cmake --fresh -S "$IBDMIX_DIR" -B "$IBDMIX_DIR/build"
  cmake --build "$IBDMIX_DIR/build" --parallel "${GU_BUILD_JOBS:-4}"
fi

if [[ ! -d $AS3_DIR/.git ]]; then
  git clone https://github.com/Shuhua-Group/ArchaicSeeker3.0.git "$AS3_DIR"
fi
# AS3's own installer can use the already-created GU environment. It installs the
# compiled causal-conv1d/mamba-ssm pair without creating a second environment.
if ! python -c 'import causal_conv1d,mamba_ssm' >/dev/null 2>&1; then
  # Compatible AS3 runtime packages are already declared in environment.yml;
  # do not let the upstream legacy requirements file downgrade numpy/pandas.
  AS3_EMPTY_REQUIREMENTS=$(mktemp "${TMPDIR:-/tmp}/gu-as3-requirements.XXXXXX")
  (cd "$AS3_DIR" && AS3_ENV_NAME="$ENV_NAME" AS3_REQUIREMENTS="$AS3_EMPTY_REQUIREMENTS" ./install.sh)
  rm -f "$AS3_EMPTY_REQUIREMENTS"
  activate
fi

check || die "installation finished but the software health check still reports missing components"
if ! check_references; then
  echo "GRCh37 reference set is incomplete; preparing 1000 Genomes and archaic references."
  prepare_references
  check_references || die "reference preparation finished but required GRCh37 files are still incomplete"
fi
yml_hash > "$STATE"
echo OK
