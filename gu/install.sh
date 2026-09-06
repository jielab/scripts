#!/usr/bin/env bash
# Conda activate/deactivate hooks in common geospatial packages are not nounset-safe.
set -eo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=f/comm.sh
source "$ROOT/f/comm.sh"
[[ -s "$ROOT/gu.env" ]] && source "$ROOT/gu.env"
ENV_NAME=${GU_ENV_NAME:-gu}
SOFT=${GU_SOFT:-/mnt/d/software/gu}
STATE=${GU_INSTALL_STATE:-$SOFT/.gu-install.state}
YML=$ROOT/environment.yml
IBDMIX_DIR=${IBDMIX_RUNTIME:-$SOFT/ibdmix}
AS3_DIR=$ROOT/f/as3
TRACE_DIR=${TRACE_RUNTIME:-$SOFT/trace}
AS3_ENV_NAME=${AS3_ENV_NAME:-as3_mamba}
AS3_PYTHON_VERSION=${AS3_PYTHON_VERSION:-3.9}
DATA_ROOT=${GU_DATA_ROOT:-/mnt/d}
REF_ROOT=${GU_REF_ROOT:-$DATA_ROOT/data.BIG/refGen}
AS3_MODEL_DIR=${AS3_MODEL_DIR:-$REF_ROOT/archaic/38/models}
TARGET_ROOT=${GU_TARGET_ROOT:-$REF_ROOT/1kg/37}
TARGET_VCF_DIR=${GU_TARGET_VCF_DIR:-$TARGET_ROOT/vcf}
SAMPLE_PANEL=${GU_SAMPLE_PANEL:-$TARGET_ROOT/samples.txt}
ARCHAIC_ROOT=${GU_ARCHAIC_ROOT:-${GU_ARCHAIC37_ROOT:-$REF_ROOT/archaic/37/vcf}}

usage(){ cat <<'HELP'
Usage: ./install.sh [--check | --check-references | --repair-as3] [options]

  --repair-as3            Repair only the isolated AS3 environment
  --env-name NAME          General GU Conda environment [gu]
  --conda FILE             Conda executable [auto-detected]
  --software-dir DIR       GU software root [/mnt/d/software/gu]
  --ibdmix-dir DIR         IBDmix checkout/build directory
  --trace-dir DIR          TRACE checkout directory
  --data-root DIR          Data root [/mnt/d]
  --reference-root DIR     Reference root [<data-root>/data.BIG/refGen]
  --target-root DIR        GRCh37 target reference root
  --target-vcf-dir DIR     Indexed target VCF directory
  --sample-panel FILE      Target sample metadata
  --archaic-root DIR       GRCh37 archaic VCF root [<reference-root>/archaic/37/vcf]
  --as3-python-version VER AS3 environment Python version [gu.env]
  -h, --help               Show this message
HELP
}

need_value(){ [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }; }

INSTALL_ACTION=install
CONDA_INPUT=""
while (( $# )); do
  case "$1" in
    --check) INSTALL_ACTION=check; shift ;;
    --check-references) INSTALL_ACTION=check-references; shift ;;
    --repair-as3) INSTALL_ACTION=repair-as3; shift ;;
    --env-name) need_value "$@"; ENV_NAME=$2; shift 2 ;;
    --conda) need_value "$@"; CONDA_INPUT=$2; shift 2 ;;
    --software-dir) need_value "$@"; SOFT=$2; STATE=$SOFT/.gu-install.state; IBDMIX_DIR=$SOFT/ibdmix; TRACE_DIR=$SOFT/trace; shift 2 ;;
    --ibdmix-dir) need_value "$@"; IBDMIX_DIR=$2; shift 2 ;;
    --trace-dir) need_value "$@"; TRACE_DIR=$2; shift 2 ;;
    --data-root) need_value "$@"; DATA_ROOT=$2; REF_ROOT=$DATA_ROOT/data.BIG/refGen; TARGET_ROOT=$REF_ROOT/1kg/37; TARGET_VCF_DIR=$TARGET_ROOT/vcf; SAMPLE_PANEL=$TARGET_ROOT/samples.txt; ARCHAIC_ROOT=$REF_ROOT/archaic/37/vcf; shift 2 ;;
    --reference-root) need_value "$@"; REF_ROOT=$2; TARGET_ROOT=$REF_ROOT/1kg/37; TARGET_VCF_DIR=$TARGET_ROOT/vcf; SAMPLE_PANEL=$TARGET_ROOT/samples.txt; ARCHAIC_ROOT=$REF_ROOT/archaic/37/vcf; shift 2 ;;
    --target-root) need_value "$@"; TARGET_ROOT=$2; TARGET_VCF_DIR=$TARGET_ROOT/vcf; SAMPLE_PANEL=$TARGET_ROOT/samples.txt; shift 2 ;;
    --target-vcf-dir) need_value "$@"; TARGET_VCF_DIR=$2; shift 2 ;;
    --sample-panel) need_value "$@"; SAMPLE_PANEL=$2; shift 2 ;;
    --archaic-root) need_value "$@"; ARCHAIC_ROOT=$2; shift 2 ;;
    --as3-python-version) need_value "$@"; AS3_PYTHON_VERSION=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

CONDA_EXE_PATH=""
for candidate in "$CONDA_INPUT" "${GU_CONDA_EXE:-}" "${GU_CONDA_BASE:+$GU_CONDA_BASE/bin/conda}" "$HOME/anaconda3/bin/conda" "$HOME/miniconda3/bin/conda" /opt/conda/bin/conda "$(command -v conda 2>/dev/null || true)"; do
  [[ -n $candidate && -x $candidate ]] || continue
  CONDA_EXE_PATH=$candidate
  break
done
[[ -n $CONDA_EXE_PATH ]] || { echo "ERROR: conda is unavailable; install Miniconda/Anaconda or use --conda FILE" >&2; exit 1; }
CONDA_BASE_PATH=$("$CONDA_EXE_PATH" info --base 2>/dev/null) || { echo "ERROR: cannot determine Conda base from $CONDA_EXE_PATH" >&2; exit 1; }
[[ -d $CONDA_BASE_PATH ]] || { echo "ERROR: invalid Conda base from $CONDA_EXE_PATH: $CONDA_BASE_PATH" >&2; exit 1; }

die(){ echo "ERROR: $*" >&2; exit 1; }
have_env(){ "$CONDA_EXE_PATH" env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$ENV_NAME"; }
activate(){ source "$CONDA_BASE_PATH/etc/profile.d/conda.sh"; conda activate "$ENV_NAME"; }
as3_python_path(){
  local candidate
  for candidate in "${AS3_PYTHON:-}" "$CONDA_BASE_PATH/envs/$AS3_ENV_NAME/bin/python"; do
    [[ -n $candidate && -x $candidate ]] || continue
    printf '%s\n' "$candidate"
    return 0
  done
  return 1
}
as3_env_ok(){
  local py
  py=$(as3_python_path) || return 1
  AS3_EXPECT_PYTHON="$AS3_PYTHON_VERSION" "$py" "$ROOT/f/as3/gu_health.py" >/dev/null 2>&1
}
as3_health_report(){
  local py
  py=$(as3_python_path) || { echo "AS3 HEALTH ERROR: environment $AS3_ENV_NAME has no Python" >&2; return 1; }
  AS3_EXPECT_PYTHON="$AS3_PYTHON_VERSION" "$py" "$ROOT/f/as3/gu_health.py"
}
as3_repair_hdf5_abi(){
  local py=$1 built runtime
  read -r built runtime < <("$py" - <<'PY'
import warnings
warnings.filterwarnings("ignore", message="h5py is running against HDF5")
import h5py
print(".".join(map(str, h5py.version.hdf5_built_version_tuple[:2])),
      ".".join(map(str, h5py.version.hdf5_version_tuple[:2])))
PY
  ) || die "could not inspect h5py/HDF5 versions in $AS3_ENV_NAME"
  [[ -n $built && -n $runtime ]] || die "h5py did not report build/runtime HDF5 versions"
  [[ $built != "$runtime" ]] || return 0
  echo "Repairing h5py/HDF5 ABI: built=$built runtime=$runtime"
  # Match the runtime library to the ABI used by the installed h5py extension.
  # Solving the two packages together prevents a later HDF5-only upgrade from
  # recreating the mismatch.
  "$CONDA_EXE_PATH" install -n "$AS3_ENV_NAME" -c conda-forge \
    "h5py>=3.10,<4" "hdf5=$built.*" -y
}
as3_compute_cap(){
  if [[ -n ${AS3_GPU_COMPUTE_CAP:-} ]]; then printf '%s\n' "$AS3_GPU_COMPUTE_CAP"; return; fi
  nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits 2>/dev/null | awk 'NR==1{gsub(/[[:space:]]/,"");print;exit}'
}
as3_blackwell(){
  local cap major
  cap=$(as3_compute_cap); major=${cap%%.*}
  [[ $major =~ ^[0-9]+$ && $major -ge 12 ]]
}
install_as3_torch_profile(){
  local py=$1 torch_version cuda_version torchvision_version index_url actual actual_torch actual_cuda actual_torchvision
  if as3_blackwell; then
    torch_version=${AS3_TORCH_VERSION:-2.8.0}; cuda_version=${AS3_CUDA_VERSION:-12.8}; torchvision_version=${AS3_TORCHVISION_VERSION:-0.23.0}
  else
    torch_version=${AS3_TORCH_VERSION:-2.4.0}; cuda_version=${AS3_CUDA_VERSION:-12.1}; torchvision_version=${AS3_TORCHVISION_VERSION:-0.19.0}
  fi
  actual=$("$py" -c 'import torch; print(torch.__version__.split("+",1)[0]+"|"+str(torch.version.cuda or ""))' 2>/dev/null || true)
  actual_torch=${actual%%|*}; actual_cuda=${actual#*|}
  if [[ $actual_torch != "$torch_version" || $actual_cuda != "$cuda_version" ]]; then
    index_url=${AS3_TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu${cuda_version/.}}
    "$py" -m pip install "torch==$torch_version" "torchvision==$torchvision_version" --index-url "$index_url"
    return
  fi
  # torchvision is not part of the AS3 requirements file because it must match
  # the hardware-specific Torch profile.  Check it independently: a valid
  # Torch install must not make a missing or mismatched torchvision look healthy.
  actual_torchvision=$("$py" -c 'import torchvision; print(torchvision.__version__.split("+",1)[0])' 2>/dev/null || true)
  if [[ $actual_torchvision != "$torchvision_version" ]]; then
    index_url=${AS3_TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu${cuda_version/.}}
    "$py" -m pip install "torchvision==$torchvision_version" --index-url "$index_url"
  fi
}
as3_runtime_ok(){
  [[ -s $AS3_DIR/ArchaicSeeker3.1-mamba ]] || return 1
  [[ -s $AS3_DIR/src/standard_postprocess.py ]] || return 1
  [[ -s $AS3_DIR/requirements-gu.txt ]] || return 1
  [[ -s $AS3_DIR/MODEL_SHA256SUMS && -d $AS3_MODEL_DIR ]] || return 1
  grep -Fq 'ancestry = int(ancestry[0])' "$AS3_DIR/src/dataloaders/genome_datasets.py" || return 1
  grep -Fq "final_bed_df_raw['sample_hap_id'] = final_bed_df_raw['sample_hap_id'].fillna('Unknown_HapID')" "$AS3_DIR/src/stepsagnostic/traintest.py" || return 1
  (cd "$AS3_MODEL_DIR" && sha256sum -c "$AS3_DIR/MODEL_SHA256SUMS" >/dev/null 2>&1)
}

validate_bundled_as3_source(){
  local data_file="$AS3_DIR/src/dataloaders/genome_datasets.py"
  local infer_file="$AS3_DIR/src/stepsagnostic/traintest.py"
  [[ -s $AS3_DIR/ArchaicSeeker3.1-mamba && -s $AS3_DIR/src/standard_postprocess.py && -s $data_file && -s $infer_file ]] || \
    die "bundled AS3 source is incomplete under $AS3_DIR"
  grep -Fq 'ancestry = int(ancestry[0])' "$data_file" || \
    die "bundled AS3 NumPy compatibility edit is missing: $data_file"
  grep -Fq "final_bed_df_raw['sample_hap_id'] = final_bed_df_raw['sample_hap_id'].fillna('Unknown_HapID')" "$infer_file" || \
    die "bundled AS3 Pandas raw-output compatibility edit is missing: $infer_file"
  [[ -s $AS3_DIR/MODEL_SHA256SUMS && -d $AS3_MODEL_DIR ]] || \
    die "AS3 model directory or checksum manifest is missing: $AS3_MODEL_DIR"
  (cd "$AS3_MODEL_DIR" && sha256sum -c "$AS3_DIR/MODEL_SHA256SUMS" >/dev/null) || \
    die "AS3 model files are missing or fail SHA-256 validation under $AS3_MODEL_DIR"
}

repair_as3_env(){
  local as3_env_python mamba_version causal_version
  validate_bundled_as3_source
  [[ -s $AS3_DIR/requirements-gu.txt ]] || die "missing GU AS3 dependency file: $AS3_DIR/requirements-gu.txt"
  if as3_env_python=$(as3_python_path 2>/dev/null); then
    "$as3_env_python" - "$AS3_PYTHON_VERSION" <<'PY' || \
      die "existing $AS3_ENV_NAME uses $($as3_env_python -c 'import platform; print(platform.python_version())'); recreate it with Python $AS3_PYTHON_VERSION before repairing it"
import sys
expected = tuple(int(x) for x in sys.argv[1].split('.')[:2])
raise SystemExit(0 if sys.version_info[:2] == expected else 1)
PY
  else
    "$CONDA_EXE_PATH" create -n "$AS3_ENV_NAME" "python=$AS3_PYTHON_VERSION" pip git packaging ninja -c conda-forge -y
    as3_env_python=$(as3_python_path) || die "could not locate Python for newly created AS3 environment: $AS3_ENV_NAME"
  fi
  install_as3_torch_profile "$as3_env_python"
  "$as3_env_python" -m pip install -U pip setuptools wheel packaging ninja
  "$as3_env_python" -m pip install -r "$AS3_DIR/requirements-gu.txt"
  as3_repair_hdf5_abi "$as3_env_python"
  causal_version=${AS3_CAUSAL_CONV1D_VERSION:-1.6.2.post1}
  mamba_version=${AS3_MAMBA_SSM_VERSION:-2.3.1}
  if ! "$as3_env_python" -c 'import causal_conv1d, mamba_ssm' >/dev/null 2>&1; then
    CUDA_HOME=$CONDA_BASE_PATH/envs/$AS3_ENV_NAME \
      "$as3_env_python" -m pip install --no-build-isolation "causal-conv1d==$causal_version" "mamba-ssm==$mamba_version"
  fi
  as3_health_report || die "AS3 repair finished but the environment is still unhealthy"
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
             'tsinfer','tsdate','bio2zarr'):
    importlib.import_module(name)
PY
  for cmd in trace-extract trace-infer trace-summarize; do command -v "$cmd" >/dev/null 2>&1 || { echo "MISSING command: $cmd" >&2; bad=1; }; done
  for cmd in data.table ape pegas writexl shiny bslib plotly DT DBI RSQLite; do
    Rscript -e "quit(status=if(requireNamespace('$cmd', quietly=TRUE)) 0 else 1)" >/dev/null 2>&1 || { echo "MISSING or unloadable R package: $cmd" >&2; bad=1; }
  done
  for cmd in ibdmix generate_gt gt_lods; do
    [[ -x "$IBDMIX_DIR/build/src/$cmd" ]] || { echo "MISSING IBDmix executable: $IBDMIX_DIR/build/src/$cmd" >&2; bad=1; }
  done
  as3_runtime_ok || { echo "MISSING AS3 runtime: $AS3_DIR" >&2; bad=1; }
  as3_health_report || {
    echo "MISSING or broken AS3 environment: $AS3_ENV_NAME" >&2; bad=1;
  }
  (( bad == 0 ))
}

find_archaic_vcf(){
  local ref=$1 chr=$2 d found
  for d in "$ARCHAIC_ROOT/$ref" "$ARCHAIC_ROOT/avcf/$ref"; do
    [[ -d $d ]] || continue
    found=$(find "$d" -maxdepth 1 -type f \( -name "*chr${chr}_*.vcf.gz" -o -name "*chr${chr}.*.vcf.gz" -o -name "*chr${chr}.vcf.gz" \) -print -quit 2>/dev/null)
    [[ -n $found ]] && { printf '%s\n' "$found"; return 0; }
  done
  return 1
}

check_references(){
  local bad=0 chr ref vcf vcf_base pbase
  [[ -s $SAMPLE_PANEL ]] || { echo "MISSING 1KG sample panel: $SAMPLE_PANEL" >&2; bad=1; }
  if [[ -s $SAMPLE_PANEL ]]; then
    awk 'BEGIN{FS="[[:space:]]+"}NR==1{for(i=1;i<=NF;i++)h[tolower($i)]=1;exit !(h["sample"]&&h["sex"])}' "$SAMPLE_PANEL" || { echo "INVALID samples.txt: sample and sex columns are required" >&2; bad=1; }
  fi
  for chr in {1..22} X; do
    vcf_base=$(gu_target_genotype_base "$TARGET_VCF_DIR/chr" "$chr")
    pbase=$(gu_target_genotype_base "$TARGET_ROOT/pfile/chr" "$chr")
    vcf=$vcf_base.vcf.gz
    if [[ ! -s $vcf || (! -s $vcf.tbi && ! -s $vcf.csi) ]]; then
      [[ -s $pbase.pgen && ( -s $pbase.pvar || -s $pbase.pvar.zst ) && -s $pbase.psam ]] || {
        echo "MISSING target chr$chr: expected indexed $vcf or complete pfile $pbase" >&2; bad=1;
      }
    fi
  done
  for ref in Altai Chagyr Vindija Denisova; do
    for chr in {1..22} X; do
      vcf=$(find_archaic_vcf "$ref" "$chr" || true)
      [[ -n $vcf && -s $vcf && ( -s $vcf.tbi || -s $vcf.csi ) ]] || { echo "MISSING indexed archaic VCF: $ref chr$chr under $ARCHAIC_ROOT" >&2; bad=1; }
    done
  done
  (( bad == 0 ))
}

mkdir -p "$SOFT"

if [[ $INSTALL_ACTION == check ]]; then
  check && echo OK
  exit
fi
if [[ $INSTALL_ACTION == check-references ]]; then
  check_references && echo OK
  exit
fi
if [[ $INSTALL_ACTION == repair-as3 ]]; then
  validate_bundled_as3_source
  if as3_env_ok; then
    as3_health_report
  else
    repair_as3_env
  fi
  echo "OK (AS3 environment)"
  exit 0
fi

# A successful installation marker makes subsequent calls a fast health check.
if [[ -s $STATE ]] && have_env && check; then
  echo OK
  exit 0
fi

if have_env; then
  "$CONDA_EXE_PATH" env update -n "$ENV_NAME" -f "$YML" --prune
else
  "$CONDA_EXE_PATH" env create -n "$ENV_NAME" -f "$YML"
fi
activate

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

# AS3 inference source is vendored under f/as3; checkpoints live with
# the GRCh38 Ref1028 data.  No external AS3 checkout is read or updated.
validate_bundled_as3_source
if ! as3_env_ok; then
  repair_as3_env
  activate
fi

check || die "installation finished but the software health check still reports missing components"
printf 'installed\n' > "$STATE"
echo "OK (software). Run ./install.sh --check-references separately for the full reference inventory."
