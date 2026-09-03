#!/usr/bin/env bash
set -euo pipefail

# EMS120 analysis and figure pipeline.

START_STEP=all
END_STEP=""
CONDA_ENV_NAME=ai
EMS120_RUN_DX_CV=1
EMS120_CHECK_GPU_ONLY=0
RSCRIPT_BIN=""
EMS120_ANALYSIS_DIR=/mnt/d/analysis/ems120

usage() {
  cat <<'EOF'
Usage: ./ems120.sh [options]

Options:
  --start-step STEP   all|fig1|fig2|fig3|fig4|fig5|fig6|figs  (default: all)
  --end-step STEP     stop after fig1|fig2|fig3|fig4|fig5|fig6|figs
  --run-cv TRUE|FALSE run MacBERT 5-fold cross-validation [TRUE]
  --check-gpu         verify CUDA inside R/reticulate, then exit
  --conda-env NAME    conda environment containing torch/transformers (default: ai)
  --analysis-dir DIR  analysis output directory (default: /mnt/d/analysis/ems120)
  --rscript FILE      Rscript executable
  -h, --help          show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-step)
      START_STEP=${2:?missing value for --start-step}; shift 2 ;;
    --end-step)
      END_STEP=${2:?missing value for --end-step}; shift 2 ;;
    --run-cv)
      case "${2,,}" in true|1|yes) EMS120_RUN_DX_CV=1;; false|0|no) EMS120_RUN_DX_CV=0;; *) echo "ERROR: --run-cv requires TRUE or FALSE" >&2; exit 2;; esac
      shift 2
      ;;
    --check-gpu)
      EMS120_CHECK_GPU_ONLY=1; shift ;;
    --conda-env)
      CONDA_ENV_NAME=${2:?missing value for --conda-env}; shift 2 ;;
    --analysis-dir)
      EMS120_ANALYSIS_DIR=${2:?missing value for --analysis-dir}; shift 2 ;;
    --rscript)
      RSCRIPT_BIN=${2:?missing value for --rscript}; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2 ;;
  esac
done

export START_STEP END_STEP CONDA_ENV_NAME EMS120_RUN_DX_CV

export EMS120_ANALYSIS_DIR

find_conda() {
  local conda_bin
  for conda_bin in \
    "${CONDA_EXE:-}" \
    "${HOME}/anaconda3/bin/conda" \
    "${HOME}/miniconda3/bin/conda" \
    "/home/huangj/anaconda3/bin/conda" \
    "/home/jiehuang001/anaconda3/bin/conda" \
    "/opt/conda/bin/conda"
  do
    if [[ -n "${conda_bin}" && -x "${conda_bin}" ]]; then
      printf '%s\n' "${conda_bin}"
      return 0
    fi
  done
  command -v conda
}

CONDA_BIN=$(find_conda) || {
  echo "Error: cannot find conda; needed for Python env '${CONDA_ENV_NAME}'." >&2
  exit 1
}
CONDA_BASE=$("${CONDA_BIN}" info --base)

# Some conda hooks reference unset variables internally; temporarily relax nounset.
set +u
if [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
  # shellcheck source=/dev/null
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
else
  eval "$("${CONDA_BIN}" shell.bash hook)"
fi
conda activate "${CONDA_ENV_NAME}"
set -u

if [[ -z "${RSCRIPT_BIN}" ]]; then
  RSCRIPT_BIN=$(command -v Rscript || true)
fi
if [[ -z "${RSCRIPT_BIN}" ]]; then
  echo "Error: cannot find Rscript." >&2
  exit 1
fi

export RETICULATE_PYTHON="${CONDA_PREFIX}/bin/python"
export PYTHONNOUSERSITE=1
unset PYTHONHOME
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
if [[ -f "${CONDA_PREFIX}/lib/libcrypto.so.3" && -f "${CONDA_PREFIX}/lib/libssl.so.3" ]]; then
  export LD_PRELOAD="${CONDA_PREFIX}/lib/libcrypto.so.3:${CONDA_PREFIX}/lib/libssl.so.3${LD_PRELOAD:+:${LD_PRELOAD}}"
fi
# Embedded Python (reticulate) can otherwise resolve Conda's CUDA stubs before
# the real WSL driver libraries.  Preload the WSL libraries explicitly so CUDA
# works inside R as well as in a standalone Python process.
if [[ -f /usr/lib/wsl/lib/libcuda.so.1 && -f /usr/lib/wsl/lib/libnvidia-ml.so.1 ]]; then
  export LD_PRELOAD="/usr/lib/wsl/lib/libcuda.so.1:/usr/lib/wsl/lib/libnvidia-ml.so.1${LD_PRELOAD:+:${LD_PRELOAD}}"
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export EMS120_SCRIPT_DIR="${SCRIPT_DIR}"
R_MAIN="${SCRIPT_DIR}/f/ems120.R"
PY_MAIN="${SCRIPT_DIR}/f/ems120.py"
R_HELPER="${SCRIPT_DIR}/f/ems120.f.R"
for f in "${R_MAIN}" "${PY_MAIN}" "${R_HELPER}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Error: required pipeline file is missing: ${f}" >&2
    exit 1
  fi
done

# Fail early on Python syntax errors before the expensive analysis starts.
"${RETICULATE_PYTHON}" -m py_compile "${PY_MAIN}"

# Test CUDA through R/reticulate, which is the environment used by the actual
# pipeline.  A standalone Python CUDA check is insufficient for this workflow.
GPU_NAME=$("${RSCRIPT_BIN}" -e '
  suppressPackageStartupMessages(library(reticulate))
  use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
  py_run_string("import torch")
  if (!isTRUE(py_eval("torch.cuda.is_available()", convert = TRUE))) {
    stop("CUDA is unavailable inside R/reticulate; refusing CPU fallback.", call. = FALSE)
  }
  cat(py_eval("torch.cuda.get_device_name(0)", convert = TRUE))
')

echo "EMS120 pipeline"
echo "START_STEP=${START_STEP}"
[[ -n "${END_STEP}" ]] && echo "END_STEP=${END_STEP}"
echo "CONDA_ENV=${CONDA_PREFIX}"
echo "RETICULATE_PYTHON=${RETICULATE_PYTHON}"
echo "CUDA_DEVICE=${GPU_NAME}"
echo "ANALYSIS_DIR=${EMS120_ANALYSIS_DIR}"
echo "EMS120_RUN_DX_CV=${EMS120_RUN_DX_CV}"

if [[ "${EMS120_CHECK_GPU_ONLY}" == "1" ]]; then
  echo "GPU_CHECK=passed"
  exit 0
fi

a=(--start_step "${START_STEP}")
if [[ -n "${END_STEP}" ]]; then
  a+=(--end_step "${END_STEP}")
fi
"${RSCRIPT_BIN}" "${R_MAIN}" "${a[@]}"
