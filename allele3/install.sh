#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ ${1:-} == -h || ${1:-} == --help ]]; then
  echo 'Usage: ./install.sh [--check [--stage render|trajectory|project|baseline|time|scan|herv|all] [--python FILE]]'
  exit 0
fi
if [[ "${1:-}" == "--check" ]]; then
  source "$ROOT/f/0.run.sh"
  shift
  stage=trajectory
  python=${ALLELE_PYTHON:-}
  explicit_python=0
  while (($#)); do
    case "$1" in
      --stage) need_value "$@"; stage=$2; shift 2 ;;
      --python) need_value "$@"; python=$2; explicit_python=1; shift 2 ;;
      *) fail "Unknown dependency-check option: $1" ;;
    esac
  done
  case "$stage" in
    render|trajectory|project)
      if [[ ${PCA3D_CONFIG+x} && ! -f $PCA3D_CONFIG ]]; then fail "Explicit PCA3D_CONFIG does not exist: $PCA3D_CONFIG"; fi
      config=${PCA3D_CONFIG:-$ROOT/pca3d.env}
      [[ ! -f $config ]] || source "$config"
      ((explicit_python)) || python=${PCA3D_PYTHON:-}
      ;;
  esac
  python=$(workflow_python "$python" /home/huangj/anaconda3/envs/allele3-evolution/bin/python3 "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
  export PYTHONDONTWRITEBYTECODE=1
  exec "$python" "$ROOT/f/check_environment.py" --stage "$stage" --plink2 "${PCA3D_PLINK2:-plink2}"
fi
if command -v mamba >/dev/null 2>&1; then
  exec mamba env create -f "$ROOT/environment.yml"
elif command -v conda >/dev/null 2>&1; then
  exec conda env create -f "$ROOT/environment.yml"
else
  echo "Install dependencies listed in environment.yml; or install conda/mamba and rerun." >&2
  echo "Python fallback runs without bcftools/pysam when FASTA is uncompressed and indexed." >&2
  exit 2
fi
