#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# This installer/checker runs directly; for a CPU installation worker, for example:
#   bsub -q 33 -n 10 -W 4:00 -M 32GB -R "span[hosts=1] rusage[mem=32GB/host]" bash install.sh
# GPU builds/checks require a suitable GPU allocation instead.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$ROOT/../0f/0os.sh"
workflow_args "$@"
set -- "${WORKFLOW_ARGS[@]}"
export PATH="$OS_GU_CONDA/bin:$DIR_SOFTWARE/bin:$PATH"
if [[ ${1:-} == -h || ${1:-} == --help ]]; then
  echo 'Usage: ./install.sh [--check [--stage render|trajectory|project|baseline|time|scan|herv|all] [--python FILE]]'
  exit 0
fi
if [[ "${1:-}" == "--check" ]]; then
  if [[ -n ${LSB_JOBID:-} ]]; then export ALLELE_BATCH_WORKER=1; fi
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
  python=$(workflow_python "$python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3 "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
  export PATH="$(dirname -- "$python"):$PATH"
  export PYTHONDONTWRITEBYTECODE=1
  exec "$python" "$ROOT/f/check_environment.py" --stage "$stage" --plink2 "${PCA3D_PLINK2:-plink2}"
fi
if [[ $submit == hpc ]]; then
  prefix=$DIR_WORK/.conda/envs/allele3-evolution
else
  prefix=$OS_GU_CONDA/envs/allele3-evolution
fi
action=create
if [[ -x $prefix/bin/python3 ]]; then
  if bash "$ROOT/install.sh" --check --stage all --python "$prefix/bin/python3"; then exit 0; fi
  action=update
fi
if command -v mamba >/dev/null 2>&1; then
  solver=mamba
elif command -v conda >/dev/null 2>&1; then
  solver=conda
else
  echo "Install dependencies listed in environment.yml; or install conda/mamba and rerun." >&2
  echo "Python fallback runs without bcftools/pysam when FASTA is uncompressed and indexed." >&2
  exit 2
fi
"$solver" env "$action" -p "$prefix" -f "$ROOT/environment.yml"
exec bash "$ROOT/install.sh" --check --stage all --python "$prefix/bin/python3"
