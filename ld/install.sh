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
workflow_args "$@"; set -- "${WORKFLOW_ARGS[@]}"
case ${1:---check} in
  -h|--help) echo 'Usage: bash install.sh [--check|--install] [--submit hpc|local]'; exit 0 ;;
  --check) ;;
  --install)
    "$OS_GU_CONDA/bin/conda" env update -p "$DIR_WORK/.conda/envs/py311" -f "$ROOT/environment.yml" ;;
  *) echo "Unknown option: $1" >&2; exit 2 ;;
esac
source "$ROOT/f/env.sh"
"$PYTHON" -c 'import sys,numpy,scipy,pandas,yaml,zstandard; assert sys.version_info >= (3,11); print("LD Python dependencies PASS")'
"$PLINK2" --version
"$PLINK2" --help r2-phased >/dev/null
scratch=$(mktemp -d)
trap 'rm -rf -- "$scratch"' EXIT
printf '#include <zstd.h>\n#include <filesystem>\nint main(){return ZSTD_versionNumber() < 10000;}\n' > "$scratch/check.cpp"
"$CXX" -std=c++17 "$scratch/check.cpp" -lzstd -o "$scratch/check"
"$scratch/check"
"$PYTHON" -m unittest discover -s "$ROOT/2test" -q
echo 'LD software check PASS (no genotype data needed)'
