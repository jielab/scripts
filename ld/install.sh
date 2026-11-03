#!/usr/bin/env bash
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
