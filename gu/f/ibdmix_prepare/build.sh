#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
prefix=${IBDMIX_HTSLIB_PREFIX:-${CONDA_PREFIX:-/home/huangj/anaconda3/envs/gu}}
out=${IBDMIX_PREPARE_BIN:-${GU_SOFT:-/mnt/d/software/gu}/ibdmix-prepare/prepare}
compiler=${CXX:-g++}
[[ -s $prefix/include/htslib/vcf.h ]] || { echo "ERROR: HTSlib headers missing under $prefix" >&2; exit 2; }
mkdir -p "$(dirname -- "$out")"
exec 8>"$out.build.lock"
flock 8
if [[ ${1:-} == --if-needed && -x $out && $out -nt $here/prepare.cc && $out -nt $here/build.sh ]]; then exit 0; fi
tmp=$out.build.$$
trap 'rm -f "$tmp"' EXIT
"$compiler" -O3 -std=c++17 -Wall -Wextra -Werror -I"$prefix/include" \
  "$here/prepare.cc" -L"$prefix/lib" -Wl,-rpath,"$prefix/lib" -lhts -o "$tmp"
"$tmp" --version
mv -f "$tmp" "$out"
