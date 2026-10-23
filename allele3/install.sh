#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ "${1:-}" == "--check" ]]; then
  python3 -c 'import sys; print(sys.version); import matplotlib,pysam,pgenlib,pyliftover,openpyxl,plotly; print("Local scan, ancient PGEN, lift-over, XLSX, Plotly dependencies available")'
  command -v bcftools
  command -v samtools
  command -v plink2
  exit 0
fi
if command -v mamba >/dev/null 2>&1; then
  exec mamba env create -f "$ROOT/environment.yml"
elif command -v conda >/dev/null 2>&1; then
  exec conda env create -f "$ROOT/environment.yml"
else
  echo "Install dependencies listed in environment.yml; or install conda/mamba and rerun." >&2
  echo "Python fallback/tests run without bcftools/pysam when FASTA is uncompressed and indexed." >&2
  exit 2
fi
