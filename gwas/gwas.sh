#!/usr/bin/env bash
# Independent modules; see --help for full usage and examples.
# Once: bash gwas.sh prep_gwas --sparse-grm FALSE --threads 16
# Five GWAS: bash gwas.sh run_regenie --grch 37 --threads 16
# Phenotype preparation: bash /mnt/d/scripts/ukb/phe.sh --steps phe4gwas
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$here/f/gwas_runner.py" "$@"
