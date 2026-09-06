#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
case "${1:-}" in
  compare) shift; exec Rscript "$here/f/gwas_compare.R" "$@" ;;
  ldsc) shift; exec python3 "$here/f/gwas_ldsc.py" "$@" ;;
  -h|--help|"")
    printf '%s\n' 'Usage: gwas_post.sh compare|ldsc [options]' 'Use compare --help or ldsc --help for module options.' 'Other existing post-GWAS operations: /mnt/d/scripts/gwas/gwas_post.sh'
    ;;
  *) exec bash "$here/../gwas/gwas_post.sh" "$@" ;;
esac
