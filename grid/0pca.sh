#!/usr/bin/env bash
# UKB preparation: reference PCA projection -> distance QC -> ancestry assignment.
# Implementation: f/pca.sh, f/combine_disco_pca.R, f/prepare_ancestry_auto.R.
set -euo pipefail
usage(){ cat <<'HELP'
PCA projection — project UKB into the reference space used by DiscoDivas

Usage examples (WSL):
  cd /mnt/d/scripts/grid
  ./0pca.sh --check
  ./0pca.sh --dir-imp /mnt/e/ukbGen/37/imp --jobs 4
  ./0pca.sh --replace TRUE --jobs 4

Modules (run in order; existing completed outputs are reused):
  projection  Score chromosomes 1-22 with reference PCA loadings; sum PC scores.
  distance    Reference distances, PCA QC workbook and figure.
  ancestry    Phenotype-assisted ancestry assignment and QC (separate from Disco).

Important parameters:
  --dir-imp DIR                 /mnt/e/ukbGen/37/imp (chr*.pgen/pvar/psam)
  --pca-file FILE               /mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz
  --pca-weight FILE             /mnt/d/files/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight
  --med-file FILE               /mnt/d/files/DiscoDivas/med.g1000.4pop.tsv
  --cov-pcs N / --distance-pcs N 20 / 10; require 5 <= distance <= covariate PCs.
  --phe-file FILE               /mnt/d/data/ukb/phe/Rdata/phe.rds
  --ancestry-file FILE          Default: ukb.ancestry.auto.tsv.gz beside PCA file.
  --ancestry-prob-min VALUE     0.90; --anchor-max-per-group 10000
  --jobs N                     At most 4 concurrent chromosomes, memory-limited.
  --threads N                  Accepted shared option; projection uses 1 thread/job.
  --output-root DIR             /mnt/d/analysis/grid (temporary).
  --replace TRUE|FALSE          Default FALSE. TRUE rebuilds projection and QC.
  --check                       Validate existing projection or genotype inputs.
  --dry-run TRUE                Same preflight/plan behavior; no analysis.

Permanent PCA, distance, ancestry and QC files: /mnt/d/data/ukb/pca_proj/
Temporary chromosome scores/commands/logs: /mnt/d/analysis/grid/pca/
Existing PCA is reused independently of temporary files. If changing reference
loadings/genotypes, explicitly use --replace TRUE or a different --pca-file.
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0;; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/pipeline.sh" pca "$@"
