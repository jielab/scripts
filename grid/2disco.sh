#!/usr/bin/env bash
# Official DiscoDivas: align PCA/PRS -> genetic-distance interpolation -> publish.
# Implementation: f/disco.sh; pinned upstream program: f/disco/DiscoDivas.R.
# Reference: https://github.com/YunfengRuan/DiscoDivas
set -euo pipefail
usage(){ cat <<'HELP'
DiscoDivas — combine the four ancestry-specific CSx scores for each UKB individual

Usage examples (WSL):
  cd /mnt/d/scripts/grid
  ./0pca.sh --check
  ./1csx.sh --traits height,ldl,t2dm --jobs 4 --threads 1
  ./2disco.sh --traits height,ldl,t2dm --check
  ./2disco.sh --traits height,ldl,t2dm
  ./2disco.sh --trait height --a-list 1,1,1,1 --regress-pca TRUE

Prerequisites:
  Completed reference PCA projection and the permanent merged csx.pgs.gz file.
  Existing projection in pca_proj can be reused; no need to rerun it per trait.
  Scored samples absent from PCA are filtered out of all four PRS inputs.
  Negative IDs and /mnt/d/files/ukb.exclude.id are excluded.
  Logs report excluded, filtered and retained sample counts.
  DiscoDivas combines individual scores; it does not directly analyze GWAS or
  produce a single universal SNP weight file (weights vary by individual).

Important parameters:
  --trait NAME / --traits LIST   Default: height,ldl,t2dm; traits run sequentially.
  --pca-file FILE               /mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz
  --med-file FILE               /mnt/d/files/DiscoDivas/med.g1000.4pop.tsv
  --distance-pcs N              10; reference and target use the same PCs (>=5).
  --a-list A1,A2,A3,A4          Default 1,1,1,1, in AFR,EAS,EUR,SAS order.
  --regress-pca TRUE|FALSE      Default TRUE (official recommendation).
  --score-dir DIR               /mnt/d/data/ukb/pgs (input and permanent output).
  --output-root DIR             /mnt/d/analysis/grid (temporary).
  --replace TRUE|FALSE          Default FALSE; reuse matching outputs.
  --check                       Check prerequisite files/packages; no computation.
  --dry-run TRUE                Same preflight/plan behavior.
  --chrs LIST                   Match the chromosome subset used in CSx, if any.

Input scores:  <score-dir>/<trait>/csx.pgs.gz
Output scores: <score-dir>/<trait>/disco.pgs.gz (IID,disco)
Coefficients:  <score-dir>/<trait>/disco.coef.tsv.gz
Commands/logs/intermediates: /mnt/d/analysis/grid/disco/

Uses official DiscoDivas distance-matrix interpolation, not inverse-square mixing.
Defaults use reference centers and equal quality factors, without phenotype tuning.
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0;; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/pipeline.sh" disco "$@"
