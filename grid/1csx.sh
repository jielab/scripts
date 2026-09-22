#!/usr/bin/env bash
# PRS-CSx: GWAS preparation -> joint MCMC -> permanent SNP weights -> UKB PRS.
# Implementation: f/csx.sh, f/csx_score.sh; official code: f/csx/.
# Reference: https://github.com/getian107/PRScsx
set -euo pipefail
usage(){ cat <<'HELP'
PRS-CSx — height / ldl / t2dm, jointly using AFR,EAS,EUR,SAS GWAS

Usage examples (WSL):
  cd /mnt/d/scripts/grid
  ./1csx.sh --traits height,ldl,t2dm --check
  ./1csx.sh --traits height,ldl,t2dm --jobs 4 --threads 1
  ./1csx.sh --traits height,ldl,t2dm --models auto,meta --jobs 4 --threads 1
  ./1csx.sh --trait height --stage weights --jobs 4 --threads 1
  ./1csx.sh --trait height --stage score --jobs 4 --threads 8

Modules:
  weights  Normalize BETA+SE GWAS -> joint four-population PRS-CSx by chromosome.
  score    Apply saved SNP weights to UKB; four individual PRS files.
  all      Both modules (default). No PCA is needed for CSx.

Important parameters:
  --trait NAME / --traits LIST   Default: height,ldl,t2dm; traits run sequentially.
  --stage all|weights|score      Default: all.
  --models LIST                Default: populations,auto,meta.
                               Use auto,meta to append only the two combined scores.
                               auto: learned phi + posterior meta; meta: fixed --phi + posterior meta.
  --dir-gwas DIR                /mnt/e/gwas/4grid/common
  --dir-gen DIR                 /mnt/e/ukbGen/37/hap (chr*.pgen/pvar/psam)
  --csx-bim-prefix PREFIX       Target variant list; default: above DIR/ukb_array
  --csx-ref-dir DIR             /mnt/e/refLD/csx (1000 Genomes LD)
  --csx-snpinfo FILE            Default: snpinfo_mult_1kg_hm3 in default LD folder
  --phi VALUE|auto              Default: 1e-2; fixed value, NOT phenotype-tuned.
  --mcmc-iter N                 4000; --mcmc-burnin 2000; --mcmc-thin 5
  --n-gwas N|AFR=N,EAS=N,...     Default: median N over retained HM3 SNPs.
                                For t2dm, verify N is effective sample size.
  --jobs N / --threads N        Parallel chromosomes / threads per job (4 / 8).
  --chrs 1-22                  Autosomes; subsets get separate filenames/directories.
  --seed N                     20260904 (chromosome number added for each chain).
  --keep FILE / --remove FILE  PLINK sample lists; default remove: /mnt/d/files/ukb.exclude.id.
  --score-dir DIR              /mnt/d/data/ukb/pgs
  --output-root DIR            /mnt/d/analysis/grid (temporary).
  --replace TRUE|FALSE         Default FALSE; reuse only matching saved results.
  --check                      Validate inputs only; no inference/scoring.
  --dry-run TRUE               Same preflight/plan behavior; no completion markers.

Permanent outputs:
  <GWAS folder>/<original-name>.csx.sumstats.gz (normalized GWAS BETA,SE,N,...)
  Matching .csx.sumstats.json enables reuse after deleting the temporary work directory.
  <GWAS folder>/<original-name>.csx.gz  (SNP,A1,BETA,CHR,BP,A2; NOT individual PRS)
  Example: /mnt/e/gwas/4grid/common/height.AFR/gwas/height.AFR.csx.gz
  /mnt/d/data/ukb/pgs/<trait>/csx.pgs.gz (eid, csx.AFR/EAS/EUR/SAS, csx.auto, csx.meta)
  Combined SNP weights: <score-dir>/<trait>/.weights/csx.{auto,meta}.gz
  Per-population scores and score signatures: /mnt/d/analysis/grid/csx/scores/<trait>/
  t2dm.AFA is mapped to AFR LD, but its GWAS output retains the name t2dm.AFA.

Temporary files/commands/logs: /mnt/d/analysis/grid/csx/
Input coordinates must be GRCh37; the preflight checks rsID-position agreement.
For another reference/target, set --csx-snpinfo and --csx-bim-prefix explicitly.
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0;; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/f/pipeline.sh" csx "$@"
