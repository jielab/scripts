#!/usr/bin/env bash
cat <<'HELP'
GRID: Genealogy-informed effect transport for multi-ancestry PRS

Usage
  ./pca.sh [options]
  ./csx.sh [options]
  ./disco.sh [options]
  ./grid.sh grid|eval [options]
  ./grid.sh pca [options]        shortcut for ./pca.sh

Recommended order
  pca.sh (once per cohort) -> csx.sh -> disco.sh
  ./grid.sh grid -> ./grid.sh eval
  grid automatically checks ARG, extracts LD, assembles features, fits the
  model, generates weights and scores participants in one run.
  ancestry is an internal stage of pca.sh.

PCA preparation
  --pca-file PATH               existing/projected PC table; defaults to
      /mnt/d/data/ukb/phe/pca/ukb.discodivas.pca.tsv.gz
  --pca-weight PATH             /mnt/d/files/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight
  --med-file PATH               /mnt/d/files/DiscoDivas/med.g1000.4pop.tsv
  --cov-pcs 20 --distance-pcs 10
  --ancestry-file PATH          output; defaults beside the PC table
  --ancestry-prob-threshold 0.90
  Existing projected PCs skip PLINK projection even if distance/ancestry files
  are missing. Each downstream stage reuses complete outputs. Use --replace
  TRUE when changing inputs/settings to rebuild all stages.

Core defaults
  --dir-gwas /mnt/d/data.BIG/gwas/4grid
      Flat files: height.AFR.gz, height.EAS.gz, ..., t2dm.SAS.gz
  --dir-gen /mnt/d/data/ukb/gen/typ
      Target chr1..chr22 PLINK bfiles or pfiles used for PRS scoring.
      If absent, GRID automatically uses /mnt/h/ukbGen/37/hap when PGENs exist.
  --dir-imp /mnt/h/ukbGen/37/imp
  --phe-file /mnt/d/data/ukb/phe/Rdata/phe.rds
  --output-root /mnt/d/analysis/grid
  --trait height                  one of height, ldl, t2dm
  --pops AFR,EAS,EUR,SAS
  --chrs 1-22
  --jobs 4 --threads 8 --replace FALSE

ARG-Needle options
  --arg-hap-dir /mnt/h/ukbGen/37/hap
      Accepts phased BGEN + Oxford .sample, or phased PGEN + .psam files.
  --arg-out/--arg-dir DIR         reusable ARG root [beside hap directory: ../arg]

  ARG construction is shared data preparation and is no longer performed by
  grid.sh. Build it first (ARG-Needle is the default method):
    bash /mnt/d/scripts/gu/arg.sh build --method needle \
      --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
      --map-dir /mnt/d/data.BIG/refGen/maps/GRCh37

GRID model options
  --ldscore-dir DIR               generated from PRS-CSx HDF5 LD panels
  --external-age FILE             optional GEVA/other age table; ARG age remains primary
  --grid-ridge-alpha 10
  --grid-max-snps-per-chr 0       0 means all matched SNPs
  --require-ld TRUE

Evaluation
  --eval-repeats 100 --eval-folds 5 --min-group 100
  --eval-zero-shot TRUE           evaluates ancestry-matched/posterior scores without phenotype tuning
  --eval-tuned TRUE               separately reports phenotype-tuned stacking as a comparator

Examples:
  cd /mnt/d/scripts/grid

  ./grid.sh pca # project PCA to 1KG, calculate genetic distance, infer ancestry

  ./csx.sh --trait height --chrs 22
  ./disco.sh --trait height

  ./grid.sh grid --trait height --chrs 22
  ./grid.sh eval --trait height

  ./csx.sh --trait ldl --chrs 22
  ./disco.sh --trait ldl
  ./grid.sh grid --trait ldl --chrs 22
  ./grid.sh eval --trait ldl

  ./csx.sh --trait t2dm --chrs 22
  ./disco.sh --trait t2dm
  ./grid.sh grid --trait t2dm --chrs 22
  ./grid.sh eval --trait t2dm
HELP
