#!/usr/bin/env bash
set -euo pipefail
case "${1:-allele3}" in
allele3) cat <<'HELP'
Usage: ./allele3.sh [all|hotspots|density|herv_enrichment|base_simulation|base_occupancy MODULE] [OPTIONS]

Modules:
    all              One CPU job: scan -> verify -> PCA -> census -> main analyses
    hotspots         Genome-wide hotspot discovery and matched comparisons [default]
    density          Continuous multi-ALT density from a completed scan
    herv_enrichment  HERV enrichment after scan and PCA
    base_occupancy   Observed nucleotide states, functional and ancient strata
    base_simulation  Four-base mutation, selection and drift simulation

Main options:
    --submit hpc|local|auto
    --out DIR --replace TRUE|FALSE [TRUE]
    --queue NAME --job-cores N --job-memory-gb N --job-hours N
    --job-dir DIR --write-cmd --wait TRUE|FALSE --plan-only
    all: --config JSON [config/workflow.json]
    hotspots: --scan-dir DIR --with-atlas auto|TRUE|FALSE --atlas-resources TSV
    -h, --help; --help-full for scientific options

Examples:
    ./allele3.sh all --submit hpc --plan-only
    ./allele3.sh all --submit hpc --job-cores 40 --job-memory-gb 160 --job-hours 120
    ./allele3.sh all --submit hpc --replace FALSE
    ./allele3.sh base_occupancy allele-time --help

Each entry creates its own command file. Do not wrap it in a second bsub call.
Use all for one sequential allocation; individual entries run only their modules.
Default overwrite; FALSE reuses complete modules and rebuilds incomplete modules.
Retained annotation/reference data are inputs. Ordinary GWAS uses ../gwas/gwas.sh.
See 0.README.md for the ordered configuration, outputs and scientific limits.
HELP
;;
herv_enrichment) cat <<'HELP'
Usage: ./allele3.sh herv_enrichment [STAGE] [OPTIONS]

Stages:
    primary -> summary -> spatial -> windows -> alleles -> time
            -> selection -> candidates -> genotypes -> ld -> report
    IBDmix is disabled by default. Inputs include a completed scan and shared PCA.

Options:
    --out DIR --scan-dir DIR --stages CSV --ibdmix TRUE|FALSE
    --metadata TSV --coordinates FILE --basis-dir DIR
    --replace TRUE|FALSE [TRUE]
    --submit hpc|local|auto --queue NAME --job-cores N --job-memory-gb N --job-hours N
    --job-dir DIR --write-cmd --wait TRUE|FALSE --plan-only --help-full

Use allele3.sh all to run scan, PCA, census and analyses in one allocation.
Use this entry after its inputs exist to run only the HERV analysis.
HELP
;;
*) printf 'Unknown help topic: %s\n' "$1" >&2; exit 2 ;;
esac
