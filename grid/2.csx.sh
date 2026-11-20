#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
for arg in "$@"; do
    case "$arg" in -h|--help)
        cat <<'HELP'
PRS-CSx / PRS-CS: independent chromosome workers, then score and fit.

Quick start (3 traits x 22 autosomes = 66 independent chromosome jobs):
  ./2.csx.sh all --grch 38 --chr 1-22 --traits height,ldl,t2dm

Preview without submitting:
  ./2.csx.sh all --grch 38 --chr 1-22 --traits height,ldl,t2dm --plan-only

Generate weights only (no cohort, score or fit):
  ./2.csx.sh infer --grch 38 --chr 1-22 --traits height,ldl,t2dm

Usage: ./2.csx.sh [all|prepare|inputs|infer|score|fit] [options]
  all                    Cohort + inputs -> shared preparation -> chromosome jobs -> merge -> score -> fit
  infer                  Inputs -> shared preparation -> chromosome jobs -> validated genome-wide weights
  prepare                Prepare shared cohort only
  inputs                 Prepare GWAS manifest only
  score / fit            Use existing completed upstream results
  --chr, --chrs LIST             Autosomes 1-22 or subset, e.g. 1,2 [config: 1-22]; X/23 unsupported
  --traits CSV           height,ldl,t2dm
  --grch 38              Target genome build [38]
  --config FILE          YAML/JSON configuration [config/grid.yml]
  --analysis-root DIR    Output root [${ANALYSIS_ROOT}/grid/38]
  --split-file FILE|none Reuse a roster, or rebuild with cohort.seed [none; seed 20261115]
  --submit hpc|local|auto Execution mode [auto]
  --queue NAME           CPU queue [auto]
  --infer-jobs N         Concurrent model chains WITHIN EACH chromosome worker [4]
  --blas-threads N       BLAS threads PER CHAIN [1]; independent of the CPU budget
  --infer-memory-gb N    Host RAM per chromosome/input/merge worker [16 GiB]
  --cores N              CPU slots per job [shared 0f/0os.sh WORKFLOW_CORES]
  --job-hours N          Maximum duration per job [120 h]
  --threads N            CPU budget for scientific stages [--cores]
  --jobs N               Concurrent local workflow workers [1]; does not serialize HPC chromosomes
  --replace TRUE|FALSE   FALSE validates/reuses completed chains, repairs missing weights [FALSE]
                         TRUE recomputes requested stages; unfinished MCMC has no mid-chain checkpoint
  --plan-only            Write reviewable commands/plan; do not submit
  --wait TRUE|FALSE      Wait for HPC completion [FALSE]; local always waits
  -h, --help             Show this help without loading Python or accessing inputs

Additional input/scientific options:
  --genotype-dir DIR     Imputed chrN pfile directory; alias: --dir-imp
  --phenotype-file FILE  Phenotypes; alias: --pheno-file
  --cohort-file FILE     Prepared cohort input
  --pca-file FILE        PCA projections
  --ancestry-file FILE   Ancestry assignments
  --group-file FILE      Group assignments
  --relationship-file FILE  Relationship graph for family-aware splitting
  --gwas-manifest FILE   GWAS manifest (otherwise generated from config)
  --train-n N            Training sample budget [0: all TRAIN]
  --label-budget N       Combined TRAIN+TUNE budget [0: all]
  --phi-grid CSV         Positive fixed phi values [0.01]; auto is also fitted

CPU use after inference:
Score runs chromosomes concurrently within the job; fit runs independent models
and family-held-out folds concurrently, bounded by --threads and host RAM.
MCMC uses its separate --infer-jobs x --blas-threads budget (default 4 x 1),
not all reserved CPUs. There are only 10 independent chains per chromosome in
the default scientific configuration; MCMC iterations within a chain are sequential.
For future runs, all 10 chains can run concurrently with a larger RAM budget:
  ./2.csx.sh all --grch 38 --chr 1-22 --traits height,ldl,t2dm --infer-jobs 10 --infer-memory-gb 32
This is up to 10 CPU cores per chromosome, not 40; do not raise BLAS threads
blindly for the small LD-block solves. Existing completed results need no rerun.

Execution and downstream workflow:
This entry runs on CPU; it has no GPU inference option. Submit one command with
--chr 1-22: the scheduler already creates independent chromosome jobs. Do not
launch overlapping runs into the same analysis root. A subset --chr 1,2 produces
only that subset, not a genome-wide result. Use a separate root/config for tests.
After all requested traits finish merge, score AND fit, run:
  ./3.disco.sh all --grch 38 --traits height,ldl,t2dm
After Disco finishes, run:
  ./grid.sh all --grch 38 --traits height,ldl,t2dm
These are separate workflows: GRID does not submit or wait for upstream jobs.
A fully finished trait can proceed separately with --traits height throughout.
Use the same analysis root, cohort and chromosome configuration in every stage.

Each trait normalizes its GWAS once. All requested chromosome workers depend
only on that preparation, never on another chromosome. LSF schedules resources;
actual simultaneous starts depend on queue capacity. Names: csx.height.chr1.infer,
csx.height.chr2.infer, ..., csx.t2dm.chr22.infer. The merge requires all requested
chromosomes to succeed before publishing weights; score/fit wait for that merge.

Default scientific scope is unchanged: phi=0.01 + auto CSx, and fixed + auto
PRS-CS for AFR/EAS/EUR/SAS (10 chains per chromosome, 220 per trait).
CPU allocation follows shared 0f/0os.sh WORKFLOW_CORES. Chromosome workers retain 4 concurrent
model chains x 1 BLAS thread and 16 GiB; reserving 40 CPUs does not set BLAS to 40.
There is no need to set --jobs 22 on HPC. The reference workflow supports 22
autosomes, not chromosome X. Deployment requires 2.csx.sh, config/grid.yml and
the f/ helpers together; replacing only the shell entry is insufficient.
HELP
        exit 0
        ;;
    esac
done
exec bash "$ROOT/f/run_entry.sh" csx "$@"
