#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
METAG: split paired-end FASTQ by the most specific supported taxon.

Usage:
  ./metag.sh run [options]
  ./metag.sh split [options]
  ./metag.sh refine --dat-list england,Wxiamen [options]

Modules:
  run                             Split and review evidence in one resumable workflow.
                                  Resume completed samples and steps from this release.
  split                           Classify all pairs against the full Kraken union,
                                  review each mate, align all pairs to a human/plant panel,
                                  compare voted eukaryotic background when it fits, then
                                  write one paired FASTQ bin per supported taxon.
  refine                          Review completed bins for human and food-source DNA,
                                  optional local BLAST and community AMR, and viral signals.
                                  Use ./metag.sh refine -h for its options.
  Publication                     Use ./metag.paper.sh for QC, truth evaluation and reports.

Options:
  --dat-list CSV                  Dataset names under metag/raw.
  --sample DATASET/SAMPLE         Select one sample; repeatable; excludes batch options.
  --batch SELECTION               Per-dataset batch numbers: 1-2, 1,3-5, 3-, all [all].
  --batch-size N                  Override automatic batching (advanced).
  --config FILE                   Project configuration override.
  --root DIR                      Root override for root-relative configuration paths.
  --replace TRUE|FALSE            Reuse matching current results [FALSE]. TRUE rebuilds selected samples.
  --metadata TSV                  Sample material and library type for evidence review.
  --seal-audit                    Audit all Weddell species and genus bin pairs with BLAST+.
                                  Build/reuse a separate local competition panel on the worker.
  --seal-db DIR                   Use a versioned panel built by refine/seal_panel.py.
  --blast-db PREFIX               Optional local nucleotide reference for source audits.
  --blast-pairs N                 Maximum pairs sampled per eligible bin [100].
  --amr-markers FILE              Optional community-level ShortBRED AMR markers.
  --amr-annotation TSV            Optional marker family annotation.
  --refine-queue NAME             Evidence job queue [configured convert_queue].
  --refine-cores N                Evidence job cores [1 or 4 with sequence searches].
  --refine-mem-gib N              Evidence job memory reservation [32 GiB].
                                  Evidence options above apply to run. See refine -h for refine.
  --submit auto|local|hpc         Execution mode [auto].
  --queue NAME                    LSF queue [configured huge-memory queue].
  --cores N                       CPU budget [shared site default].
  --threads N                     Worker threads [max(1, cores-2)].
  --mem-gib N                     Job memory reservation [2048 GiB; RAM preload].
  --hours N                       Job duration [configured].
  --plan-only                     Generate reviewable job files.
  --confidence X                  Full-library Kraken confidence [configured, normally 0.1].
  --representatives FILE          Override species reference TSV with species_taxid, name,
                                  group, path, bases. Relative FASTA paths use reference_raw.
  --panel-taxids CSV              Prioritize named green-plant species for competitive mapping.
  --panel-gib X                   Max panel sequence bases in units of 2^30 [normally 24].
  --max-refs N                    Max panel genomes per sample [normally 60].
  --min-windows N                 Min concordant 10-kb nuclear windows per exact human/plant
                                  call in the sample [3]. Count distinct 10-kb windows.
  -h, --help                      Show this help.

Output:
  metag/species/ranked/DATASET/SAMPLE/bins/GROUP/RANK/TAXID/{R1,R2}.fastq.gz
  assignments.tsv.gz, taxa.tsv, transitions.tsv, evidence.tsv, panel.tsv and summary.json

Examples:
  ./metag.sh run --sample england/ERR14374050 --plan-only
  ./metag.sh run --sample england/ERR14374050
  ./metag.sh run --dat-list england,Wxiamen --batch 1-2 --batch-size 10
  ./metag.paper.sh run --dat-list england,Wxiamen
HELP
}
case "${1:-}" in ""|-h|--help|--help-full|help) usage; exit 0;; --plan-only) set -- run "$@";; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ ${1:-} == refine && ( ${2:-} == -h || ${2:-} == --help || $# -eq 1 ) ]]; then
    cat <<'REFINE_HELP'
METAG: review source evidence in completed ranked splits.

Usage:
  ./metag.sh refine --dat-list england,Wxiamen [options]
  ./metag.sh refine --sample DATASET/SAMPLE [options]

Options:
  --dat-list CSV                  Include completed samples in these datasets.
  --sample DATASET/SAMPLE         One completed sample; repeatable.
  --metadata TSV                  sample, sample_type, library_type [unknown if absent].
  --seal-audit                    Audit all Weddell species and genus bin pairs with BLAST+.
                                  Build/reuse a separate local competition panel on the worker.
  --seal-db DIR                   Use a versioned panel built by refine/seal_panel.py.
  --blast-db PREFIX               Local broad nucleotide BLAST DB with taxids [METAG_BLAST_DB].
  --blast-pairs N                 Uniform read pairs per eligible bin [100].
  --amr-markers FILE              ShortBRED AMR protein markers [METAG_AMR_MARKERS].
  --amr-annotation TSV            Optional family annotations [METAG_AMR_ANNOTATION].
                                  AMR uses DIAMOND blastx and curated true markers.
  --config FILE, --root DIR       Same project configuration as split.
  --submit auto|local|hpc         Execution mode [auto].
  --replace TRUE|FALSE            Recompute requested evidence steps or reuse matching ones [FALSE].
  --batch-size N                  Samples per cohort job [all for tables, 4 with searches].
  --queue NAME, --cores N         Compute queue and cores [1 baseline, 4 with external tools].
  --mem-gib N, --hours N          Job allocation [32 GiB, 24 h].
  --plan-only                     Write job files without running.
  -h, --help                      Show this help.

Output:
  analysis/metag/refine/main/DATASET/SAMPLE/  Source evidence, method status and receipt.
  analysis/metag/refine/{jobs,logs}/          Submission files and logs.

Examples:
  ./metag.sh refine --dat-list Wxiamen
  ./metag.paper.sh run --dat-list Wxiamen
REFINE_HELP
    exit 0
fi
if [[ ( ${1:-} == run || ${1:-} == split ) && ( ${2:-} == -h || ${2:-} == --help ) ]]; then usage; exit 0; fi
if [[ ${1:-} == run ]]; then
    printf 'PHASE shell-environment start host=%s pid=%s\n' "${HOSTNAME:-unknown}" "$$"
fi
source "$ROOT/f/env.sh"
if [[ ${1:-} == run ]]; then
    printf 'PHASE shell-environment done; starting Python: %s\n' "$METAG_PYTHON"
fi
if [[ ${1:-} == run ]]; then
    shift
    exec "$METAG_PYTHON" -B "$ROOT/workflow.py" "$@"
fi
if [[ ${1:-} == refine ]]; then
    shift
    exec "$METAG_PYTHON" -B "$ROOT/refine/run.py" "$@"
fi
exec "$METAG_PYTHON" "$ROOT/f/rank_workflow.py" "$@"
