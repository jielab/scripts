#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
METAG: summarize completed ranked splits for a manuscript.

Usage:
  ./metag.paper.sh run --dat-list england,Wxiamen [options]

Modules:
  run                             Summarize completed samples and generate cohort tables
                                  and figures. No arguments prints this help.

Options:
  --dat-list CSV                  Include all completed ranked samples in the named datasets.
  --sample DATASET/SAMPLE         Select a specific completed sample; repeatable.
  --result-dir DIR                Read one completed ranked sample directly.
                                  Without selectors, run includes all completed datasets.
  --check-bins                    Also stream output FASTQ and assignments for integrity/QC.
  --raw-check                     Also verify original R1/R2 and assignment-to-bin order.
  --truth FILE                    Independent TSV (also .gz): ordinal,true_taxid.
                                  For multiple samples include sample (DATASET/SAMPLE).
  --out DIR                       Report root [analysis/metag/paper].
  --config FILE                   Project configuration override.
  --root DIR                      Root override for root-relative configuration paths.
  --submit auto|local|hpc         Execution mode [auto].
  --queue NAME                    LSF queue [configured].
  --cores N                       Report job CPU allocation [2].
  --mem-gib N                     Report job RAM [32 GiB].
  --hours N                       Report job duration [24].
  --plan-only                     Write job files without submitting.
  -h, --help                      Show this help.

Output:
  main/index.html                 Interactive cohort, sample, taxon and evidence explorer.
  main/ranked/explore.html        Self-contained interactive report with table/SVG export.
  main/ranked/report.html         Detailed QC, manuscript tables and figures.
  main/ranked/                    Tables (TSV), figures (PDF, PNG), sample and evidence data.
  main/ranked/seal_*.tsv          Seal audit status, per-bin outcomes and occupied windows.
  jobs/run*/                      Commands and parameters. HPC job name: metag.paper.
  logs/run*/                      Status, report log and scheduler output.

Examples:
  ./metag.paper.sh run --dat-list england,Wxiamen
  ./metag.paper.sh run --sample england/ERR14374050 --raw-check
HELP
}

case "${1:-}" in
    ""|-h|--help|--help-full|help) usage; exit 0 ;;
    run) case "${2:-}" in -h|--help) usage; exit 0 ;; esac ;;
    *) printf 'Expected run. Use ./metag.paper.sh -h for help.\n' >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/../0f/paper_name.sh"
workflow_paper_name "${BASH_SOURCE[0]}"
source "$ROOT/f/env.sh"
exec "$METAG_PYTHON" "$ROOT/f/rank_paper.py" "$@"
