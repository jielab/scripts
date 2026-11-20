#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
case "${1:-}" in
    ""|-h|--help|help)
        cat <<'HELP'
GU: offline IBDmix segment report from completed chromosome exports.

Usage:
  ./gu.paper.sh run|ibdmix --target NAME --dir-ibdmix DIR [options]

Modules:
  run | ibdmix                    Generate NAME.segments.html; no inference jobs.

Options:
  --target NAME                   Required output prefix, e.g. 1kg or ukb.
  --dir-ibdmix DIR                Required directory of whole-chromosome segment exports.
  --dir-out DIR                   Report destination [current directory].
  --grch 37|38                    Segment coordinate build [37].
  --sample-panel FILE             Actual shared tested roster (TSV or PLINK PSAM).
                                  1kg defaults to bundled Phase 3 population panel.
  --x-males                       Use only male roster members for chrX;
                                  automatic for the bundled 1KG panel.
  --window-bp N                   Genomic summary window width [100000 bp].
  --source TEXT                   Input provenance displayed in report.
  -h, --help                      Show this help.

Examples:
  ./gu.paper.sh ibdmix --target 1kg --grch 37 --dir-ibdmix /work/sph-huangj/data/gen/1kg/37/ibdmix/segments --dir-out /work/sph-huangj/analysis/gu/reports
  ./gu.paper.sh run --target ukb --grch 37 --dir-ibdmix /work/sph-huangj/data/gen/ukb/37/ibdmix/segments --sample-panel TESTED.psam --x-males --window-bp 100000
  ./gu.paper.sh -h
HELP
        exit 0 ;;
    run|ibdmix) shift; exec python3 "$root/f/ibdmix_paper.py" "$@" ;;
    *) echo "Unknown report: $1. See ./gu.paper.sh -h." >&2; exit 2 ;;
esac
