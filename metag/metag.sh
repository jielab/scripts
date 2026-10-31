#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
usage() {
cat <<'EOF'
metag: analyze completed split data.

Usage:
  ./metag.sh composition --sample /mnt/f/metag/split/<dataset>/<run>
  ./metag.sh composition --sample /path/to/sample --out /path/to/report

Reports: species/group read fractions and recorded evidence status.
Compact split retains per-read assignments; full reference hits are not collected.
The optional --screen report contains sample-level candidates only.

Preparation steps:
  ./1.library.sh        Check references and build permanent libraries
  ./2.split.sh            Classify and split samples

Split data: /mnt/f/metag/split
Settings: f/metag.conf

Legacy build/classify/split commands remain accepted for compatibility.
EOF
}
case "${1:-}" in
  ''|-h|--help|help) usage; exit 0;;
esac
command="$1"; shift
case "$command" in
  build-db|prepare) exec bash "$HERE/1.library.sh" "$@";;
  check-ref|check_ref|ref-check) exec bash "$HERE/1.library.sh" --check-only "$@";;
  build-sourmash) exec bash "$HERE/1.library.sh" --sourmash-only "$@";;
  build-metacache) exec bash "$HERE/1.library.sh" --metacache-only "$@";;
  build-kraken) exec bash "$HERE/1.library.sh" --kraken "$@";;
  split|classify) exec bash "$HERE/2.split.sh" "$command" "$@";;
  check|gpu-check) exec bash "$HERE/f/run_workflow.sh" "$command" "$@";;
  composition) exec bash "$HERE/f/run_workflow.sh" composition "$@";;
  *) printf 'Unknown analysis command: %s\n' "$command" >&2; usage; exit 2;;
esac

