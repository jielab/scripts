#!/usr/bin/env bash
source /mnt/d/scripts/0f/console.sh
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)


# 🚩 Command-line help
# Keep public ARG examples together.
usage(){ cat <<'HELP'
cd /mnt/d/scripts/gu

./install.sh --arg
./arg.sh build --dir-gen /mnt/i/refGen/1kg/37 --method threads --format trace --threads 8 --jobs 1
./arg.sh check --dir-gen /mnt/i/refGen/1kg/37 --method threads --format trace
HELP
}


command=${1:-help}
case "$command" in
  -h|--help|help) usage; exit 0;;
  --help-all) usage; exit 0;;
esac
case "$command" in
  prep_gen|build|check) shift;;
  *) echo "ERROR: unknown command: $command (use ./arg.sh -h)" >&2; exit 2;;
esac
for option in "$@"; do
  case "$option" in
    -h|--help) usage; exit 0;;
    --help-all) usage; exit 0;;
  esac
done
export REFGEN_GEN4ARG_ONLY=0
source "$ROOT/f/arg_resources.sh"
arg_run_guarded "$command" "$@"
