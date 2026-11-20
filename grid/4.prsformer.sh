#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
stage=all
case ${1:-} in all|plan|prepare|train|predict) stage=$1;shift ;; esac
exec bash "$ROOT/f/run_entry.sh" prsformer all --prsformer-stage "$stage" "$@"
