#!/usr/bin/env bash
set -euo pipefail
test_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
test_work=${1:-$(mktemp -d /tmp/bmonth-smoke.XXXXXX)}
mkdir -p "$test_work"
Rscript --vanilla "$test_root/tests/test_core.R" "$test_root"
python3 "$test_root/tests/test_cli.py"
Rscript --vanilla "$test_root/tests/make_fixture.R" "$test_work/input"
"$test_root/bmonth.sh" all --in "$test_work/input" --out "$test_work/output" --med-boot 100 --phewas-cores 2 --replace
Rscript --vanilla "$test_root/tests/verify_fixture.R" "$test_work/output"
Rscript --vanilla "$test_root/tests/test_phewas.R" "$test_root" "$test_work"
printf 'Smoke-test outputs: %s\n' "$test_work"
