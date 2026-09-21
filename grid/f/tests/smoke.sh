#!/usr/bin/env bash
# Actual tools on synthetic data only: 120 SNPs, 80 samples, 20 MCMC iterations.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
source "$ROOT/f/environment.sh"
test_root=$(mktemp -d "${TMPDIR:-/tmp}/grid-smoke.XXXXXX")
trap 'rm -rf -- "$test_root"' EXIT
fixture="$test_root/fixture"
python3 "$ROOT/f/tests/make_fixture.py" "$fixture"
plink2 --vcf "$fixture/gen/input.vcf" --make-pgen --out "$fixture/gen/chr22" > "$test_root/plink.log" 2>&1
plink2 --pfile "$fixture/gen/chr22" --make-bed --out "$fixture/gen/ukb_array" >> "$test_root/plink.log" 2>&1
args=(--models populations --trait height --chrs 22 --dir-gwas "$fixture/gwas" --dir-gen "$fixture/gen" --csx-bim-prefix "$fixture/gen/ukb_array" --csx-ref-dir "$fixture/ref" --csx-snpinfo "$fixture/ref/snpinfo_mult_1kg_hm3" --score-dir "$fixture/scores" --mcmc-iter 20 --mcmc-burnin 10 --mcmc-thin 2 --threads 1 --jobs 1)
bash "$ROOT/1csx.sh" "${args[@]}" --output-root "$test_root/work" > "$test_root/csx.log" 2>&1
# Fresh work directory proves score-only is independent of temporary inference files.
bash "$ROOT/1csx.sh" "${args[@]}" --output-root "$test_root/fresh" --stage score --replace TRUE > "$test_root/score.log" 2>&1
bash "$ROOT/2disco.sh" --trait height --chrs 22 --score-dir "$fixture/scores" --output-root "$test_root/fresh" --pca-file "$fixture/pca.tsv.gz" --med-file "$fixture/centers.tsv" > "$test_root/disco.log" 2>&1
python3 "$ROOT/f/tests/check_fixture.py" "$fixture"
# Dry-run produces no posterior completion markers or permanent weights.
bash "$ROOT/1csx.sh" "${args[@]}" --output-root "$test_root/plan" --dry-run TRUE > "$test_root/plan.log" 2>&1
if find "$test_root/plan" -name '*.done' | read -r _; then echo 'Unexpected dry-run marker' >&2; exit 1; fi
# A nonexistent chromosome and a missing score must fail before heavy computation.
if bash "$ROOT/1csx.sh" --chrs 23 --check > "$test_root/invalid.log" 2>&1; then exit 1; fi
if bash "$ROOT/2disco.sh" --trait height --score-dir "$test_root/missing" --pca-file "$fixture/pca.tsv.gz" --output-root "$test_root/negative" --check > "$test_root/missing.log" 2>&1; then exit 1; fi
printf 'PASS: inference, scoring from permanent weights, official DiscoDivas, dry-run and failure checks.\nTemporary test artifacts (removed on exit): %s\n' "$test_root"
