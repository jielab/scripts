#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../f/gwas_post_perf.f.sh"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
FINAL="$work/trait.gz"
THIN_OUT="$work/trait.thin.gz"
printf 'CHR\tPOS\tSNP\tP\n1\t100\trs1\t0.01\n' | bgzip -c > "$FINAL"
cp "$FINAL" "$THIN_OUT"
touch -d '2020-01-01' "$FINAL"
tabix -S 1 -s 1 -b 2 -e 2 "$THIN_OUT"
THIN_MODE=TRUE
gwas_post_log() { :; }
Rscript() { echo 'Unexpected thin regeneration' >&2; return 99; }
before=$(sha256sum "$THIN_OUT" "$THIN_OUT.tbi")
for DO_STEP in mplot magma,mplot; do
  for REPLACE in FALSE TRUE; do
    # Old path-dependent markers, and markers removed by interrupted old runs,
    # must both leave existing thin artifacts untouched.
    printf 'obsolete-path-dependent-signature\n' > "$THIN_OUT.done"
    gwas_post_thin
    [[ $(cat "$THIN_OUT.done") == obsolete-path-dependent-signature ]]
    rm "$THIN_OUT.done"
    gwas_post_thin
    [[ ! -e "$THIN_OUT.done" ]]
    [[ $(sha256sum "$THIN_OUT" "$THIN_OUT.tbi") == "$before" ]]
  done
done
for steps in thin thin,mplot mplot,thin format,mplot liftover,mplot all; do
  if gwas_thin_plot_only "$steps"; then exit 1; fi
done
DO_STEP=mplot
expect_failure() {
  if gwas_post_thin 2> "$work/error"; then exit 1; fi
  grep -q 'run thin first' "$work/error"
}
mv "$THIN_OUT.tbi" "$work/index"
expect_failure
printf 'broken index\n' > "$THIN_OUT.tbi"
expect_failure
mv "$work/index" "$THIN_OUT.tbi"
touch -d '2030-01-01' "$FINAL"
expect_failure
touch -d '2020-01-01' "$FINAL"
rm "$THIN_OUT"
expect_failure
echo 'PASS: plot-only reuse, replace, legacy markers, missing/corrupt/stale artifacts'
