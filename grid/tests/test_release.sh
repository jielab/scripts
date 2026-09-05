#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
for f in "$ROOT/grid.sh" "$ROOT"/f/*.sh "$ROOT"/install/*.sh "$ROOT"/tests/*.sh; do bash -n "$f"; done
for f in "$ROOT"/f/*.py "$ROOT"/optional/*.py; do [[ -e $f ]] && python3 -m py_compile "$f"; done
[[ -s "$ROOT/csx/PRScsx.py" ]]
[[ -s "$ROOT/DiscoDivas/files/g1k_hm3_maf5_woamb_wolr.pca.weight" ]]
[[ -s "$ROOT/DiscoDivas/files/med.g1000.4pop.tsv" ]]
bash "$ROOT/grid.sh" help >/dev/null
echo PASS
