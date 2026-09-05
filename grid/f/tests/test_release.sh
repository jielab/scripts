#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
for f in "$ROOT"/*.sh "$ROOT"/f/*.sh "$ROOT"/f/tests/*.sh; do bash -n "$f"; done
python3 - "$ROOT" <<'PY'
import ast, sys
from pathlib import Path
for path in (Path(sys.argv[1]) / 'f').rglob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
PY
[[ -s "$ROOT/f/csx/PRScsx.py" ]]
[[ -s "/mnt/d/files/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight" ]]
[[ -s "/mnt/d/files/DiscoDivas/med.g1000.4pop.tsv" ]]
bash "$ROOT/grid.sh" help >/dev/null
echo PASS
