#!/usr/bin/env bash
set -euo pipefail
if (( $# )); then
    case "$1" in
        --grch) export GRCH="${2:-}"; shift; if (( $# )); then shift; fi ;;
        --grch=*|grch=*) export GRCH="${1#*=}"; shift ;;
        *) echo 'Usage: bash f/check.sh [--grch 37|38]' >&2; exit 2 ;;
    esac
    if (( $# )) || [[ "$GRCH" != 37 && "$GRCH" != 38 ]]; then
        echo '--grch must be 37 or 38' >&2; exit 2
    fi
fi
source "$(dirname -- "${BASH_SOURCE[0]}")/env.sh"
python3 - <<'PY'
import os, sys
from pathlib import Path
if sys.version_info < (3, 11):
    raise SystemExit('Need Python >= 3.11')
print('Python:', sys.version.split()[0])
missing = []
sys.path.insert(0, os.environ['LD_CODE_DIR'])
from ld_pipeline import parse_chroms
chroms = parse_chroms(os.environ.get('CHRS', '1-22'))
for c in chroms:
    for suffix in ('.pgen', '.psam'):
        p = Path(os.environ['PFILE_DIR']) / f'chr{c}{suffix}'
        if not p.is_file() or not p.stat().st_size:
            missing.append(str(p))
    prefix = Path(os.environ['PFILE_DIR']) / f'chr{c}'
    if not any(Path(str(prefix) + s).is_file() for s in ('.pvar.zst', '.pvar')):
        missing.append(str(prefix) + '.pvar[.zst]')
p = Path(os.environ['PANEL'])
if not p.is_file() or not p.stat().st_size:
    missing.append(str(p))
else:
    from collections import Counter
    with p.open() as f:
        header = next(f).split()
        si, pi = header.index('sample'), header.index('super_pop')
        rows = [line.split() for line in f if line.strip()]
    counts = Counter(r[pi] for r in rows)
    assert len({r[si] for r in rows}) == len(rows), 'Duplicate panel samples'
    expected = (dict(AFR=893, EAS=585, EUR=633, SAS=601, AMR=490) if os.environ['GRCH'] == '38'
                else dict(AFR=661, EAS=504, EUR=503, SAS=489, AMR=347))
    for pop, n in expected.items():
        assert counts[pop] == n, (pop, counts[pop], n)
    assert len(rows) == sum(expected.values()), 'Unexpected panel sample count'
    expected_ids = {r[si] for r in rows}
    import csv
    for c in chroms:
        psam = Path(os.environ['PFILE_DIR']) / f'chr{c}.psam'
        if psam.is_file():
            with psam.open() as f:
                header = next(f).lstrip('#').split()
                sample_rows = [dict(zip(header, line.split())) for line in f if line.strip()]
            ids = [r['IID'] for r in sample_rows]
            assert len(ids) == len(set(ids)) and set(ids) == expected_ids, f'chr{c}: PSAM sample mismatch'
    print(f'Panel/PSAM samples: {len(rows)} total; {sum(n for p,n in expected.items() if p != "AMR")} in AFR/EAS/EUR/SAS; AMR excluded')
manifest = Path(os.environ['LD_INPUT_MANIFEST'])
if not manifest.is_file() and os.environ['GRCH'] == '38':
    raise SystemExit(f'Missing GRCh38 input manifest: {manifest}')
for line in (manifest.read_text().splitlines()[1:] if manifest.is_file() else []):
    name, size = line.split('\t')
    if name != 'samples.txt' and int(Path(name).name.split('.')[0][3:]) not in chroms:
        continue
    path = (Path(os.environ['PANEL']) if name == 'samples.txt' else
            Path(os.environ['PFILE_DIR']) / Path(name).name)
    if not path.is_file():
        missing.append(str(path))
    elif path.stat().st_size != int(size):
        missing.append(f'{path}: size {path.stat().st_size}, expected {size} (transfer incomplete or different input)')
if missing:
    print('Missing/empty inputs:\n' + '\n'.join(missing))
    raise SystemExit(2)
print(f'GRCh{os.environ["GRCH"]} chromosomes {chroms}: PGEN/PVAR/PSAM sets and panel are present.')
print('File sizes match the local input manifest; this is not a content checksum.' if manifest.is_file()
      else 'No optional input size manifest; checked file presence and sample identities only.')
PY
"$PLINK2" --version
"$PLINK2" --help r2-phased indep-pairphase >/dev/null
scratch=$(mktemp -d)
trap 'rm -rf -- "$scratch"' EXIT
printf '#include <zstd.h>\n#include <filesystem>\nint main(){return ZSTD_versionNumber() < 10000;}\n' > "$scratch/check.cpp"
"$CXX" -std=c++17 "$scratch/check.cpp" -lzstd -o "$scratch/check"
"$scratch/check"
mkdir -p "$OUT" "$LD_TMP_DIR"
test -w "$OUT" && test -w "$LD_TMP_DIR"
df -h "$OUT" "$LD_TMP_DIR"
echo 'HPC preflight PASS (no LD computation or job submission).'
