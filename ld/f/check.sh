#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/env.sh"
for tool in "$PYTHON" "$PLINK2" "$LD_CXX" "$CXX"; do
    command -v "$tool" >/dev/null || { echo "Missing/non-executable tool: $tool" >&2; exit 2; }
done
"$PYTHON" - <<'PY'
import os, sys
from pathlib import Path
sys.path.insert(0, os.environ['LD_CODE_DIR'])
from ld_pipeline import parse_chroms
if sys.version_info < (3,11): raise SystemExit('Need Python >=3.11')
chroms=parse_chroms(os.environ['CHRS'])
panel=Path(os.environ['PANEL'])
if not panel.is_file(): raise SystemExit(f'Missing panel: {panel}')
with panel.open() as stream:
    header=next(stream).split()
    if not {'sample','super_pop'} <= set(header): raise SystemExit('Panel needs sample and super_pop columns')
    rows=[dict(zip(header,line.split())) for line in stream if line.strip()]
ids=[r['sample'] for r in rows]
if len(ids)!=len(set(ids)): raise SystemExit('Duplicate panel samples')
for pop in ('AFR','EAS','EUR','SAS'):
    if sum(r['super_pop']==pop for r in rows)<2: raise SystemExit(f'Need >=2 samples for {pop}')
for chrom in chroms:
    if os.environ['INPUT_FORMAT']=='pgen':
        prefix=Path(os.environ['PFILE_DIR'])/f'chr{chrom}'
        files=[Path(str(prefix)+s) for s in ('.pgen','.psam')]
        pvars=[Path(str(prefix)+s) for s in ('.pvar.zst','.pvar')]
        files.append(next((p for p in pvars if p.is_file()), pvars[0]))
    elif os.environ['INPUT_FORMAT']=='vcf':
        prefix=Path(os.environ['VCF_DIR'])/f'chr{chrom}.vcf.gz'
        files=[prefix,Path(str(prefix)+'.tbi')]
    else: raise SystemExit('INPUT_FORMAT must be pgen/vcf')
    for path in files:
        if not path.is_file() or not path.stat().st_size: raise SystemExit(f'Missing/empty input: {path}')
    if os.environ['INPUT_FORMAT']=='pgen':
        with files[1].open() as stream:
            header=next(stream).lstrip('#').split()
            samples=[dict(zip(header,line.split())) for line in stream if line.strip()]
        actual=[r['IID'] for r in samples]
        required={r['sample'] for r in rows if r['super_pop'] in ('AFR','EAS','EUR','SAS')}
        if len(actual)!=len(set(actual)) or not required<=set(actual): raise SystemExit(f'chr{chrom}: PSAM identities mismatch')
manifest=os.environ.get('LD_INPUT_MANIFEST')
if manifest:
    for line in Path(manifest).read_text().splitlines()[1:]:
        name,size=line.split('\t')
        if name!='samples.txt' and int(Path(name).name.split('.')[0][3:]) not in chroms: continue
        path=panel if name=='samples.txt' else Path(os.environ['PFILE_DIR'])/Path(name).name
        if path.stat().st_size!=int(size): raise SystemExit(f'Input size mismatch: {path}')
print(f'Input preflight PASS: GRCh{os.environ["GRCH"]}, chromosomes {chroms}')
PY
"$PLINK2" --version
"$PLINK2" --help r2-phased indep-pairphase >/dev/null
scratch=$(mktemp -d)
trap 'rm -rf -- "$scratch"' EXIT
printf '#include <zstd.h>\n#include <filesystem>\nint main(){return ZSTD_versionNumber() < 10000;}\n' > "$scratch/check.cpp"
"$CXX" -std=c++17 -fopenmp "$scratch/check.cpp" -lzstd -o "$scratch/check"
"$scratch/check"
mkdir -p "$OUT" "$LD_TMP_ROOT"
test -w "$OUT" && test -w "$LD_TMP_ROOT"
echo 'Preflight PASS (input files, PLINK2, C++17 and zstd).'
