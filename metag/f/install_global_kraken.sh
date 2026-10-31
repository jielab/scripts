#!/usr/bin/env bash
# Install official pinned Kraken 2 separately from the existing environment.
set -euo pipefail
VERSION=2.17.2
COMMIT=8c190b1b668825935dbf6dee5f969227dc8269bb
TOOLS="${METAG_TOOLS:-$HOME/.cache/metag/tools}"
PREFIX="$TOOLS/kraken2/$VERSION"
THREAD_COUNT="${THREADS:-16}"
if [[ "${1:-}" == "--help" ]]; then
  echo "Usage: ./3.global.sh install [--threads N]"
  echo "Installs official Kraken 2 $VERSION ($COMMIT) under METAG_TOOLS or ~/.cache/metag/tools."
  exit 0
fi
if [[ "${1:-}" == "--threads" && $# == 2 ]]; then THREAD_COUNT="$2"; shift 2; fi
[[ $# == 0 && "$THREAD_COUNT" =~ ^[1-9][0-9]*$ ]] || { echo "Usage: ./3.global.sh install [--threads N]" >&2; exit 2; }
for program in git make g++ gcc perl python3 bash flock; do
  command -v "$program" >/dev/null || { echo "Required build tool missing: $program" >&2; exit 1; }
done
if [[ -e "$PREFIX" && ! -f "$PREFIX/.metag-global-owner" ]]; then
  echo "Refusing to alter unowned installation directory: $PREFIX" >&2
  echo "Choose a separate METAG_TOOLS directory." >&2
  exit 1
fi
mkdir -p "$PREFIX"
PREFIX="$(cd -- "$PREFIX" && pwd)"
touch "$PREFIX/.metag-global-owner"
exec 9>"$PREFIX/install.lock"
flock -n 9 || { echo "Another installer is running at $PREFIX" >&2; exit 1; }
if [[ -f "$PREFIX/install.json" ]]; then
  if python3 - "$PREFIX" "$COMMIT" <<'PY'
import hashlib, json, os, pathlib, subprocess, sys
root=pathlib.Path(sys.argv[1])
try:
    saved=json.loads((root/'install.json').read_text())
    assert saved['source_tag']=='2.17.2' and saved['source_commit']==sys.argv[2]
    assert {'classify','build_db','estimate_capacity','dump_table','k2mask','kraken2','kraken2-build','kraken2-inspect'} <= set(saved['files'])
    assert {p.name for p in (root/'bin').iterdir() if p.is_file()} == set(saved['files'])
    for name, expected in saved['files'].items():
        assert pathlib.Path(name).name==name
        assert hashlib.sha256((root/'bin'/name).read_bytes()).hexdigest()==expected
    for name in ('classify','build_db','estimate_capacity','dump_table','k2mask','kraken2','kraken2-build','kraken2-inspect'):
        assert os.access(root/'bin'/name,os.X_OK)
    assert subprocess.check_output([str(root/'bin/kraken2'),'--version'],text=True,timeout=15).splitlines()[0] == saved['reported_version']
except (OSError, ValueError, AssertionError, KeyError, subprocess.SubprocessError):
    raise SystemExit(1)
PY
  then
    echo "Kraken 2 source tag $VERSION is already installed and verified: $PREFIX/bin"
    exit 0
  fi
fi
if [[ ! -d "$PREFIX/src/.git" ]]; then
  git clone --branch "$VERSION" --depth 1 https://github.com/DerrickWood/kraken2.git "$PREFIX/src"
fi
ACTUAL="$(git -C "$PREFIX/src" rev-parse HEAD)"
[[ "$ACTUAL" == "$COMMIT" ]] || { echo "Source commit mismatch: $ACTUAL != $COMMIT" >&2; exit 1; }
if [[ -n "$(git -C "$PREFIX/src" status --porcelain --untracked-files=no)" ]]; then
  echo "Pinned source has local changes; refusing to build an unverified variant." >&2
  exit 1
fi
mkdir -p "$PREFIX/bin"
cd -- "$PREFIX/src"
export MAKEFLAGS="-j$THREAD_COUNT"
# Do not inherit an external request to disable OpenMP.
unset KRAKEN2_SKIP_FOPENMP || true
bash ./install_kraken2.sh "$PREFIX/bin"
for program in estimate_capacity build_db classify dump_table lookup_accession_numbers k2mask blast_to_fasta merge; do
  [[ -s "$PREFIX/bin/$program" ]] || { echo "Missing installed native: $program" >&2; exit 1; }
  chmod +x "$PREFIX/bin/$program"
done
for program in kraken2 kraken2-build kraken2-inspect; do
  chmod +x "$PREFIX/bin/$program"
  "$PREFIX/bin/$program" --version
done
"$PREFIX/bin/classify" -h >/dev/null
python3 - "$PREFIX" "$VERSION" "$COMMIT" <<'PY'
from datetime import datetime, timezone
import hashlib, json, os, pathlib, subprocess, sys
root=pathlib.Path(sys.argv[1]); folder=root/'bin'
report={
    'source_tag':sys.argv[2], 'source_commit':sys.argv[3],
    'reported_version':subprocess.check_output([str(folder/'kraken2'),'--version'],text=True).splitlines()[0],
    'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir()) if p.is_file()},
    'installed_utc':datetime.now(timezone.utc).isoformat(),
    'note':'Official source tag 2.17.2 retains VERSION=2.17.1; source commit and installed-file hashes are authoritative.'
}
temporary=root/'install.json.partial'
with temporary.open('w') as handle:
    json.dump(report,handle,ensure_ascii=False,indent=2); handle.write('\n')
    handle.flush(); os.fsync(handle.fileno())
temporary.replace(root/'install.json')
PY
echo "PASS: Kraken 2 installed at $PREFIX/bin; ./3.global.sh will select it automatically."
