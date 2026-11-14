#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Diet dependencies

Main options:
    --check          Check R/Python dependencies [default]
    --install        Install missing R packages and Python requirements
    --with-optional  Include optional R helpers and headless MRI rendering checks

DIET_RSCRIPT and DIET_PYTHON override the configured runtimes.
Python MRI dependencies are included in require.txt.

Examples:
    ./install.sh --check
    ./install.sh --install
PUBLIC_ENTRY_HELP
    exit 0
fi
# Explicit environment preparation. The default is a read-only dependency check.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODE=--check
OPTIONAL=()
for arg in "$@"; do
    case "$arg" in
        --check|--install) MODE="$arg" ;;
        --with-optional) OPTIONAL=(--with-optional) ;;
        -h|--help)
            bash "${BASH_SOURCE[0]}" --help
            exit 0 ;;
        *) printf 'Unknown option: %s\n' "$arg" >&2; exit 2 ;;
    esac
done
HELPERS=${DIET_HELPERS:-$(dirname -- "$HERE")/0f}
source "$HELPERS/0os.sh"
source "$HERE/f/runtime_env.sh"
RSCRIPT="${DIET_RSCRIPT:-Rscript}"
PYTHON="${DIET_PYTHON:-python3}"
for interpreter in "$RSCRIPT" "$PYTHON"; do
    if ! command -v "$interpreter" >/dev/null 2>&1; then
        printf 'Interpreter not found: %s. Activate the environment described in INSTALL.md.\n' "$interpreter" >&2
        exit 1
    fi
done
export DIET_SOURCE_CACHE=${DIET_SOURCE_CACHE:-$DIR_SOFTWARE/diet/source-cache}
STATUS=0
if ! "$RSCRIPT" "$HERE/f/install.R" "$MODE" "${OPTIONAL[@]}"; then STATUS=1; fi
REQ="$HERE/require.txt"
if [[ "$MODE" == --install ]]; then
    mkdir -p "$DIET_PYTHON_LIBRARY"
    enigma_archive=${DIET_ENIGMA_ARCHIVE:-$DIR_SOFTWARE/diet/source-cache/enigma.tar.gz}
    if [[ -s $enigma_archive ]]; then
        cached_req=$(mktemp)
        trap 'rm -f -- "$cached_req"' EXIT
        "$PYTHON" - "$REQ" "$cached_req" "$enigma_archive" <<'CACHE'
from pathlib import Path
import sys
source, target, archive = map(Path, sys.argv[1:])
lines = source.read_text().splitlines()
target.write_text("\n".join("enigmatoolbox @ " + archive.resolve().as_uri() if line.startswith("enigmatoolbox @") else line for line in lines) + "\n")
CACHE
        REQ=$cached_req
    fi
    pip_options=()
    if [[ -n ${DIET_WHEELHOUSE:-} ]]; then
        [[ -d $DIET_WHEELHOUSE ]] || { echo "Missing wheelhouse: $DIET_WHEELHOUSE" >&2; exit 2; }
        pip_options=(--no-index --find-links "$DIET_WHEELHOUSE" --no-build-isolation)
    fi
    if ! "$PYTHON" -m pip install --disable-pip-version-check --upgrade --target "$DIET_PYTHON_LIBRARY" "${pip_options[@]}" -r "$REQ"; then STATUS=1; fi
fi
# Metadata/specifier checks avoid executing analysis code or heavy MRI imports.
# Explicit --install remains the only path that invokes a package installer.
if ! "$PYTHON" - "$HERE" "${OPTIONAL[@]}" <<'PY'
import importlib.metadata as md
import importlib.util
import json
from pathlib import Path
import sys

try:
    from packaging.requirements import Requirement
except ImportError:
    print("Missing Python package packaging; run diet/install.sh --install.")
    raise SystemExit(1)
if sys.version_info < (3, 10):
    print("Use Python >= 3.10; env.yml selects Python 3.11.")
    raise SystemExit(1)
root = Path(sys.argv[1])
files = [root / "require.txt"]
imports = {"pillow": "PIL", "pymupdf": "fitz", "scikit-learn": "sklearn", "scikit-image": "skimage"}
failures, count = [], 0
for path in files:
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "-r ")):
            continue
        req = Requirement(line)
        if req.marker and not req.marker.evaluate():
            continue
        count += 1
        try:
            dist = md.distribution(req.name)
            module = imports.get(req.name.lower(), req.name.replace("-", "_"))
            if importlib.util.find_spec(module) is None:
                failures.append(f"{req.name}: installed distribution has no {module} module")
            elif req.specifier and not req.specifier.contains(dist.version):
                failures.append(f"{req.name} {dist.version}: outside {req.specifier}")
            if req.url:
                origin = json.loads(dist.read_text("direct_url.json") or "{}")
                actual = origin.get("vcs_info", {}).get("commit_id")
                expected = req.url.rsplit("@", 1)[-1]
                if actual != expected:
                    print(f"  {req.name}: installed revision {actual or 'unknown'}; reference {expected}")
        except md.PackageNotFoundError:
            failures.append(f"{req.name}: missing")
print(f"Python {sys.version.split()[0]}: checked {count} distributions; {len(failures)} issues.")
for failure in failures:
    print("  " + failure)
if failures:
    print("Run diet/install.sh --install with the same optional flag in the intended environment.")
raise SystemExit(bool(failures))
PY
then STATUS=1; fi
if ((${#OPTIONAL[@]})); then
    if ! "$PYTHON" "$HERE/f/tests/check_offscreen.py"; then
        echo 'Headless VTK rendering failed; prepare Mesa/EGL or OSMesa before full MRI analysis.' >&2
        STATUS=1
    fi
fi
exit "$STATUS"
