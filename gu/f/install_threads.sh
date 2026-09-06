#!/usr/bin/env bash
# Independent of the Needle and GU environments (Threads pins arg-needle-lib).
set -euo pipefail
runtime=${GU_THREADS_ENV:-$HOME/.venvs/gu-threads}
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
python=${GU_THREADS_PYTHON:-python3}
"$python" -c 'import sys; assert sys.version_info >= (3, 12), "Threads tested environment requires Python >=3.12; set GU_THREADS_PYTHON"'
"$python" -m venv "$runtime"
"$runtime/bin/python" -m pip install -r "$here/threads-tested-requirements.txt"
"$runtime/bin/threads" infer --help
