#!/usr/bin/env bash
set -eo pipefail
: "${LSB_JOBID:?Submit with bsub -q short}"
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$ROOT/f/hpc.runtime.env"
source "$ROOT/f/runtime.env"
export PIP_NO_INDEX=1 CMAKE_BUILD_PARALLEL_LEVEL=8
export CMAKE_PREFIX_PATH=$CONDA_PREFIX
export CMAKE_ARGS="-DFETCHCONTENT_SOURCE_DIR_PYBIND11=$GU_SOFT/sources/deps/pybind11-2.13.5"
python -m pip install --no-build-isolation --no-deps "$GU_SOFT/sources/arg-needle-lib-1.2.0"
export CMAKE_ARGS="-DFETCHCONTENT_SOURCE_DIR_EIGEN=$GU_SOFT/sources/deps/eigen -DFETCHCONTENT_SOURCE_DIR_FMT=$GU_SOFT/sources/deps/fmt -DFETCHCONTENT_SOURCE_DIR_RANGE-V3=$GU_SOFT/sources/deps/range-v3 -DFETCHCONTENT_SOURCE_DIR_HIGHWAY=$GU_SOFT/sources/deps/highway -DFETCHCONTENT_SOURCE_DIR_PYBIND11=$GU_SOFT/sources/deps/pybind11"
python -m pip install --no-build-isolation --no-deps "$GU_SOFT/sources/ASMC-1.4.0"
python -m pip install --no-build-isolation --no-deps "$GU_SOFT/sources/arg-needle-1.1.0"
python f/0.arg.py build-asmc-memory --source "$GU_SOFT/sources/ASMC-1.4.0" --jobs 8 \
  --cmake-arg="-DFETCHCONTENT_SOURCE_DIR_EIGEN=$GU_SOFT/sources/deps/eigen" \
  --cmake-arg="-DFETCHCONTENT_SOURCE_DIR_FMT=$GU_SOFT/sources/deps/fmt" \
  --cmake-arg="-DFETCHCONTENT_SOURCE_DIR_RANGE-V3=$GU_SOFT/sources/deps/range-v3" \
  --cmake-arg="-DFETCHCONTENT_SOURCE_DIR_HIGHWAY=$GU_SOFT/sources/deps/highway" \
  --cmake-arg="-DFETCHCONTENT_SOURCE_DIR_PYBIND11=$GU_SOFT/sources/deps/pybind11"
export CMAKE_ARGS="-DFETCHCONTENT_SOURCE_DIR_PYBIND11=$GU_SOFT/sources/deps/pybind11-2.13.5"
"$GU_THREADS_ENV/bin/python" -m pip install --no-build-isolation --no-deps "$GU_SOFT/sources/arg-needle-lib-1.1.3"
export CMAKE_ARGS="-DFETCHCONTENT_SOURCE_DIR_PYBIND11=$GU_SOFT/sources/deps/pybind11-2.13.6"
"$GU_THREADS_ENV/bin/python" -m pip install --no-build-isolation --no-deps "$GU_SOFT/sources/threads-0.2.1"
python -m pip check
"$GU_THREADS_ENV/bin/python" -m pip check
"$GU_THREADS_ENV/bin/threads" infer --help >/dev/null
python -c 'import arg_needle, arg_needle_lib, asmc'
"$GU_SINGER_ENV/bin/python" "$GU_SINGER_HOME/singer_master" --help >/dev/null
printf 'ARG software built and checked\n' > "$GU_SOFT/logs/arg.complete"
