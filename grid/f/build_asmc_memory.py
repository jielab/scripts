"""Build an isolated ASMC 1.4.0 backend with discardable HMM workspace pages.

Usage: grid-python build_asmc_memory.py --source /path/to/ASMC-1.4.0
Requires CMake, C++17, Boost and zlib. CMake fetches the pinned upstream deps;
repeat --cmake-arg=... to supply local FETCHCONTENT_SOURCE_DIR_* overrides.
The installed ASMC package is never modified.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from argneedle_memory import backend_path

PATCH = r'''
// GU: results have been copied out before this is called. These two workspaces
// are scratch, fully rewritten on the next forward/backward pass. Discard only
// whole pages strictly inside each Eigen allocation (never allocator metadata).
void HMM::releaseWorkspacePages()
{
  const auto page = static_cast<std::uintptr_t>(::sysconf(_SC_PAGESIZE));
  if (page == 0 || page > 1024 * 1024) {
    throw std::runtime_error("Cannot determine workspace page size");
  }
  for (auto* buffer : {&m_alphaBuffer, &m_betaBuffer}) {
    const auto address = reinterpret_cast<std::uintptr_t>(buffer->data());
    const auto begin = ((address + page - 1) / page) * page;
    const auto end = ((address + buffer->size() * sizeof(float)) / page) * page;
    if (end > begin && ::madvise(reinterpret_cast<void*>(begin), end - begin, MADV_DONTNEED) != 0) {
      throw std::runtime_error("Failed to discard ASMC scratch workspace pages");
    }
  }
}
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--cmake-arg', action='append', default=[])
    a = p.parse_args()
    if sys.platform != 'linux':
        p.error('This backend requires Linux madvise (WSL is supported)')
    if 'project(asmc LANGUAGES CXX VERSION 1.4.0)' not in (a.source / 'CMakeLists.txt').read_text():
        p.error('Expected upstream ASMC v1.4.0 source')
    target = backend_path()
    work = target.parent / 'build-source'
    work.mkdir(parents=True, exist_ok=True)
    shutil.copytree(a.source, work, dirs_exist_ok=True)
    def change(name, old, new):
        path = work / 'src' / name
        source = path.read_text()
        if source.count(old) != 1:
            raise RuntimeError(f'Unexpected upstream source: {name}: {old!r}')
        path.write_text(source.replace(old, new, 1))
    change('HMM.hpp', 'DecodePairsReturnStruct& getDecodePairsReturnStruct();',
           'DecodePairsReturnStruct& getDecodePairsReturnStruct();\n  void releaseWorkspacePages();')
    change('HMM.cpp', '#include "HMM.hpp"',
           '#include "HMM.hpp"\n#include <cstdint>\n#include <sys/mman.h>\n#include <unistd.h>')
    change('HMM.cpp', 'void HMM::finishDecoding()', PATCH + '\nvoid HMM::finishDecoding()')
    change('ASMC.cpp', 'mHmm.getDecodePairsReturnStruct().finaliseCalculations();',
           'mHmm.getDecodePairsReturnStruct().finaliseCalculations();\n  mHmm.releaseWorkspacePages();')
    # chr1: 2,054,102 sites * 69 states * 16 lanes exceeds signed int32.
    simd_path = work / 'src/Simd.cpp'
    simd = simd_path.read_text()
    for old, new, count in [
        ('for (int pos = from; pos < to; ++pos)',
         'for (Eigen::Index pos = from; pos < to; ++pos)', 4),
        ('const int ind = ', 'const Eigen::Index ind = ', 3),
    ]:
        if simd.count(old) != count:
            raise RuntimeError('Unexpected SIMD index implementation')
        simd = simd.replace(old, new)
    simd_path.write_text(simd)
    change('HMM.cpp', 'm_alphaBuffer.resize(sequenceLength * states * m_batchSize);',
           'm_alphaBuffer.resize(static_cast<Eigen::Index>(sequenceLength) * states * m_batchSize);')
    change('HMM.cpp', 'm_betaBuffer.resize(sequenceLength * states * m_batchSize);',
           'm_betaBuffer.resize(static_cast<Eigen::Index>(sequenceLength) * states * m_batchSize);')
    change('pybind.cpp', 'PYBIND11_MODULE(asmc_python_bindings, m)\n{',
           'PYBIND11_MODULE(asmc_python_bindings, m)\n{\n  m.attr("_gu_bounded_workspace") = 2;')
    test = Path(__file__).resolve().parent / 'tests/test_asmc_large_index.cpp'
    with (work / 'CMakeLists.txt').open('a') as cmake:
        cmake.write(f'\nadd_executable(gu_asmc_large_index "{test}")\n'
                    'target_link_libraries(gu_asmc_large_index PRIVATE ASMC)\n')
    build = target.parent / 'build'
    subprocess.run(['cmake', '-S', str(work), '-B', str(build),
                    '-DASMC_PYTHON_BINDINGS=ON', '-DASMC_TESTING=OFF',
                    '-DCMAKE_BUILD_TYPE=Release', '-DPYTHON_EXECUTABLE=' + sys.executable,
                    '-DPython_EXECUTABLE=' + sys.executable,
                    '-DCMAKE_PREFIX_PATH=' + sys.prefix, *a.cmake_arg], check=True)
    subprocess.run(['cmake', '--build', str(build), '--target', 'asmc_python_bindings', 'gu_asmc_large_index',
                    '-j', str(a.jobs)], check=True)
    subprocess.run([str(build / 'gu_asmc_large_index')], check=True)
    artifact, = build.glob('asmc_python_bindings*.so')
    shutil.copy2(artifact, target.with_suffix('.next'))
    target.with_suffix('.next').replace(target)
    (target.parent / 'build.json').write_text(json.dumps({
        'upstream': 'https://github.com/PalamaraLab/ASMC/tree/v1.4.0',
        'patch': 2, 'python': sys.version, 'artifact': str(target),
    }, indent=2) + '\n')
    print(target)


if __name__ == '__main__':
    main()
