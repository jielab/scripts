"""Resolve FLUX requirements without upgrading compatible shared packages."""
import argparse
import importlib.metadata as metadata
from pathlib import Path
import subprocess
import sys
import tempfile

try:
    from packaging.requirements import Requirement
    from packaging.utils import canonicalize_name
except ImportError:
    from pip._vendor.packaging.requirements import Requirement
    from pip._vendor.packaging.utils import canonicalize_name


def requirements():
    return [Requirement(line) for line in Path(__file__).with_name('require.txt').read_text().splitlines()
            if line.strip() and not line.startswith('#')] + [Requirement('torch>=2.7.1,<3')]


def unsatisfied(reqs):
    result = []
    for req in reqs:
        try:
            version = metadata.version(req.name)
        except metadata.PackageNotFoundError:
            version = None
        if version is None or not req.specifier.contains(version):
            result.append(str(req))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--install', action='store_true')
    p.add_argument('--wheel', required=True)
    a = p.parse_args()
    if not (3, 11) <= sys.version_info[:2] < (3, 14):
        raise RuntimeError('Select an existing Python 3.11–3.13 environment; Python will not be replaced.')
    reqs = requirements()
    missing = unsatisfied(reqs)
    for req in reqs:
        print(('NEEDS ' if str(req) in missing else 'KEEP  ') + str(req), flush=True)
    if not a.install or not missing:
        return
    # Constrain all other installed distributions, including torch's CUDA stack.
    # A resolver conflict stops installation rather than silently changing them.
    changing = {canonicalize_name(Requirement(r).name) for r in missing}
    with tempfile.TemporaryDirectory(prefix='flux-deps-') as folder:
        constraints = Path(folder) / 'constraints.txt'
        constraints.write_text('\n'.join(
            f'{d.metadata["Name"]}=={d.version}' for d in metadata.distributions()
            if d.metadata['Name'] and canonicalize_name(d.metadata['Name']) not in changing))
        cmd = [sys.executable, '-I', '-m', 'pip', 'install', '-c', str(constraints)]
        if 'torch' in changing:
            # Use PyTorch's wheel index only for torch; resolve the complete plan together.
            import json
            report = Path(folder) / 'torch.json'
            subprocess.run([sys.executable, '-I', '-m', 'pip', 'install', '--dry-run', '--ignore-installed',
                            '--no-deps', '--report', str(report), 'torch==2.7.1',
                            '--index-url', f'https://download.pytorch.org/whl/{a.wheel}'], check=True)
            url = json.loads(report.read_text())['install'][0]['download_info']['url']
            missing = [r for r in missing if Requirement(r).name != 'torch'] + [f'torch @ {url}']
        subprocess.run(cmd + missing, check=True)


if __name__ == '__main__':
    main()
