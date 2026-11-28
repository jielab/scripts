"""Dependency isolation and real device execution checks; no participant data."""
import argparse
import importlib
import importlib.util
from pathlib import Path
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--prefix',type=Path)
    p.add_argument('--expected-wheel',choices=['cpu','cu118','cu128'])
    p.add_argument('--expected-version')
    p.add_argument('--require-gpu',action='store_true')
    a=p.parse_args()
    if a.prefix and Path(sys.prefix).resolve()!=a.prefix.resolve():
        raise RuntimeError(f'Wrong Python prefix: {sys.prefix}; expected {a.prefix}')
    if not (3,11)<=sys.version_info[:2]<(3,14):
        raise RuntimeError('FLUX requires Python 3.11–3.13; select a compatible existing environment')
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from dependencies import requirements, unsatisfied
    problems=unsatisfied(requirements())
    if problems:
        raise RuntimeError('Missing/incompatible dependencies: '+', '.join(problems)+'. Run ./install.sh run first.')
    names=('numpy','pandas','scipy','sklearn','yaml','joblib','lightgbm','matplotlib','torch')
    missing=[n for n in names if importlib.util.find_spec(n) is None]
    if missing:
        raise RuntimeError('Missing packages: '+', '.join(missing)+'. Run ./install.sh run first.')
    if a.preflight:
        return
    for name in names:
        module=importlib.import_module(name)
        print(f'{name}: {getattr(module, "__version__", "available")}')
    import torch
    if a.expected_version and torch.__version__.split('+')[0]!=a.expected_version:
        raise RuntimeError('Installed PyTorch version differs from the installer pin')
    if a.expected_wheel:
        expected={'cpu':None,'cu118':'11.8','cu128':'12.8'}[a.expected_wheel]
        if torch.version.cuda!=expected:
            raise RuntimeError(f'Expected {a.expected_wheel}; got CUDA {torch.version.cuda}')
    print('Python:',sys.executable)
    print('Environment:',sys.prefix)
    print('CUDA wheel:',torch.version.cuda)
    available=torch.cuda.is_available()
    if a.require_gpu and not available:
        raise RuntimeError('GPU requested but CUDA is unavailable; check the NVIDIA driver and GPU allocation')
    devices=['cpu']+[f'cuda:{i}' for i in range(torch.cuda.device_count())] if available else ['cpu']
    for device in devices:
        x=torch.ones((16,16),device=device,requires_grad=True)
        (x@x).sum().backward()
        if not torch.isfinite(x.grad).all().item():
            raise RuntimeError(f'Nonfinite device smoke test: {device}')
        if device!='cpu':
            torch.cuda.synchronize(device)
            print(device,torch.cuda.get_device_name(device),torch.cuda.get_device_capability(device))
        print(f'{device}: forward/backward PASS')
    if torch.version.cuda and not available:
        print('GPU wheel installed; no GPU visible here. GPU execution must be checked on an allocated node.')


if __name__=='__main__':
    try:
        main()
    except (ImportError,RuntimeError,OSError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        sys.exit(2)
