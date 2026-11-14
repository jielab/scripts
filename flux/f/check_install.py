"""Dependency/device check, without accessing participant data."""
import importlib
import sys

for name in ('numpy', 'pandas', 'scipy', 'sklearn', 'yaml', 'joblib', 'lightgbm', 'torch'):
    module = importlib.import_module(name)
    print(f'{name}: {getattr(module, "__version__", "available")}')
import torch
print('Python:', sys.executable)
print('CUDA available:', torch.cuda.is_available())
print('Wheel CUDA:', torch.version.cuda)
print('Compiled GPU architectures:', torch.cuda.get_arch_list())
for i in range(torch.cuda.device_count()):
    print('GPU', i, torch.cuda.get_device_name(i), torch.cuda.get_device_capability(i))
