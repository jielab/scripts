#!/usr/bin/env python3
"""Compatibility entry point; implementation lives in scripts/0f."""
from pathlib import Path
import runpy
import sys
_shared = next(p / '0f' for p in Path(__file__).resolve().parents
               if (p / '0f' / 'process_lifetime.py').is_file())
sys.path.insert(0, str(_shared))
if __name__ == '__main__':
    runpy.run_path(str(_shared / 'process_lifetime.py'), run_name='__main__')
else:
    globals().update({k: v for k, v in runpy.run_path(str(_shared / 'process_lifetime.py')).items()
                      if not k.startswith('__')})
