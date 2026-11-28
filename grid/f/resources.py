"""Compatibility adapter for the shared QiMing resource policy."""
from pathlib import Path
import sys
_shared = next(p / '0f' for p in Path(__file__).resolve().parents
               if (p / '0f' / 'qiming_resources.py').is_file())
sys.path.insert(0, str(_shared))
from qiming_resources import *
if __name__ == '__main__':
    main()
