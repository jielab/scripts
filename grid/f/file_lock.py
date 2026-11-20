#!/usr/bin/env python3
"""Project adapter for the shared Linux OFD locking implementation."""
from __future__ import annotations

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / "0f"))
from workflow_lock import flock, file_lock, main


if __name__ == "__main__":
    raise SystemExit(main())
