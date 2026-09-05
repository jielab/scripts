#!/usr/bin/env python3
"""Compatibility helpers for GU TSV files with very wide fields."""
from __future__ import annotations

import csv
import sys


def enable_wide_csv_fields() -> int:
    """Raise the process-wide CSV field limit to the platform maximum."""
    limit = sys.maxsize
    while limit > 0:
        try:
            csv.field_size_limit(limit)
            return limit
        except OverflowError:
            # Some Python builds expose a C long smaller than sys.maxsize.
            limit //= 10
    raise RuntimeError("unable to configure a usable CSV field-size limit")

