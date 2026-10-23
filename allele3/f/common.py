"""Adapters reusing the original scanner's chromosome and I/O conventions."""
from __future__ import annotations
import math
try:
    from .core import open_text
    from .gnomad_io import canonical_chrom, chrom_key as chrom_sort_key
except ImportError:
    from core import open_text
    from gnomad_io import canonical_chrom, chrom_key as chrom_sort_key

def af_bin(x):
    try:
        v = float(x)
    except (ValueError, TypeError):
        return 'af_missing'
    if not math.isfinite(v) or not 0 <= v <= 1:
        return 'af_missing'
    for edge, name in [(1e-5,'lt_1e-5'),(1e-4,'1e-5_to_1e-4'),(1e-3,'1e-4_to_1e-3'),(.01,'0.001_to_0.01'),(.05,'0.01_to_0.05'),(.5,'0.05_to_0.5')]:
        if v < edge: return name
    return 'ge_0.5'

def an_bin(x):
    try:
        v = float(x)
        if not math.isfinite(v) or v < 0: return 'an_missing'
    except (ValueError, TypeError):
        return 'an_missing'
    for edge in (100,1000,10000,50000,100000,120000,140000,150000,200000):
        if v < edge: return 'lt_' + str(edge)
    return 'ge_200000'

def ac_tier(x):
    try: v=float(x)
    except (ValueError,TypeError): return 'ac_missing'
    if not math.isfinite(v) or v < 0: return 'ac_missing'
    if v == 0: return '0'
    if v == 1: return '1'
    if v < 5: return '2-4'
    if v < 10: return '5-9'
    if v < 100: return '10-99'
    return '100+'
