#!/usr/bin/env python3
"""Disco standalone outcome evaluation, independent of GRID learned models.

The pinned genotype/PCA interpolation remains in 3.disco.R. This module fits its
outcome-scale adapter using the shared small validation budget and reports the
full untouched test roster; it never borrows GRID build/tuning/calibration heads.
"""
from target_benchmark import main

if __name__ == "__main__":
    main(method="disco")
