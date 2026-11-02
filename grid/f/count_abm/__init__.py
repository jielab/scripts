"""Transparent GRID prediction from observed reference outcomes.

This pipeline is independent of the historical ``abm_core`` regression and
CountHeads workflow.  A genetic representation can be learned on training data;
the prediction here is an explicitly recorded count or arithmetic mean.
"""

SCHEMA = "grid.count_abm.1"
VERSION = "20261102.1"
