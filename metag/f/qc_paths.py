"""Paths used by independent QC readers; production layout stays unchanged."""
from pathlib import Path
from storage_layout import separated

__all__ = ["configured_split_root"]


def configured_split_root(cfg):
    """Result root from a user config; a frozen split manifest already stores it."""
    root = Path(cfg["output_dir"])
    return root if separated(cfg) else root / "kraken"
