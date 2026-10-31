#!/usr/bin/env python3
"""Compatibility entry point for the compact full-reference split workflow."""
from stream_pipeline import cli, main, discover_pair_files

if __name__ == "__main__":
    raise SystemExit(cli())

