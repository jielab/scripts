#!/usr/bin/env python3
"""Read declarative GRID settings; JSON remains a format for result records."""
from __future__ import annotations

from pathlib import Path
import json
import yaml


def load_document(text):
    # JSON is a YAML subset; preserve JSON numeric types in result/test records.
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return yaml.safe_load(text)


def load_settings(path):
    value = load_document(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a configuration mapping: {path}")
    return value
