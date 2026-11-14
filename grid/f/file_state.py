"""Cheap output completeness and file records; no full-file content scans."""
from pathlib import Path
import json


def file_record(path):
    p = Path(path)
    stat = p.stat()
    return "file:" + str(stat.st_size)


def complete_files(paths):
    return all(Path(p).is_file() and Path(p).stat().st_size > 0 for p in paths)
