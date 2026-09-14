#!/usr/bin/env python3
"""Bind reusable genetic outputs to explicit source files and conversion options."""
import argparse
import json
from pathlib import Path


def stamp(name):
    p = Path(name).resolve(strict=True)
    s = p.stat()
    return [str(p), s.st_size, s.st_mtime_ns]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=['check', 'record'])
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--input', action='append', default=[])
    p.add_argument('--output', action='append', default=[])
    p.add_argument('--value', action='append', default=[])
    a = p.parse_args()
    try:
        data = dict(schema=1, inputs=[stamp(x) for x in a.input],
                    outputs=[stamp(x) for x in a.output], values=a.value)
        if not all(x[1] > 0 for x in data['outputs']):
            raise ValueError('empty output')
        if a.action == 'check':
            return 0 if json.loads(a.receipt.read_text()) == data else 1
        temp = a.receipt.with_name(a.receipt.name+'.next')
        try:
            temp.write_text(json.dumps(data, indent=2)+'\n')
            temp.replace(a.receipt)
        finally:
            temp.unlink(missing_ok=True)
        return 0
    except (OSError, ValueError):
        if a.action == 'check': return 1
        raise


if __name__ == '__main__':
    raise SystemExit(main())
