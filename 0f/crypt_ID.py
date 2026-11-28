#!/usr/bin/env python3
"""Recreate custom column extracts from the encrypted, prepared phenotype table."""
import argparse
import contextlib
import gzip
import json
import os
from pathlib import Path
import tempfile


def rebuild(source, destination, definitions, overwrite=False):
    destination.mkdir(parents=True, exist_ok=True)
    staged = []
    try:
        with gzip.open(source, 'rt', newline='') as inp, contextlib.ExitStack() as stack:
            header = inp.readline().rstrip('\r\n').split('\t')
            positions = {name: i for i, name in enumerate(header)}
            outputs = []
            for filename, columns in definitions.items():
                target = destination / filename
                if target.exists() and not overwrite:
                    raise FileExistsError(f'Refusing to overwrite {target}')
                indices = [positions[name] for name in columns]
                fd, tmp = tempfile.mkstemp(prefix='.' + filename, dir=destination)
                os.close(fd)
                staged.append((Path(tmp), target))
                out = stack.enter_context(gzip.open(tmp, 'wt', newline=''))
                out.write('\t'.join(columns) + '\n')
                outputs.append((out, indices))
            count = 0
            for number, line in enumerate(inp, 2):
                row = line.rstrip('\r\n').split('\t')
                if len(row) != len(header):
                    raise ValueError(f'Column count mismatch at line {number}')
                for out, indices in outputs:
                    out.write('\t'.join(row[i] for i in indices) + '\n')
                count += 1
        for tmp, target in staged:
            if overwrite:
                # Publish each completed gzip atomically; retain old files on read/write errors.
                os.replace(tmp, target)
            else:
                os.link(tmp, target)
        print(f'Rebuilt {len(staged)} custom extracts; {count} rows each.')
    finally:
        for tmp, _ in staged:
            tmp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--file-in', type=Path, required=True,
                        help='prepared gzip TSV with encrypted IDs')
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='directory for reconstructed gzip TSV files')
    parser.add_argument('--columns', type=Path, required=True,
                        help='JSON object mapping output filenames to ordered column names')
    parser.add_argument('--overwrite', action='store_true',
                        help='replace existing outputs after all extracts are written successfully')
    args = parser.parse_args()
    definitions = json.loads(args.columns.read_text())
    if not isinstance(definitions, dict) or not definitions:
        parser.error('--columns must contain a nonempty JSON object')
    for filename, columns in definitions.items():
        if Path(filename).name != filename or not filename.endswith('.gz'):
            parser.error('output filenames must be simple .gz names without directories')
        if not isinstance(columns, list) or not columns or not all(isinstance(c, str) for c in columns):
            parser.error('each output must specify a nonempty list of column names')
        if (args.output_dir / filename).resolve() == args.file_in.resolve():
            parser.error('output must differ from input')
    rebuild(args.file_in, args.output_dir, definitions, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
