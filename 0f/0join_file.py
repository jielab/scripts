#!/usr/bin/env python3
"""Left-join two or more delimited files by configurable key columns."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import BinaryIO


SEPARATORS = {"COMMA": b",", "SPACE": b" ", "TAB": b"\t"}


def open_binary(path: Path, mode: str) -> BinaryIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def parse_spec(value: str) -> tuple[Path, bytes, int]:
    try:
        filename, separator_name, column_text = value.rsplit(",", 2)
        separator = SEPARATORS[separator_name.upper()]
        column = int(column_text)
    except (KeyError, ValueError) as error:
        raise argparse.ArgumentTypeError(
            f"invalid input specification: {value!r}; expected FILE,TAB|SPACE|COMMA,COLUMN"
        ) from error
    if column < 0:
        raise argparse.ArgumentTypeError("key column must be zero-based and non-negative")
    return Path(filename), separator, column


def split_fields(line: bytes, separator: bytes) -> list[bytes]:
    return line.rstrip(b"\r\n").split(separator)


def load_lookup(path: Path, separator: bytes, key_column: int) -> tuple[dict[bytes, bytes], bytes]:
    lookup: dict[bytes, bytes] = {}
    width: int | None = None
    with open_binary(path, "rb") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = split_fields(line, separator)
            if width is None:
                width = len(fields)
            if len(fields) != width or key_column >= len(fields) or not fields[key_column]:
                raise ValueError(f"{path}: inconsistent or missing fields at line {line_number}")
            lookup[fields[key_column]] = b" ".join(fields)
    if width is None:
        raise ValueError(f"{path}: empty input file")
    return lookup, b" ".join([b"NA"] * width)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        metavar="FILE,SEPARATOR,COLUMN",
        help="two or more input specifications; key columns are zero-based",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.specs = [parse_spec(item) for item in args.inputs]
    if len(args.specs) < 2:
        parser.error("--inputs requires at least two file specifications")
    return args


def main() -> None:
    args = parse_args()
    primary_path, primary_separator, primary_key = args.specs[0]
    lookups = [load_lookup(*spec) for spec in args.specs[1:]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_binary(primary_path, "rb") as source, open_binary(args.output, "wb") as output:
        primary_width: int | None = None
        for line_number, line in enumerate(source, start=1):
            fields = split_fields(line, primary_separator)
            if primary_width is None:
                primary_width = len(fields)
            if len(fields) != primary_width or primary_key >= len(fields) or not fields[primary_key]:
                raise ValueError(
                    f"{primary_path}: inconsistent or missing fields at line {line_number}"
                )
            key = fields[primary_key]
            joined = [b" ".join(fields)]
            joined.extend(lookup.get(key, missing) for lookup, missing in lookups)
            output.write(b" ".join(joined) + b"\n")


if __name__ == "__main__":
    main()
