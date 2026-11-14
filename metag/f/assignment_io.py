"""Small helpers for read-assignments.tsv[.gz] stored directly or in species.tar."""
from __future__ import annotations
from contextlib import contextmanager
import gzip
import io
from pathlib import Path
import tarfile


@contextmanager
def open_assignment_text(path):
    path = Path(path)
    if path.is_file() and path.suffix == ".tar":
        with tarfile.open(path, "r") as tar:
            matches = [m for m in tar.getmembers()
                       if m.isfile() and Path(m.name).name == "read-assignments.tsv.gz"]
            if len(matches) != 1:
                raise RuntimeError(f"Expected exactly one read-assignments.tsv.gz in {path}; found {len(matches)}")
            raw = tar.extractfile(matches[0])
            if raw is None:
                raise RuntimeError(f"Cannot read {matches[0].name} from {path}")
            with gzip.GzipFile(fileobj=raw, mode="rb") as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                    yield text
        return
    if not path.is_file():
        raise RuntimeError(f"Assignment source not found: {path}")
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield handle


def find_assignment_source(output_dir):
    output_dir = Path(output_dir)
    path = output_dir / "qc" / "read-assignments.tsv.gz"
    if path.is_file():
        return path
    raise RuntimeError(f"No qc/read-assignments.tsv.gz below {output_dir}")
