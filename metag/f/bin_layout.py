"""One shallow, disjoint FASTQ layout shared by producers and consumers."""
from pathlib import Path


def fastq_name(stem, mate=None):
    stem = str(stem)
    path = Path(stem)
    if path.is_absolute() or ".." in path.parts or not stem:
        raise ValueError("FASTQ bin stem must be a safe relative path")
    return f"{stem}.R{mate}.fastq.gz" if mate is not None else f"{stem}.SE.fastq.gz"


def fastq_paths(root, stem, paired=True):
    return [Path(root) / fastq_name(stem, mate) for mate in ((1, 2) if paired else (None,))]
