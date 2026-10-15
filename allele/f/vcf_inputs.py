"""Resolve one cohort's VCFs and reject known incomplete downloads before reading."""

import json
from pathlib import Path
import re

GNOMAD_PREFIX = "gnomad.genomes.v4.1.sites.chr"
AUTOSOMES = [str(i) for i in range(1, 23)]
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def chromosomes(value):
    if value == "autosomes":
        return AUTOSOMES
    if value == "all":
        return AUTOSOMES + ["X", "Y"]
    result = value.split(",")
    if not result or any(x not in AUTOSOMES + ["X", "Y"] for x in result):
        raise ValueError("--chr must be autosomes, all, or a comma list from 1..22,X,Y")
    return list(dict.fromkeys(result))


def chrom_key(c):
    x = c.removeprefix("chr")
    if x.isdigit():
        return (0, int(x), "")
    return (1, {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(x, 99), x)


def file_chrom(path):
    match = re.search(r"(?:^|\.)chr([0-9]+|X|Y|MT|M)(?:\.|$)", Path(path).name)
    return match[1] if match else Path(path).name


def resolve_inputs(vcf=None, directory=None, chrs=None):
    if vcf:
        if chrs:
            raise ValueError("--chr is for --vcf-dir; use --region with a single --vcf")
        return [Path(vcf)]
    folder = Path(directory)
    files = sorted(
        {p for pat in ("*.vcf.bgz", "*.vcf.gz", "*.vcf") for p in folder.glob(pat)}
    )
    gnomad = [p for p in files if p.name.startswith(GNOMAD_PREFIX)]
    if gnomad or folder == Path("/mnt/f/annot/gnomad"):
        # Require the entire requested set; never silently analyze only downloaded chromosomes.
        chosen = chromosomes(chrs or "autosomes")
        files = [folder / f"{GNOMAD_PREFIX}{c}.vcf.bgz" for c in chosen]
    elif chrs:
        chosen = chromosomes(chrs)
        files = [p for p in files if file_chrom(p) in chosen]
        absent = set(chosen) - {file_chrom(p) for p in files}
        if absent:
            raise ValueError(
                "Missing requested chromosomes: "
                + ",".join(sorted(absent, key=chrom_key))
            )
    if not files:
        raise ValueError(f"No VCF files found in {directory}")
    return sorted(files, key=lambda p: chrom_key(file_chrom(p)))


def validate_inputs(files, require_index=False):
    catalog = json.loads(
        Path(__file__).with_name("gnomad_v4.1_genomes.json").read_text()
    )["files"]
    errors = []
    for path in files:
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"Missing VCF: {path}")
            continue
        indexed = require_index or path.name in catalog
        indexes = [Path(str(path) + ext) for ext in (".tbi", ".csi")]
        if indexed and not any(p.is_file() and p.stat().st_size for p in indexes):
            errors.append(f"VCF requires a VCF .tbi or .csi index: {path}")
        for part in [path] + [p for p in indexes if p.exists()]:
            if Path(str(part) + ".aria2").exists():
                errors.append(f"Download incomplete (.aria2): {part}")
                continue
            expected = catalog.get(part.name)
            if expected and part.stat().st_size != expected["bytes"]:
                errors.append(
                    f"Incomplete/wrong size: {part} (expected {expected['bytes']} bytes)"
                )
                continue
            if part == path and str(path).endswith((".gz", ".bgz")):
                with path.open("rb") as handle:
                    handle.seek(max(0, path.stat().st_size - len(BGZF_EOF)))
                    if handle.read() != BGZF_EOF:
                        errors.append(
                            f"Missing BGZF end marker (truncated or not bgzip): {path}"
                        )
    if errors:
        raise ValueError("Input check failed:\n  " + "\n  ".join(errors))
