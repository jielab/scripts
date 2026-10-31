"""Resolve user-managed NCBI metadata without network access."""
import shutil
import tarfile
import tempfile
from pathlib import Path

from common import DIR_META, DIR_LIB, DIR_OUT

TAXONOMY_FILES = ("nodes.dmp", "names.dmp", "merged.dmp")


def summary_candidates(group):
    metadata = DIR_META
    return (
        metadata / "ncbi" / f"assembly_summary.{group}.txt",
        metadata / "kraken2_broad_v1" / "library" / group / "assembly_summary.txt",
        metadata / "reference_size_plan" / f"{group}.assembly_summary.txt",
        metadata / "kraken2_broad_v1" / "metag_reference_selection" / f"{group}.assembly_summary.txt",
        DIR_LIB / "reference" / f"assembly_summary.{group}.txt",
        DIR_OUT / "reference" / f"assembly_summary.{group}.txt",
        metadata / "kraken2_broad_v1" / "metag_reference_selection" / "originals" / f"{group}.assembly_summary.txt",
    )


def summary_path(group):
    candidates = summary_candidates(group)
    return next((path for path in candidates if path.is_file()), candidates[0])


def taxonomy_path():
    metadata = DIR_META
    candidates = (
        metadata / "ncbi" / "taxonomy",
        metadata / "kraken2_broad_v1" / "taxonomy",
        DIR_LIB / "metacache" / "taxonomy",
    )
    # Prefer a freshly downloaded archive over an older extracted taxonomy.
    return next((path for path in candidates if any(
        (path/name).is_file() for name in (*TAXONOMY_FILES, "taxdump.tar.gz")
    )), candidates[0])


def prepare_taxonomy(include_deleted=False, force=False):
    """Extract local taxdump when outputs are missing; never download anything."""
    source = taxonomy_path()
    required = TAXONOMY_FILES + (("delnodes.dmp",) if include_deleted else ())
    missing = [name for name in required if not (source/name).is_file()]
    if not missing and not force:
        return source
    archive = source/"taxdump.tar.gz"
    if not archive.is_file():
        return source  # The caller lists the missing local metadata files.
    if archive.with_name(archive.name + ".aria2").exists():
        raise RuntimeError(f"Taxonomy download is incomplete: {archive}; finish aria2c before building.")
    print(f"[metadata] Extract local {archive}: {', '.join(required)}", flush=True)
    try:
        # Stage all required files before publishing; a failed extraction must not
        # leave partial outputs that a later build would mistake for complete.
        with tempfile.TemporaryDirectory(prefix=".taxonomy-extract-", dir=source) as temporary:
            staging = Path(temporary)
            extracted = set()
            with archive.open("rb", buffering=8 << 20) as raw:
                with tarfile.open(fileobj=raw, mode="r|gz", bufsize=1 << 20) as bundle:
                    for member in bundle:
                        if member.name not in required:
                            continue
                        if not member.isfile() or member.name in extracted:
                            raise RuntimeError(f"Invalid taxonomy archive member: {member.name}")
                        with bundle.extractfile(member) as data, (staging/member.name).open("wb") as target:
                            shutil.copyfileobj(data, target, length=1 << 20)
                        extracted.add(member.name)
            absent = set(required) - extracted
            if absent:
                raise RuntimeError("Taxonomy archive is missing: " + ", ".join(sorted(absent)))
            # Publish missing files last: an interrupted replacement still leaves
            # a missing output, so the next build retries the complete extraction.
            order = [name for name in required if name not in missing] + missing
            for name in order:
                (staging/name).replace(source/name)
    except (tarfile.TarError, EOFError, OSError, RuntimeError) as exc:
        raise RuntimeError(f"Cannot extract local taxonomy {archive}: {exc}") from exc
    print(f"PASS: taxonomy extracted to {source}", flush=True)
    return source
