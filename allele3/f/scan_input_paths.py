"""Resolve explicit copies of scan masks without changing scientific manifests.

An override relocates an existing mask; it cannot add or replace a filter.
Verification uses the original size/edge fingerprint, not a full-file checksum.
"""
from pathlib import Path

MASK_OPTIONS = {"lcr": "architecture_lcr_bed", "segdup": "architecture_segdup_bed",
                "low_mappability": "architecture_low_mappability_bed"}


def resolve_mask_options(run, overrides=None, *, verify=False):
    """Copy runtime options; optionally check relocated masks against the scan.

    Preflight validates the saved provenance without opening mask data. The
    candidate worker verifies each explicit override before reading any sites.
    Modification times may differ after transferring the same file to HPC.
    """
    options = dict(run.get("options", run.get("identity", {}).get("options", {})))
    records = (options.get("architecture_annotations") or {}).get("files", {})
    for name, destination in (overrides or {}).items():
        if not destination:
            continue
        if name not in MASK_OPTIONS:
            raise ValueError(f"Unknown scan mask override: {name}")
        key = MASK_OPTIONS[name]
        original = options.get(key)
        expected = records.get(key)
        flag = "--" + name.replace("_", "-") + "-bed"
        if not original:
            raise ValueError(f"{flag} cannot add a mask absent from the original scan")
        if (not isinstance(expected, dict) or expected.get("path") != str(original)
                or not isinstance(expected.get("size"), int) or not expected.get("edge_sha256")):
            raise ValueError(f"{flag} requires the original scan's matching size/edge fingerprint")
        path = Path(destination).expanduser().resolve()
        if verify:
            from scan_store import fingerprint, same_content
            if not path.is_file():
                raise ValueError(f"Relocated scan mask is missing: {path}")
            if not same_content(fingerprint(path), expected):
                raise ValueError(f"Relocated scan mask differs from original size/edge fingerprint: {path}")
        options[key] = str(path)
    return options
