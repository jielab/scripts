"""Shared, read-only validation for current aggregate publications."""
import json
import re
from pathlib import Path

import pandas as pd

PAPER_FILES = (
    "paper_regional_dimensionality.tsv",
    "paper_main_performance.tsv", "paper_predictability.tsv", "paper_benefit.tsv",
    "paper_count_contrasts.tsv", "paper_key_results.tsv", "paper_method_status.tsv",
)


def aggregate(table):
    """Reject individual records before they can enter final or its catalogue."""
    from twin import IDENTIFIER_COLUMNS
    forbidden = IDENTIFIER_COLUMNS | {"y", "outcome", "prediction", "predictions",
                                     "pbase", "pcandidate", "matches", "donors"}
    bad = [str(c) for c in table if re.sub(r"[^a-z0-9]", "", str(c).lower()) in forbidden]
    if bad:
        raise ValueError("Individual-level columns cannot enter final: " + ", ".join(bad))
    return table


def read_aggregate(path):
    return aggregate(pd.read_csv(path, sep="\t"))


def current_abm_run(pointer):
    """Resolve exactly this CURRENT pointer; never search historical runs."""
    from abm_core.model import file_sha256
    pointer = Path(pointer)
    if not pointer.is_file():
        return None
    name = json.loads(pointer.read_text())["run"]
    if not isinstance(name, str) or name in ("", ".", "..") or Path(name).name != name:
        raise ValueError("ABM CURRENT run must stay inside its trait directory")
    run = pointer.parent / name
    if run.is_symlink():
        raise ValueError("ABM CURRENT cannot link to another run")
    if not all((run / name).is_file() for name in ("COMPLETE.json", "RUN_MANIFEST.json", "PHASES.json")):
        return None
    receipt = json.loads((run / "COMPLETE.json").read_text())
    if receipt.get("status") != "completed" or receipt.get("smoke"):
        return None
    if file_sha256(run / "RUN_MANIFEST.json") != receipt["manifest_sha256"]:
        raise ValueError("ABM run manifest changed")
    from kinship import active, compatible
    correction = active(pointer.parents[3])
    manifest = json.loads((run / "RUN_MANIFEST.json").read_text())
    if not compatible(manifest.get("metadata", {}), correction):
        return None
    if pointer.name == "blocknet.CURRENT.json" and manifest["contract"].get("representation") != "csx_ldblocks_fixed":
        return None
    phases = json.loads((run / "PHASES.json").read_text())["phases"]
    if not all(phases.get(phase) for phase in ("evaluate", "report")):
        return None
    for phase in ("evaluate", "report"):
        for name, expected in phases[phase].items():
            from runtime_state import artifact_matches
            if not artifact_matches(run / name, expected):
                raise ValueError("ABM report changed: " + name)
    return run
