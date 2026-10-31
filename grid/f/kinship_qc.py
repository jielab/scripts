"""Aggregate-only kinship QC export; no local IDs, paths or hashes."""
import json
from pathlib import Path
import pandas as pd
from kinship import active

COLUMNS = ["N_participants", "N_connected_components", "N_singleton_components", "N_non_singleton_components",
    "N_participants_in_non_singleton_components", "maximum_component_size", "N_test_retained",
    "N_development_with_test_relative_before_correction", "N_development_excluded", "N_development_retained",
    "N_test_with_development_relative_after_correction", "component_rule", "correction_status",
    "source_path_basename", "source_size", "source_mtime_ns"]


def public_qc(root):
    correction = active(root)
    if not correction or not correction.get("qc") or not Path(correction["qc"]).is_file():
        return pd.DataFrame(columns=COLUMNS)
    qc = json.loads(Path(correction["qc"]).read_text())
    return pd.DataFrame([{c:qc.get(c) for c in COLUMNS}], columns=COLUMNS)
