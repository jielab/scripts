"""Read-only validation of the Disco module's aggregate outcome comparison."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from final_outputs import read_aggregate
from runtime_state import artifact_matches
from kinship import active, compatible

METHOD = "DiscoDivas"


class InvalidPublication(ValueError):
    pass


def read_performance(root, trait):
    root = Path(root)
    folder = root / "disco" / trait
    try:
        receipt = json.loads((folder / "performance.json").read_text())
        if receipt.get("status") != "completed":
            return None, receipt.get("reason", "not_completed")
        if receipt.get("schema") != "grid.disco.performance.1" or receipt.get("trait") != trait:
            raise InvalidPublication("wrong schema or trait")
        correction = active(root)
        if not correction or not compatible(receipt, correction):
            raise InvalidPublication("missing or incompatible verified kinship")
        required_files = {str(folder / "performance.tsv"), str(folder / "_LOCAL_ONLY/roles.tsv.gz")}
        if not required_files <= set(receipt["files"]):
            raise InvalidPublication("missing bound aggregate/role table")
        for path, expected in receipt["files"].items():
            if not artifact_matches(path, expected):
                raise InvalidPublication("producer input/output changed")
        roster_path = root / "prsformer/cache/data.tsv.gz"
        if not artifact_matches(roster_path, receipt["roster_identity"]):
            raise InvalidPublication("frozen roster changed")
        roles = pd.read_csv(folder / "_LOCAL_ONLY/roles.tsv.gz", sep="\t", dtype=str)
        if roles.id.duplicated().any():
            raise InvalidPublication("duplicate roster IDs")
        roster = pd.read_csv(roster_path, sep="\t", usecols=["eid", "split"], dtype=str)
        if set(roles.loc[roles.role.eq("test"), "id"]) != set(roster.loc[roster.split.eq("test"), "eid"]):
            raise InvalidPublication("incomplete frozen test roster")
        groups = pd.read_csv(correction["eligibility"], sep="\t", dtype={"eid": str, "group_id": str}).set_index("eid")
        aligned = groups.loc[roles.id]
        if not aligned.eligible.all() or not np.array_equal(aligned.group_id, roles.group_id):
            raise InvalidPublication("incompatible corrected groups")
        if not roles.role.value_counts().sort_index().to_dict() == receipt["role_counts"]:
            raise InvalidPublication("role counts changed")
        table = read_aggregate(folder / "performance.tsv")
        required = {"brier", "logloss", "auc", "auprc"} if trait == "t2dm" else {"mse", "rmse", "mae", "prediction_r2"}
        if len(table) != 1 or not required <= set(table):
            raise InvalidPublication("missing aggregate metrics")
        row = table.iloc[0]
        if row.method not in (METHOD, "DiscoDivas" + " adaptation") or row.source not in ("3.disco.sh", "2." + "disco.sh") or row.trait != trait or row.stage != "test" or row.status != "completed":
            raise InvalidPublication("invalid aggregate provenance")
        if not np.isfinite(pd.to_numeric(row[list(required)], errors="coerce").to_numpy(float)).all():
            raise InvalidPublication("undefined aggregate metrics")
        table["method"] = METHOD
        table["source"] = "3.disco.sh"
        return table, "completed"
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return None, str(exc) if isinstance(exc, InvalidPublication) else "Missing or malformed Disco completion/provenance files"
