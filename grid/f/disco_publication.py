#!/usr/bin/env python3
"""Read-only validation of the Disco module's aggregate outcome comparison."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
from final_report import read_aggregate
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
        if receipt.get("schema") == "grid.target_benchmark.1":
            from shared_cohort import file_record
            if receipt.get("trait") != trait or receipt.get("method") != "disco":
                raise InvalidPublication("wrong standalone Disco trait/method")
            if not receipt.get("test_labels_used_only_after_prediction_freeze"):
                raise InvalidPublication("missing standalone prediction-freeze evidence")
            required = {str((folder / "performance.tsv").resolve()),
                        str((folder / "_LOCAL_ONLY/roles.tsv.gz").resolve()),
                        str((folder / "_LOCAL_ONLY/frozen_models.json").resolve()),
                        str((folder / "_LOCAL_ONLY/final_prs.tsv.gz").resolve())}
            if not required <= set(receipt.get("files", {})):
                raise InvalidPublication("missing bound standalone outputs")
            for path, digest in receipt["files"].items():
                if not Path(path).is_file() or file_record(path) != digest:
                    raise InvalidPublication("standalone input/output changed")
            roles = pd.read_csv(folder / "_LOCAL_ONLY/roles.tsv.gz", sep="\t", dtype=str)
            if roles.id.duplicated().any() or roles.groupby("group_id").role.nunique().max() != 1:
                raise InvalidPublication("standalone group role isolation failed")
            if int(roles.role.eq("test").sum()) != receipt["roster"]["test_N"]:
                raise InvalidPublication("standalone frozen test count changed")
            table = pd.read_csv(folder / "performance.tsv", sep="\t")
            if len(table) != 1 or not table.source.eq("3.disco.sh").all() or not table.stage.eq("test").all():
                raise InvalidPublication("standalone aggregate provenance invalid")
            return table, "completed"
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
        roster_path = root / "tmp/prsformer/cache/data.tsv.gz"
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
