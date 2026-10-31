"""Verified local kinship overlay; preserve genotype row order and every test ID."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from abm_core.model import file_sha256
from runtime_state import file_identity

VERSION = "keep-test.exclude-dev-relatives.v2"


def compatible(metadata, record):
    return (record is None or
            (record.get("version") == VERSION and bool(record.get("token"))
             and metadata.get("grouping_status") == "verified_from_source"
             and metadata.get("groups_verified") is True
             and metadata.get("kinship_version") == record["version"]
             and metadata.get("kinship_token") == record["token"]))


def publication_metadata(root):
    record = active(root)
    return (dict(grouping_status="verified_from_source", groups_verified=True,
                 kinship_version=record["version"], kinship_token=record["token"])
            if record else {})


def analysis_guard(root, *, exclusive=False):
    """Barrier held for the entire module, including pre-run input loading."""
    from contextlib import contextmanager
    import fcntl

    @contextmanager
    def held():
        folder = Path(root) / "cache/abm"
        folder.mkdir(parents=True, exist_ok=True)
        with (folder / ".kinship.activation.lock").open("a") as stream:
            try:
                fcntl.flock(stream, (fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH) | fcntl.LOCK_NB)
            except BlockingIOError:
                raise BlockingIOError("Kinship activation refused: an active GRID/PRSformer job is using the analysis root.") from None
            yield
    return held()


def active(root):
    path = Path(root) / "local/kinship/ACTIVE.json"
    if not path.is_file():
        return None
    record = json.loads(path.read_text())
    if record.get("schema") == "grid.kinship.keep_test.2":
        for filename, expected in record["files"].items():
            if not Path(filename).is_file() or file_identity(filename) != expected:
                raise ValueError("Verified kinship input changed: " + filename)
        return dict(record, signature=record["token"])
    if record.get("schema") != "grid.kinship.keep_test.1":
        raise ValueError("Unsupported local kinship correction")
    for filename, expected in record["files"].items():
        if file_sha256(filename) != expected:
            raise ValueError("Verified kinship input changed: " + filename)
    # Legacy records can be inspected, but cannot satisfy a v2 publication.
    return dict(record, signature=file_sha256(path), version="legacy.1", token=None)


def aligned(root, ids):
    record = active(root)
    if record is None:
        return None, None
    d = pd.read_csv(record["eligibility"], sep="\t", dtype={"eid": str, "group_id": str})
    if d.eid.isna().any() or d.eid.duplicated().any() or not set(map(str, ids)) <= set(d.eid):
        raise ValueError("Kinship eligibility does not cover input IDs uniquely")
    return d.set_index("eid").loc[list(map(str, ids))].reset_index(), record


def correct_benchmark(data, root, seed, *, separate_calibration=False):
    rows, record = aligned(root, data.eid)
    if record is None:
        return data, None
    data = data.copy()
    if not np.array_equal(rows.is_test.astype(bool), data.split.eq("test")):
        raise ValueError("Kinship correction cannot change the frozen test roster")
    data["group"] = rows.group_id.to_numpy()
    if "family_id" in data:
        data["family_id"] = data.group
    eligible = rows.eligible.astype(bool).to_numpy()
    if not separate_calibration:
        # Same reviewed 80/20 internal development split, now by real families.
        from data import family_assignment
        dev = eligible & ~data.split.eq("test").to_numpy()
        assignment = family_assignment(data.loc[dev, "group"], .8, seed + 1)
        data.loc[dev, "split"] = np.where(data.loc[dev, "group"].map(assignment), "train", "validation")
    data.loc[~eligible, "split"] = "unused"
    return data, record["signature"]


def correct_grid(frame, metadata, root):
    rows, record = aligned(root, frame.id)
    if record is None:
        return frame, metadata
    if not np.array_equal(rows.is_test.astype(bool), frame.split.eq("test")):
        raise ValueError("Kinship correction cannot change the frozen test roster")
    frame = frame.copy()
    frame["group_id"] = rows.group_id.to_numpy()
    eligible = rows.eligible.astype(bool).to_numpy()
    frame = frame.loc[eligible].reset_index(drop=True)
    evidence = json.loads(Path(record["group_verification"]).read_text())
    metadata = dict(metadata, group_verification=evidence,
                    kinship_correction_sha256=record["signature"],
                    kinship_version=record["version"], kinship_token=record["token"],
                    excluded_development_relatives=int((~eligible).sum()))
    if set(frame.loc[frame.split.eq("test"), "group_id"]) & set(frame.loc[~frame.split.eq("test"), "group_id"]):
        raise ValueError("Kinship correction left test relatives in development")
    return frame, metadata
