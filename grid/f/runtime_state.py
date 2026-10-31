"""Local resume identities: metadata and exact small tables, never content scans."""
from pathlib import Path
import json
import uuid

import pandas as pd


def file_identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def artifact_matches(path, expected):
    if not Path(path).is_file():
        return False
    if isinstance(expected, dict):
        return file_identity(path) == expected
    # Legacy receipts remain readable; new phase receipts use file metadata.
    from abm_core.model import file_sha256
    return file_sha256(path) == expected


def group_contract(frame, metadata, directory):
    """Bind exact ordered ID/group/role state before phase reuse is considered.

    The private table contains IDs. Only its opaque version token enters the
    public contract. A changed state gets a new token even at identical counts.
    """
    from resume import run_lock, atomic_json
    directory = Path(directory)
    rows = frame[["id", "group_id", "role"]].astype(str).reset_index(drop=True)
    if rows.id.duplicated().any() or frame[["id", "group_id", "role"]].isna().any().any():
        raise ValueError("Run grouping requires unique IDs and complete groups/roles")
    with run_lock(directory):
        marker = directory / "CURRENT.json"
        prior = json.loads(marker.read_text()) if marker.is_file() else {}
        token = prior.get("token")
        path = directory / (str(token) + ".tsv.gz")
        if not token or not path.is_file() or not rows.equals(pd.read_csv(path, sep="\t", dtype=str)):
            token = None
            for saved in directory.glob("*.tsv.gz"):
                if rows.equals(pd.read_csv(saved, sep="\t", dtype=str)):
                    token = saved.name.removesuffix(".tsv.gz")
                    break
            if token is None:
                token = uuid.uuid4().hex
                rows.to_csv(directory / (token + ".tsv.gz"), sep="\t", index=False)
            atomic_json(marker, dict(schema="grid.group_state.1", token=token))
    return dict(schema="grid.group_contract.1", state_token=token,
                n_rows=len(rows), n_groups=rows.group_id.nunique(),
                role_counts={k: int(v) for k, v in rows.role.value_counts().sort_index().items()},
                grouping_status=metadata.get("grouping_status"),
                groups_verified=metadata.get("groups_verified") is True,
                kinship_version=metadata.get("kinship_version"),
                kinship_token=metadata.get("kinship_token"))
