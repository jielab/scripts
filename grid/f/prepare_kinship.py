#!/usr/bin/env python3
"""Build and atomically activate LOCAL connected-component evidence."""
import argparse
from contextlib import ExitStack, contextmanager
import fcntl
import json
from pathlib import Path
import uuid

import pandas as pd

from kinship import VERSION, analysis_guard
from runtime_state import file_identity
from resume import atomic_json

REFUSAL = "Kinship activation refused: an active GRID/PRSformer job is using the analysis root."


def components(edges, ids):
    parents = {}
    def find(value):
        parents.setdefault(value, value)
        while parents[value] != value:
            parents[value] = parents[parents[value]]
            value = parents[value]
        return value
    for left, right in zip(edges.ID1, edges.ID2):
        a, b = find(left), find(right)
        if a != b:
            parents[max(a, b)] = min(a, b)
    return pd.Series(ids).astype(str).map(find)


@contextmanager
def activation_locks(root):
    root = Path(root)
    with ExitStack() as stack:
        stack.enter_context(analysis_guard(root, exclusive=True))
        paths = {root / name for name in (
            "prsformer/.resume.lock", "prsformer2/.resume.lock", "workflow/.resume.lock",
            "cache/abm/.dosage.lock", "cache/abm/blocks/.dosage.lock", "cache/abm/grid.wrapper.lock", "cache/abm/.pipeline.lock")}
        for folder in (root / "abm", root / "cache/abm", root / "disco"):
            if folder.is_dir():
                paths.update(folder.rglob("*.lock"))
        paths.discard(root / "cache/abm/.kinship.activation.lock")
        for path in sorted(paths):
            path.parent.mkdir(parents=True, exist_ok=True)
            handle = stack.enter_context(path.open("a"))
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise BlockingIOError(REFUSAL) from None
        yield


def verified_tables(root, mapping, eligibility):
    roster = pd.read_csv(root / "prsformer/cache/data.tsv.gz", sep="\t", dtype=str)
    groups = pd.read_csv(mapping, sep="\t", dtype={"eid": str, "group_id": str})
    rows = pd.read_csv(eligibility, sep="\t", dtype={"eid": str, "group_id": str})
    for table in (roster, groups, rows):
        if table.eid.isna().any() or table.eid.duplicated().any():
            raise ValueError("Kinship tables require unique nonmissing IDs")
    if not roster.eid.equals(rows.eid) or not groups[["eid", "group_id"]].equals(rows[["eid", "group_id"]]):
        raise ValueError("Kinship tables differ from exact ordered producer roster/groups")
    test = roster.split.eq("test")
    if not rows.is_test.equals(test) or not test.any():
        raise ValueError("Kinship correction cannot change the frozen test roster")
    expected = test | ~rows.group_id.isin(set(rows.loc[test, "group_id"]))
    if not rows.eligible.equals(expected):
        raise ValueError("Kinship eligibility is not keep-test/exclude-development-relatives")
    return rows


def activate_correction(root, out, source, mapping, eligibility, *, expected_source=None):
    """All locks precede all moves; ACTIVE is the last successful transaction write."""
    root, out = Path(root).resolve(), Path(out).resolve()
    source, mapping, eligibility = map(lambda p: Path(p).resolve(), (source, mapping, eligibility))
    with activation_locks(root):
        rows = verified_tables(root, mapping, eligibility)
        evidence = out / "group_verification.json"
        source_state = file_identity(source)
        if expected_source is not None and source_state != expected_source:
            raise ValueError("Kinship source changed before activation")
        # Random state suffix also distinguishes changed mappings with equal counts.
        token = (f"{VERSION}:{rows.group_id.nunique()}:{int((~rows.eligible).sum())}:"
                 f"{source_state['size']}:{source_state['mtime_ns']}:{uuid.uuid4().hex}")
        current = root / "local/kinship/ACTIVE.json"
        if current.is_file():
            old = json.loads(current.read_text())
            if old.get("schema") == "grid.kinship.keep_test.2":
                from kinship import active
                active(root)
                previous = pd.read_csv(old["eligibility"], sep="\t", dtype={"eid": str, "group_id": str})
                if rows.equals(previous) and old["files"].get(str(source)) == source_state:
                    return old
                raise ValueError("Active correction differs; use a separately reviewed run")
        atomic_json(evidence, dict(status="verified", source=str(mapping),
                    identity=file_identity(mapping), kinship_source=str(source),
                    kinship_source_identity=source_state,
                    component_rule="all supplied edges; no additional threshold"))
        record = dict(schema="grid.kinship.keep_test.2", version=VERSION, token=token,
                      eligibility=str(eligibility), group_verification=str(evidence),
                      qc=str(out / "qc.json"), source=str(source),
                      files={str(p): file_identity(p) for p in
                             (source, mapping, eligibility, evidence, root / "prsformer/cache/data.tsv.gz")})
        if (out / "qc.json").is_file():
            record["files"][str(out / "qc.json")] = file_identity(out / "qc.json")
        history = root / "local/kinship/history_before_correction" / uuid.uuid4().hex
        actions = []
        for pointer in (root / "abm/FINAL").glob("*/*.CURRENT.json"):
            actions.append(dict(action="archive", path=str(pointer), destination=str(history / pointer.relative_to(root))))
        folders = [root / "prsformer" / trait for trait in ("height", "ldl", "t2dm")]
        # Keep the lock inode in place while archiving fixed benchmark outputs.
        folders += [root / "prsformer2" / trait for trait in ("height", "ldl", "t2dm")]
        folders += [root / "abm_reports"]
        for folder in folders:
            if folder.exists():
                actions.append(dict(action="archive", path=str(folder), destination=str(history / folder.relative_to(root))))
        for folder in (root / "final", root / "disco"):
            if folder.is_dir():
                paths = folder.glob("*") if folder.name == "final" else folder.glob("*/performance*")
                for path in paths:
                    if path.is_file():
                        actions.append(dict(action="archive", path=str(path), destination=str(history / path.relative_to(root))))
        atomic_json(out / "activation_actions.json", dict(status="planned", actions=actions))
        moved = []
        try:
            for item in actions:
                src, dst = Path(item["path"]), Path(item["destination"])
                dst.parent.mkdir(parents=True, exist_ok=True)
                src.rename(dst)
                moved.append(item)
            if any(file_identity(p) != identity for p, identity in record["files"].items()):
                raise ValueError("Kinship input changed during activation")
            # Log success before the commit marker; a crash leaves no new ACTIVE.
            atomic_json(out / "activation_result.json", dict(status="completed", actions=moved))
            atomic_json(current, record)
        except BaseException as exc:
            rollback_errors = []
            for item in reversed(moved):
                try:
                    Path(item["destination"]).rename(item["path"])
                except OSError as error:
                    rollback_errors.append(str(error))
            atomic_json(out / "activation_failure.json", dict(error=str(exc), actions=moved, rollback_errors=rollback_errors))
            raise
        return record


def prepare(source, analysis_root, *, activate=False):
    source, root = Path(source).resolve(), Path(analysis_root).resolve()
    source_state = file_identity(source)
    roster = pd.read_csv(root / "prsformer/cache/data.tsv.gz", sep="\t", usecols=["eid", "group", "split"], dtype=str)
    if roster.isna().any().any() or roster.eid.duplicated().any():
        raise ValueError("Frozen producer roster requires unique nonmissing IDs")
    if not set(roster.split) <= {"train", "validation", "test"} or not roster.split.eq("test").any():
        raise ValueError("Frozen producer roster has invalid splits")
    edges = pd.read_csv(source, sep=r"\s+", dtype={"ID1": str, "ID2": str})
    if not {"ID1", "ID2", "Kinship"} <= set(edges) or edges[["ID1", "ID2", "Kinship"]].isna().any().any():
        raise ValueError("Kinship source requires nonmissing ID1, ID2, Kinship")
    group = components(edges, roster.eid)
    counts = group.value_counts()
    test = roster.split.eq("test")
    eligible = test | ~group.isin(set(group[test]))
    qc = dict(N_participants=len(roster), N_connected_components=len(counts),
        N_singleton_components=int(counts.eq(1).sum()), N_non_singleton_components=int(counts.gt(1).sum()),
        N_participants_in_non_singleton_components=int(counts[counts.gt(1)].sum()), maximum_component_size=int(counts.max()),
        N_test_retained=int(test.sum()), N_development_with_test_relative_before_correction=int((~eligible).sum()),
        N_development_excluded=int((~eligible).sum()), N_development_retained=int((eligible & ~test).sum()),
        N_test_with_development_relative_after_correction=int(group[test].isin(set(group[eligible & ~test])).sum()),
        component_rule="all edges supplied by ukb_rel.dat; no additional threshold",
        correction_status="active_keep_all_test_exclude_development_relatives" if activate else "prepared_not_active",
        source_path_basename=source.name, source_size=source_state["size"], source_mtime_ns=source_state["mtime_ns"])
    # Never overwrite active scientific evidence during preparation.
    out = root / "local/kinship/versions" / uuid.uuid4().hex
    out.mkdir(parents=True, mode=0o700)
    mapping, eligibility = out / "connected_components.tsv.gz", out / "eligibility.tsv.gz"
    pd.DataFrame(dict(eid=roster.eid, group_id=group)).to_csv(mapping, sep="\t", index=False)
    pd.DataFrame(dict(eid=roster.eid, group_id=group, is_test=test, eligible=eligible)).to_csv(eligibility, sep="\t", index=False)
    atomic_json(out / "qc.json", qc)
    pd.DataFrame([qc]).to_csv(out / "qc.tsv", sep="\t", index=False)
    if file_identity(source) != source_state:
        raise ValueError("Kinship source changed during preparation")
    if activate:
        activate_correction(root, out, source, mapping, eligibility, expected_source=source_state)
    print(json.dumps(dict(output=str(out), **qc), indent=2))
    return qc


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--activate", action="store_true")
    a = p.parse_args(argv)
    prepare(a.source, a.analysis_root, activate=a.activate)


if __name__ == "__main__":
    main()
