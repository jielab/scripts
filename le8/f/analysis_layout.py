"""Locations of shared reports and run support files, including legacy paths."""
from copy import deepcopy
from pathlib import Path


def cohort_directory(root):
    root = Path(root)
    current, legacy = root / "logs/freeze", root / "_protocol"
    # Existing unmigrated studies keep using their frozen files.
    if legacy.is_dir() and current.exists():
        raise RuntimeError("Both _protocol and logs/freeze exist; resolve the frozen-directory conflict first")
    return legacy if legacy.is_dir() else current


def resolve_cohort_path(path):
    """Read legacy frozen references without rewriting their signed sidecars."""
    path = Path(path)
    if path.parent.name == "_protocol" and not path.exists():
        current = path.parent.parent / "logs/freeze" / path.name
        if current.parent.is_dir():
            return current
    return path


def normalize_cohort_request(request):
    """Ignore only the supported relocation when comparing input requests."""
    request = deepcopy(request)
    for key in ("group_file", "outer_roster"):
        value = request.get("fields", {}).get(key)
        if value:
            request["fields"][key] = str(resolve_cohort_path(value))
        record = request.get("files", {}).get(key)
        if record and record.get("path"):
            record["path"] = str(resolve_cohort_path(record["path"]))
    return request


def cohort_lock_path(path):
    path = resolve_cohort_path(path)
    if path.parent.name == "freeze" and path.parent.parent.name == "logs":
        return path.parent.parent / "locks" / (path.name + ".lock")
    return Path(str(path) + ".lock")


def report_directory(root):
    return Path(root) / "reports/c1"


def summary_lock_directory(root):
    return Path(root) / "logs/locks"
