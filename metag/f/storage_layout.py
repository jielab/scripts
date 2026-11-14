"""Permanent library data and separate execution state."""
from pathlib import Path


def separated(cfg):
    return cfg.get("storage_layout") == "separated"


def library_root(cfg):
    return Path(cfg.get("parallel_reference_dir") or
                (Path(cfg["library_dir"]) if separated(cfg) else Path(cfg["library_dir"]) / "parallel"))


def library_jobs(cfg):
    return Path(cfg["work_dir"]) / "jobs" / "library" if separated(cfg) else library_root(cfg)


def category_jobs(cfg, category):
    return library_jobs(cfg) / category if separated(cfg) else library_root(cfg) / "index" / category


def category_logs(cfg, category):
    return Path(cfg["work_dir"]) / "logs" / "library" / category if separated(cfg) else category_jobs(cfg, category)


def publish_comparison(cfg, source):
    """Expose completed comparison tables outside private scheduler plans."""
    import json
    from common import atomic_json, locked
    from workflow_layout import publish_result
    source = Path(source)
    if json.loads((source / "report.json").read_text()).get("status") != "complete":
        raise ValueError("Cannot publish an incomplete comparison")
    target = Path(cfg["work_dir"]) / "compare"
    qc = target / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    with locked(qc / "publish.lock"):
        receipt = qc / "published.json"
        previous = json.loads(receipt.read_text()).get("files", []) if receipt.is_file() else []
        names = [name for name in ("benchmark.tsv", "agreement.tsv", "truth.tsv", "report.json")
                 if (source / name).is_file()]
        for name in names:
            publish_result(source / name, target / name)
        if "truth.tsv" in previous and "truth.tsv" not in names:
            (target / "truth.tsv").unlink(missing_ok=True)
        atomic_json(receipt, {"source": str(source), "files": names})
    return target


def metadata_dir(folder):
    """Canonical metadata location in the current flat sample layout."""
    return Path(folder) / "qc"


def new_run(parent):
    """Reserve run, run2, ... using mkdir, never hashes or timestamps."""
    parent = Path(parent)
    parent.mkdir(parents=True, exist_ok=True)
    number = 1
    while True:
        candidate = parent / ("run" if number == 1 else f"run{number}")
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            number += 1


def plan_directory(parent, identity):
    from common import atomic_json, locked
    import json
    parent = Path(parent)
    with locked(parent / "plan.lock"):
        for folder in sorted(parent.glob("run*")):
            record = folder / "identity.json"
            if record.is_file() and json.loads(record.read_text()).get("identity") == identity:
                return folder
        folder = new_run(parent)
        atomic_json(folder / "identity.json", {"identity": identity})
        return folder
