"""Attach completed, verified HERV caches to a fresh report directory."""
from pathlib import Path
import fcntl
import hashlib
import json
import os
import shutil


def seed_herv(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError("HERV reuse source and new output must be separate directories")
    names = ("run.json", "prepared.json", "analysis.done.json", "plan.json",
             "genome_fraction.json", "genome_fraction.tsv", "variant_summary.tsv",
             "chromosome_fraction.tsv", "boundary_audit.tsv")
    if any(not (source / name).is_file() for name in names):
        raise ValueError(f"Incomplete HERV reuse source: {source}; no annotation was started")
    run = json.loads((source / "run.json").read_text())
    done = json.loads((source / "analysis.done.json").read_text())
    if run.get("status") != "complete" or done.get("settings") != run.get("settings"):
        raise ValueError("HERV reuse requires a matching completed run and receipt")
    if any(not (source / name).is_dir() for name in ("variants", "chunks")):
        raise ValueError("Completed HERV variant/chunk caches are missing")
    destination.mkdir(parents=True, exist_ok=True)
    receipt = destination / "reuse_source.json"
    with (destination / ".reuse.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if receipt.exists():
            if json.loads(receipt.read_text()).get("source") != str(source):
                raise ValueError("Existing output uses a different HERV reuse source")
        elif any(p.name not in {"logs", ".reuse.lock"} for p in destination.iterdir()):
            raise ValueError("Attach HERV reuse only to a fresh output directory")
        state = {"source": str(source), "status": "attaching", "mode": "completed_cache_links"}
        receipt.write_text(json.dumps(state, indent=2) + "\n")
        for name in names:
            original, target = source / name, destination / name
            contents = original.read_bytes()
            if target.exists():
                if target.is_symlink() or target.read_bytes() != contents:
                    raise ValueError(f"Existing HERV metadata differs: {target}")
            else:
                shutil.copy2(original, target)
            if hashlib.sha256(target.read_bytes()).digest() != hashlib.sha256(contents).digest():
                raise ValueError(f"Copied HERV metadata changed: {target}")
        for name in ("variants", "chunks"):
            target = destination / name
            if target.is_symlink():
                if target.resolve() != source / name:
                    raise ValueError(f"HERV cache link points elsewhere: {target}")
            elif target.exists():
                raise ValueError(f"Refusing to replace existing HERV cache: {target}")
            else:
                os.symlink(source / name, target, target_is_directory=True)
        state["status"] = "attached_pending_engine_verification"
        receipt.write_text(json.dumps(state, indent=2) + "\n")
    return state
