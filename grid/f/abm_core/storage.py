"""Release completed support scratch arrays without changing frozen receipts.

Run with PYTHONPATH=f python -m abm_core.storage; add --apply to delete.
The donor bank, block inputs and shared genotypes are persistent dependencies.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import time

from resume import atomic_json
from .model import file_sha256

RELEASE_FILE = "SCRATCH_RELEASED.json"
SCRATCH = ("geometry.npy", "radius_geometry.npy")


@contextmanager
def support_lock(cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    with (cache / ".support.lock").open("a") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def busy_reason(paths):
    """Also protect legacy writers that do not hold .support.lock (including self)."""
    targets = {str(Path(p).resolve()) for p in paths}
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        for attempt in range(3):
            try:
                if proc.stat().st_uid != os.getuid():
                    break
                if (proc / "stat").read_text().rsplit(")", 1)[1].split()[0] in {"Z", "X"}:
                    break
                if (proc / "comm").read_text().strip() == "(sd-pam)":
                    break
                for line in (proc / "maps").read_text().splitlines():
                    fields = line.split(maxsplit=5)
                    if len(fields) == 6 and fields[5].removesuffix(" (deleted)") in targets:
                        return "mapped by PID " + proc.name
                for fd in (proc / "fd").iterdir():
                    try:
                        if os.readlink(fd).removesuffix(" (deleted)") in targets:
                            return "open in PID " + proc.name
                    except FileNotFoundError:
                        pass
            except (FileNotFoundError, ProcessLookupError):
                break
            except PermissionError:
                # /proc can deny maps/fd access during exit, before stat becomes
                # Z or the PID disappears. Recheck instead of treating that race
                # as a persistent, uninspectable cache user.
                try:
                    command = (proc / "cmdline").read_bytes().split(b"\0")
                except (FileNotFoundError, ProcessLookupError):
                    break
                except OSError:
                    command = []
                # The per-user service manager is nondumpable on WSL and is
                # not a GRID worker. Other persistent denials stay conservative.
                if command[:2] in ([b"/usr/lib/systemd/systemd", b"--user"],
                                   [b"/lib/systemd/systemd", b"--user"]):
                    break
                if attempt < 2:
                    time.sleep(.01)
                    continue
                return "cannot inspect PID " + proc.name
            break
    return None


def validate_support(cache):
    """Only a bound release receipt can excuse a missing scratch array."""
    cache = Path(cache)
    pointer = cache / "COMPLETE.json"
    saved = json.loads(pointer.read_text())
    released = {}
    if (cache / RELEASE_FILE).is_file():
        release = json.loads((cache / RELEASE_FILE).read_text())
        if release.get("complete_sha256") != file_sha256(pointer):
            raise ValueError("Support scratch release belongs to a different completion receipt")
        released = release["files"]
        if set(released) - set(SCRATCH):
            raise ValueError("Release receipt cannot remove persistent support files")
    required = {"support.parquet", "bank.parquet", "geometry_rows.parquet",
                "geometry.npy", "donor_geometry.npy", "neighbor_indices.npy",
                "neighbor_distances.npy", "geometry_manifest.json", "donor_bank_manifest.json",
                "geometry.npy.rows.json", "donor_geometry.npy.rows.json"}
    if not required <= saved["files"].keys():
        raise ValueError("Incomplete support receipt")
    for name, expected in saved["files"].items():
        path = cache / name
        if Path(name).name != name or path.is_symlink():
            raise ValueError("Unsafe support cache file: " + name)
        if not path.exists() and name in SCRATCH and released.get(name, {}).get("sha256") == expected:
            continue
        if not path.is_file() or file_sha256(path) != expected:
            raise ValueError("Shared support cache changed: " + name)
    return saved


def release_verified_scratch(cache, saved):
    """Caller holds support_lock and has produced or validated every saved file."""
    cache = Path(cache)
    targets = [cache / name for name in SCRATCH if (cache / name).exists()]
    if not targets:
        return 0
    if any(p.is_symlink() for p in targets):
        raise ValueError("Refusing to release symlinked support scratch")
    reason = busy_reason([*targets, cache / "donor_geometry.npy"])
    if reason:
        print(f"[GRID cleanup] RETAIN {cache}: {reason}", flush=True)
        return 0
    pointer_sha = file_sha256(cache / "COMPLETE.json")
    release_path = cache / RELEASE_FILE
    release = (json.loads(release_path.read_text()) if release_path.exists()
               else dict(complete_sha256=pointer_sha, files={}))
    if release.get("complete_sha256") != pointer_sha:
        raise ValueError("Support completion receipt changed during release")
    total = 0
    for path in targets:
        size = path.stat().st_size
        release["files"][path.name] = dict(bytes=size, sha256=saved["files"].get(path.name))
        total += size
    # Commit first: interruption before unlink is safe and can be retried.
    # COMPLETE.json stays byte-for-byte unchanged for existing model pointers.
    atomic_json(release_path, release)
    for path in targets:
        path.unlink()
    print(f"[GRID cleanup] Released {total / 2**30:.2f} GiB: {cache}", flush=True)
    return total


def clean_cache(cache, *, apply=False):
    cache = Path(cache)
    if cache.is_symlink() or (cache / "COMPLETE.json").is_symlink():
        raise ValueError("Refusing symlinked support cache")
    if not (cache / "COMPLETE.json").is_file():
        print(f"[GRID cleanup] SKIP incomplete support: {cache}", flush=True)
        return 0
    with support_lock(cache):
        targets = [cache / name for name in SCRATCH if (cache / name).exists()]
        if not targets:
            return 0
        reason = busy_reason([*targets, cache / "donor_geometry.npy"])
        if reason:
            print(f"[GRID cleanup] SKIP {cache}: {reason}", flush=True)
            return 0
        size = sum(p.stat().st_size for p in targets)
        if not apply:
            print(f"[GRID cleanup] CANDIDATE {size / 2**30:.2f} GiB (validation on --apply): {cache}", flush=True)
            return size
        print(f"[GRID cleanup] Verifying completed support: {cache}", flush=True)
        saved = validate_support(cache)
        return release_verified_scratch(cache, saved)


def clean_completed_support(analysis_root, *, traits=None, apply=False, best_effort=False):
    """Sweep old completion keys too, including when the runner reuses all phases."""
    total = 0
    for namespace in ("support", "SMOKE"):
        root = Path(analysis_root) / "cache/abm" / namespace
        for cache in sorted(root.glob("*/*/_LOCAL_ONLY")):
            if traits is not None and cache.parent.parent.name not in traits:
                continue
            try:
                total += clean_cache(cache, apply=apply)
            except BlockingIOError:
                print(f"[GRID cleanup] SKIP locked support: {cache}", flush=True)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                if not best_effort:
                    raise
                # An unrelated historical cache must not break a new analysis
                # or mask its original exception. Failed validation retains it.
                print(f"[GRID cleanup] RETAIN {cache}: {exc}", flush=True)
    return total


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args(argv)
    total = clean_completed_support(args.analysis_root, apply=args.apply)
    print(f"[GRID cleanup] {'Released' if args.apply else 'Candidates'} total: {total / 2**30:.2f} GiB", flush=True)


if __name__ == "__main__":
    main()
