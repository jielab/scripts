#!/usr/bin/env python3
"""
Download UKB-PPP GWAS files from Synapse.

Typical PPP European discovery run:
  python ppp_download_v2.py \
    --out /mnt/d/data.BIG/gwas/ppp/raw \
    --folder syn51365303 \
    --name-regex '\\.tar$' \
    --manifest European_discovery.download_manifest.tsv \
    --workers 8

Authentication:
  export SYNAPSE_AUTH_TOKEN="$(tr -d '\\r\\n' < /mnt/d/data/ukb/ppp/authToken.txt)"
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import synapseclient
except ModuleNotFoundError as exc:
    sys.stderr.write(
        "ERROR: Python module 'synapseclient' is not installed.\n"
        "Install it with one of:\n"
        "  python -m pip install --user synapseclient\n"
        "  mamba install -c conda-forge synapseclient\n"
    )
    raise SystemExit(2) from exc

# Silence Synapse getChildren deprecation warning. It is harmless for this script.
warnings.filterwarnings("ignore", category=DeprecationWarning)

_tls = threading.local()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download UKB-PPP Synapse files from explicit synIDs or a Synapse folder.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--out", default="/mnt/d/data.BIG/gwas/ppp/raw", help="Output directory")
    p.add_argument("--entities", nargs="*", default=[], help="Synapse entity IDs, e.g. syn51470659")
    p.add_argument("--entities-file", help="Text/TSV file containing Synapse IDs; comments with # are allowed")
    p.add_argument("--folder", help="Synapse Folder/Project ID to recursively scan, e.g. syn51365303")
    p.add_argument("--name-regex", default=".*", help="Regex applied to Synapse file names during folder scan")
    p.add_argument("--id-regex", default=r"^syn[0-9]+$", help="Regex used to recognize entity IDs")
    p.add_argument("--manifest", default="download_manifest.tsv", help="Manifest written under --out")
    p.add_argument("--profile", default=None, help="Synapse profile name in ~/.synapseConfig")
    p.add_argument("--limit", type=int, default=0, help="Limit number of files; 0 means no limit")
    p.add_argument("--sleep", type=float, default=0.0, help="Seconds to sleep after each finished download")
    p.add_argument("--retries", type=int, default=3, help="Retry count per entity")
    p.add_argument("--workers", type=int, default=1, help="Number of parallel download workers")
    p.add_argument("--dry-run", action="store_true", help="List matched files but do not download")
    p.add_argument("--check-login", action="store_true", help="Only test Synapse login, then exit")
    p.add_argument("--overwrite", action="store_true", help="Re-download even if output file already exists")
    return p.parse_args()


def read_entity_ids(path: Optional[str], id_regex: str) -> List[str]:
    if not path:
        return []
    rx = re.compile(id_regex)
    ids: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            for token in re.split(r"[\t, ]+", line):
                token = token.strip()
                if rx.match(token):
                    ids.append(token)
    return ids


def login(profile: Optional[str] = None):
    """Login once in the main process.

    synapseclient.login() automatically uses SYNAPSE_AUTH_TOKEN or ~/.synapseConfig.
    """
    if profile:
        return synapseclient.login(profile=profile, silent=True)
    return synapseclient.login(silent=True)


def worker_login(profile: Optional[str] = None):
    """Thread-local Synapse login.

    A separate client per worker is safer than sharing one client object across threads.
    """
    key = f"syn_{profile or 'default'}"
    syn = getattr(_tls, key, None)
    if syn is not None:
        return syn

    token = os.environ.get("SYNAPSE_AUTH_TOKEN", "").strip()
    syn = synapseclient.Synapse()
    if token:
        syn.login(authToken=token, silent=True)
    elif profile:
        syn = synapseclient.login(profile=profile, silent=True)
    else:
        syn = synapseclient.login(silent=True)
    setattr(_tls, key, syn)
    return syn


def recursively_list_files(syn, parent_id: str, name_regex: str) -> List[Dict[str, str]]:
    """Return FileEntity children under a Synapse folder/project.

    This avoids guessing consecutive Synapse IDs. It scans Folder/Project children,
    keeps FileEntity records whose names match --name-regex, and returns a stable,
    name-sorted list.
    """
    rx = re.compile(name_regex)
    stack = [parent_id]
    files: List[Dict[str, str]] = []
    seen = set()

    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        for child in syn.getChildren(pid):
            cid = child.get("id")
            name = child.get("name", "")
            ctype = child.get("type", "")
            if not cid:
                continue
            if ctype.endswith("Folder") or ctype.endswith("Project"):
                stack.append(cid)
            elif ctype.endswith("FileEntity") and rx.search(name):
                files.append({"entity_id": cid, "name": name, "parent": pid})

    files.sort(key=lambda x: (x.get("name", ""), x.get("entity_id", "")))
    return files


def read_expected_records(args: argparse.Namespace, syn) -> List[Dict[str, str]]:
    records: List[Dict[str, str]] = []

    # Records from a Synapse folder include file names, so skip-existing works well.
    if args.folder:
        records.extend(recursively_list_files(syn, args.folder, args.name_regex))

    # Records from explicit IDs may not have names until download.
    ids = list(args.entities) + read_entity_ids(args.entities_file, args.id_regex)
    for eid in ids:
        if re.match(args.id_regex, eid):
            records.append({"entity_id": eid, "name": "", "parent": ""})

    # Preserve first occurrence of each entity_id.
    out: List[Dict[str, str]] = []
    seen = set()
    for r in records:
        eid = r.get("entity_id", "")
        if not eid or eid in seen:
            continue
        seen.add(eid)
        out.append(r)

    if args.limit and args.limit > 0:
        out = out[: args.limit]
    return out


def download_one(record: Dict[str, str], outdir: Path, retries: int, overwrite: bool, profile: Optional[str]) -> Dict[str, str]:
    entity_id = record.get("entity_id", "")
    expected_name = record.get("name", "")

    # Fast skip when folder scan told us the expected filename.
    if expected_name and not overwrite:
        expected_path = outdir / expected_name
        if expected_path.exists() and expected_path.stat().st_size > 0:
            return {
                "entity_id": entity_id,
                "name": expected_name,
                "path": str(expected_path),
                "status": "SKIP_EXISTING",
                "error": "",
            }

    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            syn = worker_login(profile)
            ent = syn.get(entity_id, downloadLocation=str(outdir))
            path = getattr(ent, "path", "") or ""
            name = getattr(ent, "name", "") or expected_name or Path(path).name
            return {
                "entity_id": entity_id,
                "name": name,
                "path": path,
                "status": "OK",
                "error": "",
            }
        except Exception as e:
            last_error = str(e).replace("\n", " ")
            sys.stderr.write(f"WARN: {entity_id} attempt {attempt}/{retries} failed: {last_error}\n")
            sys.stderr.flush()
            time.sleep(min(60, 2 * attempt))

    return {
        "entity_id": entity_id,
        "name": expected_name,
        "path": "",
        "status": "FAIL",
        "error": last_error,
    }


def write_manifest(path: Path, rows: Iterable[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["entity_id", "name", "path", "status", "error"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})


def main() -> int:
    args = parse_args()
    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        syn = login(args.profile)
    except Exception as e:
        sys.stderr.write(f"ERROR: Synapse login failed: {e}\n")
        return 2

    if args.check_login:
        print("Synapse login OK")
        return 0

    records = read_expected_records(args, syn)
    if not records:
        sys.stderr.write("ERROR: no Synapse entity IDs were provided or discovered.\n")
        return 2

    sys.stderr.write(f"Matched {len(records)} Synapse file/entity IDs. Output: {outdir}\n")
    sys.stderr.flush()

    manifest_path = outdir / args.manifest

    if args.dry_run:
        rows = []
        for r in records:
            rows.append({
                "entity_id": r.get("entity_id", ""),
                "name": r.get("name", ""),
                "path": "",
                "status": "DRY_RUN",
                "error": "",
            })
        write_manifest(manifest_path, rows)
        for row in rows[:20]:
            print(row["entity_id"], row["name"])
        if len(rows) > 20:
            print(f"... {len(rows) - 20} more")
        sys.stderr.write(f"Dry-run manifest: {manifest_path}\n")
        return 0

    workers = max(1, int(args.workers))
    rows: List[Dict[str, str]] = []
    n_total = len(records)
    n_done = 0

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(download_one, r, outdir, args.retries, args.overwrite, args.profile) for r in records]
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            n_done += 1
            sys.stderr.write(f"[{n_done}/{n_total}] {row['status']} {row['entity_id']} {row.get('name', '')}\n")
            sys.stderr.flush()
            if args.sleep:
                time.sleep(args.sleep)

    # Sort manifest by name for reproducibility.
    rows.sort(key=lambda x: (x.get("name", ""), x.get("entity_id", "")))
    write_manifest(manifest_path, rows)

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_skip = sum(1 for r in rows if r["status"] == "SKIP_EXISTING")
    n_fail = sum(1 for r in rows if r["status"] == "FAIL")
    sys.stderr.write(f"Done. OK={n_ok}; SKIP_EXISTING={n_skip}; FAIL={n_fail}; manifest={manifest_path}\n")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
