#!/usr/bin/env python3
"""Adopt a completed gnomAD scan without rescanning or rewriting its manifests.

The store is portable because path mappings are relative to this directory.
Original run/COMPLETE identities remain byte-for-byte unchanged. Verification
uses the original size/edge fingerprints, not a new full-file checksum. Raw VCF
and reference FASTA inputs outside the old result directory are never copied.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import copy
import datetime as dt
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

if "allele_common" not in sys.modules:
    spec = importlib.util.spec_from_file_location("allele_common", Path(__file__).with_name("0.common.py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
from allele_common import dump_json, fingerprint, stable_hash

SCHEMA = "allele3-scan-store-1.0"
DEFAULT_SOURCE = "/mnt/d/analysis/allele3/final/gnomad_local"
DEFAULT_STORE = "/mnt/f/annot/gnomad/allele3"
OWNED = {"relocation.json", "adoption_receipt.json", "adoption_state.json", ".adoption-transfer"}
FP_KEYS = ("size", "mtime_ns", "edge_sha256")


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def absolute(path, root=None):
    p = Path(path).expanduser()
    return (p if p.is_absolute() else Path(root or Path.cwd()) / p).resolve()


def relative_target(root, relative):
    p = Path(relative)
    if p.is_absolute() or ".." in p.parts:
        raise ValueError(f"Invalid relative store path: {relative}")
    target = (root / p).resolve()
    if target != root and not target.is_relative_to(root):
        raise ValueError(f"Store path escapes its root: {relative}")
    return target


def same_content(actual, expected):
    return all(actual.get(k) == expected.get(k) for k in ("size", "edge_sha256"))


def _raw_manifest(root):
    path = root / "run.json"
    if not path.is_file():
        raise ValueError(f"Missing completed scan manifest: {path}. Adopt the existing completed scan with 3.scan.sh adopt.")
    run = json.loads(path.read_text())
    if run.get("status") not in {"complete", "scanned"}:
        raise ValueError(f"Scan is not complete at {root}; adoption cannot resume or replace a scan")
    if not run.get("signature") or not run.get("shards"):
        raise ValueError(f"Missing original scan signature or shards: {path}")
    if run.get("identity") and stable_hash(run["identity"]) != run["signature"]:
        raise ValueError(f"Original scan identity/signature mismatch: {path}; restore the original manifest")
    return run


def load_run_metadata(out):
    """Cheap path-aware reader; does not open VCF, FASTA, or shard data.

    ``identity`` and ``signature`` are preserved. Runtime ``options``, ``shards``,
    ``cache`` and ``summary`` resolve through a verified relocation sidecar.
    Use scan.load_run() for the original per-shard output integrity checks.
    """
    root = absolute(out)
    run = _raw_manifest(root)
    mapping, files, original = [], {}, str(root)
    sidecar = root / "relocation.json"
    if ((root / "adoption_state.json").exists() or (root / "adoption_receipt.json").exists()) and not sidecar.is_file():
        raise ValueError(f"Incomplete adoption at {root}; do not read original D:-drive shards through this partial copy. Resume 3.scan.sh adopt.")
    if sidecar.is_file():
        relocation = json.loads(sidecar.read_text())
        receipt_path = root / "adoption_receipt.json"
        if relocation.get("schema") != SCHEMA or relocation.get("status") != "complete" or not receipt_path.is_file():
            raise ValueError(f"Incomplete/unrecognized adopted store at {root}; resume 3.scan.sh adopt")
        receipt = json.loads(receipt_path.read_text())
        if (receipt.get("schema") != SCHEMA or receipt.get("status") != "complete"
                or relocation.get("run_signature") != run["signature"]
                or receipt.get("run_signature") != run["signature"]
                or receipt.get("path_mappings") != relocation.get("path_mappings")):
            raise ValueError(f"Relocation receipt disagrees with the original scan at {root}")
        recorded = relocation.get("run_manifest", {})
        if not same_content(fingerprint(root / "run.json"), recorded):
            raise ValueError("Adopted run.json changed; original scientific manifests must remain unchanged")
        mapping = relocation["path_mappings"]
        original = relocation["source_root"]
        for item in mapping:
            if not Path(item["source"]).is_absolute():
                raise ValueError("Relocation source mappings must be absolute")
            relative_target(root, item["destination"])
        files = {r["relative_path"]: r for r in receipt.get("files", [])}
        if len(files) != len(receipt.get("files", [])):
            raise ValueError("Duplicate file paths in adoption receipt")
    run["_store"] = {"root": str(root), "source_root": original, "mappings": mapping, "files": files}

    def project(value):
        if isinstance(value, dict):
            if {"path", "size", "mtime_ns"} <= value.keys():
                return resolve_recorded_fingerprint(run, value)
            return {k: project(v) for k, v in value.items()}
        if isinstance(value, list):
            return [project(v) for v in value]
        if isinstance(value, str) and Path(value).is_absolute():
            return str(resolve_run_path(run, value))
        return value

    run["shards"] = [str(resolve_run_path(run, p)) for p in run["shards"]]
    if len(set(run["shards"])) != len(run["shards"]):
        raise ValueError("Duplicate completed shards in the scan store")
    for name in ("cache", "summary"):
        if run.get(name): run[name] = str(resolve_run_path(run, run[name]))
    run["options"] = project(run.get("options", run.get("identity", {}).get("options", {})))
    return run


def resolve_run_path(run, path):
    """Resolve only recorded scan-root/copy mappings; unrelated inputs stay put."""
    state = run.get("_store", {})
    root = absolute(state.get("root", "."))
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = absolute(p, state.get("source_root", root))
    else:
        p = p.resolve()
    for entry in sorted(state.get("mappings", []), key=lambda r: len(Path(r["source"]).parts), reverse=True):
        source = Path(entry["source"])
        # Only recorded completed shards were copied from an external cache.
        # Its other children may be raw inputs and must retain their old paths.
        if entry.get("kind") == "cache_root" and p != source:
            continue
        if p == source or p.is_relative_to(source):
            target = relative_target(root, entry["destination"])
            return target / p.relative_to(source)
    return p


def resolve_recorded_fingerprint(run, expected):
    """Project a fingerprint without changing original run/COMPLETE identities.

    copy2 may round timestamps on another filesystem. In that case the adoption
    receipt records the verified destination timestamp; later changes still fail.
    Size/edge digest must remain equal to the original, irrespective of location.
    """
    result = copy.deepcopy(expected)
    if not isinstance(result, dict) or "path" not in result:
        return result
    p = resolve_run_path(run, result["path"])
    result["path"] = str(p)
    root = Path(run.get("_store", {}).get("root", "/nonexistent"))
    if p.is_relative_to(root):
        saved = run.get("_store", {}).get("files", {}).get(p.relative_to(root).as_posix())
        if saved:
            if not same_content(saved["destination"], expected):
                raise ValueError(f"Copied file disagrees with its original input fingerprint: {p}")
            result["mtime_ns"] = saved["destination"]["mtime_ns"]
    return result


def verify_recorded_input(run, expected):
    """Check an original bounded fingerprint at its current recorded location."""
    if not isinstance(expected, dict) or "path" not in expected:
        return
    target = resolve_recorded_fingerprint(run, expected)
    actual = fingerprint(target["path"])
    keys = [k for k in FP_KEYS if k in target]
    if any(actual[k] != target[k] for k in keys):
        raise ValueError(f"Recorded scan input changed: {target['path']}; restore it, do not rescan to repair a move")
    return actual


def verify_cached_fingerprint(path, expected):
    """Preserve cache checks while allowing an explicitly receipted store move."""
    p = absolute(path)
    actual = fingerprint(p)
    if actual == expected:
        return actual
    for root in p.parents:
        if (root / "relocation.json").is_file():
            run = load_run_metadata(root)
            target = resolve_recorded_fingerprint(run, expected)
            if target.get("path") == str(p) and all(actual[k] == target.get(k) for k in FP_KEYS):
                return actual
            break
    raise ValueError(f"Recorded annotation cache changed: {p}; restore the cache or select a new cache directory")


def verify_shard(shard, signature):
    """The original size/edge checks; stored absolute output paths are not keys."""
    root = absolute(shard)
    marker = root / "COMPLETE.json"
    if not marker.is_file():
        raise ValueError(f"Missing shard completion marker: {marker}")
    complete = json.loads(marker.read_text())
    if complete.get("signature") != signature or not complete.get("outputs"):
        raise ValueError(f"Completed shard identity differs: {root}")
    for name, saved in complete["outputs"].items():
        if Path(name).name != name:
            raise ValueError(f"Non-local output name in COMPLETE.json: {name}")
        p = root / name
        if not p.is_file() or not same_content(fingerprint(p), saved):
            raise ValueError(f"Completed cache file is missing or changed: {p}; restore it, do not rescan during adoption")
    return complete


def _tree(directory):
    """Stable inventory. Refuse links instead of silently leaving broken links."""
    directory = absolute(directory)
    if not directory.is_dir():
        raise ValueError(f"Recorded generated directory is missing: {directory}")
    directories, files = [], []
    for current, dirs, names in os.walk(directory, followlinks=False):
        dirs.sort(); names.sort()
        base = Path(current)
        for name in (*dirs, *names):
            p = base / name
            if p.is_symlink():
                raise ValueError(f"Symlink in scan artifacts: {p}. Materialize this result in a separate copy before adoption; source is untouched.")
        directories.append(base)
        for name in names:
            p = base / name
            if not p.is_file():
                raise ValueError(f"Non-regular file in scan artifacts: {p}")
            files.append(p)
    return directories, files


def build_plan(source, store):
    source, store = absolute(source), absolute(store)
    if source == store or source.is_relative_to(store) or store.is_relative_to(source):
        raise ValueError("Source and permanent store must be separate, non-nested directories")
    run = _raw_manifest(source)
    if run["status"] != "complete":
        raise ValueError("Adoption requires a completed original scan and summary; finish the original scan first")
    if any((source / name).exists() for name in OWNED):
        raise ValueError("Source already contains store sidecars; copy that portable store as a whole instead of adopting it again")
    options = run.get("options", run.get("identity", {}).get("options", {}))
    copies = [(source, Path("."), "scan_root")]
    mappings = [{"source": str(source), "destination": ".", "kind": "scan_root"}]
    original_shards = [absolute(p, source) for p in run["shards"]]
    if len(set(original_shards)) != len(original_shards):
        raise ValueError("Duplicate shards in the original run")
    external = [p for p in original_shards if not p.is_relative_to(source)]
    cache = absolute(run.get("cache", source / "shards"), source)
    for shard in external:
        relative = Path("_external/shards") / (shard.name + "." + stable_hash(str(shard))[:12])
        copies.append((shard, relative, "completed_shard"))
        mappings.append({"source": str(shard), "destination": relative.as_posix(), "kind": "completed_shard"})
    if external and not cache.is_relative_to(source):
        mappings.append({"source": str(cache), "destination": "_external/shards", "kind": "cache_root"})
    annotation = options.get("architecture_annotation_cache")
    if annotation:
        annotation = absolute(annotation, source)
        if not annotation.is_relative_to(source):
            copies.append((annotation, Path("_external/annotations"), "annotation_cache"))
            mappings.append({"source": str(annotation), "destination": "_external/annotations", "kind": "annotation_cache"})
    files, dirs, destinations = [], set(), set()
    for origin, prefix, kind in copies:
        if store == origin or store.is_relative_to(origin) or origin.is_relative_to(store):
            raise ValueError(f"Store overlaps a recorded generated input: {origin}")
        tree_dirs, tree_files = _tree(origin)
        for d in tree_dirs:
            dirs.add((prefix / d.relative_to(origin)).as_posix())
        for p in tree_files:
            relative = (prefix / p.relative_to(origin)).as_posix()
            if relative in destinations:
                raise ValueError(f"Conflicting copied paths: {relative}")
            destinations.add(relative)
            stat = p.stat()
            files.append({"source": str(p), "relative_path": relative, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    recorded_inputs = list(run.get("inputs", []))
    for key in ("fasta", "callable_bed", "erv_joined", "ld_block_dir", "architecture_segdup_bed", "architecture_lcr_bed", "architecture_low_mappability_bed"):
        if options.get(key): recorded_inputs.append(options[key])
    recorded_inputs.extend(options.get("enrichment_tracks", {}).values())
    for item in run.get("identity", {}).values():
        if isinstance(item, dict) and "path" in item:
            recorded_inputs.append(item["path"])
    for item in (options.get("architecture_annotations") or {}).get("files", {}).values():
        if isinstance(item, dict) and "path" in item:
            recorded_inputs.append(item["path"])
    for path in run.get("inputs", []):
        for suffix in (".tbi", ".csi"):
            indexed = Path(str(path) + suffix)
            if indexed.is_file(): recorded_inputs.append(str(indexed))
    dependencies = [{"path": str(p), "exists": p.exists(), "copied": False}
                    for p in sorted({absolute(p, source) for p in recorded_inputs}) if not p.is_relative_to(source)]
    return {"schema": SCHEMA, "status": "planned", "source_root": str(source), "store": str(store),
            "run_signature": run["signature"], "run_manifest": fingerprint(source / "run.json"),
            "path_mappings": mappings, "directories": sorted(dirs), "files": sorted(files, key=lambda r: r["relative_path"]),
            "file_count": len(files), "total_bytes": sum(r["size"] for r in files),
            "shards": [str(p) for p in original_shards], "external_inputs": dependencies,
            "validation": "original signatures plus file size/edge fingerprints; no new full-file checksums",
            "source_deletion": "never", "vcf_reads": 0, "rescan": False}


@contextmanager
def source_lock(source):
    path = source / ".allele3.lock"
    stream = path.open("rb") if path.is_file() else None
    try:
        if stream:
            try: fcntl.flock(stream, fcntl.LOCK_SH | fcntl.LOCK_NB)
            except BlockingIOError as exc: raise ValueError("Original scan is active; adoption was not started") from exc
        yield
    finally:
        if stream: stream.close()


def verify(store):
    root = absolute(store)
    run = load_run_metadata(root)
    for shard in run["shards"]: verify_shard(shard, run["signature"])
    files = run.get("_store", {}).get("files", {})
    for relative, saved in files.items():
        p = relative_target(root, relative)
        if not p.is_file():
            raise ValueError(f"Adopted file is missing: {p}")
        current = fingerprint(p)
        if not same_content(current, saved["source"]) or any(current[k] != saved["destination"][k] for k in FP_KEYS):
            raise ValueError(f"Adopted file changed: {p}; restore it from the untouched source")
    return {"schema": SCHEMA, "status": "verified", "store": str(root), "run_signature": run["signature"],
            "shards": len(run["shards"]), "verified_adopted_files": len(files), "vcf_reads": 0,
            "validation": "size/edge fingerprints and recorded destination timestamps; not full-file checksums"}


def _adopt_locked(source, store):
    source, store = absolute(source), absolute(store)
    with source_lock(source):
        plan = build_plan(source, store)
        identity = {k: plan[k] for k in ("schema", "source_root", "store", "run_signature", "run_manifest", "path_mappings", "files")}
        if store.exists():
            result = verify(store)
            saved = json.loads((store / "adoption_receipt.json").read_text())
            if saved.get("source_root") != str(source) or saved.get("run_signature") != plan["run_signature"]:
                raise ValueError("Existing store belongs to a different completed scan; choose a separate --store")
            journal = store / "adoption_state.json"
            if journal.is_file() and json.loads(journal.read_text()).get("identity") == identity:
                journal.unlink()  # A prior process may have stopped just after publication.
            return result | {"reused": True}
        for shard in plan["shards"]: verify_shard(shard, plan["run_signature"])
        store.parent.mkdir(parents=True, exist_ok=True)
        staging = store.with_name(store.name + ".adopting")
        if staging.exists():
            marker = staging / "adoption_state.json"
            if not marker.is_file():
                # Older versions removed this journal before the final rename.
                # Recover only a fully verified copy of precisely this source;
                # a partial/unrecognized directory remains untouched.
                try:
                    verify(staging)
                    saved = json.loads((staging / "adoption_receipt.json").read_text())
                    same_plan = all(saved.get(k) == plan[k] for k in
                                    ("schema", "source_root", "run_signature", "run_manifest", "path_mappings"))
                    saved_files = sorted((dict(source=r["source"]["path"], relative_path=r["relative_path"],
                                              size=r["source"]["size"], mtime_ns=r["source"]["mtime_ns"])
                                          for r in saved["files"]), key=lambda r: r["relative_path"])
                    if not same_plan or saved_files != plan["files"]:
                        raise ValueError("Completed copy does not match the current adoption plan")
                except (ValueError, OSError, KeyError, TypeError) as exc:
                    raise ValueError(f"Unrecognized/changed interrupted copy retained at {staging}; use another store or inspect it manually") from exc
                dump_json(marker, {"schema": SCHEMA, "status": "copying", "identity": identity,
                                   "started": now(), "recovered_verified_commit": True})
            if json.loads(marker.read_text()).get("identity") != identity:
                raise ValueError(f"Unrecognized/changed interrupted copy retained at {staging}; use another store or inspect it manually")
        else:
            staging.mkdir()
            dump_json(staging / "adoption_state.json", {"schema": SCHEMA, "status": "copying", "identity": identity, "started": now()})
        transfer = staging / ".adoption-transfer"
        transfer.mkdir(exist_ok=True)
        needed = sum(r["size"] for r in plan["files"] if not relative_target(staging, r["relative_path"]).is_file())
        if shutil.disk_usage(staging).free < needed:
            raise ValueError(f"Insufficient destination space: need up to {needed:,} additional bytes. Source remains intact; resume adoption after freeing space.")
        for relative in plan["directories"]: relative_target(staging, relative).mkdir(parents=True, exist_ok=True)
        receipt_files = []
        for index, item in enumerate(plan["files"], 1):
            original = Path(item["source"])
            before = fingerprint(original)
            if before["size"] != item["size"] or before["mtime_ns"] != item["mtime_ns"]:
                raise ValueError(f"Source changed since adoption planning: {original}; no source was removed")
            target = relative_target(staging, item["relative_path"])
            target.parent.mkdir(parents=True, exist_ok=True)
            current = fingerprint(target) if target.is_file() else None
            if current is None or not same_content(current, before):
                temporary = transfer / stable_hash(item["relative_path"])
                shutil.copy2(original, temporary)
                copied = fingerprint(temporary)
                if not same_content(copied, before):
                    raise ValueError(f"Copy verification failed: {original}. Source is retained; resume adoption, never rescan.")
                os.replace(temporary, target)
                current = fingerprint(target)
            if fingerprint(original) != before:
                raise ValueError(f"Source changed during copying: {original}; completed store was not published")
            receipt_files.append({"relative_path": item["relative_path"], "source": before,
                                  "destination": {k: current[k] for k in FP_KEYS}})
            if index == 1 or index % 100 == 0 or index == len(plan["files"]):
                print(f"[adopt] verified {index}/{len(plan['files'])} files; original scan is untouched", flush=True)
        for item in receipt_files:
            if fingerprint(item["source"]["path"]) != item["source"]:
                raise ValueError(f"Source changed before commit: {item['source']['path']}; copy retained for inspection")
        receipt = {k: plan[k] for k in ("schema", "source_root", "run_signature", "run_manifest", "path_mappings", "file_count", "total_bytes", "external_inputs", "validation")}
        receipt.update(status="complete", finished=now(), files=receipt_files, source_deleted=False, vcf_reads=0,
                       timestamp_rounding_files=sum(r["source"]["mtime_ns"] != r["destination"]["mtime_ns"] for r in receipt_files))
        dump_json(staging / "adoption_receipt.json", receipt)
        dump_json(staging / "relocation.json", {k: receipt[k] for k in ("schema", "status", "source_root", "run_signature", "run_manifest", "path_mappings")})
        # Validate the fully portable view before publishing the store directory.
        verify(staging)
        transfer.rmdir()
        os.replace(staging, store)
        # Keep the resumable journal until the atomic directory publication.
        # A stop after publication is harmless: verified reuse can remove it.
        (store / "adoption_state.json").unlink()
        return verify(store) | {"source_deleted": False, "copied_files": len(receipt_files), "reused": False}


def adopt(source, store):
    source, store = absolute(source), absolute(store)
    # Validate before creating a destination lock or copying anything. Rebuild
    # under the source read-lock below in case the source changes meanwhile.
    build_plan(source, store)
    store.parent.mkdir(parents=True, exist_ok=True)
    lock_path = store.with_name(store.name + ".adoption.lock")
    with lock_path.open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(f"Another adoption is using {store}; no second copy was started") from exc
        return _adopt_locked(source, store)


def scan_command(store, extra):
    """Original final profile, with expensive generated outputs on F: by default."""
    extra = [part for word in extra for part in (word.split('=', 1) if word.startswith('--') and '=' in word else [word])]
    def has(*names):
        return any(word.split("=", 1)[0] in names for word in extra)
    command = [sys.executable, str(Path(__file__).with_name("scan.py")), "run", "--out", str(store)]
    defaults = [("--fasta", os.environ.get("ALLELE_FASTA", "/mnt/f/gen/1hgp/38/GRCH38.fasta.gz")),
                ("--jobs", os.environ.get("ALLELE_JOBS", "2")), ("--min-an", "100"), ("--min-af", "0"),
                ("--thresholds", "1,2,5,10"), ("--backend", "auto"), ("--candidate-storage", "deferred"),
                ("--ld-block-dir", os.environ.get("ALLELE_LD_BLOCK_DIR", "/mnt/f/LD/block")),
                ("--architecture-build", "38"),
                ("--erv-joined", os.environ.get("ALLELE_ERV_JOINED", "/mnt/f/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz")),
                ("--architecture-annotation-cache", str(store / "annotations")),
                ("--architecture-low-mappability-bed", os.environ.get("ALLELE_LOW_MAPPABILITY_BED", "/mnt/d/files/allele3/umap/hg38.k100.lt1.bed.gz"))]
    if not has("--vcf", "--dir-vcf", "--vcf-dir"):
        command += ["--dir-vcf", os.environ.get("ALLELE_VCF_DIR", "/mnt/f/annot/gnomad/sites")]
    if not has("--dir-cache", "--cache-dir"):
        # Do not inherit an old D:-drive ALLELE_CACHE_DIR for a new permanent scan.
        command += ["--dir-cache", str(store / "shards")]
    for option, value in defaults:
        if not has(option): command += [option, value]
    if not has("--with-architecture"): command += ["--with-architecture"]
    if has("--out", "--dir-out"):
        raise ValueError("Set the scan output using --store; do not pass a second scanner --out")
    return command + list(extra)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, epilog="Only the explicit 'run' command can scan VCF. Extra run arguments are forwarded to scan.py.")
    parser.add_argument("command", choices=["plan", "adopt", "verify", "run"], nargs="?", default="plan")
    parser.add_argument("--source", default=os.environ.get("ALLELE_SCAN_SOURCE", DEFAULT_SOURCE))
    parser.add_argument("--store", "--dir-out", "--out", dest="store", default=os.environ.get("ALLELE_SCAN_STORE", DEFAULT_STORE))
    parser.add_argument("--replace", choices=["TRUE","FALSE"], default="TRUE")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--preflight", action="store_true")
    args, extra = parser.parse_known_args(argv)
    if extra[:1] == ["--"]: extra = extra[1:]
    if extra and args.command != "run": parser.error("Scanner arguments are accepted only after the explicit run command")
    source, store = absolute(args.source), absolute(args.store)
    if args.command == "run":
        command = scan_command(store, extra)
        if args.dry_run or args.preflight:
            print(json.dumps({"status": "planned_explicit_scan", "store": str(store), "command": shlex.join(command),
                              "existing_manifest": (store / "run.json").exists(), "vcf_reads": 0}, indent=2))
            return 0
        from fresh_run import output_run
        preserved = os.environ.get('ALLELE_SCAN_STORE', DEFAULT_STORE)
        input_options = {'--dir-vcf','--vcf-dir','--vcf','--fasta','--erv-joined','--ld-block-dir','--architecture-low-mappability-bed'}
        # The selected store is an output, not a source. A fresh run may use
        # that exact path; output_run archives any previous products there.
        protected = [] if store == absolute(preserved) else [preserved]
        protected += [str(absolute(preserved).with_name(absolute(preserved).name + '-OLD'))]
        protected += [command[i+1] for i,x in enumerate(command[:-1]) if x in input_options]
        with output_run(store, args.replace == 'TRUE', protected):
            print("[scan] Full scan requested: " + shlex.join(command), flush=True)
            return subprocess.run(command, check=False).returncode
    if args.command == "verify":
        if args.dry_run or args.preflight:
            metadata = load_run_metadata(store)
            result = {"status": "planned_verification", "store": str(store), "run_signature": metadata["signature"],
                      "shards": len(metadata["shards"]), "vcf_reads": 0}
        else:
            result = verify(store)
    elif args.command == "adopt" and not (args.dry_run or args.preflight):
        result = adopt(source, store)
    else:
        result = build_plan(source, store)
        result = {k: v for k, v in result.items() if k not in {"files", "directories"}}
        result["commands"] = [f"./3.scan.sh adopt --source {shlex.quote(str(source))} --store {shlex.quote(str(store))}",
                              f"./3.scan.sh verify --store {shlex.quote(str(store))}"]
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        raise SystemExit(f"scan-store: {exc}")
