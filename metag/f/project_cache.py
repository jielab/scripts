"""Reusable project preparation; never sample classifications or abundance estimates.

Planning caches bounded FASTQ probes. Compute workers compile the unchanged NCBI
lineage interpretation once into an immutable, read-only SQLite lookup. No cache
operation reads or hashes the large classifier indexes.
"""
from contextlib import ExitStack
import gzip
import itertools
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
import time
from urllib.parse import quote

from common import atomic_json, fingerprint, locked, utcnow
from fastq_stream import read_id, validate_pair_ids
from taxonomy import Taxonomy, safe_name

SCHEMA = 1
FASTQ_RE = re.compile(r"(.+?)(?:_R?([12])(?:_001)?|\.([12]))?\.fastq\.gz$", re.I)
PROBE_FRAGMENTS = 4


def _read(path):
    try:
        value = json.loads(Path(path).read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, TypeError):
        return {}


def _stamp(path):
    path = Path(path).resolve(strict=True)
    st = path.stat()
    return {"path": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns,
            "device": st.st_dev, "inode": st.st_ino}


def _root(cfg):
    path = Path(cfg["work_dir"]).absolute() / "cache"
    if path.resolve() != path:
        raise RuntimeError(f"Project cache must not use a symlink/alias: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _project(cfg, dataset):
    if not dataset or dataset in (".", "..") or re.search(r"[/\\\x00-\x1f]", dataset):
        raise ValueError(f"Unsafe project name: {dataset!r}")
    path = _root(cfg) / "projects" / dataset
    if path.resolve() != path:
        raise RuntimeError(f"Project cache must not use a symlink/alias: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def settings(cfg):
    """Scientific/storage protocol only; scheduling and --replace are independent."""
    keys = ("classifier_engine", "confidence", "minimum_hit_groups", "minimum_base_quality",
            "materialize_species_min_fragments", "always_materialize_taxids")
    keys += tuple(sorted(key for key in cfg if key.startswith("human_host_") and key != "human_host_threads"))
    return {"schema": SCHEMA, "fastq": "four-line paired .fastq.gz", "probe_fragments": PROBE_FRAGMENTS,
            "parameters": {key: cfg.get(key) for key in keys}}


def _scan(folder):
    directories, paths = {}, []
    for directory, children, files in os.walk(folder, followlinks=False):
        children[:] = sorted(name for name in children if not (Path(directory) / name).is_symlink())
        directories[str(Path(directory))] = _stamp(directory)
        paths.extend(str(Path(directory) / name) for name in sorted(files) if FASTQ_RE.fullmatch(name))
    return directories, paths


def inventory(cfg, dataset):
    """Stat known entries; rescan directories only when their contents changed."""
    started = time.monotonic()
    folder = (Path(cfg["reads_raw"]) / dataset).resolve(strict=True)
    cache = _project(cfg, dataset)
    with locked(cache / "prepare.lock"):
        saved = _read(cache / "prepare.json")
        valid = (saved.get("schema") == SCHEMA and saved.get("source") == str(folder)
                 and saved.get("inventory_complete") is True and isinstance(saved.get("directories"), dict)
                 and isinstance(saved.get("files"), dict))
        directory_hit = valid
        if directory_hit:
            try:
                directory_hit = bool(saved["directories"]) and all(
                    _stamp(path) == value for path, value in saved["directories"].items())
            except OSError:
                directory_hit = False
        if directory_hit:
            directories, paths = saved["directories"], sorted(saved["files"])
        else:
            directories, paths = _scan(folder)
        files = {}
        for path in paths:
            try:
                if Path(path).is_file():
                    state = _stamp(path)
                    state["downloading"] = (Path(path + ".aria2").exists() or
                                            Path(state["path"] + ".aria2").exists())
                    files[path] = state
            except FileNotFoundError:
                continue  # The ordinary discovery code reports incomplete arrivals.
        protocol = settings(cfg)
        hit = directory_hit and files == saved.get("files") and protocol == saved.get("settings")
        record = saved if valid else {"schema": SCHEMA, "source": str(folder), "probes": {}}
        record.update(inventory_complete=True, directories=directories, files=files, settings=protocol)
        if not hit:
            record["inventory_revision"] = fingerprint(files)
            record["updated_at"] = utcnow()
            atomic_json(cache / "prepare.json", record)
    info = {"inventory_hit": hit, "directory_scan": not directory_hit, "files": len(files),
            "inventory_seconds": time.monotonic() - started, "probe_hits": 0, "probed_samples": 0}
    cfg.setdefault("project_cache_discovery", {})[dataset] = info
    return [Path(path) for path in sorted(files)]


def _record(handle, path):
    header = handle.readline(1 << 20)
    if not header:
        return None
    rest = [handle.readline(1 << 20) for _ in range(3)]
    if any(len(line) == 1 << 20 for line in [header, *rest]):
        raise RuntimeError(f"FASTQ probe line exceeds 1 MiB: {path}")
    sequence, plus, quality = rest
    read_id(header)
    if not sequence or not plus.startswith(b"+") or not quality:
        raise RuntimeError(f"Invalid/truncated FASTQ prefix: {path}")
    size = len(sequence.rstrip(b"\r\n"))
    if not size or size != len(quality.rstrip(b"\r\n")):
        raise RuntimeError(f"Sequence/quality length mismatch in FASTQ prefix: {path}")
    return header, size


def _probe(sample):
    lengths = [[], []]
    with ExitStack() as stack:
        handles = [stack.enter_context(gzip.open(sample[key], "rb")) for key in ("r1", "r2")]
        for _ in range(PROBE_FRAGMENTS):
            a, b = [_record(handle, sample[key]) for handle, key in zip(handles, ("r1", "r2"))]
            if a is None or b is None:
                if a is not b:
                    raise RuntimeError(f"FASTQ mate counts differ in prefix: {sample['sample_id']}")
                break
            validate_pair_ids(a[0], b[0])
            lengths[0].append(a[1]); lengths[1].append(b[1])
    if not lengths[0]:
        raise RuntimeError(f"FASTQ contains no records: {sample['sample_id']}")
    return {"fragments_checked": len(lengths[0]), "format": "paired four-line FASTQ",
            "r1_length_range": [min(lengths[0]), max(lengths[0])],
            "r2_length_range": [min(lengths[1]), max(lengths[1])],
            "scope": "bounded prefix only; every fragment is revalidated during classification"}


def probe_samples(cfg, samples):
    """Probe selected ready pairs; unaffected pairs retain their independent cache."""
    for dataset in sorted({sample["dataset"] for sample in samples}):
        cache = _project(cfg, dataset)
        started = time.monotonic()
        hits = rebuilt = 0
        with locked(cache / "prepare.lock"):
            saved = _read(cache / "prepare.json")
            if saved.get("schema") != SCHEMA or not saved.get("inventory_complete"):
                raise RuntimeError(f"Project inventory is incomplete: {dataset}")
            if not isinstance(saved.get("probes"), dict):
                saved["probes"] = {}
            probes = saved["probes"]
            for sample in (item for item in samples if item["dataset"] == dataset):
                states = {key: _stamp(sample[key]) for key in ("r1", "r2")}
                previous = probes.get(sample["sample"], {})
                if (previous.get("inputs") == states and previous.get("probe_fragments") == PROBE_FRAGMENTS and
                        previous.get("status") == "complete" and isinstance(previous.get("result"), dict)):
                    hits += 1
                    continue
                result = _probe(sample)
                if states != {key: _stamp(sample[key]) for key in ("r1", "r2")}:
                    raise RuntimeError(f"FASTQ changed during project preparation: {sample['sample_id']}")
                probes[sample["sample"]] = {"status": "complete", "inputs": states,
                    "probe_fragments": PROBE_FRAGMENTS, "result": result}
                rebuilt += 1
            if rebuilt:
                saved["updated_at"] = utcnow()
                atomic_json(cache / "prepare.json", saved)
        info = cfg.setdefault("project_cache_discovery", {}).setdefault(dataset, {})
        info.update(probe_hits=hits, probed_samples=rebuilt, probe_seconds=time.monotonic() - started)


def _taxonomy_state(snapshot):
    folder = Path(snapshot["taxonomy_dir"])
    return {name: _stamp(folder / name) for name in ("nodes.dmp", "names.dmp", "merged.dmp")
            if (folder / name).is_file()}


class PreparedTaxonomy:
    """One indexed lookup per distinct observed taxid, then in-memory dictionaries."""
    def __init__(self, path, metadata, snapshot, source_states):
        self.connection = sqlite3.connect("file:" + quote(str(path)) + "?mode=ro&immutable=1", uri=True)
        self.connection.execute("PRAGMA query_only=ON")
        if self.connection.execute("SELECT identity FROM metadata").fetchone() != (metadata["taxonomy_cache_identity"],):
            self.connection.close()
            raise RuntimeError("Prepared taxonomy identity changed while opening its immutable snapshot")
        self._cache = {}
        self.preparation = metadata
        self._snapshot = {"taxonomy_dir": snapshot["taxonomy_dir"]}
        self._source_states = source_states

    def validate_current(self):
        if self._source_states != _taxonomy_state(self._snapshot):
            raise RuntimeError("Reference taxonomy changed after opening the prepared routing snapshot")

    def describe(self, taxid):
        taxid = int(taxid)
        if taxid not in self._cache:
            row = self.connection.execute("SELECT assigned, species, rank, grp, name, safe_name, status, human_ancestor, vertebrate_lineage FROM routing WHERE taxid=?", (taxid,)).fetchone()
            if row is None:
                raise RuntimeError(f"Classification taxid absent from database taxonomy: {taxid}")
            self._cache[taxid] = dict(taxid=row[0], species_taxid=row[1], rank=row[2], group=row[3],
                name=row[4], safe_name=row[5], status=row[6], human_ancestor=bool(row[7]), vertebrate_lineage=bool(row[8]))
        return self._cache[taxid]

    def species(self, taxid):
        return self.describe(taxid)["species_taxid"]

    def close(self):
        self.connection.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def _compile(snapshot, path, identity):
    taxonomy = Taxonomy(snapshot["taxonomy_dir"])
    with sqlite3.connect(path) as connection:
        connection.execute("PRAGMA journal_mode=OFF")  # private, unpublished build
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA user_version=1")
        connection.execute("CREATE TABLE routing (taxid INTEGER PRIMARY KEY, assigned INTEGER, species INTEGER, rank TEXT, grp TEXT, name TEXT, safe_name TEXT, status TEXT, human_ancestor INTEGER, vertebrate_lineage INTEGER)")
        connection.execute("CREATE TABLE metadata (identity TEXT, rows INTEGER)")
        pending, count = [], 0
        ids = itertools.chain((0,), taxonomy.parent, (tid for tid in taxonomy.merged if tid not in taxonomy.parent))
        for tid in ids:
            info = taxonomy.describe(tid)
            pending.append((tid, info["taxid"], info["species_taxid"], info["rank"], info["group"],
                info["name"], safe_name(info["name"]), info["status"],
                int(info.get("human_ancestor", False)), int(info.get("vertebrate_lineage", False))))
            count += 1
            if len(pending) >= 10_000:
                connection.executemany("INSERT INTO routing VALUES (?,?,?,?,?,?,?,?,?,?)", pending)
                connection.commit()
                pending.clear()
                taxonomy._cache = {0: taxonomy._cache[0]}
        if pending:
            connection.executemany("INSERT INTO routing VALUES (?,?,?,?,?,?,?,?,?,?)", pending)
        connection.execute("INSERT INTO metadata VALUES (?,?)", (identity, count))
        connection.commit()
    connection.close()
    return count


def _bind_projects(cfg, datasets, snapshot, taxonomy_identity):
    reference = {"identity": snapshot.get("identity"), "database": str(Path(snapshot["db_dir"]).resolve()),
                 "taxonomy_cache_identity": taxonomy_identity}
    # Sizes/mtimes only. Library readiness is still established by its ordinary
    # quick check; the cache cannot declare a missing/truncated index usable.
    reference["files"] = {name: _stamp(Path(snapshot["db_dir"]) / name) for name in snapshot.get("files", {})}
    result = {}
    for dataset in sorted(set(datasets)):
        cache = _project(cfg, dataset)
        with locked(cache / "prepare.lock"):
            saved = _read(cache / "prepare.json")
            if saved.get("schema") != SCHEMA:
                saved = {"schema": SCHEMA}
            protocol = settings(cfg)
            hit = (saved.get("reference") == reference and saved.get("bound_settings") == protocol and
                   saved.get("bound_inventory_revision") == saved.get("inventory_revision"))
            if not hit:
                saved.update(reference=reference, bound_settings=protocol,
                             bound_inventory_revision=saved.get("inventory_revision"), updated_at=utcnow())
                atomic_json(cache / "prepare.json", saved)
        result[dataset] = {"preparation_hit": hit,
            "discovery": cfg.get("project_cache_discovery", {}).get(dataset, {})}
    return result


def load_taxonomy(cfg, snapshot, datasets):
    """Compute-worker-only compilation; plan generation never calls this function."""
    started = time.monotonic()
    folder = _root(cfg) / "taxonomy"
    folder.mkdir(exist_ok=True)
    if folder.resolve() != folder:
        raise RuntimeError(f"Taxonomy cache must not use a symlink/alias: {folder}")
    path = folder / "routing.sqlite"
    rebuilt_seconds = 0.0
    with locked(folder / "prepare.lock"):
        # Capture after waiting: a different worker/reference update may have
        # completed while this process was blocked on the preparation lock.
        states = _taxonomy_state(snapshot)
        identity = fingerprint({"schema": SCHEMA, "taxonomy": states})
        saved = _read(folder / "complete.json")
        hit = (saved.get("schema") == SCHEMA and saved.get("identity") == identity and saved.get("status") == "complete"
               and path.is_file() and path.stat().st_size == saved.get("size") and
               isinstance(saved.get("rows"), int) and saved["rows"] > 0)
        connection = None
        if hit:
            try:
                connection = sqlite3.connect("file:" + quote(str(path)) + "?mode=ro&immutable=1", uri=True)
                hit = (connection.execute("PRAGMA user_version").fetchone()[0] == SCHEMA and
                       connection.execute("SELECT identity, rows FROM metadata").fetchone() == (identity, saved["rows"]))
            except sqlite3.DatabaseError:
                hit = False
            finally:
                if connection is not None:
                    connection.close()
        if not hit:
            build_started = time.monotonic()
            fd, name = tempfile.mkstemp(prefix=".routing-", suffix=".sqlite", dir=folder)
            os.close(fd)
            temporary = Path(name)
            try:
                count = _compile(snapshot, temporary, identity)
                if states != _taxonomy_state(snapshot):
                    raise RuntimeError("Reference taxonomy changed during project cache preparation")
                with temporary.open("rb") as handle:
                    os.fsync(handle.fileno())
                os.replace(temporary, path)
                rebuilt_seconds = time.monotonic() - build_started
                saved = {"schema": SCHEMA, "status": "complete", "identity": identity, "taxonomy": states,
                         "size": path.stat().st_size, "rows": count, "build_seconds": rebuilt_seconds, "finished_at": utcnow()}
                atomic_json(folder / "complete.json", saved)
            finally:
                temporary.unlink(missing_ok=True)
        metadata = {"taxonomy_cache_hit": hit, "taxonomy_cache_path": str(path), "taxonomy_cache_identity": identity,
                    "taxonomy_cache_rows": saved["rows"], "taxonomy_cache_build_seconds": rebuilt_seconds}
        if states != _taxonomy_state(snapshot):
            raise RuntimeError("Reference taxonomy changed while checking the preparation cache")
        taxonomy = PreparedTaxonomy(path, metadata, snapshot, states)  # bind the immutable inode while locked
    try:
        metadata["projects"] = _bind_projects(cfg, datasets, snapshot, identity)
        taxonomy.validate_current()
        metadata["project_cache_hit"] = hit and all(item["preparation_hit"] for item in metadata["projects"].values())
        metadata["project_preparation_seconds"] = time.monotonic() - started
        return taxonomy
    except BaseException:
        taxonomy.close()
        raise
