"""Pinned Kun-peng installation and lossless conversion of an uncapped Kraken index.

The official hashshard command is run only in a private directory. Conversion
preserves every hash cell and the native taxonomy/options, so method comparisons
use the same reference evidence rather than separately built databases.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import struct
import subprocess
import tarfile
import time

from hpc_common import (software_sources, CODE_ROOT, atomic_json, check_space, file_state, fingerprint,
                        locked, sha256, utcnow)
import hpc_reference

SOURCE_COMMIT = "034dabb547956be472a412a1521a207c30afef99"
VERSION = "0.7.13"
OWNER = ".metag-kunpeng-owner.json"
SCHEMA = 1


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _sync_dir(path):
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _record(path, digest=None):
    return {"state": file_state(path), "sha256": digest or sha256(path)}


def _check(path, record):
    if Path(path).is_symlink() or not Path(path).is_file() or file_state(path) != record["state"]:
        raise RuntimeError(f"Kun-peng file changed since publication: {path}")


def _extract(archive, target, prefix):
    """Extract only regular files/directories from our checksum-verified archives."""
    target.mkdir(parents=True, exist_ok=False)
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            name = Path(member.name)
            if name.is_absolute() or ".." in name.parts or not name.parts or name.parts[0] != prefix:
                raise RuntimeError(f"Unsafe bundled archive member: {member.name}")
            if not (member.isfile() or member.isdir()):
                raise RuntimeError(f"Unsupported bundled archive member: {member.name}")
            destination = target / Path(*name.parts[1:])
            if member.isdir():
                destination.mkdir(parents=True, exist_ok=True)
            else:
                destination.parent.mkdir(parents=True, exist_ok=True)
                with bundle.extractfile(member) as inp, destination.open("wb") as out:
                    shutil.copyfileobj(inp, out)
                destination.chmod(member.mode & 0o777)


def _verify_binary(path):
    if not Path(path).is_file() or not os.access(path, os.X_OK):
        raise RuntimeError(f"Missing Kun-peng executable: {path}")
    version = subprocess.check_output([str(path), "--version"], text=True, timeout=30).strip()
    if version != f"kun_peng {VERSION}":
        raise RuntimeError(f"Expected pinned Kun-peng {VERSION}; got {version!r}")
    for command in ("hashshard", "direct", "splitr", "annotate", "resolve"):
        subprocess.run([str(path), command, "--help"], check=True, timeout=30,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return version


def install(cfg):
    """Compile from vendored crates offline, or import an explicitly selected binary."""
    vendor = CODE_ROOT / "f" / "vendor"
    manifest = software_sources(CODE_ROOT)["kunpeng"]
    if manifest["source_commit"] != SOURCE_COMMIT:
        raise RuntimeError("Unsupported Kun-peng source commit")
    prefix = Path(cfg["tools_dir"]) / "kunpeng"
    if prefix.is_symlink():
        raise RuntimeError(f"Refusing symlink installation directory: {prefix}")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    mode = cfg.get("kunpeng_install_mode", "binary")
    if mode == "source":
        raise ValueError("This clean deployment uses the official binary; set kunpeng_install_mode=binary")
    if mode not in ("binary", "source"):
        raise ValueError("kunpeng_install_mode must be binary or source")
    official = not cfg.get("kunpeng_binary") and mode == "binary"
    if official:
        if platform.system() != "Linux" or platform.machine().lower() not in ("x86_64", "amd64"):
            raise RuntimeError("Bundled official Kun-peng binary is Linux x86_64; use kunpeng_install_mode=source")
        external = vendor / manifest["binary_file"]
    else:
        external = Path(cfg["kunpeng_binary"]).expanduser().absolute() if cfg.get("kunpeng_binary") else None
    external_hash = sha256(external) if external else None
    if official and external_hash != manifest["binary_sha256"]:
        raise RuntimeError("Bundled official Kun-peng binary checksum mismatch")
    if external and cfg.get("kunpeng_binary_sha256") and external_hash != cfg["kunpeng_binary_sha256"]:
        raise RuntimeError("Explicit Kun-peng binary checksum mismatch")
    if external is None and not all(manifest.get(key) for key in ("dependencies_archive", "dependencies_sha256")):
        raise RuntimeError("Offline Kun-peng source compilation needs the bundled dependency archive; "
                           "use the official binary installation or supply the complete source bundle")
    request = {"source_commit": SOURCE_COMMIT, "source_sha256": manifest["sha256"],
               "dependencies_sha256": manifest.get("dependencies_sha256") if external is None else None,
               "external_binary_sha256": external_hash,
               "install_mode": "official_release_binary" if official else ("explicit_binary" if external else "source")}
    with locked(prefix.parent / "kunpeng.install.lock"):
        owner = prefix / OWNER
        if prefix.exists() and not owner.is_file():
            raise RuntimeError(f"Unowned Kun-peng installation directory: {prefix}")
        prefix.mkdir(exist_ok=True)
        if owner.exists() and _read(owner).get("request") != request:
            raise RuntimeError("Different Kun-peng installation owns this directory; choose a new tools_dir")
        if not owner.exists():
            atomic_json(owner, {"schema": SCHEMA, "request": request})
        receipt = prefix / "install.json"
        executable = prefix / "bin/kun_peng"
        if receipt.exists():
            result = _read(receipt)
            if result.get("request") != request or sha256(executable) != result.get("binary_sha256"):
                raise RuntimeError("Published Kun-peng installation changed")
            print(f"REUSE verified Kun-peng: {executable}", flush=True)
            return result
        executable.parent.mkdir(exist_ok=True)
        log = prefix / "install.log"
        if external:
            temporary = executable.with_suffix(".partial")
            shutil.copyfile(external, temporary)
            temporary.chmod(0o755)
            # Windows uploads may discard executable mode on bundled files.
            # Verify the owned executable copy without modifying the source.
            version = _verify_binary(temporary)
            temporary.replace(executable)
            provenance = "official_release_asset" if official else "explicit_external_binary_source_not_verified"
        else:
            cargo = cfg.get("kunpeng_cargo") or shutil.which("cargo")
            if not cargo:
                raise RuntimeError("Kun-peng requires Rust/Cargo. Load the site's Rust module in env.sh, "
                                   "set kunpeng_cargo, or provide an explicitly verified kunpeng_binary. "
                                   "All source dependencies are bundled for offline compilation.")
            for program in ("cc", "make"):
                if not shutil.which(program):
                    raise RuntimeError(f"Missing {program}; Kun-peng's bundled jemalloc needs a C compiler and make")
            archive = vendor / manifest["archive"]
            dependencies = vendor / manifest["dependencies_archive"]
            if sha256(archive) != manifest["sha256"] or sha256(dependencies) != manifest["dependencies_sha256"]:
                raise RuntimeError("Bundled Kun-peng source/dependency checksum mismatch")
            work = prefix / "build"
            # The marked installation owns this entire incomplete build area.
            if work.exists():
                if work.is_symlink():
                    raise RuntimeError("Refusing symlink in Kun-peng build area")
                shutil.rmtree(work)
            work.mkdir()
            src = work / "source"
            _extract(archive, src, "kunpeng")
            _extract(dependencies, src / "vendor", "vendor")
            (src / ".cargo").mkdir(exist_ok=True)
            (src / ".cargo/config.toml").write_text(
                '[source.crates-io]\nreplace-with = "vendored-sources"\n'
                '[source.vendored-sources]\ndirectory = "vendor"\n'
                '[net]\noffline = true\n', encoding="utf-8")
            environment = dict(os.environ)
            environment["CARGO_TARGET_DIR"] = str(work / "target")
            environment["CARGO_HOME"] = str(work / "cargo-home")
            # Do not inherit architecture-specific instructions from a login node.
            environment.pop("RUSTFLAGS", None)
            environment.pop("CARGO_ENCODED_RUSTFLAGS", None)
            command = [str(cargo), "build", "--locked", "--offline", "--release", "--bin", "kun_peng",
                       "--jobs", str(int(cfg.get("setup_cores", 8)))]
            print(f"Compile pinned Kun-peng offline; log={log}", flush=True)
            with log.open("ab") as out:
                result = subprocess.run(command, cwd=src, env=environment, stdout=out, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f"Kun-peng compiler exited {result.returncode}; inspect {log}")
            compiled = work / "target/release/kun_peng"
            version = _verify_binary(compiled)
            temporary = executable.with_suffix(".partial")
            shutil.copyfile(compiled, temporary)
            temporary.chmod(0o755)
            temporary.replace(executable)
            provenance = "compiled_pinned_source_with_locked_vendored_dependencies"
        with executable.open("rb") as handle:
            os.fsync(handle.fileno())
        _sync_dir(executable.parent)
        result = {"schema": SCHEMA, "request": request, "source_commit": SOURCE_COMMIT,
                  "reported_version": version, "binary_sha256": sha256(executable),
                  "binary_path": str(executable), "provenance": provenance,
                  "source_archive_sha256": manifest["sha256"],
                  "dependencies_archive_sha256": manifest.get("dependencies_sha256") if external is None else None,
                  "binary_url": manifest.get("binary_url") if official else None,
                  "installed_at": utcnow(), "platform": platform.platform()}
        atomic_json(receipt, result)
        return result


def binary(cfg):
    path = Path(cfg["tools_dir"]) / "kunpeng/bin/kun_peng"
    receipt = _read(path.parent.parent / "install.json")
    if receipt.get("source_commit") != SOURCE_COMMIT or sha256(path) != receipt.get("binary_sha256"):
        raise RuntimeError("Kun-peng executable does not match its recorded installation")
    return path


def _destination(cfg):
    result = Path(cfg.get("kunpeng_db_dir") or Path(cfg["library_dir"]) / "kunpeng").absolute()
    source = Path(cfg.get("db_dir") or Path(cfg["library_dir"]) / "db").resolve()
    resolved = result.resolve()
    for key in ("reference_raw", "reads_raw", "tools_dir"):
        protected = Path(cfg[key]).resolve()
        if resolved == protected or resolved in protected.parents or protected in resolved.parents:
            raise RuntimeError(f"Kun-peng database overlaps protected path: {result}")
    if resolved == source or resolved in source.parents or source in resolved.parents or result.is_symlink():
        raise RuntimeError("Kun-peng database must be separate from the Kraken database")
    return result


def _capacity(cfg):
    value = cfg.get("kunpeng_hash_capacity", 1073741824)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1024 or value % 1024:
        raise ValueError("kunpeng_hash_capacity must be a positive multiple of 1024 slots")
    return value


def _native_header(source):
    path = Path(source["db_dir"]) / "hash.k2d"
    with path.open("rb") as handle:
        header = handle.read(32)
    if len(header) != 32:
        raise RuntimeError("Truncated Kraken hash header")
    capacity, used, key_bits, value_bits = struct.unpack("<QQQQ", header)
    if not capacity or used > capacity or not 0 < value_bits < 32 or key_bits + value_bits != 32:
        raise RuntimeError("Unsupported Kraken hash header")
    if path.stat().st_size != 32 + capacity * 4:
        raise RuntimeError("Kraken hash file length differs from its header")
    return header, capacity, used, value_bits


def _verify_shards(folder, source, capacity_per_shard):
    header, capacity, used, value_bits = _native_header(source)
    partitions = math.ceil(capacity / capacity_per_shard)
    expected_config = (0, partitions, capacity_per_shard, capacity, used, value_bits)
    config = folder / "hash_config.k2d"
    if config.read_bytes() != struct.pack("<QQQQQQ", *expected_config):
        raise RuntimeError("Kun-peng converted hash configuration does not match Kraken")
    expected = {f"hash_{number}.k2d" for number in range(1, partitions + 1)}
    observed = {p.name for p in folder.glob("hash_*.k2d") if p.name != "hash_config.k2d"}
    if observed != expected:
        raise RuntimeError("Kun-peng conversion has missing or unexpected hash shards")
    aggregate = hashlib.sha256(header)
    records = {"hash_config.k2d": _record(config)}
    for number in range(1, partitions + 1):
        path = folder / f"hash_{number}.k2d"
        cells = min(capacity_per_shard, capacity - (number - 1) * capacity_per_shard)
        if path.is_symlink() or path.stat().st_size != 16 + cells * 4:
            raise RuntimeError(f"Wrong converted shard length: {path}")
        checksum = hashlib.sha256()
        with path.open("rb") as handle:
            shard_header = handle.read(16)
            if shard_header != struct.pack("<QQ", number, cells):
                raise RuntimeError(f"Invalid converted shard header: {path}")
            checksum.update(shard_header)
            for block in iter(lambda: handle.read(8 << 20), b""):
                checksum.update(block)
                aggregate.update(block)
            os.fsync(handle.fileno())
        records[path.name] = _record(path, checksum.hexdigest())
    if aggregate.hexdigest() != source["files"]["hash.k2d"]["sha256"]:
        raise RuntimeError("Converted hash cells differ from the original Kraken index")
    return records, {"partitions": partitions, "capacity": capacity, "size": used,
                     "value_bits": value_bits, "hash_capacity": capacity_per_shard,
                     "hash_cells_equivalent_to_source_sha256": aggregate.hexdigest()}


def database_snapshot(cfg):
    """Validate immutable receipt and file states without rehashing TB-sized shards."""
    folder = _destination(cfg)
    receipt = _read(folder / "database.json")
    context = receipt.get("context", {})
    if receipt.get("schema") != SCHEMA or receipt.get("status") != "complete" or \
            fingerprint(context) != receipt.get("identity") or receipt.get("db_dir") != str(folder):
        raise RuntimeError("Kun-peng database has no valid complete receipt")
    for key in ("files", "taxonomy_files", "parameters", "source_database_identity"):
        if receipt.get(key) != context.get(key):
            raise RuntimeError(f"Kun-peng receipt does not match its context: {key}")
    if receipt.get("taxonomy_dir") != str(folder / "taxonomy"):
        raise RuntimeError("Kun-peng frozen taxonomy path changed")
    names = set(receipt["files"])
    partitions = receipt["parameters"]["partitions"]
    expected = {"hash_config.k2d", "taxo.k2d", "opts.k2d"} | {f"hash_{i}.k2d" for i in range(1, partitions + 1)}
    if names != expected:
        raise RuntimeError("Kun-peng receipt has missing or extra hash files")
    actual_hash = {p.name for p in folder.glob("hash*.k2d")}
    if actual_hash != {name for name in expected if name.startswith("hash")}:
        raise RuntimeError("Kun-peng directory has missing or unrecorded hash files")
    if not set(hpc_reference.TAX_FILES).issubset(receipt["taxonomy_files"]):
        raise RuntimeError("Kun-peng frozen taxonomy receipt is incomplete")
    for name, record in receipt["files"].items():
        _check(folder / name, record)
    for name, record in receipt["taxonomy_files"].items():
        if Path(name).name != name:
            raise RuntimeError("Invalid taxonomy filename")
        _check(folder / "taxonomy" / name, record)
    hpc_reference._native_opts(folder / "opts.k2d")
    return receipt


def convert_database(cfg):
    """Run official hashshard once and prove that no reference hash cells were lost."""
    folder = _destination(cfg)
    capacity_per_shard = _capacity(cfg)
    executable = binary(cfg)
    install_receipt = _read(executable.parent.parent / "install.json")
    source_dir = Path(cfg.get("db_dir") or Path(cfg["library_dir"]) / "db")
    folder.parent.mkdir(parents=True, exist_ok=True)
    with locked(Path(cfg["library_dir"]) / "reference.lock", shared=True), \
            locked(source_dir / "build.lock", shared=True), \
            locked(folder.parent / (folder.name + ".convert.lock")):
        source = hpc_reference.database_snapshot(cfg)
        request = {"source_database_identity": source["identity"], "hash_capacity": capacity_per_shard,
                   "binary_sha256": install_receipt["binary_sha256"], "source_commit": SOURCE_COMMIT}
        if (folder / "database.json").exists():
            receipt = database_snapshot(cfg)
            if receipt["context"].get("request") != request:
                raise RuntimeError("Existing Kun-peng database belongs to another conversion; choose a new kunpeng_db_dir")
            print(f"REUSE verified Kun-peng conversion: {folder}", flush=True)
            return receipt
        if folder.exists():
            raise RuntimeError(f"Unpublished or unowned final Kun-peng database directory: {folder}")
        stage = folder.with_name(folder.name + ".partial")
        if stage.exists():
            if stage.is_symlink() or not (stage / OWNER).is_file() or _read(stage / OWNER).get("request") != request:
                raise RuntimeError(f"Refusing unowned or different conversion stage: {stage}")
            shutil.rmtree(stage)
        index_bytes = source["files"]["hash.k2d"]["state"]["size"]
        check_space(folder.parent, float(cfg.get("min_free_gib", 100)) + index_bytes / 2**30 + 1)
        budget = float(cfg["convert_memory_gib"])
        if capacity_per_shard * 4 / 2**30 + 2 > budget:
            raise RuntimeError("One Kun-peng hash shard plus overhead exceeds conversion memory budget")
        stage.mkdir()
        atomic_json(stage / OWNER, {"schema": SCHEMA, "request": request})
        (stage / "hash.k2d").symlink_to(Path(source["db_dir"]) / "hash.k2d")
        for name in ("opts.k2d", "taxo.k2d"):
            shutil.copyfile(Path(source["db_dir"]) / name, stage / name)
        taxonomy = stage / "taxonomy"
        taxonomy.mkdir()
        for name in source["taxonomy_files"]:
            shutil.copyfile(Path(source["taxonomy_dir"]) / name, taxonomy / name)
        started = time.monotonic()
        command = [str(executable), "hashshard", "--db", str(stage), "--hash-capacity", f"{capacity_per_shard // 1024}K"]
        print(f"Convert the same uncapped Kraken index; shard capacity={capacity_per_shard} slots", flush=True)
        with (stage / "conversion.log").open("wb") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Kun-peng hashshard failed ({result.returncode}); see {stage / 'conversion.log'}")
        conversion_seconds = time.monotonic() - started
        verify_start = time.monotonic()
        files, parameters = _verify_shards(stage, source, capacity_per_shard)
        for name in ("opts.k2d", "taxo.k2d"):
            if sha256(stage / name) != source["files"][name]["sha256"]:
                raise RuntimeError(f"Converted {name} differs from source")
            files[name] = _record(stage / name, source["files"][name]["sha256"])
        taxonomy_files = {}
        for name, source_record in source["taxonomy_files"].items():
            if sha256(taxonomy / name) != source_record["sha256"]:
                raise RuntimeError(f"Frozen taxonomy copy differs: {name}")
            taxonomy_files[name] = _record(taxonomy / name, source_record["sha256"])
        if hpc_reference.database_snapshot(cfg)["identity"] != source["identity"]:
            raise RuntimeError("Kraken database changed during conversion")
        (stage / "hash.k2d").unlink()
        for name in ("opts.k2d", "taxo.k2d"):
            with (stage / name).open("rb") as handle:
                os.fsync(handle.fileno())
        for name in taxonomy_files:
            with (taxonomy / name).open("rb") as handle:
                os.fsync(handle.fileno())
        _sync_dir(taxonomy)
        parameters.update(hpc_reference._native_opts(stage / "opts.k2d"))
        context = {"request": request, "source_database_identity": source["identity"],
                   "source_database_dir": source["db_dir"], "files": files,
                   "taxonomy_files": taxonomy_files, "parameters": parameters,
                   "conversion": "official_hashshard_exact_hash_cell_preservation"}
        receipt = {"schema": SCHEMA, "status": "complete", "identity": fingerprint(context),
                   "context": context, "db_dir": str(folder), "taxonomy_dir": str(folder / "taxonomy"),
                   "files": files, "taxonomy_files": taxonomy_files, "parameters": parameters,
                   "source_database_identity": source["identity"],
                   "aggregate_index_bytes": sum(item["state"]["size"] for item in files.values()),
                   "conversion_seconds": conversion_seconds,
                   "validation_seconds": time.monotonic() - verify_start,
                   "completed_at": utcnow(), "command": command}
        atomic_json(stage / "database.json", receipt)
        _sync_dir(stage)
        stage.replace(folder)
        _sync_dir(folder.parent)
        return database_snapshot(cfg)
