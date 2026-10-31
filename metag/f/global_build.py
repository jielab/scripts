"""Build one Kraken2 index from the validated, active local reference inventory.

The original references and MetaCache databases are read-only.  Each bounded
FASTA chunk is added in an isolated native staging database, then published as
one complete library directory.  An interrupted add can never append twice.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time

from common import DIR_LIB, DIR_REF, DIR_RAW, ROOT, THREADS
from fasta_stream import chunks as fasta_chunks
from library import active_references, digest, file_state
from stream_runtime import durable_json, locked, sha256, space_guard, sync_directory
from global_resources import memory_status, require_memory

SCHEMA = "metag-kraken2-global-v1"
FORMAT_VERSION = "inventory-fasta-k35-v1"
GIB = 1 << 30
DEFAULT_CHUNK_BASES = 512 << 20
KMER_LEN, MINIMIZER_LEN, MINIMIZER_SPACES = 35, 31, 7
OVERLAP = KMER_LEN - 1
TAXONOMY_FILES = ("nodes.dmp", "names.dmp", "merged.dmp")
NATIVE_FILES = ("hash.k2d", "taxo.k2d", "opts.k2d")


def _positive(value, name, allow_zero=False):
    value = float(value)
    if not math.isfinite(value) or value < 0 or (not allow_zero and value == 0):
        raise ValueError(f"{name} must be {'nonnegative' if allow_zero else 'positive'} and finite")
    return value


def _db_path(db=None):
    path = Path(db if db is not None else DIR_LIB / "kraken2_global/db").resolve()
    protected = (DIR_REF, DIR_RAW, ROOT, DIR_LIB / "reference", DIR_LIB / "metacache", DIR_LIB / "sourmash")
    for item in protected:
        item = Path(item).resolve()
        if path == item or path in item.parents or item in path.parents:
            raise RuntimeError(f"Global database must be separate from protected input/code/library: {path}, {item}")
    if path == DIR_LIB.resolve() or path in DIR_LIB.resolve().parents:
        raise RuntimeError("Use a dedicated global database directory, not the library root")
    return path


def _verified_inputs():
    inventory = DIR_LIB / "reference/inventory.tsv"
    check_path = inventory.parent / "check-ref.receipt.json"
    try:
        checked = json.loads(check_path.read_text())
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Validated reference receipt missing/invalid: {check_path}") from exc
    if checked.get("status") != "complete":
        raise RuntimeError("Reference validation is not complete; do not build from an unchecked inventory")
    st = inventory.stat()
    expected_inventory = checked.get("outputs", {}).get(str(inventory))
    observed_inventory = [st.st_mode, st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_dev, st.st_ino]
    if expected_inventory != observed_inventory:
        raise RuntimeError("Active inventory differs from its validated receipt; restore/revalidate the inventory first")
    rows = active_references()
    states = {}
    needed_taxids = set()
    total_bytes = 0
    print(f"[reference check] Validating {len(rows):,} reference file states", file=sys.stderr, flush=True)
    for index, row in enumerate(rows, 1):
        taxid = str(row.get("taxid", ""))
        if not taxid.isdigit() or int(taxid) < 1:
            raise RuntimeError(f"Active reference has no positive validated taxid: {row['filename']}")
        path = Path(row["path"])
        info = path.stat()
        state = [info.st_size, info.st_mtime_ns, info.st_ctime_ns, info.st_dev, info.st_ino]
        if not stat.S_ISREG(info.st_mode) or not info.st_size or checked.get("references", {}).get(row["relative_path"]) != state:
            raise RuntimeError(f"Reference differs from its validation receipt: {path}")
        row["taxid"] = int(taxid)
        row["state"] = [info.st_size, info.st_mtime_ns, info.st_ctime_ns]
        states[row["relative_path"]] = state
        total_bytes += info.st_size
        needed_taxids.add(taxid)
        if index % 2000 == 0 or index == len(rows):
            print(f"[reference check] {index:,}/{len(rows):,}", file=sys.stderr, flush=True)
    taxonomy = DIR_LIB / "metacache/taxonomy"
    tax_states = {}
    for name in TAXONOMY_FILES:
        path = taxonomy / name
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f"Existing local taxonomy missing/empty: {path}")
        tax_states[name] = file_state(path)
    with (taxonomy / "nodes.dmp").open() as handle:
        for line in handle:
            needed_taxids.discard(line.split("|", 1)[0].strip())
    if needed_taxids:
        raise RuntimeError(f"Active taxids absent from frozen MetaCache taxonomy: {sorted(needed_taxids, key=int)[:10]}")
    source = {"inventory": str(inventory), "inventory_sha256": sha256(inventory),
              "active_references": len(rows), "source_bytes": total_bytes,
              "reference_states_sha256": digest(states), "taxonomy_path": str(taxonomy),
              "taxonomy_states": tax_states, "excluded_references_are_not_added": True}
    return rows, source


def _disk_available(path):
    path = Path(path)
    while not path.exists():
        path = path.parent
    return shutil.disk_usage(path).free


def plan_build(db_gib=64, min_free_gib=20, memory_reserve_gib=16, db=None,
               chunk_bases=DEFAULT_CHUNK_BASES):
    db_gib = _positive(db_gib, "db_gib")
    min_free_gib = _positive(min_free_gib, "min_free_gib", True)
    memory_reserve_gib = _positive(memory_reserve_gib, "memory_reserve_gib", True)
    chunk_bases = int(chunk_bases)
    if chunk_bases <= KMER_LEN:
        raise ValueError("chunk_bases must exceed k=35")
    path = _db_path(db)
    rows, source = _verified_inputs()
    memory = memory_status()
    hash_bytes = int(db_gib * GIB)
    required_memory = hash_bytes + int(memory_reserve_gib * GIB)
    free = _disk_available(path)
    initial_disk = hash_bytes + int(min_free_gib * GIB) + 4 * chunk_bases + 2 * GIB
    return {"schema": SCHEMA, "path": str(path), "reference": source,
            "parameters": {"hash_cap_bytes": hash_bytes, "db_gib": db_gib,
                           "kmer_len": KMER_LEN, "minimizer_len": MINIMIZER_LEN,
                           "minimizer_spaces": MINIMIZER_SPACES, "masking": True,
                           "chunk_bases": chunk_bases, "overlap_bases": OVERLAP},
            "memory": memory, "required_available_memory_bytes": required_memory,
            "memory_fits": memory["effective_available_bytes"] >= required_memory,
            "free_disk_bytes": free, "initial_disk_headroom_bytes": initial_disk,
            "initial_disk_fits": free >= initial_disk,
            "disk_note": "The permanent uncompressed native library is retained in addition to hash/taxonomy/opts. "
                         "Its final size is unknown until streaming completes; free space is checked continuously.",
            "budget_note": "db_gib caps the Kraken2 hash table, not peak RAM. "
                           "If minimizers exceed the cap, native Kraken2 downsamples them; all active references are supplied."}


def _environment():
    # Native wrappers otherwise inherit settings that can silently change the DB.
    env = {key: value for key, value in os.environ.items() if not key.startswith("KRAKEN2_")}
    env["KRAKEN2_MASK_LC"] = "1"
    return env


def _software(no_masking):
    found = {}
    for name in ("kraken2-build", "kraken2-inspect"):
        program = shutil.which(name)
        if not program:
            raise RuntimeError(f"{name} is required; install Kraken2 before building the global database")
        found[name] = str(Path(program).resolve())
    native = Path(found["kraken2-build"]).parent
    if Path(found["kraken2-inspect"]).parent != native:
        raise RuntimeError("kraken2-build and kraken2-inspect must come from the same installation")
    # Kraken2 2.17.2 masks nucleotide input with its bundled k2mask.  The
    # upstream helper still has an outdated dustmasker comment at its top.
    masker = native / "k2mask"
    if not no_masking and not (masker.is_file() and os.access(masker, os.X_OK)):
        fallback = shutil.which("k2mask")
        if not fallback:
            raise RuntimeError("Default low-complexity masking requires Kraken2 2.17.2 k2mask. "
                               "Use the complete installation; --no-masking is an explicit analytical choice.")
        masker = Path(fallback).resolve()
    result = subprocess.run([found["kraken2-build"], "--version"], check=True, text=True,
                            capture_output=True, env=_environment())
    found["version"] = result.stdout.strip()
    # The official 2.17.2 source tag still reports VERSION=2.17.1.  The CLI
    # verifies the pinned source commit and installation hashes separately;
    # record the actual reported version here without rewriting upstream.
    found["wrapper_sha256"] = sha256(found["kraken2-build"])
    relevant = ("kraken2-build", "kraken2-inspect", "build_db", "estimate_capacity", "dump_table",
                "add_to_library.sh", "build_kraken2_db.sh", "scan_fasta_file.pl", "cp_into_tempfile.pl")
    hashes = {}
    for name in relevant:
        filename = native / name
        if not filename.is_file():
            raise RuntimeError(f"Incomplete native Kraken2 installation: {filename}")
        hashes[name] = sha256(filename)
    if not no_masking:
        hashes["k2mask"] = sha256(masker)
        hashes["mask_low_complexity.sh"] = sha256(native / "mask_low_complexity.sh")
    found["native_sha256"] = hashes
    return found


def _run(command, log, db, min_free_bytes):
    """Supervise native subprocesses; a disk failure never becomes a receipt."""
    space_guard(db, min_free_bytes)
    started = time.monotonic()
    with Path(log).open("w") as handle:
        handle.write("COMMAND " + json.dumps(list(map(str, command))) + "\n")
        handle.flush()
        process = subprocess.Popen(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT,
                                   env=_environment(), start_new_session=True)
        try:
            while process.poll() is None:
                space_guard(db, min_free_bytes)
                time.sleep(0.5)
            if process.returncode:
                raise RuntimeError(f"Native Kraken2 command failed (exit {process.returncode}); see {log}")
        finally:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
    return time.monotonic() - started


def _piece_header(ref_index, sequence, piece, taxid):
    return f">ref{ref_index:08d}_seq{sequence:09d}_piece{piece:06d}|kraken:taxid|{taxid}\n".encode()


def _pieces(rows, work, cap, reserve):
    """Only a single record piece is buffered; large pieces spill to disk."""
    for ref_index, row in enumerate(rows, 1):
        path = Path(row["path"])
        if file_state(path) != row["state"]:
            raise RuntimeError(f"Reference changed before formatting: {path}")
        current = None
        sequence = piece = bases = new_bases = 0
        had_bases = False
        tail = b""
        try:
            for block in fasta_chunks(path):
                if block is None:
                    if sequence and not had_bases:
                        raise RuntimeError(f"Empty FASTA sequence in {path}")
                    if current is not None:
                        current.seek(0)
                        yield current, bases, new_bases, int(piece == 1)
                        current.close(); current = None
                    sequence += 1; piece = 0; bases = new_bases = 0
                    tail = b""; had_bases = False
                    continue
                had_bases = True
                offset = 0
                while offset < len(block):
                    if current is None:
                        piece += 1
                        current = tempfile.SpooledTemporaryFile(max_size=min(1 << 20, cap), dir=work)
                        current.write(_piece_header(ref_index, sequence, piece, row["taxid"]))
                        bases = new_bases = 0
                        if piece > 1:
                            current.write(tail + b"\n"); bases = len(tail)
                    take = min(cap - bases, len(block) - offset)
                    data = block[offset:offset + take]
                    space_guard(work, reserve, len(data) + (1 << 20))
                    current.write(data + b"\n")
                    bases += take; new_bases += take; offset += take
                    tail = (tail + data)[-OVERLAP:]
                    if bases == cap:
                        current.seek(0)
                        yield current, bases, new_bases, int(piece == 1)
                        current.close(); current = None
            if sequence and not had_bases:
                raise RuntimeError(f"Empty FASTA sequence in {path}")
            if current is not None:
                current.seek(0)
                yield current, bases, new_bases, int(piece == 1)
            if not sequence or file_state(path) != row["state"]:
                raise RuntimeError(f"Empty or changed reference while formatting: {path}")
        finally:
            if current is not None:
                current.close()


def _formatted_chunks(rows, work, cap, reserve):
    """Pack whole sequences where possible; split only sequences larger than cap."""
    path = work / "current-chunk.fna"
    path.unlink(missing_ok=True)
    output = None
    stats = {"bases": 0, "original_bases": 0, "sequence_records": 0, "pieces": 0}
    number = 0

    def publish():
        nonlocal number, output
        output.flush(); os.fsync(output.fileno()); output.close(); output = None
        number += 1
        return {"number": number, "path": path, "sha256": sha256(path), "bytes": path.stat().st_size, **stats}

    try:
        for piece, bases, original_bases, sequences in _pieces(rows, work, cap, reserve):
            if output is not None and stats["bases"] + bases > cap:
                yield publish()
                path.unlink(missing_ok=True)
                stats = {"bases": 0, "original_bases": 0, "sequence_records": 0, "pieces": 0}
            if output is None:
                output = path.open("wb")
            space_guard(work, reserve, bases + (1 << 20))
            shutil.copyfileobj(piece, output, 4 << 20)
            stats["bases"] += bases
            stats["original_bases"] += original_bases
            stats["sequence_records"] += sequences
            stats["pieces"] += 1
        if output is not None:
            yield publish()
            path.unlink(missing_ok=True)
    finally:
        if output is not None:
            output.close()


def _chunk_receipt(folder, identity, expected=None):
    try:
        saved = json.loads((folder / "receipt.json").read_text())
        if saved.get("identity") != identity or saved.get("status") != "complete":
            return None
        if expected is not None and saved.get("formatted_sha256") != expected:
            return None
        if set(saved.get("outputs", {})) != {"sequences.fna", "prelim_map.txt"}:
            return None
        if not all(file_state(folder / name) == state for name, state in saved["outputs"].items()):
            return None
        return saved
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _add_chunk(chunk, db, identity, programs, threads, no_masking, reserve):
    folder = db / "library" / f"chunk{chunk['number']:08d}"
    if folder.exists():
        saved = _chunk_receipt(folder, identity, chunk["sha256"])
        if saved is None:
            raise RuntimeError(f"Existing chunk is incomplete/changed; refusing duplicate add: {folder}")
        print(f"REUSE reference chunk {chunk['number']}: {saved['original_bases']:,} original bases", flush=True)
        return saved
    staging = db / ".global-work" / "adding"
    if staging.exists():
        shutil.rmtree(staging)  # Our isolated, unpublished native add only.
    staging.mkdir()
    space_guard(db, reserve, 3 * chunk["bytes"] + (1 << 20))
    command = [programs["kraken2-build"], "--add-to-library", str(chunk["path"]),
               "--db", str(staging), "--threads", str(threads)]
    if no_masking:
        command.append("--no-masking")
    _run(command, staging / "add.log", db, reserve)
    added = staging / "library/added"
    sequences = list(added.glob("*.fna"))
    maps = list(added.glob("prelim_map_*.txt"))
    if len(sequences) != 1 or len(maps) != 1 or not sequences[0].stat().st_size:
        raise RuntimeError(f"Native add did not produce exactly one sequence/map pair: {staging}")
    mapped = 0
    with maps[0].open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3 or fields[0] != "TAXID" or not fields[2].isdigit() or int(fields[2]) < 1:
                raise RuntimeError(f"Native add produced a non-explicit or invalid taxid mapping: {staging}")
            mapped += 1
    if mapped != chunk["pieces"]:
        raise RuntimeError(f"Native add mapping count mismatch: {mapped} != {chunk['pieces']}")
    ready = staging / "ready"
    ready.mkdir()
    sequences[0].rename(ready / "sequences.fna")
    maps[0].rename(ready / "prelim_map.txt")
    saved = {"schema": SCHEMA, "status": "complete", "identity": identity,
             "formatted_sha256": chunk["sha256"], **{key: chunk[key] for key in
             ("number", "bases", "original_bases", "sequence_records", "pieces")},
             "outputs": {name: file_state(ready / name) for name in ("sequences.fna", "prelim_map.txt")}}
    durable_json(ready / "receipt.json", saved)
    ready.rename(folder); sync_directory(folder.parent)
    shutil.rmtree(staging)
    print(f"ADDED reference chunk {chunk['number']}: {chunk['original_bases']:,} original bases", flush=True)
    return saved


def _freeze_taxonomy(db, source, identity):
    folder = db / "taxonomy"
    receipt = db / ".global-work/taxonomy.json"
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved.get("identity") != identity or not all(file_state(folder / key) == value for key, value in saved["outputs"].items()):
            raise RuntimeError("Frozen global taxonomy changed; use a separate database directory")
        return
    folder.mkdir(exist_ok=True)
    for name, state in source["taxonomy_states"].items():
        original = Path(source["taxonomy_path"]) / name
        if file_state(original) != state:
            raise RuntimeError("Source taxonomy changed during global preparation")
        temp = folder / (name + ".partial")
        shutil.copyfile(original, temp)
        temp.replace(folder / name)
    durable_json(receipt, {"identity": identity, "outputs": {name: file_state(folder / name) for name in TAXONOMY_FILES}})


def database_snapshot(db=None):
    path = _db_path(db)
    try:
        saved = json.loads((path / "build.json").read_text())
        required = {str(path / name) for name in NATIVE_FILES} | {str(path / "taxonomy" / name) for name in TAXONOMY_FILES}
        if saved.get("schema") != SCHEMA or saved.get("status") != "complete" or not saved.get("identity"):
            raise ValueError("not complete")
        if not required.issubset(saved.get("outputs", {})) or not saved.get("native_inspect_passed"):
            raise ValueError("missing verified native files")
        if not all(file_state(filename) == state for filename, state in saved["outputs"].items()):
            raise ValueError("output state changed")
        if saved.get("path") != str(path) or saved.get("taxonomy_path") != str(path / "taxonomy"):
            raise ValueError("database moved")
        return saved
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise RuntimeError(f"No unchanged, completed global Kraken2 database at {path}; build it first") from exc


def build_database(args):
    path = _db_path(getattr(args, "db", None))
    db_gib = _positive(getattr(args, "db_gib", 64), "db_gib")
    min_free_gib = _positive(getattr(args, "min_free_gib", 20), "min_free_gib", True)
    memory_reserve = _positive(getattr(args, "memory_reserve_gib", 16), "memory_reserve_gib", True)
    threads = int(getattr(args, "threads", THREADS))
    cap = int(getattr(args, "chunk_bases", DEFAULT_CHUNK_BASES))
    no_masking = bool(getattr(args, "no_masking", False))
    if threads < 1 or cap <= KMER_LEN:
        raise ValueError("threads must be positive and chunk_bases must exceed k=35")
    rows, source = _verified_inputs()
    programs = _software(no_masking)
    hash_cap = int(db_gib * GIB)
    if hash_cap < 1024:
        raise ValueError("Hash cap is too small; db_gib must provide at least 1024 bytes")
    parameters = {"hash_cap_bytes": hash_cap, "db_gib": db_gib, "masking": not no_masking,
                  "kmer_len": KMER_LEN, "minimizer_len": MINIMIZER_LEN,
                  "minimizer_spaces": MINIMIZER_SPACES, "chunk_bases": cap,
                  "overlap_bases": OVERLAP, "format_version": FORMAT_VERSION,
                  "native_version": programs["version"], "builder_sha256": programs["wrapper_sha256"],
                  "native_sha256": programs["native_sha256"],
                  "native_installation": getattr(args, "kraken_installation", {})}
    context = {"schema": SCHEMA, "reference": source, "parameters": parameters}
    identity = digest(context)
    if (path / "build.json").exists():
        saved = database_snapshot(path)
        if saved["identity"] != identity:
            raise RuntimeError("Completed global database has different references/settings; choose a new --db")
        print(f"REUSE completed global Kraken2 database: {path}", flush=True)
        return saved
    marker = path / "global-state.json"
    if path.exists() and any(path.iterdir()) and not marker.exists():
        raise RuntimeError(f"Refusing to adopt a nonempty, unowned database directory: {path}")
    require_memory(hash_cap, reserve_gib=memory_reserve, context="Global Kraken2 build (hash cap is not a peak-RAM cap)")
    reserve = int(min_free_gib * GIB)
    space_guard(path, reserve, hash_cap + 4 * cap + 2 * GIB)
    path.mkdir(parents=True, exist_ok=True)
    with locked(path / "global-build.lock"):
        if marker.exists():
            saved_state = json.loads(marker.read_text())
            if saved_state.get("identity") != identity:
                raise RuntimeError("Incomplete global build belongs to different inputs/settings; choose a new --db")
        else:
            if any(p.name != "global-build.lock" for p in path.iterdir()):
                raise RuntimeError("Database directory changed before ownership was established")
            durable_json(marker, {"schema": SCHEMA, "status": "preparing", "identity": identity, "context": context})
        work = path / ".global-work"
        work.mkdir(exist_ok=True)
        (path / "library").mkdir(exist_ok=True)
        _freeze_taxonomy(path, source, identity)
        ready_path = work / "library-ready.json"
        chunks = []
        if ready_path.exists():
            ready = json.loads(ready_path.read_text())
            if ready.get("identity") != identity:
                raise RuntimeError("Prepared library belongs to a different build")
            for number in range(1, ready["chunks"] + 1):
                item = _chunk_receipt(path / "library" / f"chunk{number:08d}", identity)
                if item is None:
                    raise RuntimeError("Prepared reference chunk changed or is incomplete")
                chunks.append(item)
        else:
            for chunk in _formatted_chunks(rows, work, cap, reserve + hash_cap):
                chunks.append(_add_chunk(chunk, path, identity, programs, threads, no_masking, reserve + hash_cap))
            if not chunks:
                raise RuntimeError("Active references contained no sequences")
            ready = {"identity": identity, "chunks": len(chunks)}
            durable_json(ready_path, ready)
        expected_chunks = {f"chunk{i:08d}" for i in range(1, len(chunks) + 1)}
        if {p.name for p in (path / "library").iterdir()} != expected_chunks:
            raise RuntimeError("Unexpected files/chunks in owned native library; refusing extra reference inclusion")
        # A hard kill/OOM cannot run exception cleanup. Even all three final
        # files can exist while the hash table is only partially populated.
        # Without build.json they are uncommitted: never let the native wrapper
        # skip that hash on retry. Keep the expensive committed library chunks.
        for name in NATIVE_FILES:
            (path / name).unlink(missing_ok=True)
            (path / (name + ".tmp")).unlink(missing_ok=True)
        require_memory(hash_cap, reserve_gib=memory_reserve, context="Global Kraken2 native build")
        space_guard(path, reserve, hash_cap + 2 * GIB)
        command = [programs["kraken2-build"], "--build", "--db", str(path), "--threads", str(threads),
                   "--max-db-size", str(hash_cap), "--kmer-len", str(KMER_LEN),
                   "--minimizer-len", str(MINIMIZER_LEN), "--minimizer-spaces", str(MINIMIZER_SPACES)]
        durable_json(marker, {"schema": SCHEMA, "status": "building", "identity": identity, "context": context})
        inspect_command = [programs["kraken2-inspect"], "--db", str(path), "--threads", str(threads)]
        try:
            seconds = _run(command, path / "global-build.log", path, reserve)
            if any(not (path / name).is_file() or not (path / name).stat().st_size for name in NATIVE_FILES):
                raise RuntimeError("Native build did not produce nonempty hash.k2d, taxo.k2d and opts.k2d")
            _run(inspect_command, path / "global-inspect.txt", path, reserve)
        except BaseException:
            # The native wrapper skips an existing hash.k2d.  Keep the costly
            # committed reference chunks, but never let failed final files
            # make every later retry skip the rebuild.
            for name in NATIVE_FILES:
                (path / name).unlink(missing_ok=True)
                (path / (name + ".tmp")).unlink(missing_ok=True)
            sync_directory(path)
            raise
        _, after = _verified_inputs()
        if after != source:
            raise RuntimeError("References/inventory/taxonomy changed during build; completion was not published")
        outputs = {str(path / name): file_state(path / name) for name in NATIVE_FILES}
        outputs.update({str(path / "taxonomy" / name): file_state(path / "taxonomy" / name) for name in TAXONOMY_FILES})
        saved = {"schema": SCHEMA, "status": "complete", "identity": identity,
                 "path": str(path), "taxonomy_path": str(path / "taxonomy"), "outputs": outputs,
                 "parameters": parameters, "reference": source, "native_inspect_passed": True,
                 "native_build_seconds": seconds, "native_build_command": command,
                 "native_inspect_report": str(path / "global-inspect.txt"),
                 "library": {"chunks": len(chunks), "original_bases": sum(c["original_bases"] for c in chunks),
                     "sequence_records": sum(c["sequence_records"] for c in chunks),
                     "formatted_bases": sum(c["bases"] for c in chunks),
                     "overlap_bases": sum(c["bases"] - c["original_bases"] for c in chunks),
                     "formatted_sequence_pieces": sum(c["pieces"] for c in chunks),
                     "retained_uncompressed_library_bytes": sum(c["outputs"]["sequences.fna"][0] for c in chunks)},
                 "sampling": {"hash_cap_bytes": hash_cap,
                              "policy": "Native minimizer downsampling if estimated hash capacity exceeds cap; references are not removed"},
                 "completed_at": datetime.now(timezone.utc).isoformat()}
        saved["resource_limits"] = getattr(args, "resource_limits", {})
        durable_json(path / "build.json", saved)
        durable_json(marker, {"schema": SCHEMA, "status": "complete", "identity": identity, "context": context})
        print(f"PASS global Kraken2 database: {path}; active references={len(rows):,}", flush=True)
        return saved
