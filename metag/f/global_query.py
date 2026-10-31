#!/usr/bin/env python3
"""One global Kraken 2 query per sample, with bounded ordinal bookkeeping.

The producer validates the original FASTQs and streams temporary numeric IDs to
the native classifier.  No full query FASTQ or native classification table is
materialized.  Only validated, complete ordinal taxids can be published; the
original reads are reread by global_split to retain their exact FASTQ bytes.
"""
from __future__ import annotations

from contextlib import ExitStack
import hashlib
import json
import math
import mmap
import os
from pathlib import Path
import re
import signal
import struct
import subprocess
import threading
import time

from common import ROOT, DIR_OUT
from fastq_stream import (fastq_records, input_state, validate_pair_ids,
                         new_header_digest, update_header_digest)
from library import digest, file_state
from sample_state import is_complete
from split_metacache import Taxonomy
from stream_runtime import durable_json, locked, safe_destination, sha256, space_guard, _stop

VERSION = "metag-global-query-v1"
CLASSIFIER = "Kraken 2 global minimizer classifier; not equivalent to MetaCache"
_U64 = struct.Struct("<Q")
_LENGTHS = struct.Struct("<II")
_UINT = re.compile(rb"[0-9]+\Z")


def resolve_native(executable):
    """Use the installed native program: the legacy wrapper lacks interleaving."""
    executable = Path(executable).resolve(strict=True)
    native = executable if executable.name == "classify" else executable.parent / "classify"
    if not native.is_file() or not os.access(native, os.X_OK):
        raise RuntimeError(f"Native Kraken 2 classify program missing beside {executable}. "
                           "Use the complete Kraken 2 installation directory on PATH.")
    probe = subprocess.run([str(native), "-h"], text=True, capture_output=True, timeout=15)
    help_text = probe.stdout + probe.stderr
    required = ("H", "t", "o", "p", "T", "g", "Q", "S")
    if probe.returncode or any(not re.search(r"(?m)^\s*\*?\s*-" + flag + r"(?:\s|$)", help_text)
                              for flag in required):
        raise RuntimeError(f"Unsupported native Kraken 2 classifier at {native}: "
                           "required stdin/interleaved and classification flags are absent")
    return native.resolve()


def classifier_command(native, db, args, paired):
    threads = int(getattr(args, "threads", 16))
    confidence = float(getattr(args, "confidence", 0.1))
    groups = int(getattr(args, "minimum_hit_groups", 2))
    quality = int(getattr(args, "minimum_base_quality", 0))
    if threads < 1 or groups < 1 or quality < 0 or not math.isfinite(confidence) or not 0 <= confidence <= 1:
        raise ValueError("Invalid Kraken 2 threads/confidence/minimum-hit-groups/base-quality")
    db = Path(db)
    command = [str(native), "-H", str(db/"hash.k2d"), "-t", str(db/"taxo.k2d"),
               "-o", str(db/"opts.k2d"), "-p", str(threads), "-T", str(confidence),
               "-g", str(groups), "-Q", str(quality)]
    if paired:
        command.append("-S")  # Native interleaved pairs. -P instead expects two files.
    if getattr(args, "quick", False):
        command.append("-q")
    # Native classify reads stdin when no filenames are given. Do not pass -O -:
    # that option suppresses the per-fragment table required for validation.
    return command


class OrdinalStore:
    """Grow disk-backed taxids, expected lengths and a bitmap, never a read dict."""
    def __init__(self, folder, reserve=0, chunk=1 << 20):
        self.folder, self.reserve, self.chunk = Path(folder), reserve, chunk
        self.paths = [self.folder/name for name in
                      ("taxids.partial.bin", "lengths.partial.bin", "seen.partial.bin")]
        self.files, self.maps = [], []
        self.capacity = self.produced = self.received = 0
        self.mutex = threading.Lock()
        try:
            for path in self.paths:
                self.files.append(path.open("w+b"))
            self._grow(1)
        except BaseException:
            self.close()
            raise

    def _grow(self, count):
        if count <= self.capacity:
            return
        capacity = ((count + self.chunk - 1)//self.chunk)*self.chunk
        extra = (capacity-self.capacity)*16 + (capacity+7)//8 - (self.capacity+7)//8
        space_guard(self.folder, self.reserve, extra)
        sizes = (capacity*8, capacity*8, (capacity+7)//8)
        for i, size in enumerate(sizes):
            self.files[i].truncate(size)
            self.files[i].flush()
            if i < len(self.maps):
                self.maps[i].resize(size)
            else:
                self.maps.append(mmap.mmap(self.files[i].fileno(), size))
        self.capacity = capacity

    def add_input(self, length1, length2=0):
        if not 0 < length1 < 2**32 or not 0 <= length2 < 2**32:
            raise RuntimeError("FASTQ read length is outside the supported uint32 range")
        with self.mutex:
            ordinal = self.produced + 1
            self._grow(ordinal)
            _LENGTHS.pack_into(self.maps[1], (ordinal-1)*8, length1, length2)
            self.produced = ordinal  # Publish metadata before any bytes reach native stdin.
            return ordinal

    def add_output(self, ordinal, taxid, lengths):
        with self.mutex:
            if ordinal < 1 or ordinal > self.produced:
                raise RuntimeError(f"Native query ID outside produced input range: {ordinal}/{self.produced}")
            index = ordinal-1
            byte, bit = index//8, 1 << (index%8)
            if self.maps[2][byte] & bit:
                raise RuntimeError(f"Duplicate native query ID: {ordinal}")
            expected = _LENGTHS.unpack_from(self.maps[1], index*8)
            if lengths != expected:
                raise RuntimeError(f"Native read lengths differ at fragment {ordinal}: {lengths} != {expected}")
            if not 0 <= taxid < 2**64:
                raise RuntimeError("Native taxid is outside uint64 range")
            _U64.pack_into(self.maps[0], index*8, taxid)
            self.maps[2][byte] |= bit
            self.received += 1

    def publish(self, target):
        if not self.produced or self.received != self.produced:
            raise RuntimeError(f"Native query IDs incomplete: {self.received}/{self.produced}")
        # Uniqueness plus range plus equal count proves every ordinal occurred.
        self.maps[0].flush()
        os.fsync(self.files[0].fileno())
        self.close()
        with self.paths[0].open("r+b") as handle:
            handle.truncate(self.produced*8)
            handle.flush()
            os.fsync(handle.fileno())
        self.paths[0].replace(target)
        self.paths[1].unlink()
        self.paths[2].unlink()

    def close(self):
        for mapping in self.maps:
            if not mapping.closed:
                mapping.close()
        for handle in self.files:
            handle.close()


def _number(token, label):
    if not _UINT.fullmatch(token):
        raise RuntimeError(f"Invalid native {label}: {token[:120]!r}")
    return int(token)


def parse_native_line(line, paired, tax, normalized):
    columns = line.rstrip(b"\r\n").split(b"\t")
    if len(columns) != 5 or columns[0] not in (b"C", b"U"):
        raise RuntimeError(f"Malformed Kraken 2 classification row: {line[:160]!r}")
    ordinal = _number(columns[1], "query ID")
    raw_taxid = _number(columns[2], "taxid")
    if (columns[0] == b"U") != (raw_taxid == 0):
        raise RuntimeError("Native C/U label contradicts its taxid")
    lengths = columns[3].split(b"|")
    if len(lengths) != (2 if paired else 1):
        raise RuntimeError("Native paired/single read-length format differs from input")
    first = _number(lengths[0], "read length")
    second = _number(lengths[1], "mate length") if paired else 0
    if not first or (paired and not second):
        raise RuntimeError("Native classification has an empty read")
    if raw_taxid not in normalized:
        if raw_taxid and raw_taxid not in tax.parent:
            raise RuntimeError(f"Native output taxid absent from database taxonomy: {raw_taxid}")
        # A strain/subspecies call is promoted to its species. Calls above
        # species remain above species; no arbitrary species is invented.
        normalized[raw_taxid] = tax.species(raw_taxid) or raw_taxid
    return ordinal, normalized[raw_taxid], (first, second)


def query_once(sample, snapshot, executable, args, work, *, limit=0):
    """Stream and validate one complete native run; publish only its taxid array."""
    work = Path(work)
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 0:
        raise ValueError("Query fragment limit must be a nonnegative integer")
    work.mkdir(parents=True, exist_ok=True)
    reserve = int(getattr(args, "min_free_gib", 20)*2**30)
    native = resolve_native(executable)
    command = classifier_command(native, snapshot["path"], args, sample["r2"] is not None)
    sources = [input_state(p) for p in (sample["r1"], sample["r2"]) if p is not None]
    taxonomy_states = {str(Path(snapshot["taxonomy_path"])/name):
                       file_state(Path(snapshot["taxonomy_path"])/name)
                       for name in ("nodes.dmp", "names.dmp")}
    taxonomy = Taxonomy(Path(snapshot["taxonomy_path"]))
    store = OrdinalStore(work, reserve)
    result, errors, stopped = {}, [], threading.Event()
    log_path = work/"query.log"
    process = thread = None
    started = time.monotonic()

    def produce():
        digest_headers = new_header_digest()
        total_bases = total_fastq_bytes = 0
        try:
            with ExitStack() as stack:
                left = fastq_records(sample["r1"])
                right = fastq_records(sample["r2"]) if sample["r2"] is not None else None
                stack.callback(left.close)
                if right is not None:
                    stack.callback(right.close)
                count = 0
                while not limit or count < limit:
                    if stopped.is_set():
                        raise RuntimeError("Input streaming was stopped before completion")
                    a = next(left, None)
                    b = next(right, None) if right is not None else None
                    if a is None:
                        if b is not None:
                            raise RuntimeError("FASTQ mates contain different numbers of records (R2 longer)")
                        break
                    if right is not None and b is None:
                        raise RuntimeError("FASTQ mates contain different numbers of records (R1 longer)")
                    if b is not None:
                        validate_pair_ids(a[0], b[0])
                    update_header_digest(digest_headers, a, b)
                    records = (a, b) if b is not None else (a,)
                    lengths = [len(record[1].rstrip(b"\r\n")) for record in records]
                    ordinal = store.add_input(lengths[0], lengths[1] if b is not None else 0)
                    name = b"@" + str(ordinal).encode("ascii") + b"\n"
                    for record in records:
                        # Only the transient query stream uses numeric headers.
                        # Original headers, plus lines and newlines are retained
                        # in the final species outputs by rereading raw inputs.
                        process.stdin.write(name + record[1].rstrip(b"\r\n") + b"\n+\n" +
                                            record[3].rstrip(b"\r\n") + b"\n")
                        total_fastq_bytes += sum(map(len, record))
                    total_bases += sum(lengths)
                    count += 1
                    if count % 100_000 == 0:
                        space_guard(work, reserve)
                    if count % 1_000_000 == 0:
                        print(f"[global query input] {count:,} fragments streamed", flush=True)
                if not count:
                    raise RuntimeError("FASTQ input contains no fragments")
                if sources != [input_state(p) for p in (sample["r1"], sample["r2"]) if p is not None]:
                    raise RuntimeError("FASTQ inputs changed during global query")
                process.stdin.close()
                result.update(fragments=count, reads=count*(2 if right is not None else 1),
                              paired=right is not None, header_digest=digest_headers.hexdigest(),
                              total_bases=total_bases, fastq_uncompressed_bytes=total_fastq_bytes,
                              validated_scope="prefix" if limit else "entire_input")
        except BaseException as exc:
            errors.append(exc)
            stopped.set()
            _stop(process)
            try:
                process.stdin.close()
            except (OSError, ValueError):
                pass

    try:
        with log_path.open("w") as log:
            log.write("COMMAND " + json.dumps(command) + "\n")
            log.flush()
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=log, start_new_session=True, bufsize=1 << 20)
            thread = threading.Thread(target=produce, name="metag-global-input", daemon=True)
            thread.start()
            normalized = {0: 0}
            while True:
                line = process.stdout.readline(16 << 20)
                if not line:
                    break
                if not line.endswith(b"\n"):
                    raise RuntimeError("Native row truncated or exceeds 16 MiB")
                ordinal, tid, lengths = parse_native_line(line, sample["r2"] is not None,
                                                          taxonomy, normalized)
                store.add_output(ordinal, tid, lengths)
            code = process.wait()
            thread.join(timeout=10)
            if thread.is_alive():
                raise RuntimeError("Input producer did not settle after native exit")
            if errors:
                raise RuntimeError(f"FASTQ producer failed; no result committed: {errors[0]}") from errors[0]
            if code:
                raise RuntimeError(f"Kraken 2 query exited with status {code}; see {log_path}")
            if store.received != result.get("fragments") or store.produced != result.get("fragments"):
                raise RuntimeError(f"Missing native query IDs: {store.received}/{result.get('fragments')}")
            if any(file_state(Path(path)) != saved for path, saved in snapshot["outputs"].items()):
                raise RuntimeError("Global Kraken 2 database changed during query")
            if any(file_state(Path(path)) != saved for path, saved in taxonomy_states.items()):
                raise RuntimeError("Global database taxonomy changed during query")
            target = work/"taxids.bin"
            store.publish(target)
            seconds = time.monotonic()-started
            result.update(seconds=seconds, fragments_per_second=result["fragments"]/max(seconds, 1e-9),
                          command=command, query_log=str(log_path), taxids=str(target),
                          output=file_state(target), input_states=sources, classifier=CLASSIFIER,
                          native_executable=str(native), native_sha256=sha256(native),
                          requested_limit=limit, normalized_below_species=True,
                          resource_limits=getattr(args, "resource_limits", {}),
                          compact_bytes=target.stat().st_size)
            return result
    finally:
        stopped.set()
        _stop(process)
        if thread is not None:
            thread.join(timeout=10)
        if process is not None:
            if process.stdout is not None:
                process.stdout.close()
            try:
                process.stdin.close()
            except (OSError, ValueError):
                pass
        if thread is None or not thread.is_alive():
            store.close()


def query_context(sample, snapshot, executable, args):
    native = resolve_native(executable)
    return {"version": VERSION, "name": sample["name"], "output": str(sample["target"]),
            "reads": [input_state(p) for p in (sample["r1"], sample["r2"]) if p is not None],
            "database": snapshot, "classifier": CLASSIFIER,
            "taxonomy": {name: file_state(Path(snapshot["taxonomy_path"])/name)
                         for name in ("nodes.dmp", "names.dmp")},
            "native": {"path": str(native), "sha256": sha256(native)},
            "installation": getattr(args, "kraken_installation", {}),
            "parameters": {"confidence": float(getattr(args, "confidence", 0.1)),
                           "quick": bool(getattr(args, "quick", False)),
                           "read_reports": bool(getattr(args, "read_reports", True)),
                           "minimum_hit_groups": int(getattr(args, "minimum_hit_groups", 2)),
                           "minimum_base_quality": int(getattr(args, "minimum_base_quality", 0))},
            "code": {name: sha256(ROOT/"f"/name) for name in
                     ("global_query.py", "global_split.py", "fastq_stream.py", "sample_state.py",
                      "split_metacache.py", "common.py")}}


def _verify_context(sample, snapshot, executable, args, context):
    if query_context(sample, snapshot, executable, args) != context:
        raise RuntimeError("Input/code/native identity changed during global analysis")
    if any(file_state(Path(path)) != state for path, state in snapshot["outputs"].items()):
        raise RuntimeError("Global database files changed during analysis")


def _memory_check(snapshot, args):
    from global_resources import require_memory
    return require_memory(sum(state[0] for state in snapshot["outputs"].values()),
                          reserve_gib=getattr(args, "memory_reserve_gib", 16),
                          context="Global Kraken 2 query")


def _clean_completed(work, identity):
    """Release only this run's compact query, including interrupted publication."""
    receipt = work/"run.json"
    if not receipt.is_file():
        return
    saved = json.loads(receipt.read_text())
    if saved.get("identity") != identity:
        raise RuntimeError("Completed global output has a conflicting local receipt")
    saved["status"] = "complete"
    saved["work_cleaned"] = True
    (work/"taxids.bin").unlink(missing_ok=True)
    durable_json(receipt, saved)


def run_sample(sample, snapshot, executable, args):
    """A failed query restarts once; a completed query resumes the final split."""
    current = query_context(sample, snapshot, executable, args)
    identity = digest(current)
    target = safe_destination(sample["target"])
    work = safe_destination(getattr(args, "work_root", DIR_OUT/"runs_global"))/sample["name"]/identity[:20]
    if is_complete(target, identity):
        if not getattr(args, "keep_work", False) and (work/"run.json").is_file():
            with locked(work/"run.lock"):
                _clean_completed(work, identity)
        print(f"SKIP completed global sample: {sample['name']} -> {target}", flush=True)
        return target
    if target.exists():
        raise RuntimeError(f"Global output exists with a different/incomplete identity: {target}; "
                           "choose a new output suffix. Existing output is retained.")
    work.mkdir(parents=True, exist_ok=True)
    reserve = int(getattr(args, "min_free_gib", 20)*2**30)
    with locked(work/"run.lock"), locked(target.parent/("."+target.name+".global.lock")):
        if is_complete(target, identity):
            return target
        if target.exists():
            raise RuntimeError(f"Output appeared while acquiring the global lock: {target}")
        space_guard(work, reserve)
        space_guard(target, reserve)
        durable_json(work/"run.json", {"status": "running", "identity": identity, "context": current})
        receipt_path = work/"query.json"
        query = None
        if receipt_path.is_file():
            saved = json.loads(receipt_path.read_text())
            path = work/"taxids.bin"
            if (saved.get("identity") == identity and saved.get("status") == "complete"
                    and path.is_file() and saved.get("query", {}).get("output") == file_state(path)):
                query = saved["query"]
        if query is None:
            _memory_check(snapshot, args)
            print(f"[global] {sample['name']}: one Kraken 2 query of the whole sample", flush=True)
            query = query_once(sample, snapshot, executable, args, work)
            _verify_context(sample, snapshot, executable, args, current)
            durable_json(receipt_path, {"status": "complete", "identity": identity, "query": query})
        else:
            print(f"Reuse complete global query: {query['fragments']:,} fragments", flush=True)
        from global_split import split
        result = split(sample["r1"], sample["r2"], work/"taxids.bin", target,
                       expected_fragments=query["fragments"], input_identity=identity,
                       taxonomy_folder=snapshot["taxonomy_path"], classifier=CLASSIFIER,
                       expected_header_digest=query["header_digest"], min_free_bytes=reserve,
                       read_reports=getattr(args, "read_reports", True),
                       before_publish=lambda: _verify_context(sample, snapshot, executable, args, current),
                       provenance={"identity": identity, "context": current, "query": query,
                                   "reference_scope": "all references in the global database",
                                   "reference_evidence": "not_recorded"})
        durable_json(work/"run.json", {"status": "complete", "identity": identity, "context": current,
                                       "output": str(result), "query": query})
        if not getattr(args, "keep_work", False):
            _clean_completed(work, identity)
        return result


def benchmark(sample, snapshot, executable, args):
    limit = int(getattr(args, "benchmark_fragments", 1_000_000))
    if limit < 1:
        raise ValueError("benchmark_fragments must be positive")
    work = safe_destination(getattr(args, "work_root", DIR_OUT/"runs_global"))/"benchmarks"/(
        sample["name"]+"-"+time.strftime("%Y%m%d-%H%M%S")+f"-{os.getpid()}")
    work.mkdir(parents=True)
    _memory_check(snapshot, args)
    query = query_once(sample, snapshot, executable, args, work, limit=limit)
    durable_json(work/"benchmark.json", {"version": VERSION, "sample": sample["name"],
                 "database": snapshot, "query": query, "species_output_created": False,
                 "includes_database_load": True,
                 "validated_scope": "prefix only; throughput diagnostic, not recall validation"})
    if not getattr(args, "keep_work", False):
        (work/"taxids.bin").unlink(missing_ok=True)
    print(f"PASS global benchmark: {query['fragments']:,} fragments in {query['seconds']:.2f}s; "
          f"{query['fragments_per_second']:,.0f} fragments/s including database load. {work/'benchmark.json'}",
          flush=True)
    return work
