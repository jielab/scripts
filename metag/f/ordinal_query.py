"""One mmap Kraken query with a disk-backed ordinal protocol, no read-name join."""
import errno
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

from fastq_stream import fragments, new_digest, update_digest
from hpc_common import check_space, file_state, sha256
from taxonomy import Taxonomy

_U64, _LENGTHS = struct.Struct("<Q"), struct.Struct("<II")
_UINT = re.compile(rb"[0-9]+\Z")


def stop_process(process):
    if process is None or process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=3)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=3)


def input_states(sample):
    result = []
    for name in ("r1", "r2"):
        if sample.get(name) is not None:
            path = Path(sample[name]).resolve(strict=True)
            if not path.is_file():
                raise RuntimeError(f"FASTQ input is not a regular file: {path}")
            result.append({"path": str(path), "state": file_state(path)})
    return result


def validate_snapshot(snapshot):
    if snapshot.get("status") != "complete" or not snapshot.get("identity"):
        raise RuntimeError("Only a complete database snapshot can be queried")
    db, tax = Path(snapshot["db_dir"]), Path(snapshot["taxonomy_dir"])
    for name in ("hash.k2d", "taxo.k2d", "opts.k2d"):
        item = snapshot.get("files", {}).get(name)
        if not item or file_state(db / name) != item["state"]:
            raise RuntimeError(f"Database file changed or is not recorded: {db / name}")
    taxonomy_files = snapshot.get("taxonomy_files", {})
    if not all(name in taxonomy_files for name in ("nodes.dmp", "names.dmp")):
        raise RuntimeError("Database snapshot must include nodes.dmp and names.dmp file states")
    for name, item in taxonomy_files.items():
        if Path(name).name != name or file_state(tax / name) != item["state"]:
            raise RuntimeError(f"Frozen database taxonomy changed: {name}")
    marker = db / "database.json"
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved.get("status") != "complete" or saved.get("identity") != snapshot["identity"]:
            raise RuntimeError("Database completion marker changed")


def resolve_native(cfg):
    native = Path(cfg["tools_dir"]) / "kraken/bin/classify"
    native = native.resolve(strict=True)
    if not native.is_file() or not os.access(native, os.X_OK):
        raise RuntimeError(f"Native Kraken classify is not executable: {native}")
    probe = subprocess.run([str(native), "-h"], capture_output=True, text=True, timeout=15)
    help_text = probe.stdout + probe.stderr
    # Upstream prints -c and -M on the same help line: do not anchor at line start.
    if probe.returncode or any(not re.search(r"(?:^|\s)-" + flag + r"(?:\s|$)", help_text)
                              for flag in ("H", "t", "o", "p", "T", "g", "Q", "S", "M")):
        raise RuntimeError("Unsupported native Kraken classifier: mmap/interleaved flags missing")
    return native


def classifier_command(native, snapshot, cfg, paired):
    threads = int(cfg.get("classify_threads", 8))
    confidence = float(cfg.get("confidence", 0.1))
    groups = int(cfg.get("minimum_hit_groups", 2))
    quality = int(cfg.get("minimum_base_quality", 0))
    if threads < 1 or groups < 1 or quality < 0 or not math.isfinite(confidence) or not 0 <= confidence <= 1:
        raise ValueError("Invalid classification threads/confidence/hit groups/base quality")
    db = Path(snapshot["db_dir"])
    command = [str(native), "-H", str(db / "hash.k2d"), "-t", str(db / "taxo.k2d"),
               "-o", str(db / "opts.k2d"), "-M", "-p", str(threads), "-T", str(confidence),
               "-g", str(groups), "-Q", str(quality)]
    if paired:
        command.append("-S")
    return command  # No filename: stdin. -P would instead require two input files.


class OrdinalStore:
    """16 bytes plus one bitmap bit per fragment; no Python object per read."""
    def __init__(self, folder, min_free_gib=0, chunk=1 << 20):
        self.folder, self.min_free_gib, self.chunk = Path(folder), min_free_gib, chunk
        self.paths = [self.folder / name for name in
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
        capacity = ((count + self.chunk - 1) // self.chunk) * self.chunk
        extra = (capacity - self.capacity) * 16 + (capacity + 7) // 8 - (self.capacity + 7) // 8
        check_space(self.folder, self.min_free_gib + extra / 2**30)
        for i, size in enumerate((capacity * 8, capacity * 8, (capacity + 7) // 8)):
            old = os.fstat(self.files[i].fileno()).st_size
            # Reserve actual blocks where supported, avoiding a delayed sparse-file
            # SIGBUS when an mmap page is first touched on a nearly full filesystem.
            try:
                os.posix_fallocate(self.files[i].fileno(), old, size - old)
            except (AttributeError, OSError) as exc:
                if isinstance(exc, OSError) and exc.errno not in (errno.EOPNOTSUPP, errno.ENOSYS, errno.EINVAL):
                    raise
                self.files[i].truncate(size)
            self.files[i].flush()
            if i < len(self.maps):
                self.maps[i].resize(size)
            else:
                self.maps.append(mmap.mmap(self.files[i].fileno(), size))
        self.capacity = capacity

    def add_input(self, length1, length2=0):
        if not 0 < length1 < 2**32 or not 0 <= length2 < 2**32:
            raise RuntimeError("Read length is outside uint32 range")
        with self.mutex:
            ordinal = self.produced + 1
            self._grow(ordinal)
            _LENGTHS.pack_into(self.maps[1], (ordinal - 1) * 8, length1, length2)
            self.produced = ordinal
            return ordinal

    def add_output(self, ordinal, taxid, lengths):
        with self.mutex:
            if ordinal < 1 or ordinal > self.produced:
                raise RuntimeError(f"Native query ID outside produced range: {ordinal}/{self.produced}")
            index = ordinal - 1
            byte, bit = index // 8, 1 << (index % 8)
            if self.maps[2][byte] & bit:
                raise RuntimeError(f"Duplicate native query ID: {ordinal}")
            if lengths != _LENGTHS.unpack_from(self.maps[1], index * 8):
                raise RuntimeError(f"Native read lengths differ at ordinal {ordinal}")
            if not 0 <= taxid < 2**64:
                raise RuntimeError("Native taxid outside uint64 range")
            _U64.pack_into(self.maps[0], index * 8, taxid)
            self.maps[2][byte] |= bit
            self.received += 1

    def publish(self, target):
        if not self.produced or self.received != self.produced:
            raise RuntimeError(f"Native query IDs incomplete: {self.received}/{self.produced}")
        self.maps[0].flush()
        os.fsync(self.files[0].fileno())
        self.close()
        with self.paths[0].open("r+b") as handle:
            handle.truncate(self.produced * 8)
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


def _number(value, label):
    if not _UINT.fullmatch(value):
        raise RuntimeError(f"Invalid native {label}: {value[:100]!r}")
    return int(value)


def parse_native_line(line, paired, taxonomy):
    fields = line.rstrip(b"\r\n").split(b"\t")
    if len(fields) != 5 or fields[0] not in (b"C", b"U"):
        raise RuntimeError(f"Malformed Kraken classification row: {line[:160]!r}")
    ordinal, taxid = _number(fields[1], "ordinal"), _number(fields[2], "taxid")
    if (fields[0] == b"U") != (taxid == 0):
        raise RuntimeError("Native C/U status contradicts taxid")
    values = fields[3].split(b"|")
    if len(values) != (2 if paired else 1):
        raise RuntimeError("Native paired/single length format differs from input")
    lengths = tuple(_number(x, "read length") for x in values)
    if not all(lengths):
        raise RuntimeError("Native classification contains an empty read")
    return ordinal, taxonomy.describe(taxid)["taxid"], lengths if paired else (lengths[0], 0)


def query_once(cfg, sample, snapshot, work, limit=0):
    work = Path(work)
    validate_snapshot(snapshot)
    reserve = float(cfg.get("min_free_gib", 100))
    check_space(work, reserve)
    native = resolve_native(cfg)
    command = classifier_command(native, snapshot, cfg, sample.get("r2") is not None)
    environment = os.environ.copy()
    # Match the already validated -p budget. Keep site library, locale, temporary
    # directory and scheduler settings; do not mutate the parent environment.
    environment["OMP_NUM_THREADS"] = str(int(cfg.get("classify_threads", 8)))
    sources = input_states(sample)
    taxonomy = Taxonomy(snapshot["taxonomy_dir"])
    store = OrdinalStore(work, reserve)
    result, errors, stopped = {}, [], threading.Event()
    process = thread = None
    log_path = work / "query.log"
    started = time.monotonic()

    def produce():
        digest = new_digest()
        total_bases = total_bytes = count = 0
        try:
            for a, b in fragments(sample["r1"], sample.get("r2"), limit):
                if stopped.is_set():
                    raise RuntimeError("FASTQ streaming interrupted")
                update_digest(digest, a, b)
                records = (a, b) if b is not None else (a,)
                lengths = [len(record[1].rstrip(b"\r\n")) for record in records]
                ordinal = store.add_input(lengths[0], lengths[1] if b is not None else 0)
                header = b"@" + str(ordinal).encode("ascii") + b"\n"
                for record in records:
                    process.stdin.write(header + record[1].rstrip(b"\r\n") + b"\n+\n" +
                                        record[3].rstrip(b"\r\n") + b"\n")
                    total_bytes += sum(map(len, record))
                total_bases += sum(lengths)
                count += 1
                if count % 100_000 == 0:
                    check_space(work, reserve)
                if count % 1_000_000 == 0:
                    print(f"[{sample['sample_id']}] query streamed {count:,} fragments", flush=True)
            if not count:
                raise RuntimeError("FASTQ input contains no records")
            if sources != input_states(sample):
                raise RuntimeError("FASTQ input changed during classification")
            process.stdin.close()
            result.update(fragments=count, reads=count * (2 if sample.get("r2") else 1),
                          paired=sample.get("r2") is not None, input_digest=digest.hexdigest(),
                          total_bases=total_bases, fastq_uncompressed_bytes=total_bytes,
                          validated_scope="prefix" if limit else "entire_input")
        except BaseException as exc:
            errors.append(exc)
            stopped.set()
            stop_process(process)
            try:
                process.stdin.close()
            except (OSError, ValueError):
                pass

    try:
        with log_path.open("w") as log:
            log.write("COMMAND " + json.dumps(command) + "\n")
            log.write("OMP_NUM_THREADS " + environment["OMP_NUM_THREADS"] + "\n")
            log.flush()
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=log, start_new_session=True, bufsize=1 << 20,
                                       env=environment)
            thread = threading.Thread(target=produce, name="metag-fastq-producer", daemon=True)
            thread.start()
            while True:
                line = process.stdout.readline(16 << 20)
                if not line:
                    break
                if not line.endswith(b"\n"):
                    raise RuntimeError("Native output row is truncated or exceeds 16 MiB")
                store.add_output(*parse_native_line(line, sample.get("r2") is not None, taxonomy))
            code = process.wait()
            thread.join(timeout=10)
            if thread.is_alive():
                raise RuntimeError("FASTQ producer did not stop after native exit")
            if errors:
                raise RuntimeError(f"FASTQ producer failed: {errors[0]}") from errors[0]
            if code:
                raise RuntimeError(f"Kraken classify exited {code}; see {log_path}")
            if store.received != result.get("fragments") or store.produced != result.get("fragments"):
                raise RuntimeError(f"Missing native ordinals: {store.received}/{result.get('fragments')}")
            validate_snapshot(snapshot)
            target = work / "taxids.bin"
            store.publish(target)
            seconds = time.monotonic() - started
            result.update(seconds=seconds, fragments_per_second=result["fragments"] / max(seconds, 1e-9),
                          command=command, native_executable=str(native), native_sha256=sha256(native),
                          taxids=str(target), taxids_state=file_state(target), taxids_sha256=sha256(target),
                          input_states=sources, requested_limit=limit, memory_mapping=True,
                          omp_num_threads=int(environment["OMP_NUM_THREADS"]))
            return result
    finally:
        stopped.set()
        stop_process(process)
        if thread is not None:
            thread.join(timeout=10)
        if process is not None:
            process.stdout.close()
            try:
                process.stdin.close()
            except (OSError, ValueError):
                pass
        if thread is None or not thread.is_alive():
            store.close()
