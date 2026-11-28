"""One native Kraken query with a disk-backed ordinal protocol, no read-name join."""
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
from common import check_space, file_state, sha256, kraken_memory_mapping
from taxonomy import Taxonomy

_U64, _LENGTHS = struct.Struct("<Q"), struct.Struct("<II")
_UINT = re.compile(rb"[0-9]+\Z")


def _read_process_metrics(pid, proc_root=Path("/proc")):
    """Best-effort Linux telemetry; missing/denied procfs never fails a query.

    VmPeak is virtual address space, not physical memory. VmHWM is the kernel's
    RSS high-water estimate at the sampling instant. Read bytes/page faults may
    miss the final interval after the process exits and its procfs entry closes.
    """
    folder = Path(proc_root) / str(pid)
    values = {}
    fields = {"VmRSS": "rss_peak_observed_bytes", "VmHWM": "rss_peak_observed_bytes",
              "VmSize": "virtual_peak_observed_bytes", "VmPeak": "virtual_peak_observed_bytes",
              "RssAnon": "anonymous_rss_peak_observed_bytes", "RssFile": "file_rss_peak_observed_bytes"}
    try:
        for line in (folder / "status").read_text().splitlines():
            name, separator, raw = line.partition(":")
            if separator and name in fields:
                parts = raw.split()
                if len(parts) == 2 and parts[1] == "kB":
                    key = fields[name]
                    values[key] = max(values.get(key, 0), int(parts[0]) * 1024)
    except (OSError, ValueError):
        pass
    try:
        # comm (field 2) may contain spaces or parentheses. Remaining fields
        # start at state (field 3); minflt/majflt are fields 10 and 12.
        fields = (folder / "stat").read_text().rsplit(") ", 1)[1].split()
        values["minor_faults_observed"] = int(fields[7])
        values["major_faults_observed"] = int(fields[9])
    except (OSError, ValueError, IndexError):
        pass
    try:
        for line in (folder / "io").read_text().splitlines():
            name, separator, raw = line.partition(":")
            if separator and name == "read_bytes":
                values["read_bytes_observed"] = int(raw.strip())
    except (OSError, ValueError):
        pass
    return values


class NativeTelemetry:
    """One lightweight monitor per native process, never one per read."""
    def __init__(self, pid, interval=2.0):
        self.pid, self.interval = pid, interval
        self.values, self.samples = {}, 0
        self.stopped = threading.Event()
        self.thread = threading.Thread(target=self._run, name="metag-native-memory", daemon=True)
        self._sample()
        self.thread.start()

    def _sample(self):
        values = _read_process_metrics(self.pid)
        if values:
            self.samples += 1
            for key, value in values.items():
                self.values[key] = max(self.values.get(key, 0), value)

    def _run(self):
        while not self.stopped.wait(self.interval):
            self._sample()

    def finish(self):
        self.stopped.set()
        self.thread.join()
        return dict(self.values, available=bool(self.values), samples=self.samples,
                    interval_seconds=self.interval, scope="native process, whole query batch",
                    note="Sampled Linux procfs; final interval may be missed. Virtual bytes are not RAM. "
                         "RSS excludes the Python worker and unmapped OS page cache; this is not a memory limit.")


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
    if snapshot.get("layout") == "partitioned":
        from partitioned_build import validate_snapshot as validate_union
        return validate_union(snapshot)
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
    from native_union import executable_path
    native = executable_path(cfg) if cfg.get("database_layout", "single") == "partitioned" else Path(cfg["tools_dir"]) / "kraken/bin/classify"
    native = native.resolve(strict=True)
    if not native.is_file() or not os.access(native, os.X_OK):
        raise RuntimeError(f"Native Kraken classify is not executable: {native}")
    probe = subprocess.run([str(native), "-h"], capture_output=True, text=True, timeout=15)
    help_text = probe.stdout + probe.stderr
    if cfg.get("database_layout", "single") == "partitioned" and "--metag-hashes" not in help_text:
        raise RuntimeError("Native classifier does not support the partitioned union catalog")
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
    union = snapshot.get("layout") == "partitioned"
    if union != (cfg.get("database_layout", "single") == "partitioned"):
        raise RuntimeError("Classifier layout and database snapshot disagree")
    first_hash = snapshot["hashes"][0] if union else "hash.k2d"
    command = [str(native), "-H", str(db / first_hash), "-t", str(db / "taxo.k2d"),
               "-o", str(db / "opts.k2d"), "-p", str(threads), "-T", str(confidence),
               "-g", str(groups), "-Q", str(quality)]
    if union:
        command += ["--metag-hashes", str(snapshot.get('_hashes_path') or db / "hashes.txt")]
    if kraken_memory_mapping(cfg, snapshot):
        command.append("-M")
    if paired:
        command.append("-S")
    return command  # No filename: stdin. -P would instead require two input files.


def require_preload_memory(cfg, snapshot):
    """Reject an oversized private index before launching the native process."""
    if kraken_memory_mapping(cfg, snapshot):
        return
    names = list(snapshot['hashes']) if snapshot.get('layout') == 'partitioned' else ['hash.k2d']
    index_gib = sum((Path(snapshot['db_dir']) / name).stat().st_size
                    for name in [*names, 'taxo.k2d', 'opts.k2d']) / 2**30
    allocation = float(cfg['classify_memory_gib'])
    required = index_gib + float(cfg['worker_memory_gib']) + float(cfg['system_reserve_gib'])
    if not math.isfinite(required) or not math.isfinite(allocation) or required > allocation:
        raise RuntimeError(f'RAM preload needs {required:.2f} GiB including worker/reserve; '
                           f'allocation is {allocation:.2f} GiB. Increase --mem-gib or explicitly '
                           'select mapped classification.')


def group_query_config(cfg, snapshot):
    """Compact group tables fit the reserved node RAM; full unions stay mapped."""
    result = dict(cfg, database_layout=snapshot['layout'], classifier_engine='kraken',
                  union_memory_mapping=True,
                  single_memory_mapping=cfg.get('group_memory_mapping', False))
    require_preload_memory(result, snapshot)
    return result


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


def query_once(cfg, sample, snapshot, work, limit=0, taxonomy=None, *, observe=None):
    """Classify one sample; the same engine also supports a resident-index batch."""
    return query_many(cfg, [sample], snapshot, work, limit=limit, taxonomy=taxonomy, observe=observe)[0]


def query_many(cfg, samples, snapshot, work, limit=0, taxonomy=None, *, observe=None):
    """Load the unchanged Kraken indexes once and stream a bounded sample batch.

    Native query IDs are globally unique ordinals. Per-sample boundaries are
    recorded by the producer, then compact assignments are sliced by those exact
    boundaries. Read names never serve as join keys, including across samples.
    """
    samples = list(samples)
    if not samples or len({s["sample_id"] for s in samples}) != len(samples):
        raise ValueError("A query batch needs distinct sample IDs")
    paired = samples[0].get("r2") is not None
    if any((s.get("r2") is not None) != paired for s in samples):
        raise ValueError("A native batch cannot mix paired and single-end inputs")
    sample = samples[0]
    work = Path(work)
    validate_snapshot(snapshot)
    require_preload_memory(cfg, snapshot)
    reserve = float(cfg.get("min_free_gib", 100))
    check_space(work, reserve)
    native = resolve_native(cfg)
    command = classifier_command(native, snapshot, cfg, sample.get("r2") is not None)
    environment = os.environ.copy()
    # Match the already validated -p budget. Keep site library, locale, temporary
    # directory and scheduler settings; do not mutate the parent environment.
    environment["OMP_NUM_THREADS"] = str(int(cfg.get("classify_threads", 8)))
    environment["K2_DB_READ_THREADS"] = str(max(1, min(int(cfg.get("index_read_threads", 8)), int(cfg.get("classify_threads", 8)))))
    sources = [input_states(item) for item in samples]
    taxonomy = taxonomy if taxonomy is not None else Taxonomy(snapshot["taxonomy_dir"])
    store = OrdinalStore(work, reserve)
    result, errors, stopped = {}, [], threading.Event()
    process = thread = telemetry = None
    native_metrics = None
    log_path = work / "query.log"
    started = time.monotonic()

    def produce():
        global_count = 0
        results = []
        try:
            for current, source_state in zip(samples, sources):
                digest = new_digest()
                total_bases = total_bytes = count = 0
                first = global_count
                for a, b in fragments(current["r1"], current.get("r2"), limit):
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
                    global_count += 1
                    if count % 100_000 == 0:
                        check_space(work, reserve)
                        print(f"[{current['sample_id']}] query streamed {count:,} fragments", flush=True)
                if not count:
                    raise RuntimeError(f"FASTQ input contains no records: {current['sample_id']}")
                if source_state != input_states(current):
                    raise RuntimeError("FASTQ input changed during classification")
                results.append(dict(sample_id=current["sample_id"], fragments=count,
                    reads=count * (2 if paired else 1), paired=paired,
                    input_digest=digest.hexdigest(), total_bases=total_bases,
                    fastq_uncompressed_bytes=total_bytes,
                    validated_scope="prefix" if limit else "entire_input",
                    input_states=source_state, first_ordinal=first + 1))
            process.stdin.close()
            result.update(fragments=global_count, samples=results)
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
            mode = "mmap" if "-M" in command else "ram"
            log.write("MEMORY_MODE " + mode + "\n")
            if mode == "mmap":
                note = ("Demand-paged index: the complete index is not preloaded; virtual size may exceed RAM. "
                        "Resident pages depend on OS reclaim and storage I/O; no hard RSS limit is imposed.")
                log.write("MEMORY_NOTE " + note + "\n")
                print("Mapped Kraken index. " + note, flush=True)
            else:
                print(f"Preloading Kraken index into private RAM; allocation={cfg['classify_memory_gib']} GiB; "
                      f"read_threads={environment['K2_DB_READ_THREADS']}. Classification starts after loading.", flush=True)
            log.flush()
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=log, start_new_session=True, bufsize=1 << 20,
                                       env=environment)
            telemetry = NativeTelemetry(process.pid)
            thread = threading.Thread(target=produce, name="metag-fastq-producer", daemon=True)
            thread.start()
            while True:
                line = process.stdout.readline(16 << 20)
                if not line:
                    break
                if not line.endswith(b"\n"):
                    raise RuntimeError("Native output row is truncated or exceeds 16 MiB")
                parsed = parse_native_line(line, paired, taxonomy)
                store.add_output(*parsed)
                if observe is not None:
                    observe(parsed[0], parsed[1], line)
                if store.received % 100_000 == 0:
                    print(f"[{sample['sample_id']}] query classified {store.received:,} fragments; elapsed={time.monotonic() - started:.1f}s", flush=True)
            code = process.wait()
            native_metrics = telemetry.finish()
            log.write("NATIVE_METRICS " + json.dumps(native_metrics, sort_keys=True) + "\n")
            log.flush()
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
            index_seconds = sum(float(value) for value in re.findall(r"METAG_INDEX_READY seconds=([0-9.eE+\-]+)", log_path.read_text(errors="replace")))
            queries = []
            with target.open("rb", buffering=8 << 20) as compact:
                for index, item in enumerate(result["samples"]):
                    if len(samples) == 1:
                        destination = target
                    else:
                        destination = work / f"taxids.{index + 1}.bin"
                        remaining = item["fragments"] * 8
                        with destination.open("wb", buffering=8 << 20) as out:
                            while remaining:
                                block = compact.read(min(8 << 20, remaining))
                                if not block:
                                    raise RuntimeError("Incomplete native batch assignment slice")
                                out.write(block)
                                remaining -= len(block)
                    share = seconds * item["fragments"] / result["fragments"]
                    item.update(seconds=share, fragments_per_second=result["fragments"] / max(seconds, 1e-9),
                        command=command, native_executable=str(native),
                        taxids=str(destination), taxids_state=file_state(destination),
                        requested_limit=limit, memory_mapping="-M" in command,
                        memory_mode=mode, native_metrics=native_metrics,
                        omp_num_threads=int(environment["OMP_NUM_THREADS"]),
                        log_path=str(log_path), batch_samples=len(samples),
                        batch_fragments=result["fragments"], batch_classification_seconds=seconds,
                        classification_timing="batch wall time apportioned by fragment count" if len(samples) > 1 else "sample wall time",
                        native_index_loads=1, batch_index_load_seconds=index_seconds,
                        batch_query_seconds_excluding_index_load=max(0.0, seconds - index_seconds))
                    queries.append(item)
                if len(samples) > 1 and compact.read(1):
                    raise RuntimeError("Extra native batch assignments after sample slicing")
            if len(samples) > 1:
                target.unlink()
            return queries
    finally:
        stopped.set()
        stop_process(process)
        if telemetry is not None and native_metrics is None:
            # Preserve measurements on failure too; logging must not replace the
            # actual classifier error (e.g. an out-of-memory kill).
            native_metrics = telemetry.finish()
            try:
                with log_path.open("a") as log:
                    log.write("NATIVE_METRICS " + json.dumps(native_metrics, sort_keys=True) + "\n")
            except OSError:
                pass
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

