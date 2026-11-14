"""Validated ordinal adapters for pinned Kun-peng direct and staged CPU modes.

The upstream staged resolver emits only reads with annotated hits.  A missing
row is therefore *not* accepted as an unclassified read unless the complete
splitr input map and the pre-resolve annotation files prove that it had no hit.
"""
from contextlib import ExitStack
import errno
import gzip
import json
import math
import mmap
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import time

from fastq_stream import fragments, new_digest, update_digest
from common import check_space, file_state, sha256
from ordinal_query import OrdinalStore, input_states, parse_native_line, stop_process
from taxonomy import Taxonomy

MODES = {"kunpeng_full": "direct", "kunpeng_memory": "staged"}
_U64 = struct.Struct("<Q")
_LENGTHS = struct.Struct("<II")
_ROW = struct.Struct("<III")


def parameters(cfg):
    engine = cfg.get("classifier_engine", "kraken")
    if engine not in MODES:
        raise ValueError(f"Not a supported Kun-peng mode: {engine}")
    result = {"classifier_engine": engine,
              "kunpeng_input_compression": cfg.get("kunpeng_input_compression", "gzip"),
              "kunpeng_buffer_bytes": int(cfg.get("kunpeng_buffer_bytes", 16 << 20)),
              "kunpeng_batch_size": int(cfg.get("kunpeng_batch_size", 4)),
              "kunpeng_max_work_gib": float(cfg.get("kunpeng_max_work_gib", 0))}
    if result["kunpeng_input_compression"] not in ("gzip", "none"):
        raise ValueError("kunpeng_input_compression must be gzip or none")
    if result["kunpeng_buffer_bytes"] < 16 or not 1 <= result["kunpeng_batch_size"] <= 32:
        raise ValueError("Invalid Kun-peng buffer or batch size")
    if not math.isfinite(result["kunpeng_max_work_gib"]) or result["kunpeng_max_work_gib"] < 0:
        raise ValueError("kunpeng_max_work_gib must be finite and nonnegative")
    threads = int(cfg.get("classify_threads", 8))
    confidence = float(cfg.get("confidence", 0.1))
    groups = int(cfg.get("minimum_hit_groups", 2))
    quality = int(cfg.get("minimum_base_quality", 0))
    if threads < 3:
        raise ValueError("Pinned Kun-peng requires classify_threads >= 3 (seqkmer reader constraint)")
    if groups < 1 or not 0 <= quality <= 255 or not math.isfinite(confidence) or not 0 <= confidence <= 1:
        raise ValueError("Invalid Kun-peng threads/confidence/hit groups/base quality")
    return result


def resolve_native(cfg):
    native = (Path(cfg["tools_dir"]) / "kunpeng/bin/kun_peng").resolve(strict=True)
    if not native.is_file() or not os.access(native, os.X_OK):
        raise RuntimeError(f"Kun-peng executable missing: {native}")
    probe = subprocess.run([str(native), "--help"], capture_output=True, text=True, timeout=15)
    if probe.returncode or any(name not in probe.stdout + probe.stderr
                               for name in ("direct", "splitr", "annotate", "resolve")):
        raise RuntimeError("Kun-peng executable does not expose the verified native stages")
    return native


def validate_snapshot(snapshot):
    if snapshot.get("status") != "complete" or not snapshot.get("identity"):
        raise RuntimeError("Only a complete Kun-peng database snapshot can be queried")
    if not snapshot.get("source_database_identity"):
        raise RuntimeError("Kun-peng snapshot lacks its verified Kraken source identity")
    db, taxonomy = Path(snapshot["db_dir"]), Path(snapshot["taxonomy_dir"])
    files = snapshot.get("files", {})
    if not {"hash_config.k2d", "opts.k2d", "taxo.k2d"} <= files.keys():
        raise RuntimeError("Kun-peng snapshot lacks required database files")
    if not any(re.fullmatch(r"hash_\d+\.k2d", name) for name in files):
        raise RuntimeError("Kun-peng snapshot has no recorded hash shards")
    for name, item in files.items():
        if Path(name).name != name or file_state(db / name) != item["state"]:
            raise RuntimeError(f"Frozen Kun-peng database changed: {name}")
    recorded_shards = {name for name in files if re.fullmatch(r"hash_\d+\.k2d", name)}
    actual_shards = {p.name for p in db.glob("hash_*.k2d") if re.fullmatch(r"hash_\d+\.k2d", p.name)}
    if recorded_shards != actual_shards:
        raise RuntimeError("Kun-peng hash shard inventory changed")
    taxfiles = snapshot.get("taxonomy_files", {})
    if not {"nodes.dmp", "names.dmp"} <= taxfiles.keys():
        raise RuntimeError("Kun-peng snapshot lacks frozen taxonomy")
    for name, item in taxfiles.items():
        if Path(name).name != name or file_state(taxonomy / name) != item["state"]:
            raise RuntimeError(f"Frozen Kun-peng taxonomy changed: {name}")
    marker = db / "database.json"
    if marker.is_file():
        value = json.loads(marker.read_text())
        if value.get("status") != "complete" or value.get("identity") != snapshot["identity"]:
            raise RuntimeError("Kun-peng completion marker changed")


def _work_bytes(work):
    total = 0
    for folder, dirs, files in os.walk(work):
        for name in (*dirs, *files):
            path = Path(folder) / name
            if path.is_symlink():
                raise RuntimeError(f"Symlink appeared in private Kun-peng work: {path}")
        for name in files:
            try:
                total += (Path(folder) / name).stat().st_size
            except FileNotFoundError:
                pass  # Official annotate/resolve remove their consumed files.
    return total


def _guard(cfg, work):
    check_space(work, float(cfg.get("min_free_gib", 100)))
    size = _work_bytes(work)
    maximum = float(cfg.get("kunpeng_max_work_gib", 0))
    if maximum and size > maximum * 2**30:
        raise RuntimeError(f"Kun-peng private work exceeded {maximum:g} GiB: {size / 2**30:.3f} GiB")
    return size


def _process_stats(pid):
    values = {}
    try:
        for line in (Path("/proc") / str(pid) / "status").read_text().splitlines():
            if line.startswith(("VmRSS:", "VmHWM:", "VmSize:")):
                name, count, unit = line.split()
                values[name.rstrip(":")] = int(count) * 1024
    except (OSError, ValueError):
        pass
    return values


def _run(cfg, command, work, log, stage):
    environment = os.environ.copy()
    environment["RAYON_NUM_THREADS"] = str(int(cfg.get("classify_threads", 8)))
    environment["OMP_NUM_THREADS"] = environment["RAYON_NUM_THREADS"]
    log.write("STAGE " + stage + " COMMAND " + json.dumps(command) + "\n")
    log.flush()
    process = None
    started = time.monotonic()
    peak_rss = peak_vms = 0
    peak_work = _guard(cfg, work)
    last_guard = started
    try:
        process = subprocess.Popen(command, stdout=log, stderr=log, env=environment,
                                   start_new_session=True)
        while process.poll() is None:
            stats = _process_stats(process.pid)
            peak_rss = max(peak_rss, stats.get("VmRSS", 0), stats.get("VmHWM", 0))
            peak_vms = max(peak_vms, stats.get("VmSize", 0))
            if time.monotonic() - last_guard >= 1:
                peak_work = max(peak_work, _guard(cfg, work))
                last_guard = time.monotonic()
            try:
                process.wait(timeout=0.1)
            except subprocess.TimeoutExpired:
                pass
        if process.returncode:
            raise RuntimeError(f"Kun-peng {stage} exited {process.returncode}; see {work / 'query.log'}")
        peak_work = max(peak_work, _guard(cfg, work))
        return {"stage": stage, "command": command, "seconds": time.monotonic() - started,
                "peak_rss_bytes_sampled": peak_rss, "peak_vms_bytes_sampled": peak_vms,
                "peak_work_bytes_sampled": peak_work,
                "memory_measurement": "Linux /proc status sampled at about 100 ms; not exact job accounting"}
    finally:
        stop_process(process)


def _stage_inputs(cfg, sample, store, work, limit):
    paths, digest = [], new_digest()
    bases = byte_count = count = 0
    started = time.monotonic()
    compressed = cfg.get("kunpeng_input_compression", "gzip") == "gzip"
    with ExitStack() as stack:
        handles = []
        for mate in (1, 2) if sample.get("r2") else (1,):
            path = work / (f"ordinal.R{mate}.fastq" + (".gz" if compressed else ""))
            paths.append(path)
            raw = stack.enter_context(path.open("xb", buffering=1 << 20))
            handle = stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=raw,
                                                       compresslevel=1, mtime=0)) if compressed else raw
            handles.append(handle)
        for a, b in fragments(sample["r1"], sample.get("r2"), limit):
            update_digest(digest, a, b)
            records = (a, b) if b is not None else (a,)
            lengths = [len(record[1].rstrip(b"\r\n")) for record in records]
            ordinal = store.add_input(lengths[0], lengths[1] if b is not None else 0)
            if ordinal > 0xFFFFFFFE:
                raise RuntimeError("Kun-peng's pinned internal uint32 read index limit was reached")
            header = b"@" + str(ordinal).encode("ascii") + b"\n"
            for handle, record in zip(handles, records):
                handle.write(header + record[1].rstrip(b"\r\n") + b"\n+\n" +
                             record[3].rstrip(b"\r\n") + b"\n")
                byte_count += sum(map(len, record))
            count += 1
            bases += sum(lengths)
            if count % 100_000 == 0:
                _guard(cfg, work)
            if count % 1_000_000 == 0:
                print(f"[{sample['sample_id']}] Kun-peng staged {count:,} fragments", flush=True)
    if not count:
        raise RuntimeError("FASTQ input contains no records")
    for path in paths:
        with path.open("rb") as handle:
            os.fsync(handle.fileno())
    return paths, {"fragments": count, "reads": count * (2 if sample.get("r2") else 1),
                   "paired": sample.get("r2") is not None, "input_digest": digest.hexdigest(),
                   "total_bases": bases, "fastq_uncompressed_bytes": byte_count,
                   "validated_scope": "prefix" if limit else "entire_input",
                   "normalized_input_bytes": sum(p.stat().st_size for p in paths),
                   "input_staging_seconds": time.monotonic() - started}


def _test_bit(bitmap, ordinal):
    index = ordinal - 1
    return bool(bitmap[index // 8] & (1 << (index % 8)))


def _set_bit(bitmap, ordinal):
    index = ordinal - 1
    bitmap[index // 8] |= 1 << (index % 8)


def _lengths(store, ordinal):
    return _LENGTHS.unpack_from(store.maps[1], (ordinal - 1) * 8)


def _uint(text, label):
    if not re.fullmatch(rb"[0-9]+", text):
        raise RuntimeError(f"Invalid Kun-peng {label}: {text[:80]!r}")
    return int(text)


class StagedEvidence:
    """Disk-backed correspondence between native IDs, ordinals and hit presence."""
    def __init__(self, store, work):
        self.store, self.work = store, Path(work)
        self.handles, self.maps = [], []
        self.paths = [self.work / "native-ids.bin", self.work / "native-hit-bitmap.bin"]
        self.hit_count = 0
        try:
            sizes = ((store.produced + 1) * 8, (store.produced + 7) // 8)
            check_space(self.work, store.min_free_gib + sum(sizes) / 2**30)
            for path, size in zip(self.paths, sizes):
                handle = path.open("w+b")
                self.handles.append(handle)
                try:
                    os.posix_fallocate(handle.fileno(), 0, size)
                except (AttributeError, OSError) as exc:
                    if isinstance(exc, OSError) and exc.errno not in (errno.EOPNOTSUPP, errno.ENOSYS, errno.EINVAL):
                        raise
                    handle.truncate(size)
                handle.flush()
                self.maps.append(mmap.mmap(handle.fileno(), size))
        except BaseException:
            self.close()
            raise

    def verify_map(self, chunk):
        paths = list(Path(chunk).glob("sample_id*.map"))
        if len(paths) != 1 or paths[0].name != "sample_id_1.map":
            raise RuntimeError("Kun-peng splitr did not produce exactly one input ID map")
        count = 0
        with paths[0].open("rb") as handle:
            for line in handle:
                fields = line.rstrip(b"\r\n").split(b"\t")
                if len(fields) != 4:
                    raise RuntimeError("Malformed Kun-peng splitr input map")
                internal = _uint(fields[0], "internal ID")
                ordinal = _uint(fields[1], "ordinal ID")
                if not 0 <= internal <= self.store.produced or not 1 <= ordinal <= self.store.produced:
                    raise RuntimeError("Kun-peng input map ID out of range")
                if _U64.unpack_from(self.maps[0], internal * 8)[0] or _test_bit(self.store.maps[2], ordinal):
                    raise RuntimeError("Duplicate internal ID or ordinal in Kun-peng input map")
                parts = fields[2].split(b"|")
                expected = _lengths(self.store, ordinal)
                lengths = tuple(_uint(x, "input-map length") for x in parts)
                if lengths != (expected if expected[1] else (expected[0],)):
                    raise RuntimeError("Kun-peng splitr input-map read lengths differ")
                sizes = fields[3].split(b"|")
                if len(sizes) != len(parts) or any(not re.fullmatch(rb"[0-9]+", value) for value in sizes):
                    raise RuntimeError("Kun-peng splitr minimizer count fields malformed")
                _U64.pack_into(self.maps[0], internal * 8, ordinal)
                _set_bit(self.store.maps[2], ordinal)
                count += 1
        if count != self.store.produced:
            raise RuntimeError(f"Kun-peng splitr omitted input ordinals: {count}/{self.store.produced}")
        # The store bitmap is now reserved for exactly-once output validation.
        for start in range(0, len(self.store.maps[2]), 1 << 20):
            size = min(1 << 20, len(self.store.maps[2]) - start)
            self.store.maps[2][start:start + size] = b"\0" * size
        return count

    def verify_annotations(self, chunk):
        if sys.byteorder != "little":
            raise RuntimeError("Pinned Kun-peng annotation format requires a little-endian host")
        paths = sorted(Path(chunk).glob("sample_file*.bin"))
        row_count = 0
        for path in paths:
            if not re.fullmatch(r"sample_file_1_\d+\.bin", path.name):
                raise RuntimeError("Unexpected sample file in Kun-peng annotations")
            with path.open("rb") as handle:
                while True:
                    block = handle.read(_ROW.size * 65536)
                    if not block:
                        break
                    if len(block) % _ROW.size:
                        raise RuntimeError("Truncated Kun-peng annotation row")
                    for value, internal, minimizer in _ROW.iter_unpack(block):
                        if internal > self.store.produced:
                            raise RuntimeError("Kun-peng annotation ID out of range")
                        ordinal = _U64.unpack_from(self.maps[0], internal * 8)[0]
                        if not ordinal:
                            raise RuntimeError("Kun-peng annotation ID absent from verified input map")
                        if not _test_bit(self.maps[1], ordinal):
                            _set_bit(self.maps[1], ordinal)
                            self.hit_count += 1
                        row_count += 1
        return row_count

    def consume_outputs(self, output, taxonomy):
        files = list(Path(output).glob("output_*.txt"))
        if (self.hit_count and len(files) != 1) or (files and (len(files) != 1 or files[0].name != "output_1.txt")):
            raise RuntimeError("Kun-peng resolve output inventory differs from verified hits")
        if files:
            _consume(files[0], self.store, taxonomy, expected_hits=self.maps[1])
        if self.store.received != self.hit_count:
            raise RuntimeError(f"Kun-peng resolve omitted hit-bearing ordinals: {self.store.received}/{self.hit_count}")
        for ordinal in range(1, self.store.produced + 1):
            if not _test_bit(self.maps[1], ordinal):
                self.store.add_output(ordinal, 0, _lengths(self.store, ordinal))
        return self.store.produced - self.hit_count

    def close(self):
        for mapping in self.maps:
            if not mapping.closed:
                mapping.close()
        for handle in self.handles:
            handle.close()


def _verify_query_chunks(snapshot, chunk):
    shards = sorted(int(re.fullmatch(r"hash_(\d+)\.k2d", name)[1])
                    for name in snapshot["files"] if re.fullmatch(r"hash_\d+\.k2d", name))
    if shards != list(range(1, len(shards) + 1)):
        raise RuntimeError("Kun-peng database shard numbering is not contiguous from 1")
    expected = {f"sample_{index}.k2" for index in shards}
    paths = {p.name: p for p in Path(chunk).glob("sample_*.k2")}
    if set(paths) != expected:
        raise RuntimeError("Kun-peng splitr query-chunk inventory is incomplete")
    for index in shards:
        path = paths[f"sample_{index}.k2"]
        if path.stat().st_size < 16 or path.stat().st_size % 16:
            raise RuntimeError("Kun-peng splitr query chunk is truncated")
        with path.open("rb") as handle:
            partition, capacity = struct.unpack("<QQ", handle.read(16))
        if partition != index - 1 or not capacity:
            raise RuntimeError("Kun-peng splitr query-chunk header differs from database")
        if snapshot.get("hash_capacity") and capacity != snapshot["hash_capacity"]:
            raise RuntimeError("Kun-peng query-chunk capacity differs from database")


def _consume(path, store, taxonomy, expected_hits=None):
    paired = _lengths(store, 1)[1] != 0
    with Path(path).open("rb") as handle:
        while True:
            line = handle.readline(16 << 20)
            if not line:
                break
            if not line.endswith(b"\n"):
                raise RuntimeError("Kun-peng output row truncated or exceeds 16 MiB")
            values = parse_native_line(line, paired, taxonomy)
            if expected_hits is not None:
                if not 1 <= values[0] <= store.produced or not _test_bit(expected_hits, values[0]):
                    raise RuntimeError("Kun-peng resolve emitted an ordinal without verified hits")
            store.add_output(*values)


def query_once(cfg, sample, snapshot, work, limit=0, taxonomy=None):
    params = parameters(cfg)
    validate_snapshot(snapshot)
    work = Path(work)
    _guard(cfg, work)
    native = resolve_native(cfg)
    native_sha = sha256(native)
    sources = input_states(sample)
    store = OrdinalStore(work, float(cfg.get("min_free_gib", 100)))
    evidence = None
    started = time.monotonic()
    stages = []
    try:
        inputs, result = _stage_inputs(cfg, sample, store, work, limit)
        if sources != input_states(sample):
            raise RuntimeError("FASTQ changed while staging Kun-peng inputs")
        output = work / "native-output"
        output.mkdir()
        common = ["--db", snapshot["db_dir"], "-p", str(int(cfg.get("classify_threads", 8)))]
        thresholds = ["-T", str(float(cfg.get("confidence", 0.1))),
                      "-g", str(int(cfg.get("minimum_hit_groups", 2)))]
        quality = ["-Q", str(int(cfg.get("minimum_base_quality", 0)))]
        paired = ["-P"] if sample.get("r2") else []
        taxonomy = taxonomy if taxonomy is not None else Taxonomy(snapshot["taxonomy_dir"])
        zero_hits = 0
        with (work / "query.log").open("w") as log:
            log.write("ADAPTER_PARAMETERS " + json.dumps(params, sort_keys=True) + "\n")
            log.write("NORMALIZED_INPUT " + json.dumps({k: result[k] for k in
                ("fragments", "normalized_input_bytes", "input_staging_seconds")}) + "\n")
            if params["classifier_engine"] == "kunpeng_full":
                command = [str(native), "direct", *common, *thresholds, *quality,
                           "--output-dir", str(output), *paired, *map(str, inputs)]
                stages.append(_run(cfg, command, work, log, "direct"))
                paths = list(output.glob("output_*.txt"))
                if len(paths) != 1 or paths[0].name != "output_1.txt":
                    raise RuntimeError("Kun-peng direct did not produce exactly one classification file")
                _consume(paths[0], store, taxonomy)
            else:
                chunk = work / "native-chunks"
                chunk.mkdir()
                command = [str(native), "splitr", *common, *quality, "--chunk-dir", str(chunk),
                           *paired, *map(str, inputs)]
                stages.append(_run(cfg, command, work, log, "splitr"))
                evidence = StagedEvidence(store, work)
                evidence.verify_map(chunk)
                _verify_query_chunks(snapshot, chunk)
                command = [str(native), "annotate", *common, "--chunk-dir", str(chunk),
                           "--buffer-size", str(params["kunpeng_buffer_bytes"]),
                           "--batch-size", str(params["kunpeng_batch_size"])]
                stages.append(_run(cfg, command, work, log, "annotate"))
                if list(chunk.glob("sample_*.k2")):
                    raise RuntimeError("Kun-peng annotate left unconsumed query chunks")
                annotation_rows = evidence.verify_annotations(chunk)
                command = [str(native), "resolve", *common, *thresholds,
                           "--chunk-dir", str(chunk), "--output-dir", str(output)]
                stages.append(_run(cfg, command, work, log, "resolve"))
                zero_hits = evidence.consume_outputs(output, taxonomy)
                result.update(native_annotation_rows=annotation_rows,
                              native_hit_fragments=evidence.hit_count,
                              verified_zero_hit_fragments=zero_hits,
                              zero_hit_policy="filled U only after complete splitr map and annotation bitmap validation")
        if store.received != result["fragments"]:
            raise RuntimeError(f"Kun-peng output incomplete: {store.received}/{result['fragments']}")
        if sources != input_states(sample) or sha256(native) != native_sha:
            raise RuntimeError("Kun-peng input or executable changed during classification")
        validate_snapshot(snapshot)
        target = work / "taxids.bin"
        store.publish(target)
        seconds = time.monotonic() - started
        result.update(seconds=seconds, fragments_per_second=result["fragments"] / max(seconds, 1e-9),
                      native_seconds=sum(stage["seconds"] for stage in stages), stages=stages,
                      classifier_engine=params["classifier_engine"], parameters=params,
                      native_executable=str(native), native_sha256=native_sha,
                      native_peak_rss_bytes_sampled=max(stage["peak_rss_bytes_sampled"] for stage in stages),
                      peak_work_bytes_sampled=max(_guard(cfg, work), *(stage["peak_work_bytes_sampled"] for stage in stages)),
                      taxids=str(target), taxids_state=file_state(target), taxids_sha256=sha256(target),
                      input_states=sources, requested_limit=limit, memory_mapping=False,
                      timing_note="seconds includes input normalization, native stages, and strict output validation")
        return result
    finally:
        if evidence is not None:
            evidence.close()
        store.close()
