"""One pass over original FASTQ, writing exact records into taxonomic bins."""
from collections import Counter, OrderedDict
from contextlib import ExitStack
import csv
import gzip
import json
from pathlib import Path
import struct
import time

from fastq_stream import fragments, new_digest, read_id, update_digest
from hpc_common import atomic_json, check_space, file_state
from taxonomy import Taxonomy, safe_name

class BufferedWriters:
    """Bound both pending FASTQ bytes and open gzip files independently of reads.

    Small per-bin buffers avoid opening and closing a gzip member for each read
    when hundreds of species are interleaved.  Least-recently-used buffers spill
    before the total budget is exceeded; file handles have a separate LRU bound.
    Eviction creates valid concatenated gzip members when a bin is reopened.
    """

    def __init__(self, root, *, max_open=32, buffer_bytes=64 << 20,
                 file_buffer_bytes=1 << 20):
        if min(max_open, buffer_bytes, file_buffer_bytes) < 1:
            raise ValueError("Writer limits must be positive")
        self.root = Path(root)
        self.max_open = max_open
        self.buffer_bytes = buffer_bytes
        self.file_buffer_bytes = min(file_buffer_bytes, buffer_bytes)
        self.buffers = OrderedDict()
        self.handles = OrderedDict()
        self.pending_bytes = 0
        self.paths = set()

    def _handle(self, key):
        if key in self.handles:
            handle = self.handles.pop(key)
            self.handles[key] = handle
            return handle[0]
        while len(self.handles) >= self.max_open:
            _, (compressed, raw) = self.handles.popitem(last=False)
            try:
                compressed.close()
            finally:
                raw.close()
        path = self.root / key
        path.parent.mkdir(parents=True, exist_ok=True)
        raw = path.open("ab", buffering=256 << 10)
        try:
            compressed = gzip.GzipFile(filename="", fileobj=raw, mode="ab",
                                       compresslevel=1, mtime=0)
        except BaseException:
            raw.close()
            raise
        self.handles[key] = (compressed, raw)
        return compressed

    def _flush(self, key):
        pending = self.buffers.pop(key)
        self.pending_bytes -= len(pending)
        self._handle(key).write(pending)

    def write(self, key, record):
        relative = Path(key)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Bin path must stay within the output directory")
        self.paths.add(key)
        pending = self.buffers.pop(key, bytearray())
        size = sum(map(len, record))
        for line in record:
            pending.extend(line)
        self.pending_bytes += size
        self.buffers[key] = pending
        if len(pending) >= self.file_buffer_bytes:
            self._flush(key)
        while self.pending_bytes > self.buffer_bytes:
            self._flush(next(iter(self.buffers)))

    def close(self, *, flush=True):
        try:
            if flush:
                while self.buffers:
                    self._flush(next(iter(self.buffers)))
            else:
                self.buffers.clear()
                self.pending_bytes = 0
        finally:
            while self.handles:
                _, (compressed, raw) = self.handles.popitem(last=False)
                try:
                    compressed.close()
                finally:
                    raw.close()


def taxids(path, expected_fragments):
    """Read fixed-size blocks; reject truncation or appended assignments."""
    with Path(path).open("rb", buffering=8 << 20) as handle:
        count = 0
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            if len(block) % 8:
                raise RuntimeError("Truncated compact classification record")
            for (taxid,) in struct.iter_unpack("<Q", block):
                count += 1
                if count > expected_fragments:
                    raise RuntimeError("More compact assignments than expected fragments")
                yield taxid
        if count != expected_fragments:
            raise RuntimeError(f"Compact assignment count differs: {count} != {expected_fragments}")


def split_originals(cfg, sample, query, snapshot, out):
    """Write into a new private directory; the orchestrator publishes it later."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    reserve = float(cfg.get("min_free_gib", 100))
    check_space(out, reserve)
    expected = int(query["fragments"])
    classifications = Path(query["taxids"])
    original_state = file_state(classifications)
    if original_state != query["taxids_state"] or original_state["size"] != expected * 8:
        raise RuntimeError("Ordinal classification file changed or is incomplete")
    taxonomy = Taxonomy(snapshot["taxonomy_dir"])
    counts, cache = Counter(), {}
    writers = BufferedWriters(out, max_open=int(cfg.get("max_open_bins", 32)))
    digest = new_digest()
    started = time.monotonic()
    n = bases = raw_bytes = 0
    completed = False
    iterator = fragments(sample["r1"], sample.get("r2"), int(query["requested_limit"]))
    try:
        with gzip.open(out / "read-assignments.tsv.gz", "wt", compresslevel=1,
                       encoding="utf-8", newline="") as handle:
            assignments = csv.writer(handle, delimiter="\t", lineterminator="\n")
            assignments.writerow(("ordinal", "r1_read_id", "r2_read_id", "assigned_taxid",
                                  "assigned_rank", "species_taxid", "group", "status", "bin"))
            for tid in taxids(classifications, expected):
                pair = next(iterator, None)
                if pair is None:
                    raise RuntimeError("More ordinal assignments than original FASTQ fragments")
                a, b = pair
                update_digest(digest, a, b)
                if tid not in cache:
                    info = taxonomy.describe(tid)
                    if info["taxid"] != tid:
                        raise RuntimeError("Ordinal taxid was not normalized to species")
                    group, name = info["group"], safe_name(info["name"])
                    if info["status"] == "assigned_species":
                        base = f"{group}/{tid}_{name}"
                    elif tid:
                        base = f"unresolved/{group}/{tid}_{name}"
                    else:
                        base = "unresolved/unclassified"
                    cache[tid] = (info, base)
                info, base = cache[tid]
                n += 1
                assignments.writerow((n, read_id(a[0]).decode("utf-8", "backslashreplace"),
                                      read_id(b[0]).decode("utf-8", "backslashreplace") if b else "",
                                      tid, info["rank"], info["species_taxid"], info["group"],
                                      info["status"], base))
                if b is None:
                    writers.write(f"{base}/SE.fastq.gz", a)
                else:
                    writers.write(f"{base}/PE.R1.fastq.gz", a)
                    writers.write(f"{base}/PE.R2.fastq.gz", b)
                for record in (a, b) if b is not None else (a,):
                    bases += len(record[1].rstrip(b"\r\n"))
                    raw_bytes += sum(map(len, record))
                counts[tid] += 1
                if n % 100_000 == 0:
                    check_space(out, reserve)
                if n % 1_000_000 == 0:
                    print(f"[{sample['sample_id']}] split {n:,}/{expected:,} fragments", flush=True)
            if next(iterator, None) is not None:
                raise RuntimeError("Fewer ordinal assignments than original FASTQ fragments")
            completed = True
    finally:
        try:
            writers.close(flush=completed)
        finally:
            iterator.close()
    if (n != expected or digest.hexdigest() != query["input_digest"] or
            bases != query["total_bases"] or raw_bytes != query["fastq_uncompressed_bytes"]):
        raise RuntimeError("Original FASTQ bytes/order/count differ from classification pass")
    if file_state(classifications) != original_state:
        raise RuntimeError("Ordinal classification file changed while splitting")
    rows = []
    for tid in sorted(counts):
        info, base = cache[tid]
        rows.append({"assigned_taxid": tid, "assigned_rank": info["rank"],
                     "species_taxid": info["species_taxid"], "name": info["name"],
                     "group": info["group"], "status": info["status"], "bin": base,
                     "fragments": counts[tid], "reads": counts[tid] * (2 if sample.get("r2") else 1)})
    with (out / "counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    status_counts, group_counts = Counter(), Counter()
    for row in rows:
        status_counts[row["status"]] += row["fragments"]
        group_counts[row["group"]] += row["fragments"]
    if sum(counts.values()) != expected or sum(row["reads"] for row in rows) != query["reads"]:
        raise RuntimeError("Output count conservation check failed")
    result = {"schema": "metag-hpc-split-v1", "sample_id": sample["sample_id"],
              "fragments": n, "reads": query["reads"], "paired": query["paired"],
              "total_bases": bases, "input_digest": digest.hexdigest(),
              "validated_scope": query["validated_scope"], "requested_limit": query["requested_limit"],
              "classification_seconds": query["seconds"], "split_seconds": time.monotonic() - started,
              "groups": dict(group_counts), "statuses": dict(status_counts),
              "fastq_files": len(writers.paths), "count_conservation": True,
              "pair_policy": "both original mates share one fragment classification and bin",
              "taxonomy_policy": "below-species calls promoted; higher ranks retained under unresolved",
              "reference_hits": "not collected; read assignments are taxonomic, not assembly alignments"}
    atomic_json(out / "summary.json", result)
    check_space(out, reserve)
    return result

