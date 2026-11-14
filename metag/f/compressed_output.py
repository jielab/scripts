"""Bounded parallel gzip output without one process or directory per taxon."""
from collections import OrderedDict, deque
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import gzip
import io
from pathlib import Path
import shutil
import subprocess
import tempfile


class PigzStream:
    """One persistent pigz process for a large stream; errors are never ignored."""
    def __init__(self, path, level, threads):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.output = self.path.open("wb", buffering=1 << 20)
        self.errors = tempfile.TemporaryFile()
        self.process = subprocess.Popen(
            [shutil.which("pigz"), "-n", f"-{level}", "-p", str(threads), "-c"],
            stdin=subprocess.PIPE, stdout=self.output, stderr=self.errors, bufsize=1 << 20)
        self.closed = False

    def write(self, data):
        return self.process.stdin.write(data)

    def close(self):
        if self.closed:
            return
        self.closed = True
        try:
            self.process.stdin.close()
            code = self.process.wait()
            if code:
                self.errors.seek(0)
                detail = self.errors.read(4096).decode("utf-8", "replace")
                raise RuntimeError(f"pigz exited {code} for {self.path}: {detail}")
        finally:
            if self.process.poll() is None:
                self.process.terminate()
                self.process.wait()
            self.errors.close()
            self.output.close()


@contextmanager
def gzip_text(path, *, compresslevel=3, threads=2):
    """Compress one text stream with pigz when installed, otherwise stdlib gzip."""
    if threads > 0 and shutil.which("pigz"):
        stream = PigzStream(path, compresslevel, threads)
        text = io.TextIOWrapper(stream.process.stdin, encoding="utf-8", newline="")
        try:
            yield text
        finally:
            try:
                text.flush()
                text.detach()
            finally:
                stream.close()
    else:
        with gzip.open(path, "wt", compresslevel=compresslevel, encoding="utf-8", newline="") as text:
            yield text


class BufferedWriters:
    """Multiplex records into bounded buffers and parallel gzip members.

    Each submitted chunk is compressed in C (zlib releases the Python GIL).
    Results are drained in submission order, retaining order within every bin.
    At most ``2 * threads`` chunks and ``buffer_bytes`` pending bytes are kept.
    A few known-large bins may use persistent pigz processes; tiny bins never
    spawn a compressor process. Files are valid concatenated-member gzip.
    """
    def __init__(self, root, *, max_open=32, buffer_bytes=64 << 20,
                 file_buffer_bytes=1 << 20, compresslevel=3, threads=4,
                 hot_keys=(), pigz_threads=2):
        if min(max_open, buffer_bytes, file_buffer_bytes, threads) < 1:
            raise ValueError("Writer limits and threads must be positive")
        if isinstance(compresslevel, bool) or not 1 <= compresslevel <= 9:
            raise ValueError("split_gzip_level must be an integer in 1..9")
        self.root = Path(root)
        self.max_open = int(max_open)
        self.buffer_bytes = int(buffer_bytes)
        self.file_buffer_bytes = min(int(file_buffer_bytes), self.buffer_bytes)
        self.compresslevel = int(compresslevel)
        self.buffers, self.handles = OrderedDict(), OrderedDict()
        self.pending_bytes = 0
        self.paths = set()
        self.executor = ThreadPoolExecutor(max_workers=int(threads), thread_name_prefix="metag-gzip")
        self.inflight = deque()
        self.max_pending = max(1, 2 * int(threads))
        self.hot_keys = set(hot_keys) if shutil.which("pigz") else set()
        self.pigz_threads = max(1, int(pigz_threads))
        self.hot = {}
        self.chunks = 0
        self.closed = False

    def _handle(self, key):
        if key in self.handles:
            handle = self.handles.pop(key)
            self.handles[key] = handle
            return handle
        while len(self.handles) >= self.max_open:
            _, handle = self.handles.popitem(last=False)
            handle.close()
        path = self.root / key
        path.parent.mkdir(parents=True, exist_ok=True)
        handle = path.open("ab", buffering=256 << 10)
        self.handles[key] = handle
        return handle

    def _drain(self):
        key, future = self.inflight.popleft()
        self._handle(key).write(future.result())

    def _flush(self, key):
        pending = self.buffers.pop(key)
        self.pending_bytes -= len(pending)
        if key in self.hot_keys:
            if key not in self.hot:
                self.hot[key] = PigzStream(self.root / key, self.compresslevel, self.pigz_threads)
            self.hot[key].write(pending)
        else:
            self.inflight.append((key, self.executor.submit(
                gzip.compress, bytes(pending), compresslevel=self.compresslevel, mtime=0)))
            if len(self.inflight) >= self.max_pending:
                self._drain()
        self.chunks += 1

    def write(self, key, record):
        if key not in self.paths:
            relative = Path(key)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Bin path must stay within the output directory")
        if self.closed:
            raise ValueError("Writer is closed")
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
        if self.closed:
            return
        self.closed = True
        try:
            if flush:
                while self.buffers:
                    self._flush(next(iter(self.buffers)))
                while self.inflight:
                    self._drain()
            else:
                for _, future in self.inflight:
                    future.cancel()
                self.buffers.clear()
                self.pending_bytes = 0
        finally:
            self.executor.shutdown(wait=True, cancel_futures=not flush)
            errors = []
            for handle in [*self.hot.values(), *self.handles.values()]:
                try:
                    handle.close()
                except BaseException as exc:
                    errors.append(exc)
            self.hot.clear()
            self.handles.clear()
            if errors:
                raise errors[0]
