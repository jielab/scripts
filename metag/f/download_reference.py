"""Resumable HTTP range download for a large reference, with bounded connections."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
from pathlib import Path
import shutil
import time
from urllib.request import Request, urlopen


def download(url, output, connections=8):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with Path(str(output) + ".download.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if output.is_file() and output.stat().st_size:
            return output
        with urlopen(Request(url, method="HEAD", headers={"Accept-Encoding": "identity"}), timeout=60) as response:
            total = int(response.headers["Content-Length"])
            version = response.headers.get("ETag") or response.headers.get("Last-Modified")
        if total <= 0:
            raise RuntimeError("Empty remote reference")
        parts = Path(str(output) + ".parts"); parts.mkdir(exist_ok=True)
        identity = dict(url=url, bytes=total, version=version, connections=connections)
        state = parts / "source.json"
        if state.exists() and json.loads(state.read_text()) != identity:
            raise RuntimeError(f"Remote source or connection count changed; choose a new output path: {state}")
        state.write_text(json.dumps(identity, indent=2))
        stride = (total + connections - 1) // connections
        ranges = [(i, i*stride, min(total, (i+1)*stride)-1) for i in range(connections) if i*stride < total]

        def part(item):
            i, start, end = item
            path = parts / f"{i:03d}.part"
            for attempt in range(5):
                have = path.stat().st_size if path.exists() else 0
                if have == end-start+1:
                    return path
                if have > end-start+1:
                    raise RuntimeError(f"Oversized download part: {path}")
                headers = {"Range": f"bytes={start+have}-{end}", "Accept-Encoding": "identity"}
                if version:
                    headers["If-Range"] = version
                try:
                    with urlopen(Request(url, headers=headers), timeout=90) as response:
                        expected = f"bytes {start+have}-{end}/{total}"
                        if response.status != 206 or response.headers.get("Content-Range") != expected:
                            raise ValueError("Server did not honor the exact byte range or the source changed")
                        with path.open("ab") as out:
                            remaining = end-start+1-have
                            while remaining:
                                block = response.read(min(1 << 20, remaining))
                                if not block:
                                    raise OSError("Truncated HTTP response")
                                out.write(block); remaining -= len(block)
                    print(f"Downloaded part {i+1}/{len(ranges)}", flush=True)
                    return path
                except (OSError, TimeoutError):
                    if attempt == 4:
                        raise
                    time.sleep(5)
            raise RuntimeError("Unreachable download state")

        with ThreadPoolExecutor(max_workers=connections) as pool:
            paths = list(pool.map(part, ranges))
        temporary = Path(str(output) + ".assembled")
        with temporary.open("wb") as out:
            for path in paths:
                with path.open("rb") as inp:
                    shutil.copyfileobj(inp, out, 8 << 20)
        if temporary.stat().st_size != total:
            raise RuntimeError("Assembled reference size differs from server")
        temporary.replace(output)
        # Keep the tiny source manifest; reclaim only this downloader's parts.
        for path in paths:
            path.unlink()
        print(f"Downloaded {total:,} bytes: {output}", flush=True)
        return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--url", required=True); p.add_argument("--output", required=True)
    p.add_argument("--connections", type=int, default=8)
    a = p.parse_args()
    if not 1 <= a.connections <= 16:
        p.error("--connections must be in 1..16")
    download(a.url, a.output, a.connections)
