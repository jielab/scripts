#!/usr/bin/env python3
"""Bound analysis PNGs to 1 MB; preserve plots, aspect ratio and transparency.

Existing small PNGs are untouched. Oversized files are losslessly recompressed,
then quantized to 256 colours, then resized only if necessary. Replacements are
verified and atomic. The CLI can back up originals and write an audit report.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
from pathlib import Path
import stat
import tempfile
import zipfile

MAX_BYTES = 1_000_000
PNG_DPI = 220


def limit_png(path, max_bytes=MAX_BYTES):
    from PIL import Image, PngImagePlugin, features

    path = Path(path)
    if not 128 <= max_bytes <= MAX_BYTES:
        raise ValueError("PNG byte limit must be between 128 and 1,000,000")
    before = path.stat()
    record = dict(path=str(path), before_bytes=before.st_size,
                  after_bytes=before.st_size, method="unchanged")
    if before.st_size <= max_bytes:
        return record
    original = path.read_bytes()
    record["before_sha256"] = hashlib.sha256(original).hexdigest()
    with Image.open(io.BytesIO(original)) as opened:
        if opened.format != "PNG" or getattr(opened, "n_frames", 1) != 1:
            raise ValueError(f"Expected a static PNG: {path}")
        opened.load()
        info = opened.info.copy()
        source = opened.copy()
    record["before_size"] = list(source.size)
    metadata = PngImagePlugin.PngInfo()
    for key, value in info.items():
        if isinstance(value, str):
            metadata.add_itxt(key, value)

    def encode(image):
        buffer = io.BytesIO()
        kwargs = dict(optimize=True, compress_level=9, pnginfo=metadata)
        for key in ("icc_profile", "exif"):
            if key in info:
                kwargs[key] = info[key]
        if "dpi" in info:
            kwargs["dpi"] = tuple(info["dpi"][i] * image.size[i] / source.size[i] for i in (0, 1))
        image.save(buffer, format="PNG", **kwargs)
        return buffer.getvalue()

    # Explicit alpha also handles palette PNGs carrying a tRNS chunk.
    if source.mode in ("RGB", "RGBA", "P", "L", "LA"):
        source = source.convert("RGBA")
        if source.getchannel("A").getextrema() == (255, 255):
            source = source.convert("RGB")
    else:
        # Do not silently reduce scientific 16-bit grayscale to 8-bit colour.
        raise ValueError(f"Unsupported scientific PNG mode {source.mode}: {path}")

    candidate = encode(source)
    method = "lossless"
    working = source
    while len(candidate) > max_bytes:
        quantizer = (Image.Quantize.MEDIANCUT if working.mode == "RGB" else
                     Image.Quantize.LIBIMAGEQUANT if features.check_feature("libimagequant") else
                     Image.Quantize.FASTOCTREE)
        with working.quantize(colors=256, method=quantizer, dither=Image.Dither.NONE) as palette:
            candidate = encode(palette)
        method = "palette256" if working.size == source.size else "resize+palette256"
        if len(candidate) <= max_bytes:
            break
        if working.size == (1, 1):
            raise ValueError(f"PNG metadata alone exceeds byte limit: {path}")
        scale = min(0.9, math.sqrt(max_bytes / len(candidate)) * 0.95)
        size = tuple(max(1, int(n * scale)) for n in working.size)
        if working is not source:
            working.close()
        # Always resample from the original, avoiding cumulative blur.
        working = source.resize(size, Image.Resampling.LANCZOS)

    with Image.open(io.BytesIO(candidate)) as check:
        check.verify()
    record.update(after_bytes=len(candidate), after_size=list(working.size), method=method,
                  after_sha256=hashlib.sha256(candidate).hexdigest())
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".png-limit-", suffix=".tmp", delete=False) as out:
        temporary = Path(out.name)
        out.write(candidate)
        out.flush()
        os.fsync(out.fileno())
    try:
        # Abort if an analysis process replaced the input while it was read.
        current = path.stat()
        if (current.st_size, current.st_mtime_ns) != (before.st_size, before.st_mtime_ns):
            raise RuntimeError(f"PNG changed during optimization: {path}")
        os.chmod(temporary, stat.S_IMODE(before.st_mode))
        os.utime(temporary, ns=(before.st_atime_ns, before.st_mtime_ns))
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
        if working is not source:
            working.close()
        source.close()
    return record


def save_png_figure(figure, filename, **kwargs):
    """Matplotlib savefig with bounded PNG output; vector output is unchanged."""
    is_png = str(kwargs.get("format", Path(filename).suffix.lstrip("."))).lower() == "png"
    if is_png:
        dpi = kwargs.get("dpi", figure.dpi)
        kwargs["dpi"] = min(float(figure.dpi if dpi == "figure" else dpi), PNG_DPI)
    figure.savefig(filename, **kwargs)
    if is_png:
        limit_png(filename)


def save_png_image(image, filename, **kwargs):
    image.save(filename, **kwargs)
    limit_png(filename)


def png_paths(paths):
    seen = set()
    for name in paths:
        root = Path(name)
        if not root.exists():
            raise FileNotFoundError(root)
        if root.is_file():
            candidates = [root]
        else:
            def fail(error):
                raise error
            candidates = (Path(d) / f for d, _, fs in os.walk(root, onerror=fail)
                          for f in sorted(fs) if f.lower().endswith(".png"))
        for path in candidates:
            if path.suffix.lower() == ".png" and path.resolve() not in seen:
                seen.add(path.resolve())
                yield path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--check", action="store_true", help="Report only; fail if any PNG exceeds 1 MB")
    parser.add_argument("--backup", type=Path, help="Create a ZIP of oversized originals before replacing them")
    parser.add_argument("--report", type=Path, help="JSON audit path")
    args = parser.parse_args()
    files = list(png_paths(args.paths))
    large = [p for p in files if p.stat().st_size > MAX_BYTES]
    # Complete and verify the recovery archive before changing any image.
    if args.backup and not args.check:
        args.backup.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(args.backup, "x", compression=zipfile.ZIP_STORED) as archive:
            for path in large:
                archive.write(path, path.resolve().as_posix().lstrip("/"))
        with zipfile.ZipFile(args.backup) as archive:
            if archive.testzip() is not None:
                raise RuntimeError("PNG backup verification failed")
    records = []
    errors = []
    for i, path in enumerate(large, 1):
        try:
            result = (dict(path=str(path), before_bytes=path.stat().st_size, method="oversized")
                      if args.check else limit_png(path))
            records.append(result)
            print(json.dumps(dict(index=i, total=len(large), **result)), flush=True)
        except Exception as error:
            errors.append(dict(path=str(path), error=str(error)))
            print(json.dumps(errors[-1]), flush=True)
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(json.dumps(dict(scanned=len(files), oversized=len(large),
                                                  files=records, errors=errors), indent=2) + "\n")
    summary = dict(scanned=len(files), oversized=len(large), processed=len(records), errors=errors)
    print(json.dumps(summary), flush=True)
    if args.report and not large:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(dict(**summary, files=[]), indent=2) + "\n")
    return 1 if errors or (args.check and large) else 0


if __name__ == "__main__":
    raise SystemExit(main())
