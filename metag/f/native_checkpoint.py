"""Install private native checkpoint executable without changing official Kraken.

Call install(cfg) after software_install.install(cfg). No downloads, no shell commands.
The bundled pinned source archive is verified and extracted independently, so a
changed existing build/source tree cannot quietly alter checkpoint semantics.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

from common import configured_path
from common import CODE_ROOT, atomic_json, locked, sha256, software_sources, utcnow

HERE = Path(__file__).resolve().parent
SOURCE_COMMIT = "8c190b1b668825935dbf6dee5f969227dc8269bb"
FILES = ("metag_build.cc", "checked_writer.txt", "native_checkpoint.py")
COMPILE_UNITS = ("mmap_file.cc", "compact_hash.cc", "taxonomy.cc", "seqreader.cc",
                 "mmscanner.cc", "omp_hack.cc", "utilities.cc")


def _patch_header(path):
    content = path.read_text()
    begin = "template<typename Cell>\nvoid CompactHashTable<Cell>::WriteTable(const char *filename) {"
    end = "\n\ntemplate<typename Cell>\nhvalue_t CompactHashTable<Cell>::Get(hkey_t key) const {"
    if content.count(begin) != 1 or content.count(end) != 1:
        raise RuntimeError("Pinned compact_hash.h WriteTable signature changed")
    left, tail = content.split(begin, 1)
    _, right = tail.split(end, 1)
    content = left + (HERE / "checked_writer.txt").read_text().rstrip() + end + right
    original = ("void CompactHashTable<Cell>::LoadTable(const char *filename, bool memory_mapping) {\n"
                "  locks_initialized_ = false;")
    replacement = ("void CompactHashTable<Cell>::LoadTable(const char *filename, bool memory_mapping) {\n"
                   "  // Private metag builder resumes mutation; ordinary Kraken only queries loaded tables.\n"
                   "  locks_initialized_ = !memory_mapping;\n"
                   "  if (locks_initialized_) {\n"
                   "    for (size_t i = 0; i < LOCK_ZONES; ++i) omp_init_lock(&zone_locks_[i]);\n"
                   "  }")
    if content.count(original) != 1:
        raise RuntimeError("Pinned compact_hash.h LoadTable signature changed")
    path.write_text(content.replace(original, replacement))


def _extract(archive, destination):
    destination.mkdir()
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            p = Path(member.name)
            if p.is_absolute() or ".." in p.parts or not p.parts or p.parts[0] != "kraken":
                raise RuntimeError(f"Unsafe pinned archive path: {member.name}")
            if not (member.isfile() or member.isdir()):
                raise RuntimeError("Pinned archive contains unsupported entry")
            target = destination.joinpath(*p.parts[1:])
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with bundle.extractfile(member) as inp, target.open("wb") as out:
                    shutil.copyfileobj(inp, out)
                target.chmod(member.mode & 0o777)


def install(cfg):
    """Return receipt with executable, binary SHA256 and both source identities."""
    for command in ("g++", "sha256sum"):
        if not shutil.which(command):
            raise RuntimeError(f"Missing {command}; load a GCC >= 9 compiler module and GNU coreutils")
    source = software_sources(CODE_ROOT)["kraken"]
    if source["source_commit"] != SOURCE_COMMIT:
        raise RuntimeError("Checkpoint builder requires its reviewed pinned Kraken source")
    archive = Path(configured_path("DIR_SOFTWARE", "metag")) / source["archive"]
    if sha256(archive) != source["sha256"]:
        raise RuntimeError("Pinned Kraken source checksum mismatch")
    identities = {name: sha256(HERE / name) for name in FILES}
    extension_sha256 = hashlib.sha256(json.dumps(identities, sort_keys=True).encode()).hexdigest()
    identity = {"schema": 1, "source_commit": SOURCE_COMMIT,
                "source_archive_sha256": source["sha256"],
                "extension_sha256": extension_sha256, "extension_files": identities}
    prefix = Path(cfg["tools_dir"]) / "metag_checkpoint" / "releases" / extension_sha256[:16]
    if prefix.is_symlink():
        raise RuntimeError(f"Refusing symlink installation: {prefix}")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with locked(prefix.parent / "metag_checkpoint.install.lock"):
        owner = prefix / ".metag-native-owner.json"
        if prefix.exists() and not owner.is_file():
            raise RuntimeError(f"Unowned native installation: {prefix}")
        prefix.mkdir(exist_ok=True)
        if owner.exists():
            if json.loads(owner.read_text()) != identity:
                raise RuntimeError("Checkpoint executable source changed; use a new tools_dir; published checkpoints bind the old executable")
        else:
            atomic_json(owner, identity)
        receipt_path = prefix / "install.json"
        executable = prefix / "bin" / "metag_build"
        if receipt_path.exists():
            saved = json.loads(receipt_path.read_text())
            if (saved.get("identity") != identity or saved.get("executable") != str(executable) or
                    not executable.is_file() or sha256(executable) != saved.get("binary_sha256")):
                raise RuntimeError("Published native checkpoint installation changed")
            print(f"REUSE native checkpoint builder: {executable}", flush=True)
            return saved
        build = prefix / "build.partial"
        if build.exists():
            if build.is_symlink():
                raise RuntimeError("Native build directory must not be a symlink")
            shutil.rmtree(build)
        _extract(archive, build)
        src = build / "src"
        _patch_header(src / "compact_hash.h")
        shutil.copyfile(HERE / "metag_build.cc", src / "metag_build.cc")
        (prefix / "bin").mkdir(exist_ok=True)
        candidate = prefix / "bin" / "metag_build.partial"
        if candidate.is_symlink():
            raise RuntimeError("Native candidate must not be a symlink")
        command = ["g++", "-std=c++17", "-O3", "-fopenmp", "-Wall", "-Wextra", "-DLINEAR_PROBING",
                   "-I", str(src), str(src / "metag_build.cc")]
        command += [str(src / name) for name in COMPILE_UNITS]
        command += ["-o", str(candidate)]
        log = prefix / "install.log"
        print(f"Compile separate checkpoint builder; log={log}", flush=True)
        with log.open("ab") as out:
            result = subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Native checkpoint compilation failed; see {log}")
        reported = subprocess.check_output([str(candidate), "--metag-version"], text=True).strip()
        if not reported.startswith("metag-native-checkpoint-v1 "):
            raise RuntimeError("Native version check failed")
        candidate.replace(executable)
        result = {"schema": 1, "identity": identity, "source_commit": SOURCE_COMMIT,
                  "source_archive_sha256": source["sha256"], "extension_sha256": extension_sha256,
                  "executable": str(executable), "binary_sha256": sha256(executable),
                  "reported_version": reported, "installed_at": utcnow(),
                  "compile_command": command,
                  "patches": ["checked POSIX writes and fsync", "mutable-load OpenMP locks"],
                  "official_tools_unchanged": True}
        atomic_json(receipt_path, result)
        # Retain tiny private source tree as exact provenance; no original tree is changed.
        print(f"PASS native checkpoint builder: {executable}", flush=True)
        return result
