"""Compile the bundled, pinned Kraken source without contacting the network."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

from hpc_common import software_sources, CODE_ROOT, atomic_json, locked, sha256, utcnow

REQUIRED = ("classify", "estimate_capacity", "build_db", "dump_table", "k2mask",
            "kraken2", "kraken2-build", "kraken2-inspect")


def install(cfg):
    for program in ("g++", "gcc", "make", "perl", "bash"):
        if not shutil.which(program):
            raise RuntimeError(f"Missing {program}. Load the site's compiler module in env.sh.")
    vendor = CODE_ROOT / "f" / "vendor"
    source = software_sources(CODE_ROOT)["kraken"]
    archive = vendor / source["archive"]
    if sha256(archive) != source["sha256"]:
        raise RuntimeError("Bundled Kraken source checksum mismatch")
    prefix = Path(cfg["tools_dir"]) / "kraken"
    if prefix.is_symlink():
        raise RuntimeError(f"Refusing symlink at tool installation: {prefix}")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with locked(prefix.parent / "kraken.install.lock"):
        owner = prefix / ".metag-hpc-owner.json"
        if prefix.exists() and not owner.is_file():
            raise RuntimeError(f"Unowned tool directory: {prefix}; choose a separate tools_dir")
        prefix.mkdir(exist_ok=True)
        if owner.exists():
            if json.loads(owner.read_text()).get("source_sha256") != source["sha256"]:
                raise RuntimeError("A different source owns this tool directory; use a new tools_dir")
        else:
            atomic_json(owner, {"schema": 1, "source_sha256": source["sha256"]})
        receipt = prefix / "install.json"
        if receipt.is_file():
            saved = json.loads(receipt.read_text())
            if (saved.get("source_commit") != source["source_commit"] or
                    saved.get("source_archive_sha256") != source["sha256"]):
                raise RuntimeError("Installed source differs from the bundled source")
            for name, checksum in saved["files"].items():
                if Path(name).name != name or sha256(prefix / "bin" / name) != checksum:
                    raise RuntimeError("Installed tool changed; refusing to overwrite a published installation")
            if not set(REQUIRED) <= set(saved["files"]):
                raise RuntimeError("Incomplete installation receipt")
            print(f"REUSE verified CPU tools: {prefix / 'bin'}", flush=True)
            return saved
        src = prefix / "source"
        if not src.exists():
            stage = prefix / "source.partial"
            if stage.exists():
                # The entire prefix is a private, marked installation owned above.
                shutil.rmtree(stage)
            stage.mkdir()
            with tarfile.open(archive, "r:gz") as bundle:
                for member in bundle:
                    p = Path(member.name)
                    if p.is_absolute() or ".." in p.parts or not p.parts or p.parts[0] != "kraken":
                        raise RuntimeError(f"Unsafe bundled archive member: {member.name}")
                    if not (member.isfile() or member.isdir()):
                        raise RuntimeError("Bundled archive must contain only files/directories")
                    relative = Path(*p.parts[1:])
                    target = stage / relative
                    if member.isdir():
                        target.mkdir(parents=True, exist_ok=True)
                    else:
                        target.parent.mkdir(parents=True, exist_ok=True)
                        with bundle.extractfile(member) as inp, target.open("wb") as out:
                            shutil.copyfileobj(inp, out)
                        target.chmod(member.mode & 0o777)
            stage.replace(src)
        output = prefix / "bin"
        output.mkdir(exist_ok=True)
        environment = dict(os.environ)
        environment.pop("KRAKEN2_SKIP_FOPENMP", None)
        environment["MAKEFLAGS"] = f"-j{int(cfg.get('setup_cores', 8))}"
        log = prefix / "install.log"
        print(f"Compile pinned Kraken CPU source; log={log}", flush=True)
        with log.open("ab") as handle:
            completed = subprocess.run(["bash", str(src / "install_kraken2.sh"), str(output)],
                                       cwd=src, env=environment, stdout=handle,
                                       stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f"Compiler exited {completed.returncode}; inspect {log}. "
                               "A C++ compiler, OpenMP and zlib development headers are required.")
        for name in REQUIRED:
            path = output / name
            if not path.is_file() or not os.access(path, os.X_OK):
                raise RuntimeError(f"Missing executable after compilation: {path}")
        reported = subprocess.check_output([str(output / "kraken2"), "--version"],
                                           text=True).splitlines()[0]
        subprocess.run([str(output / "classify"), "-h"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        result = {"schema": 1, "source_tag": source["source_tag"],
                  "source_commit": source["source_commit"],
                  "source_archive_sha256": source["sha256"],
                  "reported_version": reported, "installed_at": utcnow(),
                  "files": {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()},
                  "note": source["reported_version_note"]}
        atomic_json(receipt, result)
        print(f"PASS CPU installation: {reported}; source commit {source['source_commit']}", flush=True)
        return result
