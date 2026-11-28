"""Download missing SMR resources explicitly during installation; retain download failures."""
import json
import hashlib
import os
import shutil
import subprocess
import tarfile
import tempfile
import zipfile
from pathlib import Path
from resources import complete, inspect, allow_skips

root = Path(os.environ["FLT_RESOURCE_DIR"])
root.mkdir(parents=True, exist_ok=True)
records = []


def fetch(url, name):
    archive = root / name
    if archive.is_file():
        return archive
    partial = root / (name + ".part")
    subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
                    "--connect-timeout", "20", "--speed-time", "60", "--speed-limit", "1024",
                    "--output", str(partial), url], check=True)
    partial.replace(archive)
    return archive


def unpack(archive, destination):
    # Reject traversal and links; extraction is confined to a temporary directory.
    def safe(name):
        return not Path(name).is_absolute() and ".." not in Path(name).parts
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as z:
            if not all(safe(i.filename) and (i.external_attr >> 16) & 0o170000 != 0o120000 for i in z.infolist()):
                raise ValueError("Unsafe ZIP member")
            z.extractall(destination)
    else:
        with tarfile.open(archive) as t:
            if not all(safe(i.name) and (i.isfile() or i.isdir()) for i in t.getmembers()):
                raise ValueError("Unsafe TAR member")
            t.extractall(destination)


def prepare(key, url, archive_name, prefix=None):
    if key == "smr":
        if shutil.which(os.environ["FLT_SMR"]):
            records.append(dict(resource=key, status="existing", path=os.environ["FLT_SMR"]))
            return
    elif complete(prefix, (".besd", ".esi", ".epi")):
        records.append(dict(resource=key, status="existing", path=prefix))
        return
    try:
        print(f"Preparing {key}: {url}", flush=True)
        archive = fetch(url, archive_name)
        if key == "Whole_Blood":
            # Published alongside the official archive as Whole_Blood.zip.md5sum.
            digest = hashlib.md5()
            with archive.open("rb") as f:
                for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest() != "503bdcb7e3c9aa11e2759d7918bc19e7":
                raise ValueError("GTEx archive checksum differs from official MD5")
        with tempfile.TemporaryDirectory(dir=root) as tmp:
            unpack(archive, tmp)
            if key == "smr":
                candidates = [p for p in Path(tmp).rglob("smr") if p.is_file()]
                if len(candidates) != 1:
                    raise ValueError("Expected one SMR executable")
                dest = root / "smr/smr"
                dest.parent.mkdir(parents=True, exist_ok=True)
                with candidates[0].open("rb") as f:
                    header = f.read(12)
                if header[8:11] == b"AI\x02":
                    appimage = dest.with_name("smr.AppImage")
                    shutil.copy2(candidates[0], appimage)
                    appimage.chmod(0o755)
                    shutil.copy2(Path(__file__).with_name("smr_launcher.sh"), dest)
                else:
                    shutil.copy2(candidates[0], dest)
                dest.chmod(0o755)
                # SMR --help exits 1; --version exits 0 and verifies the loader.
                subprocess.run([str(dest), "--version"], check=True, stdout=subprocess.DEVNULL, timeout=30)
                os.environ["FLT_SMR"] = str(dest)
            else:
                candidates = list(Path(tmp).rglob("*.besd"))
                if len(candidates) != 1:
                    raise ValueError("Expected one BESD dataset")
                source = str(candidates[0])[:-5]
                if not complete(source, (".besd", ".esi", ".epi")):
                    raise ValueError("Archive lacks matching BESD/ESI/EPI")
                Path(prefix).parent.mkdir(parents=True, exist_ok=True)
                for suffix in (".esi", ".epi", ".besd"):
                    staging = prefix + suffix + ".partial"
                    shutil.move(source + suffix, staging)
                    os.replace(staging, prefix + suffix)
        records.append(dict(resource=key, status="prepared", url=url))
    except (OSError, ValueError, subprocess.SubprocessError, tarfile.TarError, zipfile.BadZipFile) as exc:
        if key == "smr":
            (root / "smr/smr").unlink(missing_ok=True)
        records.append(dict(resource=key, status="unavailable", url=url, reason=str(exc)))
        print(f"[UNAVAILABLE] {key}: {exc}", flush=True)


prepare("smr", "https://yanglab.westlake.edu.cn/software/smr/download/smr-1.4.3-linux-x86_64.zip", "smr-1.4.3-linux-x86_64.zip")
# User requested check-and-skip after the E: search and failed official downloads.
eqtl_ready = complete(os.environ["FLT_EQTLGEN"], (".besd", ".esi", ".epi"))
records.append(dict(resource="eQTLGen", status="existing" if eqtl_ready else "skipped",
                    path=os.environ["FLT_EQTLGEN"], policy="check_only_no_download"))
if not eqtl_ready:
    print("[SKIP] eQTLGen: no local BESD dataset; check-only policy, no download retry.", flush=True)
prepare("Whole_Blood", "https://yanglab-obs.lab.westlake.edu.cn/files/data/SMR/GTEx_V8_cis_eqtl_summary/Whole_Blood.zip",
        "Whole_Blood.zip", os.environ["FLT_GTEX"])
(root / "downloads.json").write_text(json.dumps(records, indent=2) + "\n")
for key, reason in inspect().items():
    print(f"[UNAVAILABLE] {key}: {reason}")
if not allow_skips() and inspect():
    raise SystemExit("Missing resources; see downloads.json")
