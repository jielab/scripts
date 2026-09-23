"""Download only the official TabICLv2 classification checkpoint; never store tokens."""
import hashlib
import json
import os
from pathlib import Path
import subprocess

REPO = "jingang/TabICL"
REVISION = "4dcd344ece2c00be9e831fdd35bed57b5ad83e19"
FILENAME = "tabicl-classifier-v2-20260212.ckpt"
SHA256 = "bdc7dbd5e4ff21f8f0456fcf90c6b7cdf72dbea960f2d05b19bec19f9b3d4ed0"
DEFAULT_MODEL = "/mnt/e/AI_model/TabICL/" + FILENAME


def verify_model(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size < 1_000_000:
        raise FileNotFoundError(f"Missing model weights (a Git LFS pointer is insufficient): {path}; use --download-model")
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != SHA256:
        raise ValueError(f"Model checksum mismatch: {path}; expected the official pinned TabICLv2 classifier")
    return dict(repo=REPO, revision=REVISION, file=str(path), sha256=digest, bytes=path.stat().st_size)


def download_model(path=DEFAULT_MODEL):
    path = Path(path)
    if path.name != FILENAME:
        raise ValueError(f"Use a model path ending in {FILENAME}")
    try:
        return verify_model(path)
    except (FileNotFoundError, ValueError):
        pass
    root = path.parent
    root.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, GIT_LFS_SKIP_SMUDGE="1", GIT_TERMINAL_PROMPT="0")
    try:
        if not root.exists():
            subprocess.run(["git", "clone", "--depth", "1", "https://huggingface.co/"+REPO, str(root)],
                           env=env, check=True, timeout=60)
        subprocess.run(["git", "-C", str(root), "fetch", "origin", REVISION], env=env, check=True, timeout=60)
        subprocess.run(["git", "-C", str(root), "lfs", "pull", "--include="+FILENAME,
                        "--exclude=", "origin", REVISION], env=env, check=True, timeout=300)
        info = verify_model(path)
    except (OSError, subprocess.SubprocessError, ValueError):
        # Same snapshot_download approach as /mnt/d/scripts/0f/0hf_download.py,
        # parameterized for TabICL. Authentication, if needed, stays in HF_TOKEN.
        from huggingface_hub import snapshot_download
        snapshot_download(repo_id=REPO, revision=REVISION, local_dir=str(root),
                          allow_patterns=[FILENAME, "README.md"], token=os.environ.get("HF_TOKEN", False))
        info = verify_model(path)
    (root/"panome_download.json").write_text(json.dumps(info, indent=2)+"\n")
    return info
