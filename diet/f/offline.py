"""Validate local resources before expensive downstream analysis begins."""
from pathlib import Path
import re
from importlib.util import find_spec

GMT_PATTERNS = (
    r"GO_Biological_Process_[0-9]{4}", r"Reactome_[0-9]{4}",
    r"KEGG_[0-9]{4}_Human", r"WikiPathways?_[0-9]{4}_Human",
)
STRING_FILES = tuple(f"9606.protein.{kind}.v12.0.txt.gz" for kind in ("info", "aliases", "links"))


def fsaverage_dir(env):
    if env.get("DOSE_FSAVERAGE_DIR"):
        return Path(env["DOSE_FSAVERAGE_DIR"])
    local = Path(env["DOSE_FILES_DIR"]) / "fsaverage5"
    if local.is_dir():
        return local
    # Read data shipped with the installed package, without calling a fetch API.
    spec = find_spec("nilearn")
    if spec and spec.origin:
        bundled = Path(spec.origin).parent / "datasets/data/fsaverage5"
        if bundled.is_dir():
            return bundled
    return local


def check_resources(env, steps):
    files = Path(env["DOSE_FILES_DIR"])
    missing = []
    if "pathways" in steps:
        root = Path(env.get("DOSE_GMT_DIR", str(files / "pathways")))
        names = [p.stem for p in root.glob("*.gmt") if p.is_file() and p.stat().st_size]
        missing += [f"{root}: {pattern}.gmt" for pattern in GMT_PATTERNS
                    if not any(re.fullmatch(pattern, name) for name in names)]
    if "network" in steps:
        root = Path(env.get("DOSE_STRING_DIR", str(files / "STRING_network")))
        missing += [str(root / name) for name in STRING_FILES
                    if not (root / name).is_file() or not (root / name).stat().st_size]
    if "mri_plot" in steps:
        root = fsaverage_dir(env)
        missing += [str(root / f"{kind}_{side}.gii.gz")
                    for kind in ("pial", "infl", "white") for side in ("left", "right")
                    if not (root / f"{kind}_{side}.gii.gz").is_file()]
    if missing:
        raise RuntimeError("Missing pre-downloaded offline resources; no analysis started:\n  "
                           + "\n  ".join(missing)
                           + "\nPrepare these locally before rerunning. No downloads/installations are performed.")
