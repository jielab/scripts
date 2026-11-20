"""Small shared utilities for GRID v2. No cohort or genotype IO at import time."""
from __future__ import annotations
import contextlib
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import pandas as pd

VERSION = "grid-final-20261116"
POPS = ("EUR", "AFR", "EAS", "SAS")

@contextlib.contextmanager
def text_open(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as stream: yield stream
    elif str(path).endswith(".zst"):
        import zstandard
        with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as dec:
            with io.TextIOWrapper(dec) as stream: yield stream
    else:
        with path.open() as stream: yield stream

def read_table(path, columns=None):
    """Keep IDs as strings, preserve empty TSV cells, optionally read RDS columns."""
    path = Path(path)
    if path.suffix.lower() == ".rds":
        with tempfile.TemporaryDirectory(prefix="grid-columns-") as folder:
            out = Path(folder) / "columns.tsv"
            script = Path(folder) / "read.R"
            script.write_text('''a <- commandArgs(TRUE)
x <- readRDS(a[1]); if (is.list(x) && !is.data.frame(x) && length(x)==1) x <- x[[1]]
if (!is.data.frame(x)) stop("RDS must contain a data.frame")
if (length(a)>2) x <- x[, intersect(names(x), a[-c(1,2)]), drop=FALSE]
options(digits=17)
write.table(x, a[2], sep="\\t", row.names=FALSE, quote=TRUE, na="NA")
''')
            subprocess.run([os.environ.get("GRID_RSCRIPT", "Rscript"), "--vanilla", str(script), str(path), str(out), *(columns or [])], check=True)
            return pd.read_csv(out, sep="\t", dtype=str)
    skip = 0
    with text_open(path) as stream:
        for line in stream:
            if line.startswith("##") or not line.strip(): skip += 1; continue
            sep = "\t" if "\t" in line else "," if "," in line else r"\s+"
            break
        else: raise ValueError(f"Empty table: {path}")
    with text_open(path) as stream:
        return pd.read_csv(stream, sep=sep, skiprows=skip, dtype=str,
                           usecols=(lambda c: c in columns) if columns else None)

def identified(frame, label="table"):
    frame = frame.copy()
    key = next((c for c in ("id", "eid", "IID", "#IID", "ID_2") if c in frame), None)
    if key is None: raise ValueError(f"{label} needs id/eid/IID")
    frame["id"] = frame[key].astype("string").str.strip().str.replace(r"^([0-9]+)\.0+$", r"\1", regex=True)
    if frame.id.isna().any() or frame.id.isin(["", "NA", "nan"]).any() or frame.id.duplicated().any():
        raise ValueError(f"{label} has missing or duplicate IDs")
    frame["id"] = frame.id.astype(str)
    return frame

def aligned(frame, ids, label="table", allow_missing=False):
    frame = identified(frame, label).set_index("id")
    missing = pd.Index(ids).difference(frame.index)
    if len(missing) and not allow_missing: raise ValueError(f"{label}: {len(missing)} required IDs missing")
    return frame.reindex(ids).rename_axis("id").reset_index()

def atomic_table(frame, path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.partial")
    frame.to_csv(temporary, sep="\t", index=False, na_rep="NA", compression="gzip" if str(path).endswith(".gz") else None)
    temporary.replace(path)

def json_default(obj):
    if isinstance(obj, Path): return str(obj)
    if isinstance(obj, np.ndarray): return obj.tolist()
    if isinstance(obj, np.generic): return obj.item()
    raise TypeError(type(obj).__name__)

def atomic_json(value, path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.partial")
    temporary.write_text(json.dumps(value, indent=2, default=json_default, ensure_ascii=False) + "\n")
    temporary.replace(path)

def stamp(path):
    path = Path(path).resolve(); stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}

def group_fraction(group, seed, salt="split"):
    # Hashing identifiers is deterministic allocation, not a file-integrity workflow.
    keys = pd.Series(group, dtype=str)
    values = {g: int.from_bytes(hashlib.blake2b(f"{seed}:{salt}:{g}".encode(), digest_size=8).digest(), "big") / 2**64 for g in keys.unique()}
    return keys.map(values).to_numpy(float)

def group_folds(group, n_splits=5, seed=20261115):
    groups = pd.Series(group, dtype=str)
    if groups.nunique() < n_splits: raise ValueError("Too few family groups for OOF folds")
    unique = np.array(sorted(groups.unique()), dtype=object)
    rng = np.random.default_rng(seed); rng.shuffle(unique)
    counts = groups.value_counts(); loads = np.zeros(n_splits, dtype=int); mapping = {}
    # Greedy, group-wise balance without inspecting phenotypes.
    for g in unique:
        k = int(np.argmin(loads)); mapping[g] = k; loads[k] += int(counts[g])
    return groups.map(mapping).to_numpy(int)

def check_roles(cohort):
    required = {"id", "group_id", "role", "ancestry"}
    if not required <= set(cohort): raise ValueError(f"Cohort missing {required-set(cohort)}")
    if cohort.group_id.isna().any() or cohort.groupby("group_id").role.nunique().max() != 1:
        raise ValueError("Family group crosses roles or is missing")
    if not cohort.role.isin(["train", "tune", "calibration", "test"]).all():
        raise ValueError("Roles must be train/tune/calibration/test")
    for role in ("train", "tune", "calibration", "test"):
        if not cohort.role.eq(role).any(): raise ValueError(f"Empty {role} role")

def expand(value):
    if isinstance(value, dict): return {k: expand(v) for k, v in value.items()}
    if isinstance(value, list): return [expand(v) for v in value]
    if isinstance(value, str): return os.path.expanduser(os.path.expandvars(value))
    return value

def load_config(path):
    import yaml
    with open(path) as stream: value = yaml.safe_load(stream)
    if not isinstance(value, dict): raise ValueError("Config must be a YAML/JSON mapping")
    return expand(value)

class NumericDesign:
    """Training-only median imputation and scaling. No per-test-set normalization."""
    def __init__(self, columns): self.columns = list(columns)
    def fit(self, frame):
        x = frame[self.columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        if not np.isfinite(x).any(axis=0).all(): raise ValueError("An input feature is entirely missing in training")
        x[~np.isfinite(x)] = np.nan
        self.median = np.nanmedian(x, axis=0)
        x = np.where(np.isfinite(x), x, self.median)
        self.mean = x.mean(axis=0); self.scale = x.std(axis=0)
        self.scale[self.scale < 1e-12] = 1
        return self
    def transform(self, frame):
        x = frame[self.columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        return (np.where(np.isfinite(x), x, self.median) - self.mean) / self.scale
    def fit_transform(self, frame): return self.fit(frame).transform(frame)

