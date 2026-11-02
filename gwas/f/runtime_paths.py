"""Defaults for direct Python helper use; shell entry points export these paths."""
import os
from pathlib import Path
import shutil

mode = os.environ.get("submit", os.environ.get("SUBMIT", "auto"))
if mode == "auto":
    mode = "hpc" if os.environ.get("LSB_JOBID") or shutil.which("bsub") else "local"
WORK = Path(os.environ.get("DIR_WORK", "/work/sph-huangj" if mode == "hpc" else "/mnt/d"))
DATA = Path(os.environ.get("GWAS_DATA_ROOT", os.environ.get("DIR_DATA", "/data/sph-huangj" if mode == "hpc" else "/mnt/f")))
FILES = Path(os.environ.get("GWAS_FILES_ROOT", str((DATA if mode == "hpc" else WORK) / "files")))
SOFTWARE = Path(os.environ.get("GWAS_SOFTWARE_ROOT", os.environ.get("DIR_SOFTWARE", str(WORK / "software"))))
scratch = Path(os.environ.get("DIR_SCRATCH", f"/scratch/{os.environ.get('DATE_SCRATCH', '2026-10-31')}/sph-huangj" if mode == "hpc" else str(WORK)))
ANALYSIS = Path(os.environ.get("GWAS_ANALYSIS_ROOT", str(scratch / "analysis/gwas")))
