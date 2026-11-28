"""Availability of optional genetic resources; no silent substitution of reference panels."""
import os
import shutil
from pathlib import Path


def complete(prefix, suffixes):
    return all(Path(str(prefix) + s).is_file() and Path(str(prefix) + s).stat().st_size > 0 for s in suffixes)


def inspect():
    issues = {}
    installed_smr = Path(os.environ["FLT_RESOURCE_DIR"]) / "smr/smr"
    if not shutil.which(os.environ["FLT_SMR"]) and os.access(installed_smr, os.X_OK):
        os.environ["FLT_SMR"] = str(installed_smr)
    if not shutil.which(os.environ["FLT_PLINK"]) or not complete(os.environ["FLT_MR_LD_PREFIX"], (".bed", ".bim", ".fam")):
        issues["mr_ld"] = "PLINK or complete EUR BED/BIM/FAM reference unavailable"
    if not shutil.which(os.environ["FLT_SMR"]):
        issues["smr"] = "SMR executable unavailable"
    for tissue, key in (("eQTLGen", "FLT_EQTLGEN"), ("Whole_Blood", "FLT_GTEX")):
        if not complete(os.environ[key], (".besd", ".esi", ".epi")):
            issues[tissue] = f"Complete BESD/ESI/EPI unavailable: {os.environ[key]}"
    hdl = Path(os.getenv("FLT_HDL_REF", str(Path(os.environ["FLT_LDREF"]) / "hdl/YGB_imputed_SVD_eigen99_extraction")))
    files = list(hdl.glob("*"))
    bims = [p for p in files if p.suffix == ".bim"]
    eigens = [p for p in files if "eigen" in p.name and p.suffix.lower() in {".rda", ".rdata"}]
    if not (bims and len(eigens) >= len(bims) and all(p.stat().st_size > 0 for p in bims + eigens)
            and any("_snp_counter" in p.name for p in files) and any("_snp_list" in p.name for p in files)):
        issues["hdl"] = f"HDL panel incomplete or absent: {hdl}"
    return issues


def allow_skips():
    return os.getenv("FLT_SKIP_MISSING_RESOURCES", "TRUE").upper() == "TRUE"
