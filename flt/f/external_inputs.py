"""Read-only, stage-specific checks for manually acquired FLT inputs."""
import csv
import os
import re
from pathlib import Path


def check_inputs(stages):
    errors, notes = [], []
    if not {"phewas", "fig5"}.intersection(stages):
        return errors, ["Selected modules do not consume FUMA/Atlas/Metascape inputs."]
    fuma = Path(os.environ["FLT_FUMA_DIR"])

    def require(path):
        try:
            with path.open("rb") as stream:
                if not stream.read(1):
                    raise ValueError("empty file")
            return True
        except (OSError, ValueError) as exc:
            errors.append(f"{path}: {exc}")
            return False

    def table(path, columns, delimiter):
        if not require(path):
            return None
        try:
            with path.open(encoding="utf-8-sig", newline="") as stream:
                reader = csv.DictReader(stream, delimiter=delimiter)
                missing = set(columns) - set(reader.fieldnames or [])
                if missing:
                    raise ValueError("missing columns: " + ", ".join(sorted(missing)))
                return list(reader)
        except (OSError, ValueError, csv.Error) as exc:
            errors.append(f"{path}: {exc}")
            return None

    leads = table(fuma / "leadSNPs.txt", ["rsID", "GenomicLocus", "chr", "pos", "p", "uniqID"], "\t")
    if leads is not None:
        ids = [r["rsID"] for r in leads]
        if not ids or len(ids) != len(set(ids)) or any(not re.fullmatch(r"rs\d+", x or "") for x in ids):
            errors.append(f"{fuma / 'leadSNPs.txt'}: empty, duplicate or invalid rsIDs")
        elif "phewas" in stages:
            atlas = Path(os.environ["FLT_ATLAS_DIR"])
            columns = ["atlas ID", "PMID", "Year", "Domain", "Trait", "P-value", "N", "EA", "NEA"]
            n_found = 0
            for snp in ids:
                path = atlas / f"{snp}.csv"
                if path.exists() and table(path, columns, ",") is not None:
                    n_found += 1
            notes.append(f"Atlas: {n_found}/{len(ids)} lead-SNP exports found in {atlas}")
            notes.append("OpenGWAS/Catalog use existing analysis caches or online queries; online OpenGWAS requires a valid JWT.")
    if "fig5" in stages:
        for name in ("IndSigSNPs.txt", "genes.txt", "GenomicRiskLoci.txt", "magma.genes.out",
                     "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out",
                     "magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out"):
            require(fuma / name)
        sets = fuma / "magma.gsa.out"
        require(sets if sets.exists() else fuma / "magma.sets.out")
        meta = Path(os.environ["FLT_METASCAPE_RESULTS"])
        rows = table(meta, ["Description", "Log(q-value)", "LogP", "Enrichment", "#GeneInGOAndHitList"], ",") if meta.exists() else None
        if rows is not None and (not rows or not ({"GroupID", "GROUP_ID"} & rows[0].keys())):
            errors.append(f"{meta}: empty table or missing GroupID/GROUP_ID")
        notes.append(f"Fig5: checked FUMA/MAGMA files in {fuma} and Metascape result {meta}")
    return errors, notes
