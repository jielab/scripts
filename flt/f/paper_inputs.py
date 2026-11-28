"""Validate input paths and bind archived FUMA outputs to the focal PLACO run."""
import csv
import gzip
import math
import os
from pathlib import Path

FUMA_MISSING = "FUMA lead table is absent"
FOCAL = "G732-mdd_2025_PLACO"


def placo_file(env, suffix):
    folder = Path(env["FLT_PLACO_DIR"]) if env.get("FLT_PLACO_DIR") else Path(env["FLT_OUT"]) / "03placo"
    return folder / (env.get("FLT_PLACO_PREFIX", FOCAL) + suffix)
REQUIRES = {
    "fig1": ["01cox/cox_results.rds", "01cox/demographic_characteristics.csv"],
    "fig2": ["02uvmr/{panel}/uvmr_main_table.csv"],
    "fig3": ["02uvmr/{panel}/uvmr_main_table.csv", "01cox/cox_results_table.csv",
             "01intake/spearman_selected_actual_liking.csv", "01intake/radar_project_intake.csv",
             "01lm/lm_results.csv"],
    "fig4": ["01lm/mri/mri_lm_results.csv", "01protein/limma/food_pair_mapping.csv"],
    "fig5": [],
    "genetics": [],
    "tables": [],
}


def available(path):
    return path.is_file() and path.stat().st_size > 0


def missing_inputs(stage, env):
    root = Path(env["FLT_OUT"])
    panel = env["FLT_PANEL"]
    required = [root / item.format(panel=panel) for item in REQUIRES[stage]]
    if stage == "fig5":
        required.append(placo_file(env, "_all.tsv.gz"))
    return [str(p) for p in required if not available(p)]


def review_fuma(env):
    lead = Path(env["FLT_FUMA_DIR"]) / "leadSNPs.txt"
    hits = placo_file(env, "_hits.tsv")
    if not available(lead):
        return FUMA_MISSING
    if not available(hits):
        return "PLACO hit table is absent"
    try:
        with lead.open(encoding="utf-8-sig", newline="") as src:
            leads = list(csv.DictReader(src, delimiter="\t"))
        with hits.open(encoding="utf-8-sig", newline="") as src:
            local = {row["SNP"]: row for row in csv.DictReader(src, delimiter="\t")}
        if not leads:
            return "FUMA has no lead SNPs"
        if len(leads) != len({row["rsID"] for row in leads}):
            return "FUMA lead SNP IDs are duplicated"
        mismatches = []
        for row in leads:
            snp = row["rsID"]
            hit = local.get(snp)
            if hit is None or hit["CHR1"] != row["chr"] or not math.isclose(
                    float(hit["P_PLACO"]), float(row["p"]), rel_tol=0.005):
                mismatches.append(snp)
        if mismatches:
            return f"{len(mismatches)} FUMA leads differ from native PLACO in ID, chromosome or P"
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        return f"Cannot match FUMA lead table: {exc}"
    return ""


def check_graph(jobs, env):
    from external_inputs import check_inputs
    from resources import inspect
    from workflow import records

    issues = set()
    notes = set()
    stages = {job["stage"] for job in jobs.values()}
    root = Path(env["FLT_GWAS"])
    ygb = Path(env["FLT_YGB"])

    def require(path):
        path = Path(path)
        if not available(path):
            issues.add(str(path))

    for stage in stages & REQUIRES.keys():
        issues.update(missing_inputs(stage, env))
    if stages & {"prepare", "cox", "cox-sensitivity", "uvmr-instruments", "ldsc-input", "fig1", "fig2", "fig3"}:
        require(env["FLT_DICTIONARY"])
    if "prepare" in stages:
        for name in ("all", "icd10.code0", "icd10.date0"):
            require(ygb / "Rdata" / f"{name}.rds")
        require(ygb / "common/icd10_brain.lst")
        for name in ("flt", "cmd_anthropometry_lifestyle", "behavioral"):
            require(ygb / "common" / f"{name}.lst")
            if not available(ygb / "rap" / f"{name}.tab.gz"):
                require(os.getenv("FLT_RAW", ygb / "rap/raw/pheno.tab.gz"))
    if stages & {"intake", "protein"}:
        require(ygb / "Rdata/diet.fudan.rds")
    if "protein" in stages:
        require(ygb / "Rdata/prot.rds")
    if "mri" in stages:
        require(os.getenv("FLT_MRI_RAW", ygb / "rap/img.tab.gz"))
        require(env["FLT_MRI_DICTIONARY"])
    gwas = set()
    if stages & {"ldsc-input", "uvmr-instruments"}:
        for row in records("panels.tsv"):
            if row["panel"] != "extended" or env["FLT_PANEL"] == "extended":
                gwas.add(("brain", row["id"]))
        foods = sorted((root / "flt/common").glob("*/gwas/*.gz"))
        if not foods:
            issues.add(str(root / "flt/common/ID/gwas/ID.gz"))
        else:
            notes.add(f"{len(foods)} FLT files found; R matches the 182 exact dictionary IDs")
            gwas.add(("flt", "G732"))
    if stages & {"placo", "coloc", "smr", "mvmr"}:
        gwas.update({("flt", "G732"), ("brain", "mdd_2025")})
    if "mvmr" in stages:
        gwas.update({("covs", "income"), ("covs", "physical_activity"), ("brain", "edu")})
    for group, name in sorted(gwas):
        path = root / group / "common" / name / "gwas" / f"{name}.gz"
        require(path)
        if available(path):
            try:
                with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
                    header = set(handle.readline().split())
                missing = set("SNP CHR POS EA NEA EAF N BETA SE P".split()) - header
                if missing:
                    issues.add(f"{path}: missing columns {', '.join(sorted(missing))}")
            except (OSError, UnicodeError) as exc:
                issues.add(f"{path}: {exc}")
    if "ldsc-input" in stages:
        panel = Path(env.get("FLT_LDSC_REF", Path(env["FLT_LDREF"]) / "ldsc/hm3"))
        require(panel / "w_hm3.snplist")
        for chrom in range(1, 23):
            for suffix in (".l2.ldscore.gz", ".l2.M_5_50"):
                require(panel / f"{chrom}{suffix}")
    if stages & {"uvmr", "smr", "ldsc", "mvmr"}:
        notes.update(inspect().values())
    if stages & {"phewas", "fig5"}:
        if reason := review_fuma(env):
            if reason == "PLACO hit table is absent" and "placo" in stages:
                notes.add("PLACO hit table will be generated by the prerequisite PLACO stage")
            else:
                notes.add(reason)
        else:
            external, _ = check_inputs(sorted(stages & {"phewas", "fig5"}))
            notes.update(external)
    for message in sorted(notes):
        print(f"[INPUT NOTE] {message}")
    for message in sorted(issues):
        print(f"[INPUT MISSING] {message}")
    if not issues:
        print("[INPUTS OK] Selected starting paths are available.")
    return int(bool(issues))
