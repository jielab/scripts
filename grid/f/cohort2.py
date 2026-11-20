"""Phenotypes, projected PCs and family-wise roles, without a genotype matrix cache."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from core2 import (VERSION, read_table, identified, aligned, atomic_table, atomic_json,
                   check_roles, group_fraction, stamp)

def family_components(ids, edges):
    """Union all supplied relationship edges, including links via non-retained IDs."""
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges[["ID1", "ID2"]].astype(str).itertuples(index=False, name=None):
        a, b = find(a), find(b)
        if a != b: parent[max(a,b)] = min(a,b)
    return [find(str(x)) for x in ids]

def choose_roles(frame, cfg, old_split=None):
    """Preserve provided outer test IDs; remove development relatives of that test."""
    d = frame.copy(); removed = 0; seed = int(cfg.get("seed", 20261115))
    if old_split is not None:
        # A frozen roster defines the eligible population as well as TEST.
        # Never assign newly encountered IDs to development implicitly.
        old_split = identified(old_split, "frozen split")
        d = d.loc[d.id.isin(old_split.id)].copy()
        s = aligned(old_split, d.id, "frozen split")
        key = "role" if "role" in s else "split"
        test = s[key].str.lower().eq("test").to_numpy()
        if not test.any(): raise ValueError("Provided outer split has no test participants")
        bad = d.group_id.isin(set(d.loc[test, "group_id"])).to_numpy() & ~test
        removed = int(bad.sum()); d = d.loc[~bad].copy(); test = test[~bad]
    else:
        test = group_fraction(d.group_id, seed, "outer-test") < float(cfg.get("test_fraction", .5))
    d["role"] = "test"
    dev = ~test; u = group_fraction(d.loc[dev, "group_id"], seed, "development-role")
    proportions = cfg.get("development_fractions", [.70, .15, .15])
    if len(proportions) != 3 or min(proportions) <= 0 or abs(sum(proportions)-1) > 1e-8:
        raise ValueError("development_fractions must be positive train,tune,calibration fractions summing to1")
    d.loc[dev, "role"] = np.where(u < proportions[0], "train", np.where(u < sum(proportions[:2]), "tune", "calibration"))
    # Optional training budget uses whole families and outcome-blind order. Unused
    # families are removed from ALL methods, including ABM donors.
    def cap_role(d,role,budget):
        selected = d[d.role.eq(role)]
        counts = selected.group_id.value_counts(); groups = counts.index.to_numpy()
        order = np.argsort(group_fraction(groups, seed, role+"-budget"))
        chosen = []; n = 0
        for g in groups[order]:
            if n >= budget: break
            chosen.append(g); n += int(counts[g])
        return d.loc[~d.role.eq(role) | d.group_id.isin(chosen)].copy()
    label_budget=int(cfg.get("label_budget",0) or 0)
    if label_budget:
        n_train=int(label_budget*proportions[0]/sum(proportions[:2]))
        if n_train<2 or label_budget-n_train<2:raise ValueError("Label budget is too small for separate TRAIN and TUNE")
        d=cap_role(d,"train",n_train);d=cap_role(d,"tune",label_budget-n_train)
    budget = int(cfg.get("train_n", 0) or 0)
    if budget:d=cap_role(d,"train",budget)
    check_roles(d)
    return d.reset_index(drop=True), removed

def phenotype_values(raw, traits, cfg):
    definitions = {}; result = pd.DataFrame({"id": raw.id})
    columns = cfg.get("phenotype_columns", {})
    for trait in traits:
        column = columns.get(trait, "auto" if trait == "t2dm" else trait)
        if trait == "height" and column == "height" and column not in raw and "bb_Ht" in raw: column = "bb_Ht"
        if trait == "ldl" and column == "ldl" and column not in raw and "bb_LDL" in raw: column = "bb_LDL"
        if trait == "t2dm" and column == "auto":
            if not {"t2dm.Yr2e", "t2dm.Yt2e"} <= set(raw):
                # The existing shared UKB endpoint definition is reused verbatim.
                from data import derive_t2dm_endpoints
                endpoint, definition = derive_t2dm_endpoints(raw)
                raw = raw.copy(); raw[["t2dm.Yr2e", "t2dm.Yt2e"]] = endpoint
            else: definition = {"source": "t2dm.Yr2e,t2dm.Yt2e"}
            pre = pd.to_numeric(raw["t2dm.Yr2e"], errors="coerce")
            inc = pd.to_numeric(raw["t2dm.Yt2e"], errors="coerce")
            y = np.where(pre == 1, 1., np.where((pre == 0) | inc.isin([0.,1.]), 0., np.nan))
            definitions[trait] = {"endpoint": "baseline prevalent 0/1; incident valid followup can establish baseline control", **definition}
        else:
            if column not in raw: raise ValueError(f"Missing phenotype {trait}: column {column}")
            y = pd.to_numeric(raw[column], errors="coerce").to_numpy(float)
            definitions[trait] = {"source_column": column, "transformation": "identity; no unit conversion or medication correction"}
        if trait == "t2dm" and not np.isin(np.asarray(y)[np.isfinite(y)], [0.,1.]).all():
            raise ValueError("T2DM must be an explicit baseline 0/1 phenotype")
        result[trait] = y
    return result, definitions

def extend_traits(cfg, root, traits):
    """Add outcomes only after the driver verifies the frozen cohort receipt."""
    out=Path(root)/"tmp/cohort"
    d=identified(read_table(out/"participants.tsv.gz"),"frozen cohort")
    check_roles(d)
    inputs=cfg["inputs"];opts=cfg["cohort"]
    source=inputs.get("cohort_file") or inputs["phenotype_file"]
    selected={"id","eid","IID","#IID","ID_2","height","bb_Ht","ldl","bb_LDL",
              "t2dm.Yr2e","t2dm.Yt2e","date_attend","date_lost","date_death","birth_date",
              "fod_icd10_t2dm","cnt_icd10_t2dm",*opts.get("phenotype_columns",{}).values()}
    raw=identified(read_table(source,sorted(selected)),"phenotype/cohort")
    outcomes,definitions=phenotype_values(raw,traits,opts)
    dev=~d.role.eq("test")
    pending=[]
    for trait in traits:
        y=aligned(outcomes[["id",trait]].rename(columns={trait:"y"}),d.id,"outcome")
        for table,path in ((y.loc[dev],out/f"{trait}.development.tsv.gz"),
                           (y.loc[~dev],out/"_private"/f"{trait}.test.tsv.gz")):
            if path.exists():
                existing=read_table(path)
                existing["y"]=pd.to_numeric(existing["y"],errors="raise")
                try:
                    pd.testing.assert_frame_equal(existing.reset_index(drop=True),table.reset_index(drop=True),check_dtype=False)
                except AssertionError as exc:
                    raise ValueError(f"Existing phenotype labels differ: {path}") from exc
            else:pending.append((table,path))
    # Validate every existing label before publishing any missing output.
    for table,path in pending:atomic_table(table,path)
    import json
    metadata=json.loads((out/"cohort.json").read_text())
    metadata.setdefault("phenotype_definitions",{}).update(definitions)
    atomic_json(metadata,out/"cohort.json")


def prepare(cfg, root, traits):
    root = Path(root); out = root / "tmp/cohort"; inputs = cfg["inputs"]; opts = cfg["cohort"]
    cov = cfg["benchmark"]["covariates"]
    selected = {"id","eid","IID","#IID","ID_2","age","sex","bb_age","bb_sex", "height", "bb_Ht", "ldl","bb_LDL",
                "t2dm.Yr2e","t2dm.Yt2e","date_attend","date_lost","date_death","birth_date","fod_icd10_t2dm","cnt_icd10_t2dm",
                "group_id","group","family_id","ancestry","genetic_ancestry","target","role","split", *cov,
                *opts.get("phenotype_columns", {}).values(), *opts.get("covariate_columns", {}).values(),
                *[f"PC{k}" for k in range(1,21)]}
    source = inputs.get("cohort_file") or inputs["phenotype_file"]
    raw = identified(read_table(source, sorted(selected)), "phenotype/cohort")
    raw = raw[~raw.id.str.startswith("-")].copy()
    outcomes, definitions = phenotype_values(raw, traits, opts)
    d = pd.DataFrame({"id": raw.id})
    for name in cov:
        column = opts.get("covariate_columns", {}).get(name, name)
        if column not in raw and name in {"age","sex"} and f"bb_{name}" in raw: column = f"bb_{name}"
        if column in raw: d[name] = pd.to_numeric(raw[column], errors="coerce")
    pc_cols = [f"PC{k}" for k in range(1, max(20, int(cfg["disco"].get("distance_pcs", 10))) + 1)]
    # Explicit projected PCs override phenotype PCs, which may be UKB-internal PCs.
    pca_path = inputs.get("pca_file")
    if pca_path:
        pcs = aligned(read_table(pca_path), d.id, "1KG-projected PCs", allow_missing=True)
        for col in pc_cols:
            if col not in pcs: raise ValueError(f"Projected PCA missing {col}")
            d[col] = pd.to_numeric(pcs[col], errors="coerce").to_numpy()
    else:
        for col in pc_cols:
            if col in raw: d[col] = pd.to_numeric(raw[col], errors="coerce")
    if inputs.get("ancestry_file"):
        ancestry = aligned(read_table(inputs["ancestry_file"]), d.id, "ancestry", allow_missing=True)
    else: ancestry = raw
    acol = next((c for c in ("ancestry","genetic_ancestry","target") if c in ancestry), None)
    if acol is None: raise ValueError("An ancestry table is required for anchors and stratified evaluation")
    d["ancestry"] = ancestry[acol].fillna("OTH").astype(str).to_numpy()
    if inputs.get("group_file"):
        families = aligned(read_table(inputs["group_file"]), d.id, "family groups")
        gcol = next((c for c in ("group_id","group","family_id") if c in families), None)
        if gcol is None: raise ValueError("group_file needs group_id/group/family_id")
        d["group_id"] = families[gcol].to_numpy()
        family_source = inputs["group_file"]
    elif any(c in raw for c in ("group_id","group","family_id")):
        gcol = next(c for c in ("group_id","group","family_id") if c in raw)
        d["group_id"] = raw[gcol].to_numpy(); family_source = source
    else:
        relationship = inputs.get("relationship_file")
        if not relationship: raise ValueError("Supply group_file or relationship_file; family IDs are mandatory")
        edges = read_table(relationship)
        if not {"ID1","ID2"} <= set(edges) or edges[["ID1","ID2"]].isna().any().any():
            raise ValueError("Relationship file needs nonmissing ID1/ID2")
        d["group_id"] = family_components(d.id, edges); family_source = relationship
    # Align metadata across every chromosome without decoding any dosage matrix.
    common = set(d.id); first = None
    for chromosome in cfg.get("chromosomes", list(range(1,23))):
        prefix = Path(inputs["genotype_dir"]) / f"chr{chromosome}"
        if not prefix.with_suffix(".pgen").is_file(): raise FileNotFoundError(str(prefix)+".pgen")
        if not prefix.with_suffix(".pvar").is_file() and not Path(str(prefix)+".pvar.zst").is_file(): raise FileNotFoundError(str(prefix)+".pvar[.zst]")
        samples = identified(read_table(str(prefix)+".psam"), "PGEN PSAM")
        if first is None: first = samples
        common.intersection_update(samples.id)
    d = d[d.id.isin(common)].copy()
    missing = [c for c in [*cov,*pc_cols] if c not in d]
    if missing: raise ValueError(f"Missing covariates/projected PCs: {missing}")
    # PC distances cannot tolerate missing PC entries. Covariates use TRAIN medians.
    eligible = d[pc_cols].notna().all(axis=1)
    d = d.loc[eligible].copy()
    if d.empty:
        raise ValueError("No eligible samples after genotype/PCA alignment; check participant ID formats and nonmissing projected PCs")
    d["iid"] = d.id
    fcol = next((c for c in ("#FID","FID","fid") if c in first), None)
    d["fid"] = d.id.map(first.set_index("id")[fcol]) if fcol else "0"
    d = d.reset_index(drop=True)
    old = identified(read_table(inputs["split_file"]), "outer split") if inputs.get("split_file") else None
    if old is None and "split" in raw: old = raw[["id","split"]]
    outside_roster = int((~d.id.isin(old.id)).sum()) if old is not None else 0
    d, removed = choose_roles(d, opts, old)
    if d.id.empty: raise ValueError("No eligible samples")
    atomic_table(d, out/"participants.tsv.gz")
    for trait in traits:
        y = aligned(outcomes[["id",trait]].rename(columns={trait:"y"}), d.id, "outcome")
        dev = ~d.role.eq("test")
        atomic_table(y.loc[dev], out/f"{trait}.development.tsv.gz")
        atomic_table(y.loc[~dev], out/"_private"/f"{trait}.test.tsv.gz")
    atomic_table(d.groupby(["role","ancestry"]).agg(n=("id","size"),families=("group_id","nunique")).reset_index(), root/"qc/cohort.tsv")
    atomic_json({"schema":VERSION,"sources":{"cohort":stamp(source),"family":stamp(family_source)},
        "phenotype_definitions":definitions,"role_counts":d.role.value_counts().to_dict(),
        "removed_development_relatives_of_old_test":removed,
        "excluded_ids_outside_frozen_split":outside_roster,
        "split_uses_outcome_values":False,"outer_test_preserved":old is not None,
        "validation_status":opts.get("validation_status","exploratory_previously_inspected_test"),
        "ancestry_note":"Discrete analysis labels; OTH is not a complete definition of admixture",
        "genome_build":38,"genotype_dir":inputs["genotype_dir"]}, out/"cohort.json")
    return d

if __name__ == "__main__":
    from core2 import load_config
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--config",required=True);p.add_argument("--root",required=True)
    a=p.parse_args(); c=load_config(a.config);prepare(c,a.root,c.get("traits",["height","ldl","t2dm"]))
