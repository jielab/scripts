"""PRSformer-owned cohort and HM3 mapping, with no fitted-model inputs."""
from pathlib import Path
import json
import fcntl

import pandas as pd
from core2 import atomic_json, stamp
from output_layout import module_root


def input_dir(root, traits):
    return module_root(root, "prsformer") / "tmp" / "-".join(traits) / "inputs"


def cohort_dir(root, traits):
    return input_dir(root, traits) / "cohort"


def target_manifest(cfg, root, traits):
    base = cfg.get("prsformer", {})
    section = base.get("prepare", {})
    options = {**base, **(section if isinstance(section, dict) else {})}
    return Path(options.get("target_manifest") or input_dir(root, traits)/"target/manifest.tsv").expanduser().resolve()


def reference_path(cfg):
    explicit = cfg.get("prsformer", {}).get("reference_snps")
    if explicit:
        path = Path(explicit).expanduser()
    else:
        # Fix the default universe to 1KG HM3 even when both panels are installed.
        path = Path(cfg["inputs"]["ld_ref"])/"snpinfo_mult_1kg_hm3"
    if not path.is_file():
        raise FileNotFoundError(f"PRSformer reference SNP metadata missing: {path}; set prsformer.reference_snps if needed")
    return path



def input_contract(cfg, traits):
    inp = cfg["inputs"]
    source = inp.get("cohort_file") or inp["phenotype_file"]
    files = [Path(source)]
    for key in ("pca_file", "ancestry_file", "group_file", "relationship_file", "split_file", "rsid_map"):
        if inp.get(key):
            files.append(Path(inp[key]))
    options = cfg.get("prsformer", {})
    section = options.get("prepare", {})
    custom = (section.get("target_manifest") if isinstance(section, dict) else None) or options.get("target_manifest")
    if custom:
        files.append(Path(custom))
    else:
        files.append(reference_path(cfg))
    for chrom in cfg["chromosomes"]:
        prefix = Path(inp["genotype_dir"])/f"chr{chrom}"
        pvar = Path(str(prefix)+".pvar")
        files.extend([Path(str(prefix)+".pgen"), Path(str(prefix)+".psam"),
                      pvar if pvar.is_file() else Path(str(prefix)+".pvar.zst")])
    return dict(schema=1, traits=list(traits), cohort=cfg["cohort"],
                covariates=cfg["benchmark"]["covariates"],
                distance_pcs=cfg.get("disco", {}).get("distance_pcs", 10),
                chromosomes=cfg["chromosomes"], custom_target=str(custom) if custom else None,
                files=[stamp(p) if p.is_file() else {"path":str(p.resolve()),"missing":True} for p in files])


def ensure_inputs(cfg, root, traits):
    """Freeze inputs per task; validate reuse without touching another workflow."""
    folder = input_dir(root, traits)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/".inputs.lock").open("a+") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        contract = input_contract(cfg, traits)
        receipt = folder/"inputs.json"
        if receipt.is_file():
            previous = json.loads(receipt.read_text())
            if previous["inputs"] != contract:
                raise ValueError("PRSformer raw inputs changed; use a new analysis root")
            if not all(Path(p).is_file() and stamp(Path(p)) == record
                       for p, record in previous["outputs"].items()):
                raise ValueError("PRSformer independent inputs are incomplete or changed; use a new analysis root")
            return
        from cohort2 import prepare
        prepare(cfg, root, traits, out=cohort_dir(root, traits))
        if not contract["custom_target"]:
            # Only metadata is required, never the LD matrices or CSx posteriors.
            reference = pd.read_csv(reference_path(cfg), sep=r"\s+", dtype={"SNP":str})
            if not {"CHR", "SNP", "A1", "A2"} <= set(reference):
                raise ValueError("HM3 metadata requires CHR, SNP, A1 and A2")
            reference = reference.loc[reference.CHR.isin(cfg["chromosomes"])].copy()
            reference["SNP"] = reference.SNP.str.lower()
            for col in ("A1", "A2"):
                reference[col] = reference[col].str.upper()
            if reference.empty or reference.SNP.duplicated().any():
                raise ValueError("HM3 metadata must have nonempty, unique SNP IDs")
            from infer import build_target_reference
            build_target_reference(reference, cfg["inputs"]["genotype_dir"], folder,
                                   cfg["chromosomes"], cfg["inputs"].get("rsid_map"))
        outputs = [p for p in folder.rglob("*") if p.is_file() and p.name != ".inputs.lock"]
        atomic_json(dict(inputs=contract, outputs={str(p):stamp(p) for p in outputs}), receipt)
