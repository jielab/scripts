"""Independent prediction-entry task graphs; none invokes another entry."""
from output_layout import project_path, trait_path, module_root
from pathlib import Path


def build_jobs(cfg, args, add):
    entry = args.entry
    if entry == "csx" and args.module in {"all", "infer"}:
        cohort = add("prepare") if args.module == "all" else None
        inputs = add("inputs")
        for trait in cfg["traits"]:
            prepared = add("infer_prepare", trait, [inputs])
            partitions = [add("infer_chr", trait, [prepared], chromosome=c)
                          for c in cfg["chromosomes"]]
            inference = add("infer", trait, partitions)
            if args.module == "all":
                scoring = add("score", trait, [cohort, inference])
                add("fit", trait, [scoring])
    elif entry=="prsformer" and args.prsformer_stage=="all":
        from prsformer2 import _section
        opts=cfg["prsformer"]
        def enabled(stage):
            return _section(opts,stage).get("enabled",True)
        for trait in cfg["traits"]:
            prepare=add("prsformer",trait,action="prepare") if enabled("prepare") else None
            action="train" if enabled("train") else "predict" if enabled("predict") else None
            if action:
                add("prsformer",trait,[prepare] if prepare else [],action=action)
    elif args.module in {"pca", "prepare", "inputs"}:
        add(args.module)
    else:
        for trait in cfg["traits"]:
            add(args.module, trait)


def validate_inputs(cfg, args):
    root = Path(cfg["analysis_root"])
    paths = []
    if args.entry == "csx" and args.module in {"all", "prepare"}:
        paths = [Path(cfg["inputs"][key]) for key in ("pca_file", "ancestry_file", "split_file") if cfg["inputs"].get(key)]
    elif args.entry == "disco":
        paths = [project_path(root, "tmp/cohort/participants.tsv.gz")]
        for trait in cfg["traits"]:
            paths.extend([project_path(root, f"tmp/cohort/{trait}.development.tsv.gz"),
                          project_path(root, f"tmp/{trait}/scores.tsv.gz"), project_path(root, f"tmp/{trait}/baseline.predictions.tsv.gz"),
                          project_path(root, f"tmp/{trait}/benchmark.json"), project_path(root, f"tmp/{trait}/inference/weights/manifest.tsv")])
    elif args.entry == "prsformer":
        # Preparation owns its cohort and target map. Validate the same raw
        # input contract as the worker, never artifacts from the CSx entry.
        from prsformer_inputs import input_contract
        contract = input_contract(cfg, cfg["traits"])
        paths = [Path(record["path"]) for record in contract["files"]]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise ValueError("Independent workflow requires existing input artifacts; no other workflow "
                         "will be launched. No jobs submitted.\n" + "\n".join(missing))
    if args.entry == "csx" and args.module in {"all", "infer"}:
        from genomics import pfile
        for chrom in cfg["chromosomes"]:
            pfile(Path(cfg["inputs"]["genotype_dir"])/f"chr{chrom}")
        root = Path(cfg["inputs"]["ld_ref"])
        from infer import select_ld_reference
        select_ld_reference(root,cfg["chromosomes"])
