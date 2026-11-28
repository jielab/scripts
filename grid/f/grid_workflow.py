"""GRID's own task graph. External baseline artifacts are read-only inputs."""
from output_layout import project_path, trait_path, module_root
from pathlib import Path


def build_jobs(cfg, args, add):
    if args.module == "paper-run":
        evaluations = [add("evaluate", trait) for trait in cfg["traits"]]
        add("paper", depends=evaluations)
        return
    if args.module != "run":
        if args.module == "paper":
            add("paper")
        else:
            for trait in cfg["traits"]:
                add(args.module, trait)
        return
    for trait in cfg["traits"]:
        regions = None
        if (cfg["abm"]["geometry"] != "global_prs" or cfg["block_models"]["enabled"]
                or cfg["block_deep"]["enabled"]):
            regions = add("regions", trait)
            add("architecture", trait, [regions])
        abm = add("abm", trait, [regions] if regions else [])
        add("twins", trait, [abm])
        if cfg["block_models"]["enabled"]:
            add("block_models", trait, [regions])
        if cfg["block_deep"]["enabled"]:
            add("abm_deep", trait, [regions])
            if cfg.get("direct_blocknet"):
                direct=add("direct_blocknet", trait, [regions])
                add("blocknet_ablation",trait,[regions])
                add("blocknet_twins",trait,[direct,abm])


def required_inputs(cfg, module):
    root = Path(cfg["analysis_root"])
    if module == "paper-run":
        module = "evaluate"
    if module in {"direct_blocknet","blocknet_ablation"}:
        module = "block_models"
    files = []
    if module in {"run","check","regions"}:
        from ancestry_regions import read_maps
        read_maps(cfg)
    for trait in cfg["traits"]:
        if module in {"run","regions","check"}:
            from driver2 import regional_weights
            regional_weights(cfg, root, trait)
        if module == "check":
            files.extend([trait_path(root,trait,"inference/target/manifest.tsv"),
                          trait_path(root,trait,"disco/model.json")])
        if module in {"run", "check", "regions", "block_models", "abm", "abm_deep", "evaluate", "architecture", "robustness"}:
            files.append(project_path(root, "tmp/cohort/participants.tsv.gz"))
        if module in {"run", "check", "block_models", "abm", "abm_deep", "robustness"}:
            files.extend([project_path(root, f"tmp/cohort/{trait}.development.tsv.gz"),
                          trait_path(root, trait, "baseline.predictions.tsv.gz"), trait_path(root, trait, "predictions.tsv.gz")])
        if module in {"run", "check", "regions"}:
            files.extend([trait_path(root, trait, "benchmark.json"), trait_path(root, trait, "inference/weights/manifest.tsv")])
        if module in {"check", "evaluate", "robustness"}:
            files.extend([project_path(root, f"tmp/cohort/_private/{trait}.test.tsv.gz"), trait_path(root, trait, "predictions.tsv.gz")])
        if module in {"block_models", "abm_deep", "architecture", "robustness"} or (module == "abm" and cfg["abm"]["geometry"] != "global_prs"):
            files.extend(trait_path(root, trait, "regions") / f for f in ("scores.npy", "ids.tsv", "features.tsv"))
        if module in {"twins", "evaluate"}:
            files.append(trait_path(root, trait, "abm/abm.predictions.tsv.gz"))
        if module == "paper":
            files.extend([project_path(root, f"grid.paper/{trait}.performance.tsv"), project_path(root, f"grid.paper/{trait}.ancestry.tsv"),
                          trait_path(root, trait, "benchmark.json"), trait_path(root, trait, "abm/model.json")])
        if module=="blocknet_twins":
            files.extend([project_path(root,"tmp/cohort/participants.tsv.gz"),project_path(root,f"tmp/cohort/{trait}.development.tsv.gz"),
                          trait_path(root,trait,"direct_blocknet/embedding.npy"),trait_path(root,trait,"direct_blocknet/predictions.tsv.gz"),
                          trait_path(root,trait,"abm/abm.predictions.tsv.gz")])
    return list(dict.fromkeys(files))


def validate_inputs(cfg, module):
    missing = [str(p) for p in required_inputs(cfg, module) if not p.is_file()]
    if missing:
        raise ValueError("Required input artifacts are missing. GRID only consumes upstream results; "
                         "it never generates them automatically. No jobs submitted.\n" + "\n".join(missing))
