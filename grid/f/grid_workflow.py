"""GRID's own task graph. External baseline artifacts are read-only inputs."""
from pathlib import Path


def build_jobs(cfg, args, add):
    if args.module == "paper-run":
        evaluations = [add("evaluate", trait) for trait in cfg["traits"]]
        add("final", depends=evaluations)
        return
    if args.module != "run":
        if args.module == "final":
            add("final")
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


def required_inputs(cfg, module):
    root = Path(cfg["analysis_root"])
    if module == "paper-run":
        module = "evaluate"
    files = []
    for trait in cfg["traits"]:
        work = root / "tmp" / trait
        if module in {"run", "check", "regions", "block_models", "abm", "abm_deep", "evaluate", "architecture", "robustness"}:
            files.append(root / "tmp/cohort/participants.tsv.gz")
        if module in {"run", "check", "block_models", "abm", "abm_deep", "robustness"}:
            files.extend([root / f"tmp/cohort/{trait}.development.tsv.gz",
                          work / "baseline.predictions.tsv.gz", work / "predictions.tsv.gz"])
        if module in {"run", "check", "regions"}:
            files.extend([work / "benchmark.json", work / "inference/weights/manifest.tsv"])
        if module in {"check", "evaluate", "robustness"}:
            files.extend([root / f"tmp/cohort/_private/{trait}.test.tsv.gz", work / "predictions.tsv.gz"])
        if module in {"block_models", "abm_deep", "architecture", "robustness"} or (module == "abm" and cfg["abm"]["geometry"] != "global_prs"):
            files.extend(work / "regions" / f for f in ("scores.npy", "ids.tsv", "features.tsv"))
        if module in {"twins", "evaluate"}:
            files.append(work / "abm/abm.predictions.tsv.gz")
        if module == "final":
            files.extend([root / f"final/{trait}.performance.tsv", root / f"final/{trait}.ancestry.tsv",
                          work / "benchmark.json", work / "abm/model.json"])
    return list(dict.fromkeys(files))


def validate_inputs(cfg, module):
    missing = [str(p) for p in required_inputs(cfg, module) if not p.is_file()]
    if missing:
        raise ValueError("Required input artifacts are missing. GRID only consumes upstream results; "
                         "it never generates them automatically. No jobs submitted.\n" + "\n".join(missing))
