"""Independent prediction-entry task graphs; none invokes another entry."""
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
    elif args.module in {"pca", "prepare", "inputs", "prsformer"}:
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
        paths = [root / "tmp/cohort/participants.tsv.gz"]
        for trait in cfg["traits"]:
            paths.extend([root / f"tmp/cohort/{trait}.development.tsv.gz",
                          root / f"tmp/{trait}/scores.tsv.gz", root / f"tmp/{trait}/baseline.predictions.tsv.gz",
                          root / f"tmp/{trait}/benchmark.json", root / f"tmp/{trait}/inference/weights/manifest.tsv"])
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise ValueError("Independent workflow requires existing input artifacts; no other workflow "
                         "will be launched. No jobs submitted.\n" + "\n".join(missing))
    if args.entry == "csx" and args.module in {"all", "infer"}:
        from genomics import pfile
        for chrom in cfg["chromosomes"]:
            pfile(Path(cfg["inputs"]["genotype_dir"])/f"chr{chrom}")
        root = Path(cfg["inputs"]["ld_ref"])
        panels = []
        for kind in ("1kg", "ukbb"):
            complete = (root/f"snpinfo_mult_{kind}_hm3").is_file()
            for pop in ("AFR", "EAS", "EUR", "SAS"):
                available = [root/f"ldblk_{kind}_{suffix}" for suffix in (pop.lower(), pop)]
                complete &= sum(all((p/f"ldblk_{kind}_chr{c}.hdf5").is_file()
                                    for c in cfg["chromosomes"]) for p in available) == 1
            if complete:panels.append(kind)
        if len(panels) != 1:
            raise ValueError("CSx requires exactly one complete four-population LD panel for the requested "
                             "autosomes 1-22. No jobs submitted.")
