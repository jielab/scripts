#!/usr/bin/env python3
"""Lightweight stage/trait DAG planning; no numerical imports or input preparation."""
from __future__ import annotations
import argparse
import json
import math
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "0f"))
from workflow_jobs import Job, resources, run_plan, write_plan

ENTRIES = {"grid": "grid.sh", "pca": "1.pca.sh", "csx": "2.csx.sh", "disco": "3.disco.sh", "prsformer": "4.prsformer.sh", "install": "install.sh"}
IDENTITY = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")
COMMON_ENV = {"submit", "SUBMIT", "DATE_SCRATCH", "DIR_WORK", "DIR_DATA", "DIR_SCRATCH", "DIR_SOFTWARE",
              "ANALYSIS_ROOT", "SCRIPTS_ROOT", "OS_FILES_ROOT", "OS_SNP_LIST", "SHARED_PYTHON", "WORKFLOW_PYTHON"}


def _env():
    return {k: v for k, v in os.environ.items()
            if (k in COMMON_ENV or k.startswith(("GRID_", "PRSFORMER_")))
            and not re.search(r"TOKEN|PASSWORD|SECRET|CREDENTIAL", k, re.I)}


def option(tokens, names):
    for i, token in enumerate(tokens):
        if token in names:
            return tokens[i + 1] if i + 1 < len(tokens) else None
        for name in names:
            if token.startswith(name + "="):
                return token.split("=", 1)[1]
    return None


def subset(tokens, options, flags=()):
    """Project-specific prerequisite argv, never pass unrelated model options."""
    result = []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        name, sep, value = token.partition("=")
        if name in flags:
            result.append(token)
        elif name in options:
            if sep:
                result += [name, value]
            elif i + 1 < len(tokens):
                result += [name, tokens[i + 1]]
                i += 1
            else:
                raise ValueError(f"Missing value: {name}")
        i += 1
    return result


def without(tokens, options):
    result, i = [], 0
    while i < len(tokens):
        token = tokens[i]
        if token.split("=", 1)[0] in options:
            i += 1 if "=" in token else 2
        else:
            result.append(token)
            i += 1
    return result


def score_subset(tokens):
    """Use the native scorer's chromosome order and permanent subdirectory."""
    text = option(tokens, {"--chrs", "--chr"}) or os.getenv("GRID_CHRS", "1-22")
    chromosomes = []
    for part in re.split(r"[,;\s]+", text):
        if not part:
            continue
        if re.fullmatch(r"\d+-\d+", part):
            first, last = map(int, part.split("-"))
            chromosomes.extend(range(first, last + 1))
        else:
            part = part.replace("chr", "").replace("CHR", "")
            if not part.isdigit():
                raise ValueError("Invalid autosome list: " + text)
            chromosomes.append(int(part))
    chromosomes = list(dict.fromkeys(chromosomes))
    if not chromosomes or any(c < 1 or c > 22 for c in chromosomes):
        raise ValueError("Invalid autosome list: " + text)
    return "" if chromosomes == list(range(1, 23)) else "chr" + ",".join(map(str, chromosomes))


def _profile(args):
    cfg = json.loads(Path(args.hpc_config).read_text())
    if cfg.get("schema_version") != 1:
        raise ValueError("Unsupported HPC resource profile")
    return cfg


def _resources(args, config, stage):
    r = {**config.get("defaults", {}), **config.get("stages", {}).get(stage, {})}
    local = os.getenv("submit", "local") == "local"
    if local:
        r.update(config.get("local", {"cores": 8, "memory_gib": 16, "gpus": 0}))
        if os.getenv("GRID_MEMORY_CAP_GB"):
            r["memory_gib"] = float(os.environ["GRID_MEMORY_CAP_GB"])
    gpu_capable = stage in {"weights", "score", "blockprepare", "blocknet", "prsformer", "prsformer2", "countprepare", "twins", "install"}
    wants_gpu = gpu_capable and (args.device == "cuda" or (args.device == "auto" and r.get("gpus", 0) > 0))
    if wants_gpu:
        if not local:
            r.update(config.get("gpu", {}))
        else:
            r["gpus"] = 1
        if args.gpu_queue:
            r["queue"] = args.gpu_queue
        if args.gpu_count is not None:
            if args.gpu_count < 1:
                raise ValueError("A CUDA stage requires --gpu-count >= 1")
            r["gpus"] = args.gpu_count
    else:
        r["gpus"] = 0
        if args.queue:
            r["queue"] = args.queue
    for attr, key in (("job_cores", "cores"), ("job_memory_gb", "memory_gib"), ("job_hours", "hours")):
        value = getattr(args, attr)
        if value is not None:
            r[key] = value
    if str(r.get("queue")) == "short" and int(r.get("cores", 0)) < 40:
        r["cores"] = 40
    if args.threads and args.threads > int(r["cores"]) and stage not in {"architecture", "final", "evaluate", "targetfit"}:
        raise ValueError(f"--threads exceeds allocated cores for {stage}")
    if not local and not (args.gpu_queue if wants_gpu else args.queue) and os.environ.get('HPC_CPU_QUEUE'):
        if wants_gpu:
            r['queue'] = os.environ['HPC_GPU_QUEUE']
        else:
            memory = float(r['memory_gib']); cores = int(r['cores'])
            tier = 'HPC_CPU_QUEUE' if memory <= 180 and cores <= 40 else 'HPC_MEDIUM_QUEUE' if memory <= 480 and cores <= 64 else 'HPC_HIGHMEM_QUEUE' if memory <= 960 and cores <= 120 else 'HPC_HUGEMEM_QUEUE'
            r['queue'] = os.environ[tier]
    return resources(r)


RAW_OPTIONS = {"--pheno-file", "--ancestry-file", "--group-col", "--group-file", "--split-group-file", "--relationship-file", "--dir-gen", "--snpinfo", "--chrs", "--remove", "--height-col", "--ldl-col", "--t2dm-col"}
PREP_OPTIONS = {"--cohort", "--prepared", "--covariates", "--config", "--block-config"} | RAW_OPTIONS
FIT_OPTIONS = {"--cohort", "--prepared", "--family-col", "--scores", "--score-file", "--score-dir", "--validation-n", "--covariates", "--seed", "--prs-standardization", "--models", "--block-config"}
BLOCK_OPTIONS = {"--cohort", "--prepared", "--family-col", "--config", "--block-config", "--ld-block-dir", "--ld-block-files", "--blocks", "--score-dir", "--covariates", "--validation-n", "--seed", "--epochs", "--batch-size"}


def count_block_config(tokens):
    """Resolve the same outcome-free producer config used by count workers."""
    path = Path(option(tokens, {"--config"}) or ROOT / "config/count_abm.json").expanduser()
    if path.is_dir():
        path /= "count_abm.json"
    config = json.loads(path.read_text())
    if config.get("schema_version") != "grid.count_abm.config.1":
        raise ValueError("GRID count requires count_abm.json, not a historical fitted-model config")
    explicit = option(tokens, {"--block-config"})
    if explicit:
        return str(Path(explicit).expanduser().resolve())
    block = Path(config["block_config"]).expanduser()
    return str((block if block.is_absolute() else path.parent / block).resolve())


def build_jobs(args, extra):
    cfg = _profile(args)
    traits = args.traits.split(",")
    if not traits or any(not IDENTITY.fullmatch(t) for t in traits) or len(set(traits)) != len(traits):
        raise ValueError("--traits requires unique safe trait names")
    root = str(Path(args.analysis_root).expanduser().resolve())
    if args.entry == "grid" and args.module in {"prepare", "all"} and not option(extra, {"--input-manifest"}):
        extra = without(extra, {"--block-config"}) + ["--block-config", count_block_config(extra)]
    base_env = _env()
    base_env.update(GRID_OUTPUT_ROOT=root, GRID_ANALYSIS_ROOT=root, WORKFLOW_WORKER="1")
    for key, names in {
        "GRID_PHE_FILE": {"--pheno-file", "--phe-file", "--phe"},
        "GRID_TARGET_DIR": {"--dir-gen"}, "GRID_CSX_SNPINFO": {"--snpinfo", "--csx-snpinfo"},
        "GRID_ANCESTRY_FILE": {"--ancestry-file"}, "GRID_GROUP_FILE": {"--group-file", "--split-group-file"},
        "GRID_REL_FILE": {"--relationship-file"}, "GRID_SCORE_DIR": {"--score-dir", "--dir-score"},
        "GRID_GWAS_DIR": {"--dir-gwas"}, "GRID_REMOVE": {"--remove", "--exclude-ids"},
    }.items():
        value = option(extra, names)
        if value is not None:
            base_env[key] = value
    local = os.getenv("submit", "local") == "local"
    jobs = []
    def add(key, stage, entry, module, selected, depends=(), passed=None, device=False):
        r = _resources(args, cfg, stage)
        stage_threads = min(args.threads or 8, r["cores"])
        env = dict(base_env)
        env.update(GRID_MEMORY_CAP_GB=str(math.ceil(r["memory_gib"])), GRID_THREADS=str(stage_threads),
                   GRID_JOBS=str(max(1, min(4, r["cores"] // stage_threads))))
        worker_device = args.device if local else ("cuda" if r["gpus"] else "cpu")
        env["GRID_CSX_BACKEND"] = worker_device
        env["WORKFLOW_DEVICE"] = worker_device
        command = ["bash", str(ROOT / ENTRIES[entry])]
        if entry not in {"pca", "install"}:
            command.append(module)
        if entry != "install" and module != "self-check":
            command += ["--traits", ",".join(selected)]
            output_flag = "--output-root" if entry in {"csx", "disco", "pca"} and module in {"all", "weights", "score"} else "--analysis-root"
            command += [output_flag, root]
        command += list(extra if passed is None else passed)
        if entry == "csx" and stage in {"weights", "score"}:
            children = option(command, {"--jobs"})
            if children and int(children) * int(env["GRID_THREADS"]) > r["cores"]:
                raise ValueError("CSx --jobs * --threads exceeds --job-cores; reduce within-job parallelism or request more cores")
        if device:
            command += ["--device", worker_device,
                        "--threads", str(stage_threads)]
        # Preserve the public caller's meaning of relative scientific paths.
        jobs.append(Job(key, command, list(depends), r, env, str(Path.cwd())))
        return key

    module = args.module
    checking = "--check" in extra and args.entry != "install"
    if checking:
        supported = ((args.entry == "pca" and module == "all")
                     or (args.entry == "csx" and module in {"all", "weights", "score"})
                     or (args.entry == "disco" and module in {"all", "score"})
                     or (args.entry == "prsformer" and module == "prsformer2"))
        if not supported:
            raise ValueError(f"{args.entry} {module} does not support --check; use --plan-only to inspect its DAG")
        # Native check paths validate inputs without preparing the shared cohort
        # or fitting target models. Do not append ordinary downstream workers.
        if args.entry in {"csx", "disco"}:
            stage = "weights" if args.entry == "csx" and module != "score" else "score"
            for trait in traits:
                add(f"{trait}.check", stage, args.entry, module, [trait])
            return jobs
    if args.entry in {"pca", "install"}:
        if module != "all":
            raise ValueError(f"{args.entry} accepts options without a positional module")
        add(args.entry, args.entry, args.entry, "all", traits)
    elif args.entry == "grid":
        if module in {"self-check", "cleanup"}:
            add(module, "final", "grid", module, traits)
            return jobs
        if module in {"assess", "abm"}:
            module = "evaluate"
        if module not in {"prepare", "twins", "evaluate", "architecture", "final", "all"}:
            raise ValueError("GRID fitted block models use 2.csx.sh blockridge/blocknet")
        if module == "architecture":
            add("architecture", "architecture", "grid", "architecture", traits,
                passed=without(extra, RAW_OPTIONS))
            return jobs
        preparation = []
        if module in {"prepare", "all"} and not option(extra, {"--input-manifest"}):
            preparation = [add("shared.prepare", "prepare", "csx", "prepare", traits,
                               passed=subset(extra, PREP_OPTIONS - {"--config"}))]
            # The outcome-free producer owns a shared genotype/dosage write
            # lock. Finish its trait writes before parallel read-only models;
            # never send independent writers to fail on one nonblocking lock.
            for trait in traits:
                preparation = [add(f"{trait}.blockprepare", "blockprepare", "csx", "blockprepare",
                    [trait], preparation, passed=subset(extra, BLOCK_OPTIONS - {"--config"}), device=True)]
        end = []
        for trait in traits:
            parent = preparation
            for phase in (["prepare", "twins", "evaluate"] if module == "all" else [module]):
                if phase == "final":
                    break
                resource_stage = "countprepare" if phase == "prepare" else phase
                key = add(f"{trait}.{phase}", resource_stage, "grid", phase, [trait], parent,
                          passed=without(extra, RAW_OPTIONS), device=phase in {"prepare", "twins"})
                parent = [key]
            end.extend(parent)
        if module == "all":
            end = [add("architecture", "architecture", "grid", "architecture", traits, end,
                       passed=without(extra, RAW_OPTIONS))]
        if module in {"all", "final"}:
            add("final", "final", "grid", "final", traits, end, passed=without(extra, RAW_OPTIONS))
    elif args.entry == "csx":
        if module not in {"all", "weights", "score", "prepare", "targetfit", "blockprepare", "blockridge", "blocknet"}:
            raise ValueError(f"Unknown CSx module: {module}")
        if module in {"all", "weights", "score"}:
            end = []
            legacy = [t for t in extra]
            # target-only options belong to targetfit, never to upstream MCMC.
            clean = []
            i = 0
            while i < len(legacy):
                name = legacy[i].split("=", 1)[0]
                if name in {"--cohort", "--prepared", "--family-col", "--validation-n", "--covariates", "--prs-standardization", "--scores", "--score-file", "--block-config"}:
                    i += 1 if "=" in legacy[i] else 2
                else:
                    clean.append(legacy[i]); i += 1
            for trait in traits:
                parent = []
                if module in {"all", "weights"}:
                    parent = [add(f"{trait}.weights", "weights", "csx", "weights", [trait], passed=clean)]
                if module in {"all", "score"}:
                    parent = [add(f"{trait}.score", "score", "csx", "score", [trait], parent, passed=clean)]
                end.extend(parent)
            if module == "all":
                parents = end
                if not option(extra, {"--cohort", "--prepared"}):
                    parents = [add("shared.prepare", "prepare", "csx", "prepare", traits, end, passed=subset(extra, PREP_OPTIONS - {"--score-dir"}))]
                fitted_options = subset(extra, FIT_OPTIONS)
                chrom_subset = score_subset(extra)
                if chrom_subset:
                    fitted_options += ["--score-subset", chrom_subset]
                for trait in traits:
                    add(f"{trait}.targetfit", "targetfit", "csx", "targetfit", [trait], parents, passed=fitted_options)
        elif module == "prepare":
            add("shared.prepare", "prepare", "csx", "prepare", traits, passed=extra)
        elif module == "targetfit":
            for trait in traits:
                add(f"{trait}.targetfit", "targetfit", "csx", "targetfit", [trait], passed=extra)
        else:
            deps = []
            if not option(extra, {"--cohort", "--prepared"}):
                deps = [add("shared.prepare", "prepare", "csx", "prepare", traits, passed=subset(extra, PREP_OPTIONS))]
            for trait in traits:
                deps = [add(f"{trait}.blockprepare", "blockprepare", "csx", "blockprepare", [trait], deps,
                            passed=subset(extra, BLOCK_OPTIONS, {"--replace"}), device=True)]
            if module != "blockprepare":
                for trait in traits:
                    add(f"{trait}.{module}", module, "csx", module, [trait], deps,
                        passed=without(extra, RAW_OPTIONS), device=True)
    elif args.entry == "disco":
        if module not in {"all", "score", "evaluate", "prepare"}:
            raise ValueError("Disco modules: all, score, evaluate, prepare")
        if module == "prepare":
            add("shared.prepare", "prepare", "disco", "prepare", traits, passed=subset(extra, PREP_OPTIONS))
        else:
            deps = []
            if module == "all":
                if not option(extra, {"--cohort", "--prepared"}):
                    deps = [add("shared.prepare", "prepare", "disco", "prepare", traits, passed=subset(extra, PREP_OPTIONS - {"--score-dir"}))]
            for trait in traits:
                add(f"{trait}.{module}", "evaluate" if module == "evaluate" else "score", "disco", module, [trait], deps)
    elif args.entry == "prsformer":
        if module not in {"prepare", "prsformer", "prsformer2", "all"}:
            raise ValueError("Unknown PRSformer module")
        parent = []
        if module in {"all", "prepare"}:
            parent = [add("shared.prepare", "prepare", "prsformer", "prepare", traits, passed=subset(extra, PREP_OPTIONS))]
        for method in (["prsformer", "prsformer2"] if module == "all" else [module]):
            if method == "prepare":
                continue
            for trait in traits:
                add(f"{trait}.{method}", method, "prsformer", method, [trait], parent, device=True)
    return jobs


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument("--entry", choices=ENTRIES, required=True)
    p.add_argument("--workflow-module", dest="module", default="all", help=argparse.SUPPRESS)
    p.add_argument("--traits", "--trait", default="height,ldl,t2dm")
    p.add_argument("--analysis-root", "--output-root", "--dir-analysis-root", "--dir-output-root",
                   default=os.getenv("GRID_OUTPUT_ROOT", "/mnt/d/analysis/grid"))
    p.add_argument("--hpc-config", default=os.getenv("GRID_HPC_CONFIG", str(ROOT / "config/hpc.json")))
    p.add_argument("--job-dir")
    p.add_argument("--resume-plan", type=Path)
    p.add_argument("--queue")
    p.add_argument("--gpu-queue")
    p.add_argument("--job-cores", type=int)
    p.add_argument("--job-memory-gb", type=float)
    p.add_argument("--job-hours", type=float)
    p.add_argument("--gpu-count", type=int)
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--threads", type=int)
    p.add_argument("--jobs-N", type=int, default=int(os.getenv("jobs_N", "8")))
    p.add_argument("--wait", nargs="?", const="TRUE", choices=["TRUE", "FALSE"], default=os.getenv("SUBMIT_WAIT", "FALSE"))
    p.add_argument("--plan-only", action="store_true")
    p.add_argument("--dry-run", nargs="?", const="TRUE", choices=["TRUE", "FALSE"], default="FALSE")
    return p


def main(argv=None):
    p = parser()
    tokens = list(sys.argv[1:] if argv is None else argv)
    # The public module is the first project token after --entry ENTRY. Do not
    # let argparse mistake an unknown scientific option's VALUE for a module.
    start = 2 if len(tokens) >= 2 and tokens[0] == "--entry" else 0
    if len(tokens) > start and not tokens[start].startswith("-"):
        module = tokens.pop(start)
        tokens += ["--workflow-module", module]
    args, extra = p.parse_known_args(tokens)
    try:
        mode = os.getenv("submit", "local")
        if args.entry == "csx":
            legacy_stage = option(extra, {"--stage"})
            if legacy_stage:
                if args.module != "all" and args.module != legacy_stage:
                    raise ValueError("Positional CSx module and --stage disagree")
                args.module = legacy_stage
                extra = without(extra, {"--stage"})
            backend = option(extra, {"--backend"})
            if backend:
                if backend not in {"auto", "cpu", "cuda"} or (args.device != "auto" and args.device != backend):
                    raise ValueError("--backend must be auto|cpu|cuda and agree with --device")
                args.device = backend
                extra = without(extra, {"--backend"})
        if args.resume_plan:
            path = args.resume_plan.resolve()
        else:
            jobs = build_jobs(args, extra)
            path = write_plan(jobs, args.job_dir or Path(args.analysis_root) / "jobs", mode=mode, name=f"grid-{args.entry}")
        preview = args.plan_only or args.dry_run == "TRUE"
        if not preview:
            import shutil
            if not shutil.which(os.getenv("WORKFLOW_PYTHON", "python3")):
                raise ValueError("WORKFLOW_PYTHON is unavailable on this submission host; set the shared worker interpreter before submitting")
        result = run_plan(path, submit=mode, jobs=args.jobs_N, wait=args.wait == "TRUE", dry_run=preview)
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        p.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
