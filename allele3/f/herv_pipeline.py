#!/usr/bin/env python3
"""Run current HERV analyses with explicit input and result availability.

Source data are read-only. Missing data are reported without inventing estimates.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import datetime as dt
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from analysis_integrity import identity as content_identity, validate_manifest

VERSION = "6.0.1"
SOURCE = Path(__file__).resolve().parent
DEFAULT_STAGES = ("primary", "archaic", "summary", "spatial", "windows", "alleles", "time", "selection", "candidates", "genotypes", "ld", "report")
STAGES = DEFAULT_STAGES
PREFIX = "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes"
SCAN_STORE = "/mnt/f/annot/gnomad/allele3"


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def state(path):
    return content_identity(path)


def save(path, value):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    temporary = p.with_name(p.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n")
    os.replace(temporary, p)


def first_file(paths):
    return next((Path(p) for p in paths if Path(p).is_file()), None)


def resolve_output_dir(analysis_root, output=None):
    root = Path(analysis_root).expanduser().resolve()
    out = Path(output).expanduser().resolve() if output else root
    if out == root:
        out = root / "final/herv_focus_20261028"
    if out in root.parents:
        raise ValueError("Output cannot be a parent of the analysis root")
    return out


def resolve(args):
    root = Path(args.dir_analysis_root).expanduser().resolve()
    out = resolve_output_dir(root, args.dir_out)
    modern = Path(args.dir_modern).expanduser().resolve() if args.dir_modern else Path(SCAN_STORE)
    if out == modern or out in modern.parents or modern in out.parents:
        raise ValueError("Output and source scan must use separate, non-nested directories")
    ancient = Path(args.dir_ancient).expanduser().resolve()
    prefix = Path(args.pgen_prefix).expanduser().resolve() if args.pgen_prefix else ancient / PREFIX
    pvar = first_file([str(prefix) + ".pvar", str(prefix) + ".pvar.zst"]) or Path(str(prefix) + ".pvar")
    old_metadata = first_file([root / "final/evolution/metadata/sample_metadata.tsv", root / "evolution/metadata/sample_metadata.tsv",
                               root / "pca3d/layered_1000yr/metadata/sample_metadata.tsv"])
    metadata_mode = "explicit" if args.metadata else "existing_metadata" if old_metadata else "prepare_from_original_supplement"
    metadata = Path(args.metadata).expanduser().resolve() if args.metadata else old_metadata or out / "metadata/sample_metadata.tsv"
    coordinate_candidates = []
    for directory in (root / "pca3d/layered_1000yr/pca", root / "final/evolution/pca", root / "pca3d/pca", root / "evolution/pca"):
        coordinate_candidates += [directory / "all.sscore", directory / "all.sscore.gz", directory / "shared.eigenvec", directory / "pca.eigenvec"]
    coordinates = Path(args.coordinates).expanduser().resolve() if args.coordinates else first_file(coordinate_candidates)
    return dict(root=root, out=out, modern=modern, ancient=ancient, prefix=prefix, pvar=pvar,
        targets=Path(args.targets).expanduser().resolve() if args.targets else out / "candidates/targets.tsv",
        catalog=Path(args.allele_catalog).expanduser().resolve() if args.allele_catalog else out / "alleles/site_catalog.tsv.gz",
        selection=Path(args.selection_stats).expanduser().resolve() if args.selection_stats else ancient / "Selection_Summary_Statistics_01OCT2025.tsv.gz",
        metadata=metadata, metadata_mode=metadata_mode, coordinates=coordinates,
        hq=Path(args.hq_variants).expanduser().resolve() if args.hq_variants else ancient / "HighQualityVariants.tsv",
        supplement=Path(args.supplement_xlsx).expanduser().resolve() if args.supplement_xlsx else ancient / "41586_2026_10358_MOESM3_ESM.xlsx",
        chrom_sizes=Path(args.chrom_sizes).expanduser().resolve() if args.chrom_sizes else Path(args.fasta38 + ".fai"))


def stage_selection(args):
    chosen = ([s.strip() for s in args.stages.split(",") if s.strip()] if args.stages
              else list(DEFAULT_STAGES) if args.command == "pipeline" else [args.command])
    unknown = set(chosen) - set(STAGES)
    if unknown or not chosen or len(chosen) != len(set(chosen)):
        raise ValueError(f"Select unique stages from {','.join(STAGES)}; invalid={sorted(unknown)}")
    return [stage for stage in STAGES if stage in chosen]


def modern_files(modern, kind):
    """Audit the files actually read, including merged-counter/shard fallback.

    Keep this metadata-only: preflight must report missing scientific packages
    rather than importing them before the dependency audit can run.
    """
    from scan_store import load_run_metadata, resolve_run_path
    modern = Path(modern)
    sidecars = [modern / name for name in ("relocation.json", "adoption_receipt.json") if (modern/name).is_file()]
    run_path = modern / "run.json"
    inputs = [run_path, *sidecars]
    if not run_path.is_file():
        return inputs
    run = load_run_metadata(modern)
    options = run.get("options", {})
    prefix = options.get("architecture_prefix", "architecture")
    targeted = kind in {"windows", "spatial"}
    counter_kinds = {"primary": ("global", "exposure"), "summary": ("global", "exposure"), "windows": ("global",),
                     "spatial": ("global", "erv") + (("exposure",) if options.get("ld_block_dir") else ())}.get(kind, ())
    # Match herv_modern.architecture_paths: each kind independently prefers its
    # merged file. An absent merged file is legitimate when complete shards exist.
    fallback = []
    for counter in counter_kinds:
        merged = modern / "summary" / f"architecture.{counter}_strata.tsv.gz"
        if merged.is_file():
            inputs.append(merged)
        else:
            fallback.append(counter)
    need_shards = targeted or kind == "candidates" or bool(fallback)
    for directory in run.get("shards", []):
        directory = Path(directory)
        if not directory.is_absolute():
            directory = modern / directory
        if need_shards:
            inputs += [directory / "COMPLETE.json", directory / "scan.json"]
        if kind == "candidates":
            inputs.append(directory / "multiallelic_sites.tsv.gz")
        inputs += [directory / f"{prefix}.{counter}_strata.tsv.gz" for counter in fallback]
        if (targeted or bool(fallback)) and (directory / "COMPLETE.json").is_file():
            # Existing scan.load_run verifies all recorded outputs; include these
            # in the cheap missing-input audit instead of reporting a false error.
            receipt = json.loads((directory / "COMPLETE.json").read_text())
            inputs += [directory / x for x in receipt.get("outputs", {})]
        if targeted and (directory / "scan.json").is_file():
            info = json.loads((directory / "scan.json").read_text())
            path = info.get("input")
            if isinstance(path, dict):
                path = path.get("path")
            if path:
                path = resolve_run_path(run, path)
                inputs.append(Path(path))
                index = first_file([str(path) + ".tbi", str(path) + ".csi"])
                inputs.append(index or Path(str(path) + ".tbi"))
    for key in ("architecture_segdup_bed", "architecture_lcr_bed", "architecture_low_mappability_bed"):
        if (kind == "candidates" or targeted) and options.get(key):
            inputs.append(Path(options[key]))
    if targeted:
        # These are verified/read by targeted_inputs and the regional scanners,
        # even though their old fingerprints are also embedded in run.json.
        for key in ("callable_bed", "fasta", "erv_joined"):
            if options.get(key): inputs.append(Path(options[key]))
        if options.get("fasta"):
            inputs.append(Path(options["fasta"] + ".fai"))
            if options["fasta"].endswith((".gz", ".bgz")):
                inputs.append(Path(options["fasta"] + ".gzi"))
        if options.get("ld_block_dir"):
            inputs.append(Path(options["ld_block_dir"]) / f"EUR.{options.get('architecture_build', '38')}.bed")
        identity = run.get("identity", {})
        annotations = (options.get("architecture_annotations") or {}).get("files", {})
        recorded = [identity.get(k) for k in ("fasta", "fasta_index", "fasta_gzi", "callable_bed")]
        recorded += [annotations.get(k) for k in ("erv_joined", "architecture_segdup_bed", "architecture_lcr_bed", "architecture_low_mappability_bed", "EUR")]
        inputs += [resolve_run_path(run, r["path"]) for r in recorded if isinstance(r, dict) and r.get("path")]
    if (targeted or kind == "candidates") and options.get("architecture_annotation_cache"):
        cache = Path(options["architecture_annotation_cache"])
        marker = cache / "manifest.json"
        if targeted and cache.is_dir():
            inputs.append(marker)
        if marker.is_file():
            inputs.append(marker)
            # prepare_joined_cache validates every committed chromosome file.
            # Candidate selection can instead rebuild an invalid optional split
            # cache, so missing optional files are omitted in that branch only.
            cached = json.loads(marker.read_text())
            inputs += [cache / name for name in cached.get("files", {})
                       if targeted or (cache / name).is_file()]
    return list(dict.fromkeys(inputs))


@dataclass
class Stage:
    name: str
    command: list
    inputs: list
    out: Path
    modules: tuple = ()
    packages: tuple = ()
    executables: tuple = ()
    manifest: str | None = None
    details: dict | None = None


def plan(args, paths, names):
    py = sys.executable
    out, modern = paths["out"], paths["modern"]
    report = str(SOURCE / "herv_report.py")
    commands = {}
    common_report = ["--analysis-root", str(paths["root"]), "--modern-dir", str(modern)]
    commands["primary"] = Stage("primary", [py, str(SOURCE / "herv_primary.py"),
        "--modern-dir", str(modern), "--out", str(out / "primary"), "--bootstrap", str(args.bootstrap), "--seed", str(args.seed)],
        modern_files(modern, "primary") if "primary" in names else [], out / "primary",
        ("herv_primary.py", "herv_modern.py", "scan.py", "scan_store.py", "0.common.py"),
        ("numpy", "matplotlib"), manifest="primary_manifest.json")
    commands["summary"] = Stage("summary", [py, str(SOURCE / "herv_modern.py"), "summary", "--modern-dir", str(modern),
        "--out", str(out / "summary"), "--bootstrap", str(args.bootstrap), "--seed", str(args.seed), "--min-training-sites", str(args.min_training_sites), "--min-support", str(args.min_support)],
        modern_files(modern, "summary") if "summary" in names else [], out / "summary", ("herv_modern.py", "scan.py", "scan_store.py", "0.common.py"), ("numpy", "scipy", "matplotlib"), manifest="summary_manifest.json")
    if "archaic" in names:
        from herv_archaic import discover, modern_inputs
        datasets = discover(args.archaic_root, args.archaic_datasets, args.archaic_module, args.archaic_chromosomes)
        ready = [d for d in datasets if d['status'] == 'ready']
        needs_lift = any(d['build'] != '38' for d in ready)
        cache = args.archaic_cache or str(modern.with_name(modern.name + '-archaic'))
        archaic_cmd = [py, str(SOURCE / 'herv_archaic.py'), '--gu-root', args.archaic_root,
            '--datasets', args.archaic_datasets, '--module', args.archaic_module,
            '--chromosomes', args.archaic_chromosomes,
            '--modern-dir', str(modern), '--out', str(out/'archaic'), '--cache-dir', cache,
            '--chain', args.chain37to38, '--liftover', args.archaic_liftover, '--bcftools', args.bcftools,
            '--min-lod', str(args.archaic_min_lod), '--min-length', str(args.archaic_min_length),
            '--carrier-thresholds', args.archaic_carrier_thresholds, '--bootstrap', str(args.bootstrap),
            '--seed', str(args.archaic_seed), '--min-training-sites', str(args.archaic_min_training_sites)]
        archaic_inputs = [*modern_inputs(modern), *(Path(p) for d in ready for p in d['inputs'])]
        if needs_lift: archaic_inputs.append(Path(args.chain37to38))
        commands['archaic'] = Stage('archaic', archaic_cmd, archaic_inputs, out/'archaic',
            ('herv_archaic.py', 'herv_modern.py', 'scan.py', 'scan_store.py', '0.common.py'),
            ('numpy', 'matplotlib', 'pysam'), ('sort', args.bcftools) + ((args.archaic_liftover,) if needs_lift else ()),
            manifest='manifest.json', details={'datasets': [{k: v for k, v in d.items() if k not in {'files', 'inputs', 'sample_panel'}} for d in datasets],
                'dataset_cache': cache, 'gnomad_access': 'indexed tract-union queries; completed background remains read-only'})
    commands["candidates"] = Stage("candidates", [py, str(SOURCE / "herv_candidates.py"), "--modern-dir", str(modern),
        "--erv-joined", args.erv_joined, "--out", str(out / "candidates"), "--max-sites", str(args.max_sites),
        "--controls-per-case", str(args.controls_per_case), "--seed", str(args.seed),
        "--end-bp", str(args.end_bp), "--body-bins", str(args.body_bins)],
        [*(modern_files(modern, "candidates") if "candidates" in names else []), Path(args.erv_joined)], out / "candidates", ("herv_candidates.py", "herv_modern.py", "scan.py", "scan_store.py", "0.common.py"), ("numpy",), manifest="candidate_manifest.json")
    commands["alleles"] = Stage("alleles", [py, str(SOURCE / "herv_candidates.py"), "--modern-dir", str(modern),
        "--erv-joined", args.erv_joined, "--out", str(out / "alleles"), "--max-sites", str(args.allele_max_sites),
        "--controls-per-case", str(args.allele_controls_per_case), "--all-alleles-catalog", "--seed", str(args.seed),
        "--end-bp", str(args.end_bp), "--body-bins", str(args.body_bins)],
        [*(modern_files(modern, "candidates") if "alleles" in names else []), Path(args.erv_joined)], out / "alleles",
        ("herv_candidates.py", "herv_modern.py", "scan.py", "scan_store.py", "0.common.py"), ("numpy",), manifest="candidate_manifest.json")
    kg_inputs = [paths["targets"], Path(args.kg_sample_panel)]
    # Only chromosomes represented in frozen targets are required for the targeted fetch.
    chromosomes = []
    if "genotypes" in names and paths["targets"].is_file():
        import csv
        with paths["targets"].open() as handle:
            chromosomes = sorted({r["chrom"].removeprefix("chr") for r in csv.DictReader(handle, delimiter="\t")}, key=lambda x: int(x))
    else:
        chromosomes = [str(x) for x in range(1, 23)]
    for chrom in chromosomes:
        vcf = Path(args.kg_vcf_dir) / args.kg_vcf_pattern.format(chrom=chrom)
        kg_inputs += [vcf, first_file([str(vcf) + ".tbi", str(vcf) + ".csi"]) or Path(str(vcf) + ".tbi")]
    commands["genotypes"] = Stage("genotypes", [py, str(SOURCE / "herv_genotypes.py"), "--targets", str(paths["targets"]),
        "--vcf-dir", args.kg_vcf_dir, "--sample-panel", args.kg_sample_panel, "--vcf-pattern", args.kg_vcf_pattern,
        "--out", str(out / "genotypes"), "--seed", str(args.seed), "--min-ac", str(args.kg_min_ac), "--min-call-rate", str(args.min_call_rate)],
        kg_inputs, out / "genotypes", ("herv_genotypes.py",), ("pysam", "numpy"), manifest="genotypes_manifest.json")
    window_cmd = [py, str(SOURCE / "herv_modern.py"), "windows", "--modern-dir", str(modern), "--erv-joined", args.erv_joined,
        "--fasta", args.fasta38, "--out", str(out / "windows"), "--window", str(args.window), "--step", str(args.step), "--flank", str(args.flank),
        "--bootstrap", str(args.bootstrap), "--seed", str(args.seed), "--min-training-sites", str(args.min_training_sites),
        "--min-support", str(args.min_support), "--bcftools", args.bcftools]
    if args.chromosomes:
        window_cmd += ["--chromosomes", args.chromosomes]
    target_inputs = [Path(args.erv_joined), Path(args.fasta38), Path(args.fasta38 + ".fai")]
    if args.fasta38.endswith((".gz", ".bgz")):
        target_inputs.append(Path(args.fasta38 + ".gzi"))
    for key in ("segdup", "lcr", "low_mappability"):
        path = getattr(args, key + "_bed")
        if path:
            window_cmd += ["--" + key.replace("_", "-") + "-bed", path]
            target_inputs.append(Path(path))
    window_inputs = [*(modern_files(modern, "windows") if "windows" in names else []), *target_inputs]
    commands["windows"] = Stage("windows", window_cmd, window_inputs, out / "windows", ("herv_modern.py", "scan.py", "scan_store.py", "0.common.py"),
                                 ("numpy", "scipy", "matplotlib") + (("pysam",) if args.fasta38.endswith((".gz", ".bgz")) else ()), (args.bcftools,), manifest="windows_manifest.json")
    spatial_cmd = [py, str(SOURCE / "herv_modern.py"), "spatial", "--modern-dir", str(modern), "--erv-joined", args.erv_joined,
        "--fasta", args.fasta38, "--out", str(out / "spatial"), "--end-bp", str(args.end_bp), "--body-bins", str(args.body_bins),
        "--bootstrap", str(args.bootstrap), "--seed", str(args.seed), "--min-training-sites", str(args.min_training_sites),
        "--min-support", str(args.min_support), "--bcftools", args.bcftools]
    if args.chromosomes:
        spatial_cmd += ["--chromosomes", args.chromosomes]
    for key in ("segdup", "lcr", "low_mappability"):
        if getattr(args, key + "_bed"):
            spatial_cmd += ["--" + key.replace("_", "-") + "-bed", getattr(args, key + "_bed")]
    spatial_inputs = [*modern_files(modern, "spatial"), *target_inputs] if "spatial" in names else []
    # Shared regional inputs plus the ERV/exposure counters used for reconciliation.
    commands["spatial"] = Stage("spatial", spatial_cmd, list(dict.fromkeys(spatial_inputs)), out / "spatial",
        ("herv_modern.py", "scan.py", "scan_store.py", "0.common.py"), ("numpy", "scipy", "matplotlib") + (("pysam",) if args.fasta38.endswith((".gz", ".bgz")) else ()),
        (args.bcftools,), manifest="spatial_manifest.json")
    time_cmd = [py, str(SOURCE / "herv_allele_time.py"), "run", "--catalog", str(paths["catalog"]), "--pgen-prefix", str(paths["prefix"]),
        "--metadata", str(paths["metadata"]), "--hq-variants", str(paths["hq"]), "--fasta38", args.fasta38, "--fasta37", args.fasta37,
        "--chain38to37", args.chain38to37, "--chain37to38", args.chain37to38, "--out", str(out / "time"),
        "--bin-years", str(args.bin_years), "--min-samples", str(args.min_samples), "--min-call-rate", str(args.min_call_rate),
        "--max-age-bp", str(args.max_age_bp), "--chunk-size", str(args.time_chunk_size), "--pcs", str(args.time_pcs),
        "--bootstrap", str(args.bootstrap), "--seed", str(args.seed)]
    time_inputs = [paths["catalog"], Path(str(paths["prefix"]) + ".pgen"), Path(str(paths["prefix"]) + ".psam"), paths["pvar"],
        paths["hq"], Path(args.fasta38), Path(args.fasta38 + ".fai"), Path(args.fasta37), Path(args.fasta37 + ".fai"), Path(args.chain38to37), Path(args.chain37to38)]
    time_inputs += [Path(path + ".gzi") for path in (args.fasta38, args.fasta37) if path.endswith((".gz", ".bgz"))]
    catalog_manifest = paths["catalog"].parent / "candidate_manifest.json"
    if catalog_manifest.is_file():
        time_inputs.append(catalog_manifest)
    if paths["metadata_mode"] == "prepare_from_original_supplement":
        time_inputs += [paths["supplement"]]
    else:
        time_inputs += [paths["metadata"]]
    if paths["coordinates"]:
        time_cmd += ["--coordinates", str(paths["coordinates"])]
        time_inputs.append(paths["coordinates"])
    for option in ('basis_dir', 'coordinate_receipt', 'provenance_manifest'):
        value = getattr(args, option, None)
        if value:
            time_cmd += ['--' + option.replace('_', '-'), value]
            if option != 'basis_dir': time_inputs.append(Path(value))
    commands["time"] = Stage("time", time_cmd, time_inputs, out / "time", ("herv_allele_time.py", "herv_time.py", "evolve.py", "0.common.py"),
        ("numpy", "scipy", "matplotlib", "pgenlib", "pyliftover", "pysam") + (("openpyxl",) if paths["metadata_mode"] == "prepare_from_original_supplement" else ())
        + (("zstandard",) if any(str(p).endswith(".zst") for p in time_inputs) else ()), manifest="manifest.json")
    selection_inputs = [paths["catalog"], paths["selection"], Path(args.fasta38), Path(args.fasta38 + ".fai"),
        Path(args.fasta37), Path(args.fasta37 + ".fai"), Path(args.chain38to37), Path(args.chain37to38)]
    selection_inputs += [Path(p+".gzi") for p in (args.fasta38,args.fasta37) if p.endswith((".gz",".bgz"))]
    commands["selection"] = Stage("selection", [py, str(SOURCE / "herv_selection.py"), "--catalog", str(paths["catalog"]),
        "--selection-stats", str(paths["selection"]), "--fasta38", args.fasta38, "--fasta37", args.fasta37,
        "--chain38to37", args.chain38to37, "--chain37to38", args.chain37to38,
        "--effect-allele", args.selection_effect_allele, "--out", str(out/"selection")],
        selection_inputs, out/"selection", ("herv_selection.py", "herv_time.py", "0.common.py"), ("pyliftover", "pysam")
        + (("zstandard",) if any(str(p).endswith(".zst") for p in selection_inputs) else ()), manifest="manifest.json")
    ld_cmd = [py, str(SOURCE / "herv_ld.py"), "--ld-block-dir", args.ld_block_dir, "--erv-joined", args.erv_joined,
              "--build", "38", "--out", str(out / "ld"), "--chrom-sizes", str(paths["chrom_sizes"]), "--bootstrap", str(args.ld_bootstrap), "--seed", str(args.seed)]
    ld_inputs = [Path(args.ld_block_dir) / f"{ancestry}.38.bed" for ancestry in ("AFR", "EAS", "EUR", "SAS")]
    ld_inputs += [Path(args.erv_joined), paths["chrom_sizes"]]
    if args.ld_covariates:
        ld_cmd += ["--covariates", args.ld_covariates]
        ld_inputs.append(Path(args.ld_covariates))
    commands["ld"] = Stage("ld", ld_cmd, ld_inputs, out / "ld", ("herv_ld.py", "scan.py", "0.common.py"), ("numpy", "scipy", "matplotlib"), manifest="manifest.json")
    # Optional regional association. A request never silently starts VCF counting.
    commands["report"] = Stage("report", [py, report, "report", *common_report, "--out", str(out)], [], out, ("herv_report.py", "herv_modern.py", "scan.py", "0.common.py"), ("numpy", "matplotlib"))
    if getattr(args, 'pca_temporal_dir', None):
        commands['report'].command += ['--pca-temporal-dir', args.pca_temporal_dir]
    return [commands[name] for name in names]


def availability(stage):
    missing = [str(p) for p in stage.inputs if not Path(p).is_file()]
    missing += [str(SOURCE / m) for m in stage.modules if not (SOURCE / m).is_file()]
    packages = [name for name in stage.packages if importlib.util.find_spec(name) is None]
    programs = [name for name in stage.executables if shutil.which(name) is None]
    return dict(missing_inputs=list(dict.fromkeys(missing)), missing_packages=packages, missing_programs=programs)


def snapshot(stage):
    inputs = list(dict.fromkeys(stage.inputs + [SOURCE / p for p in stage.modules] + [SOURCE / p for p in ("analysis_integrity.py", "pca_contract.py", "pca_coordinate_receipts.py", "temporal_stats.py", "herv_joint.py", "allele_measurement.py", "herv_primary.py")]))
    return dict(schema=VERSION, command=stage.command, inputs=[state(p) for p in inputs])


def outputs(directory):
    return {str(p):state(p) for p in directory.rglob("*") if p.is_file() and not p.name.endswith((".partial", ".lock")) and "logs" not in p.parts}


def completed_result(stage):
    if not stage.manifest:
        if stage.name == "report" and not (stage.out / "report.html").is_file():
            raise RuntimeError("Report command exited without report.html")
        return {}
    path = stage.out / stage.manifest
    if not path.is_file():
        if not path or not path.is_file():
            raise RuntimeError(f"{stage.name} exited successfully without required manifest: {stage.manifest}")
    result = validate_manifest(stage.out, stage.manifest)
    allowed = {
        'primary': {'complete'}, 'summary': {'complete'}, 'spatial': {'complete'}, 'windows': {'complete'},
        'alleles': {'complete'}, 'candidates': {'complete'}, 'selection': {'complete'}, 'ld': {'complete'},
        'genotypes': {'complete', 'not_run'},
        'time': {'complete_descriptive_imputed_panel', 'complete_descriptive_all_allele_panel',
                 'complete_adjustment_unavailable', 'no_usable_exact_high_quality_contrast', 'no_usable_contrast'},
    }
    if result.get('status') not in allowed.get(stage.name, {'complete'}):
        raise RuntimeError(f'{stage.name} has invalid terminal status: {result.get("status")}')
    return {"result_status":result.get("status", "unknown"), "result_manifest":str(path)}


def run_stage(stage, receipts, always=False):
    available = availability(stage)
    if any(available.values()):
        print(f"[not_available] {stage.name}: inputs={len(available['missing_inputs'])}, packages={available['missing_packages']}, programs={available['missing_programs']}", flush=True)
        return dict(stage=stage.name, status="not_available", **available)
    identity = snapshot(stage)
    receipt = receipts / (stage.name + ".json")
    if not always and receipt.is_file():
        old = json.loads(receipt.read_text())
        if old.get("identity") == identity and old.get("outputs") and all(state(path) == saved for path, saved in old["outputs"].items()):
            completed_result(stage)
            print(f"[reuse] {stage.name}: input/parameter/output content identities match", flush=True)
            return dict(stage=stage.name, status="reused", **old.get("result", {}))
    print(f"[stage] {stage.name}: {shlex.join(stage.command)}", flush=True)
    start = now()
    try:
        subprocess.run(stage.command, check=True)
        result = completed_result(stage)
        if snapshot(stage) != identity:
            raise RuntimeError(f"{stage.name}: inputs changed during execution")
        if stage.name != "report":
            produced = outputs(stage.out)
            if not produced:
                raise RuntimeError(f"{stage.name}: no outputs produced")
            save(receipt, dict(identity=identity, outputs=produced, result=result, started=start, completed=now()))
        return dict(stage=stage.name, status="complete", started=start, completed=now(), **result)
    except subprocess.CalledProcessError as exc:
        return dict(stage=stage.name, status="failed", exit_code=exc.returncode, reason="stage command failed; see complete log", started=start, completed=now())
    except (OSError, ValueError, RuntimeError) as exc:
        return dict(stage=stage.name, status="failed", reason=str(exc), started=start, completed=now())


def prepare_metadata_if_needed(paths, receipts):
    if paths["metadata_mode"] != "prepare_from_original_supplement":
        return
    identity = dict(schema=VERSION, input=[state(str(paths["prefix"]) + ".psam"), state(paths["supplement"]), state(SOURCE / "evolve.py")])
    receipt = receipts / "metadata.json"
    if receipt.is_file() and paths["metadata"].is_file():
        old = json.loads(receipt.read_text())
        if old.get("identity") == identity and old.get("output") == state(paths["metadata"]):
            return
    from evolve import prepare_metadata
    print("[metadata] Preparing canonical sample table from the original supplement into the new output", flush=True)
    prepare_metadata(str(paths["prefix"]) + ".psam", paths["supplement"], paths["metadata"].parent)
    save(receipt, dict(identity=identity, output=state(paths["metadata"])))


def overall(records):
    if any(r["status"] == "failed" for r in records):
        return "failed"
    if any(r["status"] == "not_available" for r in records):
        return "partial_missing_inputs"
    return "complete_requested_stages"


def execute(args):
    paths = resolve(args)
    names = stage_selection(args)
    configuration = {k:str(v) if isinstance(v, Path) else v for k,v in paths.items()}
    configuration.update(stages=names, version=VERSION, reuse="size/mtime and parameters; no added large-file checksum")
    if args.dry_run or args.preflight:
        stages = plan(args, paths, names)
        details = []
        for stage in stages:
            details.append(dict(stage=stage.name, command=shlex.join(stage.command), **availability(stage),
                                **({'details': stage.details} if stage.details is not None else {})))
        print(json.dumps(dict(status="plan_only", configuration=configuration, stages=details,
            notes=["No files written or analyses started", "Run the integrated pipeline after the scan and PCA inputs have completed; primary reads the selected current scan counters", "Full allele catalog and small validation targets generated earlier in this plan may be absent at preflight", "Time uses the full allele catalog and all eligible ancient samples, not the small genotype-validation panel", "Existing shared PCA is reused; absent coordinates do not trigger PCA fitting", "Actual execution through allele3.sh retains the systemd/cgroup memory cap"]), indent=2))
        return 0
    out = paths["out"]
    out.mkdir(parents=True, exist_ok=True)
    receipts = out / ".receipts"
    receipts.mkdir(exist_ok=True)
    import fcntl
    with (out / ".pipeline.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another focused pipeline is using {out}") from exc
        records = []
        status_path = out / "pipeline_status.json"
        for name in names:
            # Re-resolve generated candidate files and their exact chromosomes.
            stage = plan(args, paths, [name])[0]
            candidate_status = next((r["status"] for r in records if r["stage"] == "candidates"), None)
            if name == "genotypes" and not args.targets and candidate_status not in {None, "complete", "reused"}:
                records.append(dict(stage=name, status="not_available", reason="candidates stage did not complete; stale candidate output is not reused"))
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                continue
            catalog_status = next((r["status"] for r in records if r["stage"] == "alleles"), None)
            if name in {"time", "selection"} and not args.allele_catalog and catalog_status not in {None, "complete", "reused"}:
                records.append(dict(stage=name, status="not_available", reason="full allele catalog stage did not complete; stale catalog is not reused"))
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                continue
            if name == "time" and not any(availability(stage).values()):
                try:
                    prepare_metadata_if_needed(paths, receipts)
                    if paths["metadata"] not in stage.inputs:
                        stage.inputs.append(paths["metadata"])
                except (OSError, ValueError, RuntimeError) as exc:
                    records.append(dict(stage="time", status="failed", reason="metadata adapter: " + str(exc)))
                    save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                    continue
            if name == "report":
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
            # Archaic checks its persistent per-cohort caches itself, including
            # their inventories, even when the aggregate table is unchanged.
            result = run_stage(stage, receipts, always=name in {"report", "archaic"})
            records.append(result)
            save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
        final = overall(records)
        print(f"[herv-focus] {final}: {status_path}", flush=True)
        return 2 if final == "failed" else 3 if final == "partial_missing_inputs" else 0


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", nargs="?", choices=("pipeline", *STAGES), default="pipeline")
    p.add_argument("--stages", help="Comma-separated stage subset; only requested stages run")
    p.add_argument("--dir-analysis-root", default=os.environ.get("ALLELE_ANALYSIS_ROOT", "/mnt/d/analysis/allele3"))
    p.add_argument("--dir-modern", "--modern-dir", "--scan-dir", default=os.environ.get("ALLELE_SCAN_STORE", SCAN_STORE))
    p.add_argument("--dir-out", "--out", default=os.environ.get("ALLELE_HERV_FOCUS_OUT"),
                   help="Output directory; the analysis root expands to ROOT/final/herv_focus_20261028")
    p.add_argument("--dir-ancient", default=os.environ.get("ALLELE_ANCIENT_DIR", "/mnt/f/gen/ancient"))
    p.add_argument("--erv-joined", default=os.environ.get("ALLELE_ERV_JOINED", "/mnt/f/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz"))
    p.add_argument("--fasta38", default=os.environ.get("ALLELE_FASTA", "/mnt/f/gen/1hgp/38/GRCH38.fasta.gz"))
    p.add_argument("--fasta37", default=os.environ.get("ALLELE_FASTA37", "/mnt/f/gen/1hgp/37/GRCH37.fasta.gz"))
    p.add_argument("--chain38to37", default=os.environ.get("ALLELE_CHAIN38TO37", "/mnt/f/gen/chains/hg38ToHg19.over.chain.gz"))
    p.add_argument("--chain37to38", default=os.environ.get("ALLELE_CHAIN37TO38", "/mnt/f/gen/chains/hg19ToHg38.over.chain.gz"))
    p.add_argument("--ld-block-dir", default=os.environ.get("ALLELE_LD_BLOCK_DIR", "/mnt/f/LD/block"))
    for name in ("pgen-prefix", "hq-variants", "metadata", "coordinates", "supplement-xlsx", "targets", "allele-catalog", "selection-stats", "ld-covariates", "chrom-sizes", "chromosomes"):
        p.add_argument("--" + name)
    for name in ("basis-dir", "coordinate-receipt", "provenance-manifest", "pca-temporal-dir"):
        p.add_argument("--" + name)
    p.add_argument("--selection-effect-allele", choices=("raw_only", "ALT"), default="raw_only", help="Keep published signs raw unless the source's ALT effect contract has been verified")
    for name in ("segdup", "lcr", "low-mappability"):
        p.add_argument("--" + name + "-bed")
    p.add_argument("--kg-vcf-dir", default=os.environ.get("ALLELE_1KG_VCF_DIR", "/mnt/f/gen/1kg/38/vcf"))
    p.add_argument("--kg-sample-panel", default=os.environ.get("ALLELE_1KG_PANEL", "/mnt/f/gen/1kg/38/samples.txt"))
    p.add_argument("--kg-vcf-pattern", default="chr{chrom}.vcf.gz")
    p.add_argument("--kg-min-ac", type=int, default=2)
    p.add_argument('--archaic-root', default=os.environ.get('ALLELE_ARCHAIC_ROOT', '/mnt/f/analysis/gu'))
    p.add_argument('--archaic-datasets', default=os.environ.get('ALLELE_ARCHAIC_DATASETS', '1kg'), help='Comma-separated names/globs [default: %(default)s]; auto selects all complete 1kg/ukbNNN datasets')
    p.add_argument('--archaic-module', choices=('auto', 'ibdmix', 'ibdmix-GPU'), default='auto')
    p.add_argument('--archaic-chromosomes', default='auto', help='auto includes available X/Y separately; or autosomes,X,Y / X,Y / X / Y')
    p.add_argument('--archaic-cache', default=os.environ.get('ALLELE_ARCHAIC_CACHE'), help='Persistent per-cohort products; default SCAN-DIR-archaic')
    p.add_argument('--archaic-liftover', default='/mnt/d/software/bin/liftOver')
    p.add_argument('--archaic-min-lod', type=float, default=4.)
    p.add_argument('--archaic-min-length', type=int, default=50000)
    p.add_argument('--archaic-carrier-thresholds', default='0.01,0.05,0.10,0.20')
    p.add_argument('--archaic-min-training-sites', type=int, default=1)
    p.add_argument('--archaic-seed', type=int, default=20261024)
    for name, default in (("window",200),("step",50),("flank",2000),("bootstrap",2000),("min-training-sites",100),
                          ("max-sites",500),("controls-per-case",4),("bin-years",1000),("min-samples",20),
                          ("allele-max-sites",0),("allele-controls-per-case",1),("max-age-bp",10000),
                          ("time-chunk-size",64),("time-pcs",10),("ld-bootstrap",1000),("end-bp",200),("body-bins",20),("seed",20261027)):
        p.add_argument("--" + name, type=int, default=default)
    p.add_argument("--min-support", type=float, default=.99)
    p.add_argument("--min-call-rate", type=float, default=.9)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--preflight", action="store_true")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        return execute(args)
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"herv-pipeline: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
