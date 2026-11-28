"""Shared preparation and validated publication for independent chromosome jobs.

Only preparation writes normalized GWAS/target metadata. Each chromosome owns
its partition directory; only merge publishes the genome-wide weight manifest.
"""
from __future__ import annotations

from copy import copy
import json
import filecmp
import shutil
import os
from pathlib import Path

import pandas as pd


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def prep_settings(args, source_ld_ref=None):
    return dict(trait=args.trait,
                manifest=str(Path(args.manifest).resolve()),
                genotype_dir=str(Path(args.genotype_dir).resolve()),
                ld_ref=str(Path(source_ld_ref or args.ld_ref).resolve()),
                rsid_map=str(Path(args.rsid_map).resolve()) if args.rsid_map else None,
                allow_ambiguous=args.allow_ambiguous,
                strand_aligned_provenance=args.strand_aligned_provenance)


def method_settings(args):
    from infer import parse_phi_grid
    return dict(trait=args.trait, phi_grid=parse_phi_grid(args.phi_grid),
                cs_fixed=args.cs_fixed, n_iter=args.n_iter,
                n_burnin=args.n_burnin, thin=args.thin, seed=args.seed,
                backend=args.backend, threads=args.threads)


def verify_files(records, label):
    from infer import file_stat
    for recorded in records:
        try:
            current = file_stat(recorded["path"])
        except OSError as exc:
            raise ValueError(f"{label}: missing file {recorded['path']}") from exc
        if current != recorded:
            raise ValueError(f"{label}: changed file {recorded['path']}")


def prepared_outputs(out, chromosomes):
    from infer import POPULATIONS
    out = Path(out)
    paths = [out/"prepared.json", out/"target/manifest.tsv", out/"target/mapping_counts.tsv"]
    paths += [out/"target"/f"chr{c}.{ext}" for c in chromosomes for ext in ("tsv", "bim")]
    paths += [out/"normalized"/p/f"chr{c}.tsv" for p in POPULATIONS for c in chromosomes]
    # The preparation receipt binds canonical reference links and SNP metadata.
    paths += sorted((out/"singlepop_reference").glob("*/snpinfo_*"))
    return paths


def write_prepared(args, chromosomes, source_ld_ref, records, kind,
                   reference_file, normalized, normalization_audit, singlepop):
    from infer import HERE, POPULATIONS, file_stat
    from genomics import pfile
    out = Path(args.output).resolve()
    inputs = [Path(args.manifest), reference_file, HERE/"infer.py", HERE/"genomics.py", Path(__file__)]
    if args.rsid_map:
        inputs.append(Path(args.rsid_map))
    for pop in POPULATIONS:
        inputs.extend(Path(p) for p in records[pop]["files"])
        inputs.extend(Path(args.ld_ref)/f"ldblk_{kind}_{pop.lower()}"/f"ldblk_{kind}_chr{c}.hdf5"
                      for c in chromosomes)
    inputs.extend(pfile(Path(args.genotype_dir)/f"chr{c}")[1] for c in chromosomes)
    outputs = prepared_outputs(out, chromosomes)[1:]
    # Include symlink targets so the receipt detects changed LD references.
    outputs.extend(p for p in (out/"singlepop_reference").glob("*/*") if p.is_file())
    result = dict(schema="grid-csx-prepared-v1", chromosomes=chromosomes,
                  settings=prep_settings(args, source_ld_ref),
                  inputs=[file_stat(p) for p in dict.fromkeys(inputs)],
                  outputs=[file_stat(p) for p in dict.fromkeys(outputs)],
                  records=records, kind=kind, ld_root=str(args.ld_ref),
                  reference_file=str(reference_file),
                  normalized={pop: {str(c): str(p) for c, p in paths.items()}
                              for pop, paths in normalized.items()},
                  normalization_audit=normalization_audit,
                  singlepop={p: str(path) for p, path in singlepop.items()})
    atomic_json(out/"prepared.json", result)
    return out/"prepared.json"


def load_prepared(args, chromosomes, *, load_maps=True):
    directory = Path(args.prepared_inputs)
    path = directory/"prepared.json"
    if not path.is_file():
        raise ValueError(f"Missing shared preparation: {path}; run infer preparation first")
    data = json.loads(path.read_text())
    if data.get("schema") != "grid-csx-prepared-v1" or data["settings"] != prep_settings(args):
        raise ValueError("Shared inference preparation settings changed; rerun preparation")
    if not set(chromosomes) <= set(data["chromosomes"]):
        raise ValueError("Requested chromosome is absent from the shared preparation")
    verify_files(data["inputs"], "Shared preparation inputs")
    verify_files(data["outputs"], "Shared preparation outputs")
    data["normalized"] = {p: {int(c): Path(v) for c, v in paths.items()}
                          for p, paths in data["normalized"].items()}
    data["maps"] = ({c: pd.read_csv(directory/"target"/f"chr{c}.tsv", sep="\t")
                     for c in chromosomes} if load_maps else {})
    return data


def write_partition_receipt(args, chromosomes):
    from infer import file_stat
    if len(chromosomes) != 1:
        raise ValueError("An independent chromosome worker must process exactly one chromosome")
    out = Path(args.output)
    files = [out/"provenance.json", out/"weights/manifest.tsv", out/"weights/snp_counts.tsv"]
    files += sorted((out/"weights").glob("*.chr*.tsv"))
    markers = sorted((out/"raw").glob("*/complete.json"))
    expected = json.loads((out/"provenance.json").read_text())["mcmc_jobs_total"]
    if len(markers) != expected:
        raise ValueError("Missing chain receipts; refusing completed chromosome receipt")
    files += markers
    for marker in markers:
        files.extend(Path(p["path"]) for p in json.loads(marker.read_text())["raw_outputs"])
    atomic_json(out/"partition.complete.json", dict(
        schema="grid-csx-partition-v1", chromosomes=chromosomes,
        method=method_settings(args),
        preparation=file_stat(Path(args.prepared_inputs)/"prepared.json"),
        outputs=[file_stat(p) for p in files]))


def reuse_partition(args, chromosomes):
    """Preserve valid exports and their timestamps on a no-change resume."""
    from infer import file_stat
    path = Path(args.output)/"partition.complete.json"
    if not path.is_file():
        return False
    try:
        load_prepared(args, chromosomes, load_maps=False)
        data = json.loads(path.read_text())
        if (data.get("schema") != "grid-csx-partition-v1" or data["chromosomes"] != chromosomes
                or data["method"] != method_settings(args)
                or data["preparation"] != file_stat(Path(args.prepared_inputs)/"prepared.json")):
            return False
        verify_files(data["outputs"], "Chromosome outputs")
        for record in data["outputs"]:
            if Path(record["path"]).name == "complete.json":
                marker = json.loads(Path(record["path"]).read_text())
                verify_files(marker["signature"]["original_inputs"], "Sampler inputs/code")
        return True
    except (ValueError, OSError, KeyError):
        return False


def merge_partitions(args):
    publication = Path(args.output)/"weights/manifest.tsv"
    try:
        return _merge_partitions(args)
    except Exception:
        publication.unlink(missing_ok=True)
        raise


def _merge_partitions(args):
    from infer import file_stat, model_manifest, parse_phi_grid
    from genomics import parse_chromosomes
    out = Path(args.output).resolve()
    weights = out/"weights"
    weights.mkdir(parents=True, exist_ok=True)
    publication = weights/"manifest.tsv"
    chromosomes = parse_chromosomes(args.chrs)
    shared_args = copy(args)
    shared_args.prepared_inputs = str(out)
    shared = load_prepared(shared_args, chromosomes, load_maps=False)
    if shared["chromosomes"] != chromosomes:
        raise ValueError("Merge must cover exactly the prepared chromosome set")
    models = model_manifest(out, parse_phi_grid(args.phi_grid), args.cs_fixed)
    preparation = file_stat(out/"prepared.json")
    counts, provenances, sources = [], [], []
    chains_per_chromosome = len(parse_phi_grid(args.phi_grid)) + 1
    chains_per_chromosome += 4 * (len(parse_phi_grid(args.phi_grid)) + 1 if args.cs_fixed else 1)
    for chrom in chromosomes:
        part = out/"chromosomes"/f"chr{chrom}"
        receipt = part/"partition.complete.json"
        if not receipt.is_file():
            raise ValueError(f"chr{chrom}: missing successful partition receipt; no genome-wide publication")
        data = json.loads(receipt.read_text())
        if (data.get("schema") != "grid-csx-partition-v1" or data["chromosomes"] != [chrom]
                or data["method"] != method_settings(args) or data["preparation"] != preparation):
            raise ValueError(f"chr{chrom}: partition belongs to different inputs/settings")
        verify_files(data["outputs"], f"chr{chrom} partition")
        for record in data["outputs"]:
            if Path(record["path"]).name == "complete.json":
                marker = json.loads(Path(record["path"]).read_text())
                verify_files(marker["signature"]["original_inputs"], f"chr{chrom} sampler inputs/code")
        recorded = {p["path"] for p in data["outputs"]}
        provenance = json.loads((part/"provenance.json").read_text())
        if provenance["chromosomes"] != [chrom] or provenance["mcmc_jobs_total"] != chains_per_chromosome:
            raise ValueError(f"chr{chrom}: incomplete model set")
        manifest = pd.read_csv(part/"weights/manifest.tsv", sep="\t", dtype=str)
        expected = models.drop(columns="weight_pattern").astype(str)
        if not manifest.drop(columns="weight_pattern").equals(expected):
            raise ValueError(f"chr{chrom}: model manifest differs from requested models")
        frame = pd.read_csv(part/"weights/snp_counts.tsv", sep="\t")
        if (len(frame) != len(models) or set(frame.score_name) != set(models.score_name)
                or not frame.chromosome.eq(chrom).all() or not frame.n_variants.gt(0).all()):
            raise ValueError(f"chr{chrom}: incomplete variant counts")
        for name in models.score_name:
            path = part/"weights"/f"{name}.chr{chrom}.tsv"
            if str(path.resolve()) not in recorded or not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"chr{chrom}: missing validated weights for {name}")
            sources.append(path)
        counts.append(frame)
        provenances.append(provenance)
    provenance = dict(provenances[0], chromosomes=chromosomes,
                      execution_layout="independent-chromosome-jobs",
                      chromosome_workers=len(chromosomes),
                      mcmc_jobs_total=sum(p["mcmc_jobs_total"] for p in provenances),
                      mcmc_jobs_reused=sum(p["mcmc_jobs_reused"] for p in provenances),
                      preparation=preparation)
    count_text = pd.concat(counts, ignore_index=True).sort_values(["score_name", "chromosome"]).to_csv(sep="\t", index=False)
    model_text = models.to_csv(sep="\t", index=False)
    count_path, provenance_path = weights/"snp_counts.tsv", out/"provenance.json"
    if (publication.is_file() and publication.read_text() == model_text
            and count_path.is_file() and count_path.read_text() == count_text
            and provenance_path.is_file() and json.loads(provenance_path.read_text()) == provenance
            and {p.name for p in weights.glob("*.chr*.tsv")} == {p.name for p in sources}
            and all((weights/p.name).is_file() and not (weights/p.name).is_symlink()
                    and filecmp.cmp(weights/p.name, p, shallow=False) for p in sources)):
        return publication
    # Withdraw obsolete wildcard matches, then publish only the validated set.
    publication.unlink(missing_ok=True)
    for path in weights.glob("*.chr*.tsv"):
        path.unlink()
    for source in sources:
        shutil.copy2(source, weights/source.name)
    count_path.write_text(count_text)
    atomic_json(out/"provenance.json", provenance)
    temporary = publication.with_suffix(".tmp")
    temporary.write_text(model_text)
    temporary.replace(publication)
    return publication
