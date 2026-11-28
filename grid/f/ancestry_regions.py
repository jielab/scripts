"""Build population-specific regional scores on one shared participant roster."""
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import hashlib
import json
import shutil

import numpy as np
import pandas as pd

from consensus3 import read_consensus
from core2 import POPS, atomic_json


def configure_beds(cfg, args):
    """Resolve BED files once at submission. Bind four files by ancestry name, never list order."""
    import os
    import re
    inputs, regions = cfg["inputs"], cfg["regions"]
    build = str(args.ld_bed_build or regions.get("bed_build",38))
    layout = args.ld_layout or regions.get("layout","consensus")
    layout = "consensus" if layout=="uniform" else layout
    names = [x.strip() for x in args.ld_bed.split(",")] if args.ld_bed is not None else None
    if names is not None:
        if len(names) not in (1,4) or any(not x for x in names):
            raise ValueError("--ld-bed needs one filename or four ancestry filenames")
        inferred = "consensus" if len(names)==1 else "ancestry"
        if args.ld_layout and layout!=inferred:
            raise ValueError("--ld-layout conflicts with the number of BED files")
        layout = inferred
    elif args.ld_block_dir and not args.ld_layout:
        layout = "ancestry"
    rep = args.representation or regions["representation"]
    if layout not in {"consensus","ancestry"} or rep not in {"ancestry_channels","meta"}:
        raise ValueError("Unknown regional layout or representation")
    if (layout=="consensus" and build!="38") or (layout=="ancestry" and rep!="ancestry_channels"):
        raise ValueError("Shared BEDs use GRCh38; population BEDs require four CSx channels")
    directory = inputs["consensus_dir"] if layout=="consensus" else inputs["ld_block_dir"]
    if directory=="auto":
        directory = Path(inputs["consensus_dir"]).parent.parent/"block"
    directory = Path(args.dir_ld_bed or directory).expanduser()
    if not directory.is_absolute() and directory.parts[:2]==("data","LD") and os.environ.get("DIR_LD"):
        directory = Path(os.environ["DIR_LD"]).joinpath(*directory.parts[2:])
    directory = directory.resolve()
    if layout=="consensus":
        # An explicitly supplied catalog can also be a directory of chr*.bed.
        source = directory if args.consensus_dir and names is None and not args.dir_ld_bed else directory/(names or ["uni.bed"])[0]
        files = {"shared": str(source.resolve())}
    else:
        files = {}
        for filename in names or [f"{p}.{build}.bed" for p in POPS]:
            match = re.fullmatch(r"(EUR|AFR|EAS|SAS)\.(37|38)\.bed",Path(filename).name,re.I)
            if not match or match[2]!=build or match[1].upper() in files:
                raise ValueError(f"Expected one distinct <ancestry>.{build}.bed for each of EUR, AFR, EAS, SAS: {filename}")
            files[match[1].upper()] = str((directory/filename).resolve())
        files = {pop: files[pop] for pop in POPS}
    regions.update(layout=layout, representation=rep, bed_build=int(build), bed_files=files)


def bed_paths(cfg):
    if "bed_files" in cfg["regions"]:
        return {k: Path(v) for k,v in cfg["regions"]["bed_files"].items()}
    inputs = cfg["inputs"]
    if cfg["regions"].get("layout", "consensus") == "consensus":
        return {"shared": inputs["consensus_dir"]}
    directory = Path(inputs["ld_block_dir"])
    return {pop: directory / f"{pop}.{cfg['regions'].get('bed_build',38)}.bed" for pop in POPS}


def read_maps(cfg):
    return {name: read_consensus(path, expected_chromosomes=cfg["chromosomes"],
                                 build="GRCh38" if name == "shared" else cfg["regions"].get("bed_build",38),
                                 require_sorted=name == "shared")
            for name, path in bed_paths(cfg).items()}


def score(cfg, roster_path, weights, output, phi):
    from genomics import read_roster
    from scoring_parallel import regional_scores

    regions, inputs, runtime = cfg["regions"], cfg["inputs"], cfg["runtime"]
    maps = read_maps(cfg)
    paths = bed_paths(cfg)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    options = dict(plink=runtime["plink"], threads=runtime["threads"],
                   memory_mb=runtime["plink_memory_mb"], batch_size=cfg["scoring"]["batch_size"])
    common = "shared" in maps
    build = "38" if common else str(regions.get("bed_build",38))
    options["block_build"] = "GRCh"+build
    if regions["representation"]=="ancestry_channels" and set(weights) != set(POPS):
        raise ValueError("Four-channel regional scoring requires exactly EUR, AFR, EAS and SAS posterior channels")
    for name, (_, audit) in maps.items():
        for file in audit["files"]:
            file["sha256"] = hashlib.sha256(Path(file["path"]).read_bytes()).hexdigest()
    provenance = dict(layout=regions.get("layout", "consensus"), representation=regions["representation"],
                      geometry_phi=phi, block_build="GRCh"+build,
                      assignment_position="REF_BP37" if build=="37" else "BP38", shared_boundaries=common, block_files={
                          name: audit["files"] for name, (_, audit) in maps.items()})
    if common:
        regional_scores(roster_path, weights, paths["shared"], inputs["genotype_dir"], output, **options)
        meta = json.loads((output / "metadata.json").read_text())
        atomic_json({**meta, **provenance}, output / "metadata.json")
        return

    roster = read_roster(roster_path)
    width = sum(len(blocks) for blocks, _ in maps.values())
    # Publish metadata last. Each ancestry scratch matrix is discarded after copying.
    (output / "metadata.json").unlink(missing_ok=True)
    matrix = np.lib.format.open_memmap(output / "scores.npy", mode="w+", dtype="float32",
                                        shape=(len(roster), width))
    totals = roster[["id"]].copy()
    cols, catalogs, counts, checks = [], [], [], []
    offset = 0
    with gzip.open(output / "unassigned_variants.tsv.gz", "wt") as missing:
        for pop, (blocks, _) in maps.items():
            print(f"[regions] {pop}  {len(blocks)} blocks  {paths[pop]}", flush=True)
            with TemporaryDirectory(prefix=f"score-{pop}-", dir=output) as tmp:
                part = Path(tmp)
                normalized = part / "blocks.bed"
                blocks[[s+build for s in ("chrom","start","end")]].to_csv(normalized, sep="\t", header=False, index=False)
                regional_scores(roster_path, {pop: weights[pop]}, normalized,
                                inputs["genotype_dir"], part, **options)
                def table(name):
                    return pd.read_csv(part / name, sep="\t", dtype={"id": str})
                ids = table("ids.tsv")
                total = table("covered_global.tsv.gz")
                if not ids.id.equals(roster.id) or not total.id.equals(roster.id):
                    raise ValueError(f"{pop} regional scores disagree with the frozen roster")
                features = table("features.tsv")
                if not features.block_id.equals(blocks.block_id) or not features["map"].eq(pop).all():
                    raise ValueError(f"{pop} regional feature catalog mismatch")
                features["source_file"] = blocks.source_file
                features["source_line"] = blocks.source_line
                features["feature_index"] = np.arange(offset, offset+len(blocks))
                cols.append(features)
                catalogs.append(blocks.assign(map=pop, feature_index=features.feature_index))
                counts.append(table("region_counts.tsv"))
                checks.append(table("sum_check.tsv"))
                totals[pop] = total[pop]
                part_matrix = np.load(part / "scores.npy", mmap_mode="r")
                if part_matrix.shape != (len(roster), len(blocks)):
                    raise ValueError(f"Unexpected regional matrix shape for {pop}")
                for start in range(0, len(roster), 4096):
                    matrix[start:start+4096, offset:offset+len(blocks)] = part_matrix[start:start+4096]
                del part_matrix
                with gzip.open(part / "unassigned_variants.tsv.gz", "rt") as stream:
                    header = next(stream)
                    if offset == 0:
                        missing.write(header)
                    shutil.copyfileobj(stream, missing)
                offset += len(blocks)
    matrix.flush()
    del matrix
    from core2 import atomic_table
    for frames, filename in ((cols, "features.tsv"), (catalogs, "catalog.tsv"),
                             (counts, "region_counts.tsv"), (checks, "sum_check.tsv")):
        atomic_table(pd.concat(frames, ignore_index=True), output / filename)
    atomic_table(roster[["id"]], output / "ids.tsv")
    atomic_table(totals, output / "covered_global.tsv.gz")
    atomic_json(dict(**provenance, target_build="GRCh38", bed_coordinates="0-based-half-open",
                     channels=list(maps), shape=[len(roster), width], dtype="float32",
                     genotype_source=str(inputs["genotype_dir"]), frequency_fit_role="train",
                     score_type="posterior_dosage_sum", n_blocks_by_map={p: len(b) for p,(b,_) in maps.items()},
                     n_empty_features=sum(int(f.empty_block.sum()) for f in cols),
                     n_unassigned_variants=sum(int(c.n_unassigned.sum()) for c in counts),
                     additivity="checked per channel and chromosome on its covered SNPs"), output / "metadata.json")
