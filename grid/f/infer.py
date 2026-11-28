#!/usr/bin/env python3
"""Infer joint PRS-CSx and official single-population PRS-CS posterior weights.

Reference/GWAS rsIDs join actual GRCh38 PVAR metadata. A metadata-only BIM passes
alleles to the upstream samplers; its positions are actual BP38, while the
samplers' reference BP37 is retained separately in each exported weight file.
No individual target phenotype or target genotype dosage enters inference.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import ExitStack
import glob
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

from genomics import (align_weights, allele_alignment, map_target_variants,
                      parse_chromosomes, pfile, text_input)

POPULATIONS = ("AFR", "EAS", "EUR", "SAS")
HERE = Path(__file__).resolve().parent
ALIASES = {
    "SNP": ("RSID", "RS_ID", "SNP", "MARKERNAME", "ID", "VARIANT_ID"),
    "A1": ("EA", "A1", "EFFECT_ALLELE", "ALLELE1", "ALT"),
    "A2": ("NEA", "A2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "ALLELE0", "ALLELE2", "REF"),
    "BETA": ("BETA", "EFFECT", "ESTIMATE", "LOG_ODDS", "LOGOR"),
    "OR": ("OR", "ODDS_RATIO", "ODDSRATIO"),
    "SE": ("SE", "SEBETA", "SE_BETA", "STANDARD_ERROR", "STDERR"),
    "P": ("P", "PVAL", "P_VALUE", "PVALUE", "P_BOLT_LMM", "PVAL_NOMINAL"),
}


def phi_key(value):
    number = float(value)
    if not np.isfinite(number) or number <= 0:
        raise ValueError("Fixed phi values must be finite and positive")
    mantissa, exponent = f"{number:.12e}".split("e")
    return mantissa.rstrip("0").rstrip(".")+"e"+str(int(exponent))


def parse_phi_grid(value):
    result = [phi_key(x.strip()) for x in value.split(",") if x.strip()]
    if not result or len(result) != len(set(result)):
        raise ValueError("Provide unique positive fixed phi values")
    return result


def effective_n(row):
    value = row.get("n_gwas", "")
    if pd.notna(value) and str(value).strip():
        number = float(value)
        source = "manifest_n_gwas"
    else:
        cases, controls = row.get("n_case", row.get("ncase", "")), row.get("n_control", row.get("ncontrol", ""))
        if not str(cases).strip() or not str(controls).strip():
            raise ValueError("n_gwas is required; binary traits may supply explicit n_case and n_control")
        cases, controls = float(cases), float(controls)
        if not np.isfinite([cases, controls]).all() or min(cases, controls) <= 0:
            raise ValueError("Both case and control counts must be positive")
        number, source = 4/(1/cases+1/controls), "4/(1/n_case+1/n_control)"
    if not np.isfinite(number) or number < 2:
        raise ValueError("Discovery GWAS N must be explicit, finite, and at least 2")
    # Both official command-line parsers accept integers. Keep the original value
    # and its <=0.5 rounding difference in provenance.
    return number, int(round(number)), source


def read_manifest(path, trait):
    path = Path(path).resolve()
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame.columns = [x.lower().strip() for x in frame.columns]
    if not {"trait", "population", "path"} <= set(frame):
        raise ValueError("Manifest requires trait,population,path,n_gwas,sample_overlap")
    frame = frame.loc[frame.trait.eq(trait)].copy()
    frame["population"] = frame.population.str.upper()
    if set(frame.population) != set(POPULATIONS) or frame.population.duplicated().any():
        raise ValueError(f"{trait}: exactly one manifest row is required for each of {POPULATIONS}")
    records = {}
    for raw in frame.to_dict("records"):
        files = []
        for pattern in raw["path"].split(";"):
            candidate = Path(os.path.expandvars(os.path.expanduser(pattern.strip())))
            if not candidate.is_absolute():
                candidate = path.parent/candidate
            matched = sorted(glob.glob(str(candidate)))
            if not matched:
                raise FileNotFoundError(candidate)
            files.extend(str(Path(x).resolve()) for x in matched)
        if len(files) != len(set(files)):
            raise ValueError("A GWAS shard appears more than once in the manifest")
        n_original, n_used, n_source = effective_n(raw)
        overlap = raw.get("sample_overlap", "unknown").lower().strip() or "unknown"
        if overlap not in {"no", "yes", "unknown"}:
            raise ValueError("sample_overlap must be no, yes, or unknown")
        records[raw["population"]] = dict(population=raw["population"], files=files,
            n_gwas=n_used, n_gwas_original=n_original, n_source=n_source,
            sample_overlap=overlap, source_build=raw.get("build", "unknown"))
    return records


def select_ld_reference(directory, chromosomes):
    """Identify the unique complete HM3 panel for all requested chromosomes."""
    directory = Path(directory).resolve()
    panels = []
    for kind in ("1kg", "ygbb"):
        metadata = directory / f"snpinfo_mult_{kind}_hm3"
        blocks = {}
        for pop in POPULATIONS:
            candidates = [directory / f"ldblk_{kind}_{suffix}" for suffix in (pop.lower(), pop.upper())]
            valid = [p for p in candidates if all((p / f"ldblk_{kind}_chr{c}.hdf5").is_file() for c in chromosomes)]
            if len(valid) == 1:
                blocks[pop] = valid[0]
        if metadata.is_file() and len(blocks) == len(POPULATIONS):
            panels.append((kind, metadata, blocks))
    if len(panels) != 1:
        raise ValueError("LD root must contain exactly one complete panel for the requested populations and chromosomes")
    return panels[0]


def canonical_ld_reference(directory, output, chromosomes):
    """Expose the selected panel under the lowercase paths required by PRS-CSx."""
    kind, metadata, blocks = select_ld_reference(directory, chromosomes)
    output = Path(output).resolve() / kind
    output.mkdir(parents=True, exist_ok=True)
    for name, source in [(metadata.name, metadata), *[(f"ldblk_{kind}_{pop.lower()}", p) for pop, p in blocks.items()]]:
        link = output / name
        if link.is_symlink() and link.resolve() == source.resolve():
            continue
        if link.exists() or link.is_symlink():
            raise ValueError(f"Conflicting LD reference path: {link}")
        link.symlink_to(source, target_is_directory=source.is_dir())
    return str(output)


def reference_table(directory, chromosomes):
    directory = Path(directory).resolve()
    candidates = [(kind, directory/f"snpinfo_mult_{kind}_hm3") for kind in ("1kg", "ygbb")]
    found = [(kind, path) for kind, path in candidates if path.is_file()]
    if len(found) != 1:
        raise ValueError("LD root must contain exactly one snpinfo_mult_1kg_hm3 or snpinfo_mult_ygbb_hm3")
    kind, path = found[0]
    d = pd.read_csv(path, sep=r"\s+", dtype={"SNP":str})
    required = {"CHR", "SNP", "BP", "A1", "A2"}
    required |= {key+"_"+pop for pop in POPULATIONS for key in ("FRQ", "FLP")}
    if not required <= set(d):
        raise ValueError("Multi-population LD SNP metadata is incomplete")
    d = d.loc[d.CHR.isin(chromosomes)].copy()
    d["SNP"] = d.SNP.str.lower()
    if d.empty or d.SNP.duplicated().any():
        raise ValueError("Reference SNP IDs must be nonempty and unique")
    d["A1"], d["A2"] = d.A1.str.upper(), d.A2.str.upper()
    # Published HM3 LD metadata is GRCh37. These coordinates are not used to
    # locate target genotypes or assign GRCh38 regional PRS features.
    d = d.rename(columns={"BP":"REF_BP37"})
    for pop in POPULATIONS:
        valid = d["FRQ_"+pop].between(0, 1, inclusive="neither")
        if not d.loc[valid, "FLP_"+pop].isin([-1, 1]).all():
            raise ValueError(f"Invalid native LD orientation for {pop}")
        for chrom in chromosomes:
            ld = directory/f"ldblk_{kind}_{pop.lower()}"/f"ldblk_{kind}_chr{chrom}.hdf5"
            if not ld.is_file():
                raise FileNotFoundError(ld)
    return kind, path, d


def build_target_reference(reference, genotype_dir, output_dir, chromosomes,
                           rsid_map=None, allow_ambiguous=False):
    target_dir = Path(output_dir)/"target"; target_dir.mkdir(parents=True, exist_ok=True)
    tables, audits, manifest = {}, [], []
    for chrom in chromosomes:
        _, pvar = pfile(Path(genotype_dir)/f"chr{chrom}")
        ref = reference.loc[reference.CHR.eq(chrom)].set_index("SNP")
        mapped, counts = map_target_variants(pvar, ref, chrom, rsid_map=rsid_map, allow_ambiguous=allow_ambiguous)
        if mapped.empty:
            raise ValueError(f"No eligible reference rsIDs map to chr{chrom} PVAR; provide --rsid-map for coordinate-only IDs")
        path = target_dir/f"chr{chrom}.tsv"
        mapped.to_csv(path, sep="\t", index=False)
        # Metadata only: neither .bed nor a re-encoded target genotype is made.
        bim = mapped.assign(CM=0)[["CHR", "SNP", "CM", "BP38", "ALT38", "REF38"]]
        bim.to_csv(target_dir/f"chr{chrom}.bim", sep="\t", index=False, header=False)
        tables[chrom] = mapped
        audits.append(dict(chromosome=chrom, **counts))
        manifest.append(dict(chromosome=chrom, path=str(path.resolve()), target_build="GRCh38", n_variants=len(mapped)))
    pd.DataFrame(audits).to_csv(target_dir/"mapping_counts.tsv", sep="\t", index=False)
    pd.DataFrame(manifest).to_csv(target_dir/"manifest.tsv", sep="\t", index=False)
    return tables


def header_aliases(columns):
    normalized = {str(c).lstrip("#").strip().upper().replace("-", "_"):c for c in columns}
    out = {}
    for name, aliases in ALIASES.items():
        out[name] = next((normalized[x] for x in aliases if x in normalized), None)
    if any(out[x] is None for x in ("SNP", "A1", "A2")) or not (out["BETA"] or out["OR"]) or not (out["SE"] or out["P"]):
        raise ValueError("GWAS requires rsID, effect/other alleles, BETA (or OR), and SE (or P); columns: "+str(list(columns)))
    return out


def gwas_chunks(path, chunksize=150000):
    with text_input(path) as stream:
        first = stream.readline()
        while first.startswith("##") or not first.strip():
            if not first:
                raise ValueError(f"Empty GWAS shard {path}")
            first = stream.readline()
        delimiter = "," if "," in first and "\t" not in first and len(first.split()) == 1 else r"\s+"
        names = first.strip().split(",") if delimiter == "," else first.split()
        # Explicit names avoid requiring a seekable compressed stream.
        yield from pd.read_csv(stream, sep=delimiter, names=names, header=None,
                               dtype=str, chunksize=chunksize)


def normalize_sumstats(record, reference, mappings, output_dir, allow_ambiguous=False):
    """One pass through raw shards; use exact z evidence without a P threshold."""
    pop = record["population"]
    mapped = pd.concat(mappings.values(), ignore_index=True).set_index("SNP")
    valid_reference = reference.loc[reference["FRQ_"+pop].between(0,1,inclusive="neither"), "SNP"]
    allowed = set(valid_reference) & set(mapped.index)
    lookup = mapped.loc[mapped.index.isin(allowed), ["CHR", "REF38", "ALT38"]].to_dict("index")
    output = Path(output_dir)/"normalized"/pop; output.mkdir(parents=True, exist_ok=True)
    paths = {chrom:output/f"chr{chrom}.tsv" for chrom in mappings}
    seen, counts, by_chrom = {}, Counter(), Counter()
    with ExitStack() as stack:
        handles = {chrom:stack.enter_context(path.open("w")) for chrom,path in paths.items()}
        for handle in handles.values():
            handle.write("SNP\tA1\tA2\tBETA\tSE\n")
        for path in record["files"]:
            for chunk in gwas_chunks(path):
                col = header_aliases(chunk.columns)
                counts["raw_rows"] += len(chunk)
                ids = chunk[col["SNP"]].str.lower()
                retain = ids.isin(allowed)
                chunk = chunk.loc[retain].copy(); chunk["_id"] = ids.loc[retain]
                counts["candidate_rows"] += len(chunk)
                if chunk.empty:
                    continue
                a1 = chunk[col["A1"]].str.upper().to_numpy()
                a2 = chunk[col["A2"]].str.upper().to_numpy()
                raw_effect = pd.to_numeric(chunk[col["BETA"] or col["OR"]], errors="coerce").to_numpy(float)
                with np.errstate(invalid="ignore", divide="ignore"):
                    beta = raw_effect if col["BETA"] else np.log(raw_effect)
                z = np.full(len(chunk), np.nan)
                if col["SE"]:
                    se = pd.to_numeric(chunk[col["SE"]], errors="coerce").to_numpy(float)
                    ok = np.isfinite(beta)&np.isfinite(se)&(se>0)
                    z[ok] = beta[ok]/se[ok]
                if col["P"]:
                    p = pd.to_numeric(chunk[col["P"]], errors="coerce").to_numpy(float)
                    ok = ~np.isfinite(z)&np.isfinite(beta)&np.isfinite(p)&(p>=0)&(p<=1)
                    # Official P-only convention clips zero P to 1e-323.
                    z[ok] = np.sign(beta[ok])*np.abs(norm.ppf(np.maximum(p[ok], 1e-323)/2))
                for snp, effect, other, value in zip(chunk._id, a1, a2, z):
                    target = lookup[snp]
                    if not np.isfinite(value):
                        counts["invalid_effect_or_uncertainty"] += 1; continue
                    orientation = allele_alignment(effect, other, target["REF38"], target["ALT38"], allow_ambiguous)
                    if orientation is None:
                        counts["allele_incompatible_or_ambiguous"] += 1; continue
                    # Keep the standardized marginal evidence on the target ALT
                    # strand. Both samplers use BETA/SE/sqrt(N), so z/1 is exact.
                    value = -float(value) if orientation[0] == target["REF38"] else float(value)
                    if snp in seen:
                        if not math.isclose(seen[snp], value, rel_tol=1e-10, abs_tol=1e-12):
                            raise ValueError(f"Conflicting GWAS shard records for {pop} {snp}; resolve duplicates explicitly")
                        counts["identical_duplicate_rows"] += 1; continue
                    seen[snp] = value
                    chrom = int(target["CHR"])
                    handles[chrom].write(f"{snp}\t{target['ALT38']}\t{target['REF38']}\t{value:.17g}\t1\n")
                    by_chrom[chrom] += 1
                    counts[orientation[1]] += 1
    for chrom in mappings:
        if by_chrom[chrom] == 0:
            raise ValueError(f"{pop} has no eligible GWAS SNPs on chromosome {chrom}")
    return paths, dict(population=pop, retained=len(seen), by_chromosome=dict(by_chrom), **counts)


def prepare_singlepop_reference(reference, ld_root, reference_kind, output_dir, chromosomes):
    """Adapt the *same* native LD HDF5 to official PRS-CS's six-column schema.

FRQ_POP is the frequency of the shared metadata A1; FLP_POP identifies its sign
relative to the population's native LD. A negative FLP swaps A1/A2 and frequency
in the six-column reference. HDF5 matrices are linked unchanged.
"""
    directories = {}
    for pop in POPULATIONS:
        path = Path(output_dir)/"singlepop_reference"/f"ldblk_{reference_kind}_{pop.lower()}"
        path.mkdir(parents=True, exist_ok=True)
        d = reference.loc[reference["FRQ_"+pop].between(0,1,inclusive="neither")].copy()
        flip = d["FLP_"+pop].eq(-1)
        native = pd.DataFrame(dict(CHR=d.CHR, SNP=d.SNP, BP=d.REF_BP37,
            A1=np.where(flip,d.A2,d.A1), A2=np.where(flip,d.A1,d.A2),
            MAF=np.where(flip,1-d["FRQ_"+pop],d["FRQ_"+pop])))
        native.to_csv(path/f"snpinfo_{reference_kind}_hm3", sep="\t", index=False)
        for chrom in chromosomes:
            name = f"ldblk_{reference_kind}_chr{chrom}.hdf5"
            source = (Path(ld_root)/f"ldblk_{reference_kind}_{pop.lower()}"/name).resolve()
            link = path/name
            if link.is_symlink():
                if link.resolve() != source:
                    link.unlink()
            elif link.exists():
                raise ValueError(f"Expected a reference symlink, found an existing file: {link}")
            if not link.exists():
                link.symlink_to(source)
        directories[pop] = path.resolve()
    return directories


def model_manifest(output, phis, cs_fixed):
    rows = []
    def add(name, family, pop, phi):
        rows.append(dict(score_name=name, family=family, population=pop, phi=phi,
            weight_pattern=str((Path(output)/"weights"/f"{name}.chr*.tsv").resolve()),
            reference_build="GRCh37", target_build="GRCh38"))
    for phi in [*phis, "auto"]:
        for pop in POPULATIONS:
            add(f"csx.auto.{pop}" if phi=="auto" else f"csx.phi_{phi}.{pop}", "csx", pop, phi)
        add("csx.auto" if phi=="auto" else f"csx.phi_{phi}.meta", "csx_meta", "META", phi)
    for pop in POPULATIONS:
        add(f"cs.{pop}", "cs", pop, "auto")
        if cs_fixed:
            for phi in phis:
                add(f"cs.phi_{phi}.{pop}", "cs", pop, phi)
    return pd.DataFrame(rows)


def inference_command(family, chrom, phi, output, target_prefix, sumstats, n_gwas,
                      ld_root, trait, population=None, n_iter=4000, n_burnin=2000, thin=5, seed=2026):
    if family == "csx":
        command = [sys.executable, str(HERE.parent/"csx"/"csx.py"), f"--ref_dir={ld_root}",
            f"--bim_prefix={target_prefix}", "--sst_file="+",".join(str(sumstats[p]) for p in POPULATIONS),
            "--n_gwas="+",".join(str(n_gwas[p]) for p in POPULATIONS), "--pop="+",".join(POPULATIONS),
            f"--out_dir={output}", f"--out_name={trait}", "--meta=TRUE"]
    elif family == "cs":
        command = [sys.executable, str(HERE.parent/"cs"/"cs.py"), f"--ref_dir={ld_root}",
            f"--bim_prefix={target_prefix}", f"--sst_file={sumstats[population]}",
            f"--n_gwas={n_gwas[population]}", f"--out_dir={output}", "--beta_std=FALSE", "--write_psi=FALSE"]
    else:
        raise ValueError("Unknown inference family")
    command += [f"--chrom={chrom}", "--a=1", "--b=0.5", f"--n_iter={n_iter}", f"--n_burnin={n_burnin}",
                f"--thin={thin}", "--write_pst=FALSE", f"--seed={seed}"]
    if phi != "auto":
        command.append(f"--phi={float(phi):.12g}")
    return command


def export_posterior(raw, target, output, allow_ambiguous=False):
    d = pd.read_csv(raw, sep=r"\s+", names=["CHR", "SNP", "REF_BP37", "A1", "A2", "BETA"])
    if d.empty:
        raise ValueError(f"Empty posterior: {raw}")
    exported = align_weights(d, target, allow_ambiguous)
    if len(exported) != len(d) or exported.empty:
        raise ValueError(f"Posterior does not completely map to the prepared target metadata: {raw}")
    exported.to_csv(output, sep="\t", index=False)
    return len(exported)



def file_stat(path):
    """Small, transparent resume record; no file-content hashes or genotype read."""
    path = Path(path).resolve()
    state = path.stat()
    return dict(path=str(path), size=state.st_size, mtime_ns=state.st_mtime_ns)


def can_resume(marker, signature, raw_outputs):
    try:
        recorded = json.loads(Path(marker).read_text())
        return (recorded.get("signature") == signature and
                recorded.get("raw_outputs") == [file_stat(p) for p in raw_outputs] and
                all(Path(p).stat().st_size > 0 for p in raw_outputs))
    except (OSError, ValueError, TypeError):
        return False


def write_job_marker(path, signature, raw_outputs, counts):
    path = Path(path)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(dict(signature=signature,
        raw_outputs=[file_stat(p) for p in raw_outputs], counts=counts), indent=2)+"\n")
    temporary.replace(path)


def run_inference(args):
    chromosomes = parse_chromosomes(args.chrs)
    phis = parse_phi_grid(args.phi_grid)
    if args.n_iter <= args.n_burnin or args.n_burnin < 0 or args.thin < 1 or args.n_iter%args.thin or args.n_burnin%args.thin or (args.n_iter-args.n_burnin)//max(args.thin,1) < 2:
        raise ValueError("MCMC iterations and burnin must be multiples of thin, with retained posterior draws")
    if args.jobs < 1 or args.threads < 1:
        raise ValueError("--jobs and --threads must be positive")
    if args.allow_ambiguous and not args.strand_aligned_provenance:
        raise ValueError("--allow-ambiguous requires --strand-aligned-provenance explaining independent allele/strand alignment")
    if getattr(args, "merge_only", False):
        from csx_parallel import merge_partitions
        return merge_partitions(args)
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=True)
    weights_dir = out/"weights"; weights_dir.mkdir(exist_ok=True)
    completed = weights_dir/"manifest.tsv"
    # A global manifest alone is not a reusable computation receipt. Validate
    # every per-job marker below, even when a previous global publication exists.
    # Withdraw publication before replay so a failed repair cannot leave a
    # complete manifest pointing at partially reconstructed derived weights.
    prepared = getattr(args, "prepared_inputs", None)
    if prepared and not args.replace:
        from csx_parallel import reuse_partition
        if reuse_partition(args, chromosomes):
            print(f"[reuse] chromosome {chromosomes[0]}: verified posterior and weights", flush=True)
            return completed
    if getattr(args, "prepare_only", False):
        from copy import copy
        from csx_parallel import load_prepared
        if not args.replace and (out/"prepared.json").is_file():
            previous = copy(args)
            previous.prepared_inputs = str(out)
            try:
                data = load_prepared(previous, chromosomes, load_maps=False)
                if data["chromosomes"] == chromosomes:
                    print("[reuse] shared inference preparation", flush=True)
                    return out/"prepared.json"
            except (ValueError, OSError, KeyError):
                pass
        (out/"prepared.json").unlink(missing_ok=True)
    completed.unlink(missing_ok=True)
    if prepared:
        from csx_parallel import load_prepared
        context = load_prepared(args, chromosomes)
        records, kind, reference_file = context["records"], context["kind"], context["reference_file"]
        maps, normalized = context["maps"], context["normalized"]
        normalization_audit, singlepop = context["normalization_audit"], context["singlepop"]
        args.ld_ref = context["ld_root"]
        target_root = Path(prepared)/"target"
    else:
        records = read_manifest(args.manifest, args.trait)
        source_ld_ref = args.ld_ref
        args.ld_ref = canonical_ld_reference(args.ld_ref, out/"ld_reference", chromosomes)
        kind, reference_file, reference = reference_table(args.ld_ref, chromosomes)
        maps = build_target_reference(reference, args.genotype_dir, out, chromosomes, args.rsid_map, args.allow_ambiguous)
        normalized, normalization_audit = {}, []
        for pop in POPULATIONS:
            normalized[pop], audit = normalize_sumstats(records[pop], reference, maps, out, args.allow_ambiguous)
            normalization_audit.append(audit)
        singlepop = prepare_singlepop_reference(reference, args.ld_ref, kind, out, chromosomes)
        del reference
        target_root = out/"target"
        if getattr(args, "prepare_only", False):
            from csx_parallel import write_prepared
            return write_prepared(args, chromosomes, source_ld_ref, records, kind,
                                  reference_file, normalized, normalization_audit, singlepop)
    models = model_manifest(out, phis, args.cs_fixed)
    environment = os.environ.copy()
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[var] = str(args.threads)
    environment["GRID_CSX_BACKEND"] = args.backend
    n = {p:records[p]["n_gwas"] for p in POPULATIONS}
    jobs = [("csx",c,phi,None) for phi in [*phis,"auto"] for c in chromosomes]
    jobs += [("cs",c,phi,p) for c in chromosomes for p in POPULATIONS
             for phi in (["auto",*phis] if args.cs_fixed else ["auto"])]
    def task(job):
        family, chrom, phi, pop = job
        name = f"{family}.{phi}.chr{chrom}"+(f".{pop}" if pop else "")
        directory = out/"raw"/name; directory.mkdir(parents=True, exist_ok=True)
        prefix = directory/"posterior" if family=="cs" else directory
        command = inference_command(family,chrom,phi,prefix,target_root/f"chr{chrom}",
            {p:normalized[p][chrom] for p in POPULATIONS},n,args.ld_ref if family=="csx" else singlepop[pop],
            args.trait,pop,args.n_iter,args.n_burnin,args.thin,args.seed+chrom)
        if family == "csx":
            selected = models.loc[models.family.isin(["csx","csx_meta"]) & models.phi.eq(phi)]
            raw_outputs = [directory/f"{args.trait}_META_pst_eff_chr{chrom}.txt" if row.population=="META"
                           else directory/row.population/f"{row.population}.chr{chrom}.pst_eff.txt"
                           for row in selected.itertuples(index=False)]
        else:
            selected = models.loc[models.family.eq("cs") & models.population.eq(pop) & models.phi.eq(phi)]
            token = "auto" if phi == "auto" else f"{float(phi):.0e}"
            raw_outputs = [directory/f"posterior_pst_eff_a1_b0.5_phi{token}_chr{chrom}.txt"]
        _, target_pvar = pfile(Path(args.genotype_dir)/f"chr{chrom}")
        populations = POPULATIONS if family == "csx" else (pop,)
        originals = [reference_file,target_pvar]
        if args.rsid_map:
            originals.append(args.rsid_map)
        for population in populations:
            originals += records[population]["files"]
            originals.append(Path(args.ld_ref)/f"ldblk_{kind}_{population.lower()}"/f"ldblk_{kind}_chr{chrom}.hdf5")
        # Code-stat changes invalidate a previous posterior instead of reusing
        # output from a different local sampler/normalization implementation.
        originals += [Path(__file__),HERE/"genomics.py",*sorted((HERE/family).glob("*.py"))]
        signature = dict(schema="grid-inference-job-v1", family=family, chromosome=chrom, phi=phi,
            population=pop, n_gwas={p:n[p] for p in populations}, n_iter=args.n_iter,
            n_burnin=args.n_burnin, thin=args.thin, seed=args.seed+chrom,
            backend=args.backend if family=="csx" else "cpu", threads=args.threads,
            allow_ambiguous=args.allow_ambiguous,strand_aligned_provenance=args.strand_aligned_provenance,
            original_inputs=[file_stat(p) for p in originals])
        marker = directory/"complete.json"
        weight_outputs = [weights_dir/f"{row.score_name}.chr{chrom}.tsv" for row in selected.itertuples(index=False)]
        complete_outputs = [*raw_outputs, *weight_outputs]
        reused = not args.replace and can_resume(marker,signature,complete_outputs)
        if marker.exists() and not args.replace and not reused:
            raise ValueError(f"{name}: checkpoint or outputs differ; rerun with --replace TRUE")
        if reused:
            return json.loads(marker.read_text())["counts"], True
        if not reused:
            marker.unlink(missing_ok=True)
            with (directory/"command.log").open("w") as handle:
                handle.write(" ".join(command)+"\n"); handle.flush()
                subprocess.run(command, env=environment, stdout=handle, stderr=subprocess.STDOUT, check=True)
        counts = []
        # Publish weights together with the posterior and its completion receipt.
        for row,raw in zip(selected.itertuples(index=False),raw_outputs):
            target = weights_dir/f"{row.score_name}.chr{chrom}.tsv"
            count = export_posterior(raw,maps[chrom],target,args.allow_ambiguous)
            counts.append(dict(score_name=row.score_name,chromosome=chrom,n_variants=count))
        if not reused:
            write_job_marker(marker,signature,complete_outputs,counts)
        return counts,reused
    counts = []; reused_jobs = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(task,job):job for job in jobs}
        for future in as_completed(futures):
            try:
                job_counts, reused = future.result()
                counts.extend(job_counts); reused_jobs += int(reused)
            except Exception:
                for pending in futures:
                    pending.cancel()
                raise
    pd.DataFrame(counts).sort_values(["score_name","chromosome"]).to_csv(weights_dir/"snp_counts.tsv", sep="\t", index=False)
    (out/"provenance.json").write_text(json.dumps(dict(trait=args.trait, chromosomes=chromosomes,
        reference_file=str(reference_file), reference_build="GRCh37", target_build="GRCh38",
        genotype_dir=str(Path(args.genotype_dir).resolve()), target_role="allele and availability metadata only",
        raw_gwas=records, normalization=normalization_audit, normalized_effect="target-ALT z with SE=1",
        ambiguous_variants="allowed with explicit provenance" if args.allow_ambiguous else "excluded",
        strand_aligned_provenance=args.strand_aligned_provenance,
        phi_grid=phis, csx_auto_meta="posterior inverse-variance meta, not equal-weight mean",
        singlepop_engine="getian107/PRScs fbc9ec017652f0b40a6ab47f63747d97fa7d033b",
        n_iter=args.n_iter,n_burnin=args.n_burnin,thin=args.thin,seed=args.seed,
        jobs=args.jobs,threads_per_job=args.threads,csx_backend=args.backend,
        mcmc_jobs_total=len(jobs),mcmc_jobs_reused=reused_jobs,
        resume_policy="matching method settings and original input/code/output file size+mtime_ns; --replace forces rerun",
        common_variant_policy="same HM3-target eligibility; model GWAS availability may differ"), indent=2)+"\n")
    models.to_csv(completed, sep="\t", index=False)
    if prepared:
        from csx_parallel import write_partition_receipt
        write_partition_receipt(args, chromosomes)
    return completed


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("manifest", "trait", "ld-ref", "genotype-dir", "output"):
        p.add_argument("--"+name, required=True)
    p.add_argument("--chrs", default="1-22")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--prepare-only", action="store_true", help="normalize all chromosomes once, without MCMC")
    mode.add_argument("--prepared-inputs", help="read-only shared inputs from a completed preparation")
    mode.add_argument("--merge-only", action="store_true", help="validate every chromosome before publishing weights")
    p.add_argument("--phi-grid", default="1e-2", help="comma-separated positive values; auto always runs as well")
    p.add_argument("--cs-fixed", action="store_true", help="also match every fixed phi in the four singlepop controls")
    p.add_argument("--rsid-map", help="optional target_id,rsid GRCh38 bridge for coordinate-only PVAR IDs")
    p.add_argument("--allow-ambiguous", action="store_true")
    p.add_argument("--strand-aligned-provenance", default="")
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--backend", choices=["cpu","cuda"], default="cpu")
    p.add_argument("--n-iter", type=int, default=4000)
    p.add_argument("--n-burnin", type=int, default=2000)
    p.add_argument("--thin", type=int, default=5)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--replace", action="store_true")
    args = p.parse_args(argv)
    print(run_inference(args))


if __name__ == "__main__":
    main()
