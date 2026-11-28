#!/usr/bin/env python3
"""GRCh38 PGEN metadata harmonization and bounded-memory dosage scoring.

No genotype matrix is exported. SNP identity joins reference/GWAS rsIDs to actual
GRCh38 PVAR records; reference BP37 is never used as a target coordinate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
import csv
import glob
import gzip
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd

COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASES = set("ACGT")


@contextmanager
def text_input(path):
    path = Path(path)
    if path.suffix == ".zst":
        process = subprocess.Popen(["zstdcat", str(path)], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True)
        try:
            yield process.stdout
            error = process.stderr.read()
            if process.wait():
                raise ValueError(f"Cannot read {path}: {error.strip()}")
        finally:
            if process.poll() is None:
                process.terminate()
            process.stdout.close()
            process.stderr.close()
    else:
        with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as stream:
            yield stream


def chromosome(value):
    value = str(value).removeprefix("chr").removeprefix("CHR")
    number = int(value)
    if not 1 <= number <= 22:
        raise ValueError("This scoring workflow uses autosomes 1–22")
    return number


def parse_chromosomes(value):
    result = []
    for part in str(value).replace(" ", ",").split(","):
        if not part:
            continue
        if "-" in part:
            lo, hi = map(int, part.split("-"))
            result.extend(range(lo, hi+1))
        else:
            result.append(chromosome(part))
    if not result or len(set(result)) != len(result) or not set(result) <= set(range(1, 23)):
        raise ValueError("Use unique autosomes, for example 1-22 or 21,22")
    return sorted(result)


def pfile(prefix):
    prefix = Path(prefix)
    for suffix in (".pgen", ".psam"):
        if not Path(str(prefix)+suffix).is_file():
            raise FileNotFoundError(str(prefix)+suffix)
    plain, packed = Path(str(prefix)+".pvar"), Path(str(prefix)+".pvar.zst")
    if plain.is_file():
        return ["--pfile", str(prefix)], plain
    if packed.is_file():
        return ["--pfile", str(prefix), "vzs"], packed
    raise FileNotFoundError(f"Missing {prefix}.pvar or .pvar.zst")


def rsids(value, info=""):
    found = set(re.findall(r"(?<![A-Za-z0-9])rs[0-9]+(?![A-Za-z0-9])", str(value), flags=re.I))
    result = {s.lower() for s in found}
    for token in str(info).split(";"):
        if token.startswith(("RS=", "RSID=", "rsid=")):
            for item in token.split("=", 1)[1].split(","):
                if item.isdigit():
                    result.add("rs"+item)
                elif re.fullmatch("rs[0-9]+", item, flags=re.I):
                    result.add(item.lower())
    return result


def pvar_rows(path):
    with text_input(path) as stream:
        columns = None
        for line in stream:
            if line.startswith("##") or not line.strip():
                continue
            fields = line.rstrip("\n").split()
            if line.startswith("#"):
                columns = [x.lstrip("#").upper() for x in fields]
                if not {"CHROM", "POS", "ID", "REF", "ALT"} <= set(columns):
                    raise ValueError(f"Invalid PVAR header: {path}")
                continue
            if columns is None or len(fields) < len(columns):
                raise ValueError(f"PVAR row without complete header/fields: {path}")
            row = dict(zip(columns, fields))
            yield dict(CHR=chromosome(row["CHROM"]), BP38=int(row["POS"]),
                       TARGET_ID=row["ID"], REF38=row["REF"].upper(), ALT38=row["ALT"].upper(),
                       RSIDS=rsids(row["ID"], row.get("INFO", "")))


def bridge_for_chromosome(path, chrom, candidates):
    """Optional GRCh38 bridge: target_id,rsid plus chrom38/bp38 when available."""
    result = {}
    if not path:
        return result
    with text_input(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for raw in reader:
            row = {k.lower(): v for k, v in raw.items()}
            if not {"target_id", "rsid"} <= set(row):
                raise ValueError("rsID map requires target_id and rsid columns")
            ids = rsids(row["rsid"])
            ids &= candidates
            if not ids:
                continue
            c = row.get("chrom38", row.get("chr", ""))
            if c and chromosome(c) != chrom:
                continue
            build = row.get("build", row.get("target_build", "GRCh38")).upper()
            if build not in {"GRCH38", "HG38", "38"}:
                raise ValueError("rsID bridge must declare GRCh38 target coordinates")
            result.setdefault(row["target_id"], []).append((ids, row.get("bp38", "")))
    return result


def allele_alignment(a1, a2, ref, alt, allow_ambiguous=False):
    """Return target-strand effect allele and orientation, or None.

No BETA sign change is needed for a strand complement. REF-dosage conversion
to ALT dosage is a separate arithmetic step retaining its diploid offset.
"""
    a1, a2, ref, alt = map(str.upper, (str(a1), str(a2), str(ref), str(alt)))
    if any(x not in BASES for x in (a1, a2, ref, alt)) or a1 == a2 or ref == alt:
        return None
    if {a1, a2} in ({"A", "T"}, {"C", "G"}) and not allow_ambiguous:
        return None
    if {a1, a2} == {ref, alt}:
        return a1, "direct"
    b1, b2 = a1.translate(COMPLEMENT), a2.translate(COMPLEMENT)
    if {b1, b2} == {ref, alt}:
        return b1, "strand_complement"
    return None


def map_target_variants(pvar, reference, chrom, *, rsid_map=None, allow_ambiguous=False):
    """Map one chromosome; only variants with unique compatible identity survive.

reference is an rsID-indexed DataFrame with A1,A2 and optional CHR/BP37. The
candidate dictionary bounds memory by the HM3/reference list, not all PVAR rows.
"""
    if reference.index.duplicated().any():
        raise ValueError("Reference rsIDs must be unique")
    wanted = set(reference.index.astype(str))
    bridge = bridge_for_chromosome(rsid_map, chrom, wanted)
    found, bad, seen_ids = {}, set(), Counter()
    counts = Counter()
    for row in pvar_rows(pvar):
        if row["CHR"] != chrom:
            raise ValueError(f"PVAR contains another chromosome: {pvar}")
        ids = row["RSIDS"] & wanted
        for aliases, bp in bridge.get(row["TARGET_ID"], []):
            if bp and int(bp) != row["BP38"]:
                raise ValueError("rsID bridge coordinate differs from actual GRCh38 PVAR")
            ids |= aliases
        if not ids:
            continue
        counts["candidate_records"] += 1
        seen_ids[row["TARGET_ID"]] += 1
        if row["TARGET_ID"] in {".", "0"} or row["BP38"] < 1:
            counts["invalid_identity"] += 1
            continue
        for snp in ids:
            r = reference.loc[snp]
            alignment = allele_alignment(r.A1, r.A2, row["REF38"], row["ALT38"], allow_ambiguous)
            if alignment is None:
                counts["allele_incompatible_or_ambiguous"] += 1
                continue
            if snp in found:
                bad.add(snp)
            else:
                found[snp] = {k:v for k,v in row.items() if k != "RSIDS"}
                found[snp].update(SNP=snp, reference_alignment=alignment[1])
    repeated_ids = {k for k,v in seen_ids.items() if v > 1}
    # One target ID cannot represent two independent scoring SNPs.
    target_counts = Counter(v["TARGET_ID"] for v in found.values())
    repeated_ids |= {k for k,v in target_counts.items() if v > 1}
    retained = [v for k,v in found.items() if k not in bad and v["TARGET_ID"] not in repeated_ids]
    cols = ["SNP", "CHR", "BP38", "TARGET_ID", "REF38", "ALT38", "reference_alignment"]
    table = pd.DataFrame(retained, columns=cols).sort_values(["BP38", "SNP"])
    counts.update(reference_candidates=len(reference), retained=len(table),
                  ambiguous_rsid=len(bad), duplicate_target_ids=len(repeated_ids))
    counts["not_retained"] = len(reference)-len(table)
    return table.reset_index(drop=True), dict(counts)


def align_weights(weights, target, allow_ambiguous=False):
    """Export ALT-oriented per-allele BETA and its exact autosomal offset."""
    required = {"SNP", "A1", "A2", "BETA"}
    if not required <= set(weights) or weights.SNP.duplicated().any() or target.SNP.duplicated().any():
        raise ValueError("Unique SNP,A1,A2,BETA weights and unique target mapping are required")
    joined = weights.merge(target, on="SNP", how="inner", validate="one_to_one", suffixes=("_reference", ""))
    records = []
    for row in joined.itertuples(index=False):
        aligned = allele_alignment(row.A1, row.A2, row.REF38, row.ALT38, allow_ambiguous)
        if aligned is None or not np.isfinite(float(row.BETA)):
            continue
        effect, strand = aligned
        flip = effect == row.REF38
        beta = float(row.BETA)
        records.append(dict(SNP=row.SNP, A1=row.ALT38, A2=row.REF38,
            BETA=-beta if flip else beta, OFFSET=2*beta if flip else 0.,
            ORIGINAL_A1=row.A1, ORIGINAL_BETA=beta, STRAND_ALIGNMENT=strand,
            CHR=row.CHR, REF_BP37=getattr(row, "REF_BP37", np.nan),
            TARGET_ID=row.TARGET_ID, BP38=row.BP38, REF38=row.REF38, ALT38=row.ALT38))
    return pd.DataFrame(records, columns=["SNP","A1","A2","BETA","OFFSET","ORIGINAL_A1","ORIGINAL_BETA",
        "STRAND_ALIGNMENT","CHR","REF_BP37","TARGET_ID","BP38","REF38","ALT38"])


def read_roster(path):
    d = pd.read_csv(path, sep="\t", dtype=str)
    if not {"id", "fid", "iid", "group_id", "role"} <= set(d):
        raise ValueError("Roster requires id,fid,iid,group_id,role")
    if d[["id","fid","iid","group_id","role"]].isna().any().any() or d.id.duplicated().any():
        raise ValueError("Roster identities and roles must be complete and unique")
    if d.duplicated(["fid","iid"]).any() or not d.role.eq("train").any():
        raise ValueError("Unique FID/IID and at least one train sample are required")
    if d.groupby("group_id").role.nunique().max() != 1:
        raise ValueError("A family group crosses roster roles")
    return d.reset_index(drop=True)


def write_keep(roster, path, train_only=False):
    d = roster.loc[roster.role.eq("train")] if train_only else roster
    d[["fid","iid"]].rename(columns={"fid":"#FID","iid":"IID"}).to_csv(path, sep="\t", index=False)


def plink_score_command(prefix, weights, keep, frequencies, out, n_scores, *, plink="plink2", threads=4, memory_mb=4096, extract=None):
    args, _ = pfile(prefix)
    if n_scores < 1 or threads < 1 or memory_mb < 640:
        raise ValueError("Invalid scoring resource/column settings")
    command = [plink, *args, "--keep", str(keep), "--read-freq", str(frequencies),
        "--score", str(weights), "1", "2", "header-read", "list-variants", "cols=+scoresums",
        "--score-col-nums", "3" if n_scores == 1 else f"3-{n_scores+2}",
        "--threads", str(threads), "--memory", str(memory_mb), "--out", str(out)]
    if extract is not None:
        command += ["--extract", str(extract)]
    return command


def run_command(command, log):
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    with Path(log).open("w") as stream:
        stream.write(" ".join(map(str, command))+"\n"); stream.flush()
        subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT, check=True)


def ensure_frequencies(prefix, train_keep, extract, out, *, plink, threads, memory_mb):
    # Every requested run derives AF using training genotypes only. It never
    # estimates test-batch frequencies or converts PGEN dosages to hard calls.
    args, _ = pfile(prefix)
    command = [plink, *args, "--keep", str(train_keep), "--extract", str(extract),
        "--freq", "--nonfounders", "--threads", str(threads), "--memory", str(memory_mb), "--out", str(out)]
    run_command(command, str(out)+".command.log")
    path = Path(str(out)+".afreq")
    if not path.is_file() or not path.stat().st_size:
        raise ValueError("PLINK did not produce training allele frequencies")
    return path


def read_scores(path, roster, names):
    d = pd.read_csv(path, sep=r"\s+", dtype={"#FID":str,"FID":str,"#IID":str,"IID":str})
    d = d.rename(columns={"#FID":"FID","#IID":"IID"})
    if "FID" not in d:
        d["FID"] = "0"
    if "IID" not in d or d.duplicated(["FID","IID"]).any():
        raise ValueError("Malformed PLINK score identities")
    expected = pd.MultiIndex.from_frame(roster[["fid","iid"]])
    d = d.set_index(["FID","IID"])
    if set(d.index) != set(expected):
        raise ValueError("PLINK scores do not cover the exact roster FID/IID pairs")
    columns = [name+"_SUM" for name in names]
    if not set(columns) <= set(d):
        raise ValueError("PLINK must report score SUM columns; averages are not accepted")
    values = d.loc[expected, columns].to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError("PLINK scores are incomplete or nonfinite")
    return values


def named_paths(values):
    out = {}
    for value in values:
        if "=" not in value:
            raise ValueError("Use NAME=path-pattern entries")
        key, pattern = value.split("=", 1)
        if key in out or not re.fullmatch(r"[A-Za-z0-9_.-]+", key):
            raise ValueError("Scoring names must be unique simple identifiers")
        paths = [Path(x) for x in sorted(glob.glob(pattern))]
        if not paths:
            raise FileNotFoundError(pattern)
        out[key] = paths
    return out


def resolve_weight_files(value):
    if isinstance(value, (str, Path)):
        files = [Path(p) for p in sorted(glob.glob(str(value))) ]
        if not files:
            raise FileNotFoundError(value)
        return files
    return [Path(p) for p in value]


def index_weight_files(paths):
    """Index one-chromosome exports using only a header/first row per file."""
    catalog = {}
    for name, files in paths.items():
        catalog[name] = {}
        for path in resolve_weight_files(files):
            head = pd.read_csv(path, sep="\t", nrows=1, usecols=["CHR"])
            if head.empty:
                raise ValueError(f"Empty posterior export: {path}")
            chrom = chromosome(head.CHR.iloc[0])
            catalog[name].setdefault(chrom, []).append(path)
    if not catalog:
        raise ValueError("At least one posterior weight model is required")
    return catalog


def load_weight_sets(paths, common_variants=False):
    data = {}
    for name, files in paths.items():
        files = resolve_weight_files(files)
        required = {"SNP","A1","A2","BETA","CHR","TARGET_ID","BP38","REF38","ALT38"}
        d = (pd.concat([pd.read_csv(p, sep="\t", dtype={"SNP":str,"TARGET_ID":str}) for p in files], ignore_index=True)
             if files else pd.DataFrame(columns=sorted(required|{"OFFSET"})))
        if not required <= set(d) or d.SNP.duplicated().any() or d.TARGET_ID.duplicated().any():
            raise ValueError(f"{name}: requires unique mapped GRCh38 weight exports from infer.py")
        if not np.isfinite(d.BETA.to_numpy(float)).all() or not d.A1.eq(d.ALT38).all() or not d.A2.eq(d.REF38).all():
            raise ValueError(f"{name}: weights must use the exported target ALT orientation")
        if "OFFSET" not in d:
            d["OFFSET"] = 0.
        if not np.isfinite(d.OFFSET.to_numpy(float)).all():
            raise ValueError(f"{name}: invalid allele-conversion offset")
        data[name] = d
    if common_variants and data:
        keep = set.intersection(*(set(d.SNP) for d in data.values()))
        data = {k:d.loc[d.SNP.isin(keep)].copy() for k,d in data.items()}
    return data


def chromosome_weight_sets(catalog, chrom, common_variants=False):
    models = load_weight_sets({name:by_chr.get(chrom, []) for name,by_chr in catalog.items()}, common_variants)
    for name,d in models.items():
        if not d.empty and not d.CHR.eq(chrom).all():
            raise ValueError(f"{name}: use the one-chromosome weight exports from infer.py")
    return models


def score_batch(prefix, roster, weights, columns, keep, train_keep, work, *, plink, threads, memory_mb, frequencies=None):
    """One bounded matrix of SNP weights and one PLINK dosage pass."""
    work = Path(work); work.mkdir(parents=True, exist_ok=True)
    ids = weights.TARGET_ID.to_numpy(str)
    extract = work/"variants.txt"; pd.Series(ids).to_csv(extract, index=False, header=False)
    score_path = work/"weights.tsv"
    weights[["TARGET_ID","A1",*columns]].to_csv(score_path, sep="\t", index=False)
    if frequencies is None:
        frequencies = ensure_frequencies(prefix, train_keep, extract, work/"train", plink=plink, threads=threads, memory_mb=memory_mb)
    out = work/"score"
    run_command(plink_score_command(prefix, score_path, keep, frequencies, out, len(columns),
        plink=plink, threads=threads, memory_mb=memory_mb, extract=extract), work/"score.command.log")
    used_file = Path(str(out)+".sscore.vars")
    if not used_file.is_file() or set(used_file.read_text().split()) != set(ids):
        raise ValueError("PLINK silently omitted requested scoring variants; inspect the score log")
    result = read_scores(str(out)+".sscore", roster, columns)
    # The .sscore table can be large. Retain only small command/PLINK logs.
    Path(str(out)+".sscore").unlink(); used_file.unlink()
    score_path.unlink(); extract.unlink()
    return result


def score_models(roster_path, weights_paths, genotype_dir, output, *, common_variants=False,
                 plink="plink2", threads=4, memory_mb=4096, batch_size=32):
    roster = read_roster(roster_path)
    catalog = index_weight_files(weights_paths)
    names = list(catalog)
    output = Path(output); output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent/"tmp"/"scoring"; work.mkdir(parents=True, exist_ok=True)
    keep, train_keep = work/"roster.keep", work/"train.keep"
    write_keep(roster, keep); write_keep(roster, train_keep, True)
    totals = np.lib.format.open_memmap(work/"score_sums.npy", mode="w+", dtype="float64", shape=(len(roster),len(names)))
    totals[:] = 0
    audit = []
    chromosomes = sorted(set().union(*(set(by_chr) for by_chr in catalog.values())))
    for chrom in chromosomes:
        chrom = chromosome(chrom); prefix = Path(genotype_dir)/f"chr{chrom}"
        subsets = chromosome_weight_sets(catalog, chrom, common_variants)
        if not any(len(d) for d in subsets.values()):
            continue
        metadata = pd.concat([d[["TARGET_ID","A1","BP38","REF38","ALT38"]] for d in subsets.values()])
        if (metadata.groupby("TARGET_ID")[["A1","BP38","REF38","ALT38"]].nunique()>1).any().any():
            raise ValueError("Model metadata disagrees on target identity/strand/build")
        metadata = metadata.drop_duplicates("TARGET_ID").set_index("TARGET_ID")
        cwork = work/f"chr{chrom}"; cwork.mkdir(exist_ok=True)
        extract = cwork/"frequency_variants.txt"
        pd.Series(metadata.index).to_csv(extract, index=False, header=False)
        frequency = ensure_frequencies(prefix, train_keep, extract, cwork/"train", plink=plink, threads=threads, memory_mb=memory_mb)
        for start in range(0, len(names), batch_size):
            batch = names[start:start+batch_size]
            matrix = metadata[["A1"]].copy()
            offsets = []
            for name in batch:
                d = subsets[name]
                matrix[name] = d.set_index("TARGET_ID").BETA.reindex(matrix.index).fillna(0.)
                offsets.append(float(d.OFFSET.sum()))
                audit.append(dict(score_name=name, chromosome=chrom, n_scored_variants=len(d),
                                  common_variant_mode=common_variants, target_build="GRCh38", score_type="dosage_sum"))
            values = score_batch(prefix, roster, matrix.reset_index(), batch, keep, train_keep,
                cwork/f"batch{start//batch_size}", plink=plink, threads=threads, memory_mb=memory_mb, frequencies=frequency)
            totals[:,start:start+len(batch)] += values + np.asarray(offsets)[None,:]
        totals.flush()
    if not audit:
        raise ValueError("No scoring variants remain; check the common-variant intersection")
    # Stream the final table without materializing the complete sample/score matrix.
    with (gzip.open(output, "wt") if output.suffix == ".gz" else output.open("w")) as stream:
        for start in range(0,len(roster),50000):
            d = pd.DataFrame(np.asarray(totals[start:start+50000]), columns=names)
            d.insert(0,"id",roster.id.iloc[start:start+50000].to_numpy())
            d.to_csv(stream, sep="\t", index=False, header=start==0)
    pd.DataFrame(audit).to_csv(output.parent/"score_counts.tsv",sep="\t",index=False)
    del totals
    (work/"score_sums.npy").unlink()
    return output


def read_blocks(path, population):
    """Read real GRCh38 BED, zero-based half-open; reject overlap per map."""
    records = []
    with text_input(path) as stream:
        for number,line in enumerate(stream,1):
            if not line.strip() or line.startswith(("#","track","browser")):
                continue
            row = line.split()
            if row[0].lower() in {"chr","chrom","chromosome"}:
                continue
            if len(row)<3:
                raise ValueError(f"Invalid BED row {path}:{number}")
            c,start,end = chromosome(row[0]),int(row[1]),int(row[2])
            if not 0 <= start < end:
                raise ValueError("Blocks must be 0-based half-open GRCh38 BED intervals")
            name = row[3] if len(row)>3 else f"chr{c}:{start}-{end}"
            records.append(dict(map=population,population=population,chrom38=c,start38=start,end38=end,name=population+":"+name))
    d = pd.DataFrame(records).sort_values(["chrom38","start38"]).reset_index(drop=True)
    if d.empty or d.name.duplicated().any():
        raise ValueError("Empty blocks or duplicate feature names")
    for _,g in d.groupby("chrom38"):
        if np.any(g.start38.to_numpy()[1:]<g.end38.to_numpy()[:-1]):
            raise ValueError("Overlapping blocks within one ancestry map would double-count SNPs")
    return d


def assign_blocks(bp38, blocks):
    """Target PVAR positions are one-based; BED intervals are zero-based."""
    positions = np.asarray(bp38,int)-1
    suffix = "37" if "start37" in blocks else "38"
    starts,ends = blocks["start"+suffix].to_numpy(int),blocks["end"+suffix].to_numpy(int)
    index = np.searchsorted(starts,positions,side="right")-1
    valid = index>=0
    valid[valid] &= positions[valid]<ends[index[valid]]
    return np.where(valid,index,-1)


def regional_scores(roster_path, weights_paths, block_paths, genotype_dir, output_dir, *,
                    plink="plink2",threads=4,memory_mb=4096,batch_size=32,block_build="GRCh38"):
    """One shared consensus catalog; META gives one regional PRS per LD block.

    Additional ancestry channels may be requested, but every channel uses the
    identical ordered catalog. Gaps are preserved and unassigned SNPs reported.
    An independent covered-total score in each PLINK pass verifies additivity
    without an extra full-genotype pass.
    """
    from consensus3 import read_consensus,unassigned_reasons
    if batch_size<1:
        raise ValueError("batch_size must be positive")
    roster=read_roster(roster_path)
    catalog=index_weight_files(weights_paths)
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+",name) for name in catalog):
        raise ValueError("Score channel names must be simple identifiers")
    chromosomes=sorted(set().union(*(set(value) for value in catalog.values())))
    if any(set(value)!=set(chromosomes) for value in catalog.values()):
        raise ValueError("Every requested score channel needs weights for the same chromosomes")
    blocks,consensus=read_consensus(block_paths,expected_chromosomes=chromosomes,build=block_build)
    suffix = "37" if "chrom37" in blocks else "38"
    position = "REF_BP37" if suffix == "37" else "BP38"
    lost_columns = ["SNP","TARGET_ID","BP38"] + (["REF_BP37"] if suffix=="37" else [])
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    (out/"metadata.json").unlink(missing_ok=True)
    work=out/"tmp";work.mkdir(exist_ok=True)
    keep,train_keep=work/"roster.keep",work/"train.keep"
    write_keep(roster,keep);write_keep(roster,train_keep,True)
    features=pd.concat([blocks.assign(map=name,population=name,name=name+":"+blocks.block_id)
                        for name in catalog],ignore_index=True)
    features.insert(0,"feature_index",np.arange(len(features)))
    features["n_snps"]=0
    matrix=np.lib.format.open_memmap(out/"scores.npy",mode="w+",dtype="float32",shape=(len(roster),len(features)))
    matrix[:]=0
    audit=[];checks=[]
    totals=np.zeros((len(roster),len(catalog)),dtype=np.float64)
    unassigned_path=out/"unassigned_variants.tsv.gz"
    with gzip.open(unassigned_path,"wt") as missing_stream:
        missing_stream.write("\t".join(["map","chromosome",*lost_columns,"reason"])+"\n")
        for model_index,name in enumerate(catalog):
            assigned_total=0
            f=features.loc[features["map"].eq(name)]
            for chrom,block in f.groupby("chrom"+suffix,sort=True):
                w=chromosome_weight_sets({name:catalog[name]},chrom)[name].copy()
                if w.empty:
                    raise ValueError(f"Missing posterior chromosome {chrom} for {name}")
                if position not in w:
                    raise ValueError(f"GRCh{suffix} regional assignment needs {position} in the posterior export")
                bp=pd.to_numeric(w[position],errors="raise").to_numpy(float)
                if not np.isfinite(bp).all() or (bp<1).any() or not np.equal(bp,np.floor(bp)).all():
                    raise ValueError(f"GRCh{suffix} BED assignment needs valid positive integer {position} positions")
                local=assign_blocks(bp,block)
                absolute=block.feature_index.to_numpy(int)
                w["feature_index"]=np.where(local>=0,absolute[np.maximum(local,0)],-1)
                lost=w.loc[local<0,lost_columns].copy()
                if not lost.empty:
                    lost.insert(0,"chromosome",chrom);lost.insert(0,"map",name)
                    lost["reason"]=unassigned_reasons(lost[position],block)
                    lost.to_csv(missing_stream,sep="\t",index=False,header=False)
                assigned=w.loc[w.feature_index.ge(0)]
                assigned_total+=len(assigned)
                counts=assigned.groupby("feature_index").size()
                for idx,n in counts.items():features.loc[idx,"n_snps"]=int(n)
                row=dict(map=name,chromosome=int(chrom),n_weight_snps=len(w),n_assigned=len(assigned),
                    n_unassigned=len(lost),assigned_fraction=len(assigned)/len(w),
                    n_internal_gap=int(lost.reason.eq("internal_cut_gap").sum()) if len(lost) else 0,
                    n_before_catalog=int(lost.reason.eq("before_catalog").sum()) if len(lost) else 0,
                    n_after_catalog=int(lost.reason.eq("after_catalog").sum()) if len(lost) else 0,
                    n_catalog_blocks=len(block),n_empty_blocks=int(len(block)-len(counts)),
                    target_build="GRCh38",block_build="GRCh"+suffix,assignment_position=position,block_coordinate_system="BED_0_based_half_open")
                audit.append(row)
                if assigned.empty:
                    checks.append(dict(map=name,chromosome=int(chrom),n_scored_variants=0,
                        max_abs_error=0.,max_allowed_abs_error=0.,status="empty_chromosome"))
                    continue
                cwork=work/name/f"chr{chrom}";cwork.mkdir(parents=True,exist_ok=True)
                extract=cwork/"frequency_variants.txt"
                assigned.TARGET_ID.to_csv(extract,index=False,header=False)
                prefix=Path(genotype_dir)/f"chr{chrom}"
                frequency=ensure_frequencies(prefix,train_keep,extract,cwork/"train",plink=plink,threads=threads,memory_mb=memory_mb)
                expected_global=np.zeros(len(roster),float)
                text_scale=np.zeros(len(roster),float)
                for start in range(0,len(block),batch_size):
                    chosen=block.iloc[start:start+batch_size].feature_index.to_numpy(int)
                    batch=assigned.loc[assigned.feature_index.isin(chosen)]
                    if batch.empty:continue
                    score_matrix=batch[["TARGET_ID","A1"]].copy()
                    names=[f"block_{i}" for i in chosen]
                    offsets=[]
                    for idx,column in zip(chosen,names):
                        selected=batch.feature_index.eq(idx)
                        score_matrix[column]=np.where(selected,batch.BETA,0.)
                        offsets.append(float(batch.loc[selected,"OFFSET"].sum()))
                    # Independent same-universe global sum: every variant once.
                    score_matrix["covered_total"]=batch.BETA.to_numpy(float)
                    values=score_batch(prefix,roster,score_matrix,[*names,"covered_total"],keep,train_keep,
                        cwork/f"batch{start//batch_size}",plink=plink,threads=threads,memory_mb=memory_mb,frequencies=frequency)
                    offsets=np.asarray(offsets,float)
                    matrix[:,chosen]=values[:,:-1]+offsets[None,:]
                    expected_global+=values[:,-1]+offsets.sum()
                    text_scale+=np.abs(values[:,:-1]).sum(axis=1)+np.abs(offsets).sum()
                reconstructed=np.zeros(len(roster),float)
                for row_start in range(0,len(roster),8192):
                    reconstructed[row_start:row_start+8192]=np.asarray(
                        matrix[row_start:row_start+8192,absolute],dtype=np.float64).sum(axis=1)
                error=np.abs(reconstructed-expected_global)
                # PLINK text serialization and float32 storage have finite precision.
                tolerance=1e-6+1e-5*text_scale
                if np.any(error>tolerance):
                    raise ValueError(f"Regional sum disagrees with same-universe global score for {name} chr{chrom}")
                checks.append(dict(map=name,chromosome=int(chrom),n_scored_variants=len(assigned),
                    max_abs_error=float(error.max()),max_allowed_abs_error=float(tolerance.max()),status="passed"))
                totals[:,model_index]+=expected_global
            if not assigned_total:
                raise ValueError(f"No SNPs from {name} fall in the supplied consensus catalog")
    features["empty_block"]=features.n_snps.eq(0)
    features.to_csv(out/"features.tsv",sep="\t",index=False)
    blocks.to_csv(out/"catalog.tsv",sep="\t",index=False)
    roster[["id"]].to_csv(out/"ids.tsv",sep="\t",index=False)
    pd.DataFrame(audit).to_csv(out/"region_counts.tsv",sep="\t",index=False)
    pd.DataFrame(checks).to_csv(out/"sum_check.tsv",sep="\t",index=False)
    with gzip.open(out/"covered_global.tsv.gz","wt") as stream:
        for start in range(0,len(roster),50000):
            frame=pd.DataFrame(totals[start:start+50000],columns=list(catalog))
            frame.insert(0,"id",roster.id.iloc[start:start+50000].to_numpy())
            frame.to_csv(stream,sep="\t",index=False,header=start==0)
    matrix.flush()
    (out/"metadata.json").write_text(json.dumps(dict(target_build="GRCh38",shape=list(matrix.shape),
        dtype="float32",genotype_source=str(genotype_dir),frequency_fit_role="train",score_type="posterior_dosage_sum",
        geometry_scaling="not performed here; downstream TRAIN-only transformation",
        bed_coordinates="0-based-half-open",block_build="GRCh"+suffix,assignment_position=position,consensus=consensus,channels=list(catalog),
        n_common_blocks=len(blocks),n_empty_features=int(features.empty_block.sum()),
        n_unassigned_variants=int(sum(row["n_unassigned"] for row in audit)),
        additivity="sum of blocks equals independently scored global PRS on the same covered SNP universe",
        unassigned_policy="report_and_exclude; do not fill consensus cut gaps"),indent=2)+"\n")
    del matrix
    return out


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="module",required=True)
    for command in ("score","regions"):
        p=sub.add_parser(command)
        p.add_argument("--roster",required=True)
        p.add_argument("--weights",nargs="+",required=True,help="NAME=exported-weight-glob")
        p.add_argument("--genotype-dir",required=True,help="GRCh38 imputed chr*.pgen directory")
        p.add_argument("--plink",default="plink2")
        p.add_argument("--threads",type=int,default=4)
        p.add_argument("--memory-mb",type=int,default=4096)
        p.add_argument("--batch-size",type=int,default=32)
        if command=="score":
            p.add_argument("--output",required=True)
            p.add_argument("--common-variants",action="store_true")
        else:
            p.add_argument("--blocks",nargs="+",required=True,help="Shared GRCh38 consensus directory, quoted chr*.bed glob, or BED files")
            p.add_argument("--output-dir",required=True)
    a=parser.parse_args(argv)
    if a.batch_size<1:
        parser.error("--batch-size must be positive")
    common=dict(plink=a.plink,threads=a.threads,memory_mb=a.memory_mb,batch_size=a.batch_size)
    if a.module=="score":
        score_models(a.roster,named_paths(a.weights),a.genotype_dir,a.output,common_variants=a.common_variants,**common)
    else:
        regional_scores(a.roster,named_paths(a.weights),a.blocks,a.genotype_dir,a.output_dir,**common)


if __name__=="__main__":
    main()
