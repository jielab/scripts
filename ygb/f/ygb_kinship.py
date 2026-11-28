#!/usr/bin/env python3
"""Generate and validate the shared YGB family map; no LE8 dependency."""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import csv
import math
import statistics
from collections import Counter

KINSHIP_PROTOCOL = "2026-11-20.ygb-king-components-v2"
GROUPING_SOURCE = "YGB KING kinship connected components"
KINSHIP_ALGORITHM = "full released graph connected components, then project onto all.rds; singleton unrelated IDs"
KINSHIP_SEMANTICS = "all released KING edges; no kinship-threshold filtering"

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def eid_hash(ids):
    return hashlib.sha256("".join(x + "\n" for x in sorted(set(map(str, ids)))).encode()).hexdigest()

def _validate_ids(values, context, unique=False):
    import re
    ids = [str(v) for v in values]
    if not ids or any(not re.fullmatch(r"[A-Z]{5}", v) for v in ids):
        raise ValueError(f"{context}: participant IDs must be five uppercase letters")
    if unique and len(set(ids)) != len(ids): raise ValueError(f"Duplicate participant IDs in {context}")
    return ids

def read_pairs(path):
    with open(path) as f:
        names = next(f).split()
        required = {"ID1", "ID2", "HetHet", "IBS0", "Kinship"}
        if not required <= set(names): raise ValueError("YGB KING file requires ID1 ID2 HetHet IBS0 Kinship")
        rows = []
        for line in f:
            values = line.split()
            if not values: continue
            if len(values) != len(names): raise ValueError("Malformed YGB KING row")
            rows.append(dict(zip(names, values)))
        return rows

def read_groups(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != ["eid", "group"]: raise ValueError("Family map must contain eid,group only")
        return list(reader)

def build_kinship_groups(pairs, analysis_ids):
    """Use every released edge, including paths through people outside all.rds."""
    if not pairs: raise ValueError("YGB relatedness source is empty")
    ids = _validate_ids(analysis_ids, "all.rds", unique=True)
    source_ids = _validate_ids([v for row in pairs for v in (row["ID1"], row["ID2"])], "KING source")
    parent, size, seen, kinship = {}, {}, set(), []
    def find(x):
        if x not in parent: parent[x], size[x] = x, 1
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for row in pairs:
        a, b = row["ID1"], row["ID2"]
        edge = tuple(sorted((a, b)))
        if a == b or edge in seen: raise ValueError("Self/duplicate canonical YGB relatedness pair")
        seen.add(edge)
        numbers = [float(row[k]) for k in ("HetHet", "IBS0", "Kinship")]
        if not all(math.isfinite(v) for v in numbers): raise ValueError("Nonfinite KING coefficient")
        if numbers[2] < .0442: raise ValueError("Unexpected KING coefficient below 0.0442; review source, do not discard edges")
        kinship.append(numbers[2])
        a, b = find(a), find(b)
        if a != b:
            if (size[a], a) < (size[b], b): a, b = b, a
            parent[b] = a; size[a] += size[b]
    minima = {}
    for x in parent:
        root = find(x); minima[root] = min(x, minima.get(root, x))
    groups = [{"eid": x, "group": "KIN_" + minima[find(x)] if x in parent else "IND_" + x} for x in sorted(ids)]
    counts = Counter(row["group"] for row in groups)
    related = sum(x in parent for x in ids)
    stats = dict(pair_rows=len(pairs), unique_related_ids_in_source=len(set(source_ids)),
                 kinship_min=min(kinship), kinship_median=statistics.median(kinship), kinship_max=max(kinship),
                 kinship_below_0_0442_N=0, all_rds_N=len(ids), all_rds_eid_hash=eid_hash(ids),
                 all_rds_related_N=related, all_rds_singleton_N=len(ids)-related,
                 component_N=len(counts), largest_analysis_component_size=max(counts.values()),
                 connected_components_with_at_least_2_analysis_participants=sum(n >= 2 for n in counts.values()),
                 component_size_1_N=sum(n == 1 for n in counts.values()),
                 component_size_2_N=sum(n == 2 for n in counts.values()),
                 component_size_3_N=sum(n == 3 for n in counts.values()),
                 component_size_4plus_N=sum(n >= 4 for n in counts.values()))
    return groups, stats

def _rds_export(rds,out,r_bin,group_col=""):
    code = r'''
args<-commandArgs(trailingOnly=TRUE)
x<-readRDS(args[1])
if(!is.data.frame(x)||!"eid"%in%names(x))stop("all.rds must be a data.frame with eid")
g<-args[3]
z<-data.frame(eid=as.character(x$eid),stringsAsFactors=FALSE)
if(nzchar(g)){
 if(!g%in%names(x))stop("Group column missing: ",g)
 z$group<-as.character(x[[g]])
}
write.table(z,args[2],sep="\t",quote=FALSE,na="NA",row.names=FALSE)
'''
    subprocess.run([r_bin,"--vanilla","-e",code,str(rds),str(out),group_col],check=True)


def cohort_ids(allrds, r_bin):
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "ids.tsv"
        _rds_export(allrds, out, r_bin)
        with open(out, newline="") as f:
            return _validate_ids([row["eid"] for row in csv.DictReader(f, delimiter="\t")], "all.rds", unique=True)

def _read_validate(output, ygb_phe=None, r_bin="Rscript"):
    output = Path(output).resolve()
    meta = Path(str(output) + ".protocol.json")
    if not output.is_file() or not meta.is_file():
        raise FileNotFoundError("Shared YGB family map/metadata missing; run ygb/phe.sh phe --step kinship")
    record = json.loads(meta.read_text())
    if record.get("kinship_protocol") != KINSHIP_PROTOCOL:
        raise ValueError("Shared YGB family map protocol mismatch; regenerate with ygb/phe.sh --step kinship")
    rel = (output.parent / record["source_file"]).resolve()
    allrds = Path(ygb_phe).resolve() / "Rdata/all.rds" if ygb_phe else (output.parent / record["all_rds_file"]).resolve()
    rel_sha, all_sha = sha(rel), sha(allrds)
    if sha(output) != record["group_map_sha256"] or rel_sha != record["source_sha256"]:
        raise ValueError("Shared YGB kinship map/source checksum mismatch")
    saved = read_groups(output)
    saved_ids = _validate_ids([row["eid"] for row in saved], "family map", unique=True)
    # Phenotype values may change, but membership must remain exactly the same.
    if all_sha != record["all_rds_sha256"]:
        ids = cohort_ids(allrds, r_bin)
        if eid_hash(ids) != record["all_rds_eid_hash"]:
            raise ValueError("YGB cohort IDs changed; explicitly regenerate the shared family map")
    pairs = read_pairs(rel)
    rebuilt, stats = build_kinship_groups(pairs, saved_ids)
    if saved != rebuilt or any(record.get(k) != v for k, v in stats.items()):
        raise ValueError("Family map/audit does not match the full YGB relatedness graph")
    if sha(rel) != rel_sha or sha(allrds) != all_sha:
        raise RuntimeError("YGB input changed during validation")
    # Callers receive host-resolved paths; the on-disk sidecar stays portable.
    return {**record, "source_file": str(rel), "all_rds_file": str(allrds)}

def validate_kinship_groups(output, ygb_phe=None, r_bin="Rscript"):
    """Read-only consumer validation: never create or update shared files."""
    return _read_validate(output, ygb_phe, r_bin)

def generate(ygb_phe, r_bin="Rscript", replace=False):
    ygb_phe = Path(ygb_phe).resolve()
    output = ygb_phe / "rap/raw/ygb_kinship_groups.tsv"
    rel, allrds = output.parent / "ygb_rel.dat", ygb_phe / "Rdata/all.rds"
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output) + ".lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        meta = Path(str(output) + ".protocol.json")
        if (output.exists() or meta.exists()) and not replace:
            _read_validate(output, ygb_phe, r_bin)
            print("Validated existing YGB family map:", output)
            return output
        rel_sha, all_sha = sha(rel), sha(allrds)
        ids = cohort_ids(allrds, r_bin)
        groups, stats = build_kinship_groups(read_pairs(rel), ids)
        record = dict(kinship_protocol=KINSHIP_PROTOCOL, source_file=os.path.relpath(rel, output.parent),
                      source_sha256=rel_sha, all_rds_file=os.path.relpath(allrds, output.parent),
                      all_rds_sha256=all_sha, algorithm=KINSHIP_ALGORITHM, source_semantics=KINSHIP_SEMANTICS,
                      grouping_unit="family", grouping_source=GROUPING_SOURCE, **stats)
        if sha(rel) != rel_sha or sha(allrds) != all_sha: raise RuntimeError("YGB inputs changed during generation")
        with tempfile.TemporaryDirectory(dir=output.parent) as td:
            staged = Path(td) / output.name
            with open(staged, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["eid", "group"], delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerows(groups)
            staged.chmod(0o600)
            record["group_map_sha256"] = sha(staged)
            sidecar = Path(td) / meta.name
            sidecar.write_text(json.dumps(record, indent=2) + "\n"); sidecar.chmod(0o600)
            os.replace(staged, output); os.replace(sidecar, meta)
        print(f"Generated {output}: {len(groups)} participants, {stats['component_N']} family components")
        return output

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ygb-phe", required=True)
    p.add_argument("--r-bin", default="Rscript")
    p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE")
    a = p.parse_args()
    generate(a.ygb_phe, a.r_bin, a.replace == "TRUE")

if __name__ == "__main__": main()
