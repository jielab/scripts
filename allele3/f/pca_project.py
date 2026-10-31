#!/usr/bin/env python3
"""Project chromosome PGENs onto an immutable GRCh37 PCA basis.

Only the exact shared, pre-pruned marker panel is scored. Reference and target
people receive the same reference-frequency transform; target people never fit
the axes. Chromosome SUMs and DENOMs are added before calculating a single AVG.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import contextlib
import csv
import fcntl
import glob
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

VERSION = "pca-project-1.1"
BASES = set("ACGT")
NA = {"", ".", "NA", "nan", "NaN"}
DEFAULT_GLOB = "/mnt/f/gen/ukb/37/imp/chr[0-9]*.pgen"
LIMITATIONS = [
    "The target cohort never participates in PCA fitting, frequency estimation or LD pruning.",
    "Both cohorts are re-scored with the identical retained marker subset, weights and reference frequencies; these subset coordinates must not be mixed with old full-panel coordinates.",
    "Mean-imputed missing genotypes contribute zero after reference-frequency standardization. Low call fraction can still shrink a projection; sample QC is reported, not silently interpreted as ancestry.",
    "PLINK uses stored dosages when present and otherwise hardcalls. A dosage is an expected copy count, not a directly observed chromosome count; this operation discovers no additional allele.",
    "File-level dosage presence does not prove that every selected genotype has a stored posterior dosage or that imputation quality passed.",
    "An external INFO table is a user-supplied quality annotation for this target release; no INFO value is reconstructed from hardcalls or assumed from a PGEN filename.",
    "Shared markers do not imply shared imputation pipelines. Original UKB HRC/UK10K imputation and the Reich re-imputation can have different reference, batch and uncertainty effects.",
    "Out-of-reference ancestry and projection shrinkage remain possible. This is a genetic-coordinate projection, not an ancestry proportion, historical ancestor assignment or selection test.",
    "Large input identities use size, modification time and bounded file-edge fingerprints, not whole-file hashes. Original inputs must remain immutable during execution.",
]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for data in iter(lambda: handle.read(1048576), b""): h.update(data)
    return h.hexdigest()


def state(path):
    """Bounded input identity: never hash a whole UKB PGEN/PVAR/INFO table."""
    path = Path(path).expanduser().resolve(); st = path.stat(); h = hashlib.sha256()
    with path.open("rb") as handle:
        h.update(handle.read(65536))
        if st.st_size > 65536:
            handle.seek(max(65536, st.st_size-65536)); h.update(handle.read(65536))
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns, edge_sha256=h.hexdigest())


def atomic_text(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name+".partial")
    part.write_text(value, encoding="utf-8"); os.replace(part, path)


def save_json(path, value):
    atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n")


@contextlib.contextmanager
def open_text(path):
    path = Path(path)
    if str(path).endswith(".zst"):
        try: import zstandard
        except ImportError as exc: raise RuntimeError("Reading .pvar.zst requires zstandard") from exc
        import io
        with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as reader, io.TextIOWrapper(reader) as handle:
            yield handle
    elif str(path).endswith(".gz"):
        with gzip.open(path, "rt") as handle: yield handle
    else:
        # PVARs can contain 100M rows on a mounted drive. Buffer binary lines
        # before decoding to avoid one filesystem RPC per 8 KiB text chunk.
        with path.open("rb", buffering=8 << 20) as handle:
            yield (line.decode("utf-8") for line in handle)


def rows(path, required=()):
    """Whitespace-delimited PLINK tables, preserving all original row indexes."""
    with open_text(path) as handle:
        fields = None
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("##"): continue
            values = line.rstrip().split()
            if fields is None:
                fields = [x.lstrip("#") for x in values]
                if len(set(fields)) != len(fields) or not set(required) <= set(fields):
                    raise ValueError(f"{path}: required unique columns {sorted(required)}")
                continue
            if len(values) != len(fields): raise ValueError(f"{path}:{line_no}: wrong column count")
            yield dict(zip(fields, values))
        if fields is None: raise ValueError(f"Empty table: {path}")


def write_table(path, records, fields):
    path = Path(path); part = path.with_name(path.name+".partial")
    path.parent.mkdir(parents=True, exist_ok=True)
    with (gzip.open(part, "wt", newline="") if str(path).endswith(".gz") else part.open("w", newline="")) as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in records: writer.writerow({k:"NA" if row.get(k) is None else row.get(k) for k in fields})
    os.replace(part, path)


def chrom(value):
    value = value.removeprefix("chr")
    return str(int(value)) if value.isdigit() and 1 <= int(value) <= 22 else None


def numeric(value, label):
    try: value = float(value)
    except (TypeError, ValueError) as exc: raise ValueError(f"{label}: numeric value required") from exc
    if not math.isfinite(value): raise ValueError(f"{label}: finite value required")
    return value


def pgen_files(prefix):
    prefix = Path(prefix).expanduser().resolve()
    pvars = [Path(str(prefix)+suffix) for suffix in (".pvar", ".pvar.zst") if Path(str(prefix)+suffix).is_file()]
    if len(pvars) != 1: raise ValueError(f"Need exactly one .pvar or .pvar.zst for {prefix}")
    files = dict(pgen=str(prefix)+".pgen", pvar=str(pvars[0]), psam=str(prefix)+".psam")
    for path in files.values():
        if not Path(path).is_file(): raise ValueError(f"Missing PGEN input: {path}")
    for suffix in (".pgen.pgi", ".pgi"):
        if Path(str(prefix)+suffix).is_file(): files[suffix] = str(prefix)+suffix
    return files


def sample_ids(path):
    result, seen = [], set()
    for r in rows(path, {"IID"}):
        iid = r["IID"]; fid = r.get("FID", "0")
        if iid in NA or iid in seen: raise ValueError(f"PSAM requires unique nonmissing IID: {path}: {iid}")
        if fid in NA: raise ValueError(f"Invalid PSAM FID: {path}")
        seen.add(iid); result.append((fid, iid))
    if not result: raise ValueError(f"No PSAM samples: {path}")
    return result


def _completed_outputs(basis_dir):
    """Verify the historical completion receipt without imposing new scoring QC."""
    root = Path(basis_dir).expanduser().resolve()
    plan = json.loads((root/"pca_plan.json").read_text())
    if plan.get("build") != "GRCh37": raise ValueError("Projection requires an explicitly GRCh37 frozen PCA basis")
    if (root/"PCA_COMPLETE").read_text().strip() != plan.get("signature"):
        raise ValueError("PCA_COMPLETE does not match the frozen basis signature")
    expected = {"all.sscore", "ref.acount", "ref.eigenvec.allele", "ref.eigenval"}
    checked = set()
    for line in (root/"completed.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1); name = Path(name.lstrip("*")).name
        if name not in expected or name in checked or not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
            raise ValueError("Unexpected/duplicate completed PCA checksum entry")
        if sha(root/name) != digest.lower(): raise ValueError(f"Frozen PCA output changed: {name}")
        checked.add(name)
    if checked != expected: raise ValueError("Frozen PCA completed.sha256 lacks required outputs")
    return root, plan, expected


def validate_completed_coordinates(basis_dir):
    """Read-only validation for re-rendering; does not certify target projection."""
    root, plan, expected = _completed_outputs(basis_dir)
    samples = sample_ids(root/"all_pruned.psam")
    old_ids, pcs = [], None
    for r in rows(root/"all.sscore", {"IID"}):
        if pcs is None:
            pcs = [key for key in r if re.fullmatch(r"PC\d+_AVG", key)]
            if len(pcs) < 2 or pcs != [f"PC{i}_AVG" for i in range(1,len(pcs)+1)]:
                raise ValueError("Completed scores require consecutive PC1_AVG..PCn_AVG columns")
        old_ids.append((r.get("FID","0"),r["IID"]))
        for p in pcs: numeric(r[p],p)
    if old_ids != samples:
        raise ValueError("Frozen all.sscore sample order differs from all_pruned PSAM")
    eigenvalues = [numeric(x,"eigenvalue") for x in (root/"ref.eigenval").read_text().split()]
    if len(eigenvalues) != len(pcs) or any(x <= 0 for x in eigenvalues):
        raise ValueError("Frozen eigenvalues disagree with completed score dimensions")
    return dict(status="verified_completed_coordinates", basis_dir=str(root),
                signature=plan["signature"], n_samples=len(samples), pcs=pcs,
                verified_outputs=sorted(expected), projection_eligibility="not_evaluated")


def _basis(basis_dir, *, allow_duplicate_coordinates=False):
    root, plan, expected = _completed_outputs(basis_dir)
    files = pgen_files(root/"all_pruned")
    positions, markers = {}, {}
    for index, r in enumerate(rows(files["pvar"], {"CHROM", "POS", "ID", "REF", "ALT"})):
        c, pos, ref, alt = chrom(r["CHROM"]), int(r["POS"]), r["REF"], r["ALT"]
        if c is None or pos < 1 or ref not in BASES or alt not in BASES or ref == alt:
            raise ValueError("Frozen all_pruned PVAR must contain autosomal biallelic A/C/G/T SNPs")
        if r["ID"] in NA or r["ID"] in markers or ((c,pos) in positions and not allow_duplicate_coordinates):
            raise ValueError("Frozen all_pruned requires unique variant IDs and coordinates")
        markers[r["ID"]] = dict(chrom=c, pos=pos, id=r["ID"], ref=ref, alt=alt, index=index, weights={})
        positions[c,pos] = r["ID"]
    pcs = None
    for r in rows(root/"ref.eigenvec.allele", {"CHROM", "ID", "REF", "ALT", "A1"}):
        names = [key for key in r if re.fullmatch(r"PC\d+", key)]
        if pcs is None:
            pcs = names
            if len(pcs) < 2 or pcs != [f"PC{i}" for i in range(1,len(pcs)+1)]: raise ValueError("Weights require consecutive PC1..PCn columns")
        marker = markers.get(r["ID"])
        if marker is None or (chrom(r["CHROM"]),r["REF"],r["ALT"]) != (marker["chrom"],marker["ref"],marker["alt"]):
            raise ValueError("Frozen weights disagree with original pruned PVAR identities")
        allele = r["A1"]
        if allele not in {marker["ref"],marker["alt"]} or allele in marker["weights"]:
            raise ValueError("Duplicate or invalid scored allele in frozen weights")
        marker["weights"][allele] = [numeric(r[p],p) for p in pcs]
    for r in rows(root/"ref.acount", {"CHROM", "ID", "REF", "ALT", "ALT_CTS", "OBS_CT"}):
        marker = markers.get(r["ID"])
        if marker is None: continue
        if "alt_count" in marker or (chrom(r["CHROM"]),r["REF"],r["ALT"]) != (marker["chrom"],marker["ref"],marker["alt"]):
            raise ValueError("Duplicate or discordant reference frequency identity")
        ac, an = numeric(r["ALT_CTS"],"ALT_CTS"), numeric(r["OBS_CT"],"OBS_CT")
        if not 0 < ac < an: raise ValueError("PCA variance standardization requires polymorphic reference frequencies")
        marker.update(alt_count=ac, obs_count=an)
    if not markers or any(not m["weights"] or "alt_count" not in m for m in markers.values()):
        raise ValueError("Every pruned marker needs exact allele weights and reference counts")
    eigenvalues = [numeric(x,"eigenvalue") for x in (root/"ref.eigenval").read_text().split()]
    if len(eigenvalues) != len(pcs) or any(x <= 0 for x in eigenvalues):
        raise ValueError("Frozen eigenvalues disagree with allele-weight dimensions")
    samples = sample_ids(files["psam"])
    old_ids = []
    for r in rows(root/"all.sscore", {"IID",*[p+"_AVG" for p in pcs]}):
        old_ids.append((r.get("FID","0"),r["IID"]))
        for p in pcs: numeric(r[p+"_AVG"],p+"_AVG")
    if old_ids != samples:
        raise ValueError("Frozen all.sscore sample order differs from all_pruned PSAM")
    inputs = {**files, **{name:str(root/name) for name in ["pca_plan.json", "PCA_COMPLETE", "completed.sha256", *sorted(expected)]}}
    return dict(root=str(root), plan=plan, files=files, markers=markers, pcs=pcs, samples=samples,
                inputs={key:state(path) for key,path in inputs.items()})


def validate_basis(basis_dir):
    """Public read-only basis check; no genotype decoding, writes, or target access."""
    basis = _basis(basis_dir)
    return dict(status="verified_frozen_basis", basis_dir=basis["root"], build="GRCh37",
                signature=basis["plan"]["signature"], settings=basis["plan"].get("settings",{}),
                n_markers=len(basis["markers"]), n_samples=len(basis["samples"]),
                n_training_samples=basis["plan"].get("counts",{}).get("samples_basis"),pcs=basis["pcs"], inputs=basis["inputs"])


def quality_values(path, positions):
    result = {}
    if not path: return result
    for r in rows(path, {"chrom","pos","ref","alt","info"}):
        c, pos = chrom(r["chrom"]), int(r["pos"])
        if (c,pos) not in positions: continue
        key = c,pos,r["ref"],r["alt"]
        if key in result: raise ValueError(f"Duplicate INFO contrast: {key}")
        value = None if r["info"] in NA else numeric(r["info"],"INFO")
        if value is not None and not 0 <= value <= 1:
            raise ValueError(f"INFO must be in [0,1] for target contrast: {key}")
        result[key] = value
    return result


def intersect_targets(basis, targets, info_path, min_info, allow_unknown, allow_palindromic):
    positions = {(m["chrom"],m["pos"]):m for m in basis["markers"].values()}
    found = defaultdict(list); raw_counts = Counter()
    common_samples = None
    for number, target in enumerate(targets):
        ids = sample_ids(target["psam"])
        if common_samples is None: common_samples = ids
        elif common_samples != ids: raise ValueError("Target chromosome PSAM FID/IID sample order differs; never merge by row number across different people")
        for index,r in enumerate(rows(target["pvar"], {"CHROM","POS","ID","REF","ALT"})):
            raw_counts[str(number)] += 1
            c, pos = chrom(r["CHROM"]), int(r["POS"])
            if (c,pos) in positions:
                found[c,pos].append(dict(target_index=number, pvar_index=index, target_id=r["ID"], target_ref=r["REF"], target_alt=r["ALT"]))
        print(f"[projection plan] inspected PVAR {number+1}/{len(targets)}; {raw_counts[str(number)]:,} rows; no target genotypes decoded",flush=True)
    info = quality_values(info_path,positions)
    audit, selected = [], []
    for marker in basis["markers"].values():
        record = {key:marker[key] for key in ("chrom","pos","id","ref","alt")}
        record.update(status="not_present", quality_status="not_evaluated", info=None)
        hits = found.get((marker["chrom"],marker["pos"]),[])
        record["target_coordinate_records"] = len(hits)
        if not allow_palindromic and {marker["ref"],marker["alt"]} in ({"A","T"},{"C","G"}):
            record["status"] = "palindromic_excluded"
        elif len(hits) > 1: record["status"] = "duplicate_target_coordinate"
        elif hits:
            record.update(hits[0]); pair = {record["target_ref"],record["target_alt"]}
            if pair != {marker["ref"],marker["alt"]}: record["status"] = "allele_pair_mismatch_no_strand_guess"
            elif record["target_id"] in NA: record["status"] = "invalid_target_id"
            else:
                value = info.get((marker["chrom"],marker["pos"],record["target_ref"],record["target_alt"]))
                record["info"] = value
                if value is None:
                    record["quality_status"] = "unknown_explicitly_allowed" if allow_unknown else "unknown_excluded"
                    passed = allow_unknown
                else:
                    passed = value >= min_info
                    record["quality_status"] = "external_info_pass" if passed else "external_info_below_threshold"
                record["status"] = ("matched_ref_alt_swapped" if record["target_ref"] != marker["ref"] else "matched_exact") if passed else "quality_unavailable_or_failed"
                if passed: selected.append(record)
        audit.append(record)
    # Native rsIDs can repeat at coordinates outside the reference panel. A
    # second streaming metadata pass checks only selected IDs (bounded memory).
    selected_ids = defaultdict(set)
    for r in selected: selected_ids[r["target_index"]].add(r["target_id"])
    duplicates = set()
    for i,ids in selected_ids.items():
        counts = Counter(r["ID"] for r in rows(targets[i]["pvar"], {"ID"}) if r["ID"] in ids)
        duplicates.update((i,vid) for vid,count in counts.items() if count != 1)
    for r in selected:
        if (r["target_index"],r["target_id"]) in duplicates: r["status"] = "duplicate_native_target_id"
    return audit, [r for r in selected if r["status"].startswith("matched_")], common_samples, dict(raw_counts)


def score_files(directory, records, basis, target=False):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    weights, counts, ids = [], [], []
    for r in records:
        m = basis["markers"][r["id"]]
        vid, ref, alt = (r["target_id"],r["target_ref"],r["target_alt"]) if target else (m["id"],m["ref"],m["alt"])
        ids.append(vid)
        for allele, values in m["weights"].items(): weights.append(dict(ID=vid,A1=allele,**dict(zip(basis["pcs"],values))))
        ac = m["alt_count"] if alt == m["alt"] else m["obs_count"]-m["alt_count"]
        counts.append(dict(CHROM=m["chrom"],ID=vid,REF=ref,ALT=alt,ALT_CTS=ac,OBS_CT=m["obs_count"]))
    weight_path, count_path, id_path = directory/"weights.tsv", directory/"reference.acount", directory/"markers.ids"
    write_table(weight_path,weights,["ID","A1",*basis["pcs"]])
    write_table(count_path,counts,["CHROM","ID","REF","ALT","ALT_CTS","OBS_CT"])
    # PLINK treats the leading # as optional for imported frequency tables.
    atomic_text(count_path,"#"+count_path.read_text())
    atomic_text(id_path,"".join(vid+"\n" for vid in ids))
    return dict(weights=str(weight_path),frequencies=str(count_path),extract=str(id_path))


def score_command(args, files, prepared, output):
    return [str(args.plink2),"--threads",str(args.threads),"--memory",str(args.memory_mb),
            "--pgen",files["pgen"],"--pvar",files["pvar"],"--psam",files["psam"],
            "--extract",prepared["extract"],"--chr","1-22","--min-alleles","2","--max-alleles","2",
            "--snps-only","just-acgt","--read-freq",prepared["frequencies"],"--score",prepared["weights"],
            "1","2","header-read","variance-standardize","list-variants",
            "cols=fid,nallele,denom,dosagesum,scoreavgs,scoresums","--score-col-nums",f"3-{2+args.components}","--out",str(output)]


def prepare_plan(args):
    if not 0 <= args.min_info <= 1 or not 0 < args.min_overlap <= 1 or not 0 < args.min_call_rate <= 1:
        raise ValueError("Require INFO in [0,1], overlap/call-rate in (0,1]")
    if args.threads < 1 or args.memory_mb < 640: raise ValueError("Positive threads and at least 640 MiB PLINK memory required")
    if not args.target_info and not args.allow_unknown_info:
        raise ValueError("Supply --target-info (chrom,pos,ref,alt,info) or explicitly --allow-unknown-info; PGEN alone does not establish imputation QC")
    basis = _basis(args.basis_dir); out = Path(args.out).expanduser().resolve()
    paths = sorted({Path(p).resolve() for p in glob.glob(os.path.expanduser(args.target_glob))})
    if not paths or any(p.suffix != ".pgen" for p in paths): raise ValueError("--target-glob must match one or more original .pgen files")
    targets = [pgen_files(str(path)[:-5]) for path in paths]
    from pca_coordinate_receipts import metadata_request
    coordinate_metadata = metadata_request(args, basis["plan"])
    inputs = dict(basis=basis["inputs"],targets=[{key:state(path) for key,path in t.items()} for t in targets],
                  info=state(args.target_info) if args.target_info else None, code=state(__file__),
                  coordinate_receipts_code=state(Path(__file__).with_name("pca_coordinate_receipts.py")),
                  coordinate_metadata=coordinate_metadata)
    for path in [Path(basis["root"]),*paths,*([] if not args.target_info else [Path(args.target_info).resolve()])]:
        if out == path or out.is_relative_to(path) or path.is_relative_to(out): raise ValueError("Projection output must be separate from input/basis directories and files")
    args.components = len(basis["pcs"])
    settings = {key:getattr(args,key) for key in ("min_info","min_overlap","min_call_rate","allow_unknown_info","allow_palindromic","plink2","threads","memory_mb")}
    signature = hashlib.sha256(json.dumps(dict(version=VERSION,inputs=inputs,settings=settings),sort_keys=True).encode()).hexdigest()
    plan_path = out/"projection_plan.json"
    if plan_path.exists():
        prior = json.loads(plan_path.read_text())
        if prior.get("signature") != signature: raise ValueError("Projection inputs/settings/code changed; choose a new output directory")
        for path,digest in prior["prepared_sha256"].items():
            if not Path(path).is_file() or sha(path) != digest: raise ValueError(f"Prepared projection file changed: {path}")
        return prior
    if out.exists() and any(p.name != ".projection.lock" for p in out.iterdir()):
        raise ValueError("New projection plan requires an empty output directory")
    audit, selected, samples, raw_counts = intersect_targets(basis,targets,args.target_info,args.min_info,args.allow_unknown_info,args.allow_palindromic)
    eligible = sum(r["status"] != "palindromic_excluded" for r in audit)
    out.mkdir(parents=True,exist_ok=True)
    audit_path = out/"marker_audit.tsv.gz"
    fields = ["chrom","pos","id","ref","alt","status","target_index","pvar_index","target_id","target_ref","target_alt","target_coordinate_records","info","quality_status"]
    write_table(audit_path,audit,fields)
    if not selected or not eligible or len(selected)/eligible < args.min_overlap:
        save_json(out/"projection_manifest.json",dict(status="insufficient_marker_overlap",n_original_markers=len(audit),n_eligible=eligible,n_shared=len(selected),required_fraction=args.min_overlap,marker_audit=str(audit_path)))
        raise ValueError(f"Only {len(selected)}/{eligible} eligible markers shared; below --min-overlap {args.min_overlap}. Inspect marker_audit.tsv.gz; choose a new output for revised inputs")
    stages, generated = [], [audit_path]
    batches = [("reference",basis["files"],selected,False)]
    for i,target in enumerate(targets):
        subset = [r for r in selected if r["target_index"] == i]
        if subset: batches.append((f"target_{i+1:02d}",target,subset,True))
    for label,files,subset,is_target in batches:
        prepared = score_files(out/"prepared"/label,subset,basis,is_target)
        output = out/"scores"/label; output.parent.mkdir(exist_ok=True)
        stages.append(dict(label=label,target=is_target,files=files,prepared=prepared,output=str(output),n_markers=len(subset),
                           command=score_command(args,files,prepared,output)))
        generated.extend(Path(p) for p in prepared.values())
    script = out/"commands.sh"
    # The executable entry verifies receipts, markers, tool identity and merges
    # outputs; this script is an inspectable command list, not the resume engine.
    atomic_text(script,"#!/usr/bin/env bash\nset -euo pipefail\n# Inspect only: execute pca_project.py run for verification and final aggregation.\n"+
                "\n".join(shlex.join(s["command"]) for s in stages)+"\n")
    generated.append(script)
    plan = dict(schema=VERSION,status="planned",signature=signature,inputs=inputs,settings=settings,basis_dir=basis["root"],out=str(out),
        coordinate_metadata=coordinate_metadata,
        build="GRCh37",basis_signature=basis["plan"]["signature"],basis_settings=basis["plan"].get("settings",{}),pcs=basis["pcs"],
        n_reference_samples=len(basis["samples"]),n_target_samples=len(samples),n_training_samples=basis["plan"].get("counts",{}).get("samples_basis"),
        n_original_markers=len(audit),n_eligible_markers=eligible,n_shared_markers=len(selected),
        fraction_of_eligible=len(selected)/eligible,fraction_of_original=len(selected)/len(audit),marker_status_counts=dict(Counter(r["status"] for r in audit)),
        marker_quality_counts=dict(Counter(r["quality_status"] for r in selected)),pvar_rows_inspected=raw_counts,stages=stages,
        coordinate_kind="projected_average_same_common_panel",missing_rule="reference_frequency_mean_imputation; no target frequency estimation",
        allele_matching="exact GRCh37 coordinate and unordered two-base pair; explicit REF/ALT swaps, no strand complementation or rsID guess",
        mean_imputation_denominator="twice the common biallelic marker count, identical across cohorts; chromosome SUM and DENOM aggregated",
        prepared_sha256={str(p):sha(p) for p in generated},limitations=LIMITATIONS)
    save_json(plan_path,plan); save_json(out/"projection_manifest.json",dict(status="planned",signature=signature))
    return plan


def check_output_state(saved):
    for path,expected in saved.items():
        if not Path(path).is_file() or state(path) != expected: raise ValueError(f"Completed projection output changed: {path}")


def execute_stage(stage, plan, executable, tool):
    output = Path(stage["output"]); receipt = Path(str(output)+".receipt.json")
    identity = dict(signature=plan["signature"],label=stage["label"],tool=tool)
    if receipt.exists():
        prior = json.loads(receipt.read_text())
        if prior["identity"] != identity: raise ValueError("Completed projection tool/input identity changed; use a new output directory")
        check_output_state(prior["outputs"]); print(f"[projection resume] {stage['label']}",flush=True); return prior
    command = [executable,*stage["command"][1:]]
    with Path(str(output)+".console.log").open("w") as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    expected = set(Path(stage["prepared"]["extract"]).read_text().splitlines())
    used = Path(str(output)+".sscore.vars").read_text().splitlines()
    if set(used) != expected or len(used) != len(expected): raise ValueError(f"PLINK did not score exactly the agreed marker set: {stage['label']}")
    # Metadata-only .pgen inspection: distinguishes stored dosage support from
    # hardcall-only input without treating it as per-genotype INFO evidence.
    command = [executable,"--pgen",stage["files"]["pgen"],"--pgen-info","--threads","1","--memory",str(plan["settings"]["memory_mb"]),"--out",str(output)+".storage"]
    result = subprocess.run(command,capture_output=True,text=True,check=True)
    storage_path = Path(str(output)+".storage.txt"); atomic_text(storage_path,result.stdout+result.stderr)
    storage = "hardcalls_only" if "No dosages present" in result.stdout else "stored_dosages_present" if re.search(r"dosages? present",result.stdout,re.I) else "storage_status_not_recognized"
    paths = [Path(str(output)+suffix) for suffix in (".sscore",".sscore.vars",".log")]+[storage_path]
    receipt_data = dict(status="complete",identity=identity,storage_evidence=storage,outputs={str(p):state(p) for p in paths})
    save_json(receipt,receipt_data); print(f"[projection complete] {stage['label']}: {len(used):,} common markers",flush=True)
    return receipt_data


def aggregate_scores(paths, samples, pcs, n_markers, out_path, qc_path, min_call_rate):
    """Streaming aggregation: O(number of chromosomes * PCs) numeric memory."""
    expected_denom = 2*n_markers
    generators = [rows(path,{"FID","IID","ALLELE_CT","DENOM","NAMED_ALLELE_DOSAGE_SUM",*[p+"_SUM" for p in pcs]}) for path in paths]
    output, qc = Path(out_path),Path(qc_path)
    score_fields = ["#FID","IID","ALLELE_CT","DENOM","NAMED_ALLELE_DOSAGE_SUM",*[p+"_AVG" for p in pcs],*[p+"_SUM" for p in pcs]]
    n, low = 0,0
    with output.with_name(output.name+".partial").open("w") as dest, qc.with_name(qc.name+".partial").open("w") as quality:
        writer = csv.writer(dest,delimiter="\t",lineterminator="\n"); writer.writerow(score_fields)
        qwriter = csv.writer(quality,delimiter="\t",lineterminator="\n"); qwriter.writerow(["FID","IID","n_markers","observed_alleles","denominator","call_fraction","status"])
        for sample,parts in itertools.zip_longest(samples,itertools.zip_longest(*generators)):
            if sample is None or parts is None or any(r is None for r in parts): raise ValueError("Score row count differs from original PSAM")
            if any((r["FID"],r["IID"]) != sample for r in parts): raise ValueError("Score sample order differs from original PSAM")
            denom = sum(numeric(r["DENOM"],"DENOM") for r in parts)
            observed = sum(numeric(r["ALLELE_CT"],"ALLELE_CT") for r in parts)
            named = sum(numeric(r["NAMED_ALLELE_DOSAGE_SUM"],"NAMED_ALLELE_DOSAGE_SUM") for r in parts)
            if denom != expected_denom or not 0 <= observed <= denom: raise ValueError("Unexpected score denominator; cohort marker sets or mean-imputation rule differ")
            sums = [sum(numeric(r[p+"_SUM"],p+"_SUM") for r in parts) for p in pcs]
            writer.writerow([*sample,int(observed),int(denom),format(named,".12g"),*[format(v/denom,".12g") for v in sums],*[format(v,".12g") for v in sums]])
            status = "pass_call_fraction" if observed/denom >= min_call_rate else "low_call_fraction_projection_shrinkage_risk"
            qwriter.writerow([*sample,n_markers,int(observed),int(denom),format(observed/denom,".8g"),status])
            n += 1; low += int(observed/denom < min_call_rate)
    os.replace(output.with_name(output.name+".partial"),output); os.replace(qc.with_name(qc.name+".partial"),qc)
    return dict(n_samples=n,n_low_call_fraction=low,denominator=expected_denom)


def geometry_change(old_path, new_path, pcs, out):
    import numpy as np
    old = {r["IID"]:r for r in rows(old_path,{"IID",*[p+"_AVG" for p in pcs]})}
    pairs = [(old[r["IID"]],r) for r in rows(new_path,{"IID"}) if r["IID"] in old]
    result = []
    for p in pcs:
        x = np.asarray([float(a[p+"_AVG"]) for a,_ in pairs]); y = np.asarray([float(b[p+"_AVG"]) for _,b in pairs])
        centered = x-x.mean(); spread = float(centered@centered)
        slope = float(centered@(y-y.mean())/spread) if spread else None
        result.append(dict(pc=p,n=len(pairs),pearson_r=float(np.corrcoef(x,y)[0,1]) if spread and np.var(y)>0 else None,
                           slope_new_on_old=slope,intercept=float(y.mean()-slope*x.mean()) if slope is not None else None,
                           rms_difference=float(np.sqrt(np.mean((y-x)**2))),note="Descriptive geometry audit; no rescaling or correction applied"))
    write_table(out,result,["pc","n","pearson_r","slope_new_on_old","intercept","rms_difference","note"])
    return result


def run(args, *, prepared_plan=None):
    plan = prepare_plan(args) if prepared_plan is None else prepared_plan
    out = Path(plan["out"])
    executable = shutil.which(str(args.plink2))
    if not executable: raise ValueError(f"PLINK2 executable unavailable: {args.plink2}")
    version = subprocess.run([executable,"--version"],capture_output=True,text=True,check=True).stdout.strip()
    tool = dict(version=version,executable=state(executable))
    completion = out/"PROJECTION_COMPLETE.json"
    if completion.exists():
        saved = json.loads(completion.read_text())
        if saved["signature"] != plan["signature"] or saved["tool"] != tool: raise ValueError("Completed projection identity changed; use a new output directory")
        check_output_state(saved["outputs"])
        from pca_coordinate_receipts import verify_projection_receipt
        verify_projection_receipt(out/"reference.sscore", out/"reference.coordinate_receipt.json",
                                  [p+"_AVG" for p in plan["pcs"]], plan["basis_dir"])
        return json.loads((out/"projection_manifest.json").read_text())
    save_json(out/"projection_manifest.json",dict(status="running",signature=plan["signature"]))
    receipts = {s["label"]:execute_stage(s,plan,executable,tool) for s in plan["stages"]}
    refs = [s for s in plan["stages"] if not s["target"]]; targets = [s for s in plan["stages"] if s["target"]]
    summaries = {}
    for label,stages in (("reference",refs),("target",targets)):
        samples = sample_ids(stages[0]["files"]["psam"])
        summaries[label] = aggregate_scores([s["output"]+".sscore" for s in stages],samples,plan["pcs"],plan["n_shared_markers"],
            out/f"{label}.sscore",out/f"{label}_sample_qc.tsv",args.min_call_rate)
    geometry = geometry_change(Path(plan["basis_dir"])/"all.sscore",out/"reference.sscore",plan["pcs"],out/"reference_geometry_change.tsv")
    from pca_coordinate_receipts import publish_projection_receipts
    receipt_outputs = publish_projection_receipts(plan, receipts)
    coordinate_frame_id = json.loads((out/"coordinate_frame.json").read_text())["coordinate_frame_id"]
    manifest = dict(status="complete_frozen_basis_projection",signature=plan["signature"],basis_signature=plan["basis_signature"],
        coordinate_frame_id=coordinate_frame_id,
        coordinate_receipts={role:str(out/f"{role}.coordinate_receipt.json") for role in ("reference","target")},
        build="GRCh37",tool=tool,counts={k:plan[k] for k in ("n_original_markers","n_eligible_markers","n_shared_markers","n_reference_samples","n_target_samples","n_training_samples")},
        scores=summaries,coordinate_kind=plan["coordinate_kind"],marker_quality_counts=plan["marker_quality_counts"],
        dosage_storage={k:r["storage_evidence"] for k,r in receipts.items()},geometry_change=geometry,
        reference_coordinates=str(out/"reference.sscore"),target_coordinates=str(out/"target.sscore"),
        fit_cohort="frozen reference basis only; no UKB fit",imputation_pipeline_equivalence="not established",limitations=LIMITATIONS)
    save_json(out/"projection_manifest.json",manifest)
    outputs = [out/name for name in ("reference.sscore","target.sscore","reference_sample_qc.tsv","target_sample_qc.tsv","reference_geometry_change.tsv","projection_manifest.json")]
    outputs.extend(receipt_outputs)
    outputs.extend(Path(s["output"]+".receipt.json") for s in plan["stages"])
    states = {str(p):state(p) for p in outputs}
    for receipt in receipts.values(): states.update(receipt["outputs"])
    save_json(completion,dict(signature=plan["signature"],tool=tool,outputs=states))
    return manifest


def parser():
    p = argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
    p.add_argument("command",choices=("plan","run")); p.add_argument("--basis-dir",required=True,type=Path)
    p.add_argument("--target-glob",default=DEFAULT_GLOB); p.add_argument("--out",required=True,type=Path)
    p.add_argument("--reference-metadata",type=Path,help="Reference metadata; retain frozen historical provenance checks")
    p.add_argument("--target-metadata",type=Path,help="Optional explicit target metadata snapshot with sample_id")
    p.add_argument("--provenance-manifest",type=Path,help="Existing historical reference-metadata provenance receipt")
    p.add_argument("--target-info",type=Path,help="Exact target GRCh37 chrom,pos,ref,alt,info TSV; gzip/zstd accepted")
    p.add_argument("--min-info",type=float,default=.8); p.add_argument("--allow-unknown-info",action="store_true")
    p.add_argument("--allow-palindromic",action="store_true",help="Only when both inputs have independently verified common reference-strand labels")
    p.add_argument("--min-overlap",type=float,default=.9); p.add_argument("--min-call-rate",type=float,default=.9)
    p.add_argument("--plink2",default="plink2"); p.add_argument("--threads",type=int,default=4); p.add_argument("--memory-mb",type=int,default=8192)
    return p


def main(argv=None):
    args = parser().parse_args(argv); out = args.out.expanduser().resolve()
    try:
        basis_root = args.basis_dir.expanduser().resolve()
        if out == basis_root or out.is_relative_to(basis_root) or basis_root.is_relative_to(out):
            raise ValueError("Projection output must be separate from the frozen basis")
        out.mkdir(parents=True,exist_ok=True)
        with (out/".projection.lock").open("a") as lock:
            try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError: raise ValueError("Another projection job is using this output directory") from None
            # Only this lock owner may update status, and only after its inputs
            # have been matched to the saved plan. A rejected concurrent or
            # mismatched invocation must not overwrite another run's manifest.
            plan = prepare_plan(args)
            try:
                result = plan if args.command == "plan" else run(args, prepared_plan=plan)
            except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
                if args.command == "run" and not (out/"PROJECTION_COMPLETE.json").exists():
                    save_json(out/"projection_manifest.json",dict(status="failed",signature=plan["signature"],reason=str(exc)))
                raise
        print(json.dumps({k:result[k] for k in ("status","signature")},indent=2)); return 0
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
        print(f"pca-project: {exc}",file=sys.stderr); return 2


if __name__ == "__main__": raise SystemExit(main())
