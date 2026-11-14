#!/usr/bin/env python3
"""Consolidated pca workflow; internal stages share one explicit module."""
from __future__ import annotations


# ========================================================================
# pca project
# ========================================================================
_description_pca_project = 'Project chromosome PGENs onto an immutable GRCh37 PCA basis.\n\nOnly the exact shared, pre-pruned marker panel is scored. Reference and target\npeople receive the same reference-frequency transform; target people never fit\nthe axes. Chromosome SUMs and DENOMs are added before calculating a single AVG.'
#!/usr/bin/env python3
"""Project chromosome PGENs onto an immutable GRCh37 PCA basis.

Only the exact shared, pre-pruned marker panel is scored. Reference and target
people receive the same reference-frequency transform; target people never fit
the axes. Chromosome SUMs and DENOMs are added before calculating a single AVG.
"""

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
project_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(project_paths_root / "0f"))
from workflow_paths import configured_path

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

PROJECT_VERSION = "pca-project-1.1"
PROJECT_BASES = set("ACGT")
PROJECT_NA = {"", ".", "NA", "nan", "NaN"}
PROJECT_DEFAULT_GLOB = configured_path('DIR_GEN', 'ukb/37/imp/chr[0-9]*.pgen')
PROJECT_LIMITATIONS = [
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


def project_sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for data in iter(lambda: handle.read(1048576), b""): h.update(data)
    return h.hexdigest()


def project_state(path):
    """Bounded input identity: never hash a whole UKB PGEN/PVAR/INFO table."""
    path = Path(path).expanduser().resolve(); st = path.stat(); h = hashlib.sha256()
    with path.open("rb") as handle:
        h.update(handle.read(65536))
        if st.st_size > 65536:
            handle.seek(max(65536, st.st_size-65536)); h.update(handle.read(65536))
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns, edge_sha256=h.hexdigest())


def project_atomic_text(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name+".partial")
    part.write_text(value, encoding="utf-8"); os.replace(part, path)


def project_save_json(path, value):
    project_atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n")


@contextlib.contextmanager
def project_open_text(path):
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


def project_rows(path, required=()):
    """Whitespace-delimited PLINK tables, preserving all original row indexes."""
    with project_open_text(path) as handle:
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


def project_write_table(path, records, fields):
    path = Path(path); part = path.with_name(path.name+".partial")
    path.parent.mkdir(parents=True, exist_ok=True)
    with (gzip.open(part, "wt", newline="") if str(path).endswith(".gz") else part.open("w", newline="")) as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in records: writer.writerow({k:"NA" if row.get(k) is None else row.get(k) for k in fields})
    os.replace(part, path)


def project_chrom(value):
    value = value.removeprefix("chr")
    return str(int(value)) if value.isdigit() and 1 <= int(value) <= 22 else None


def project_numeric(value, label):
    try: value = float(value)
    except (TypeError, ValueError) as exc: raise ValueError(f"{label}: numeric value required") from exc
    if not math.isfinite(value): raise ValueError(f"{label}: finite value required")
    return value


def project_pgen_files(prefix):
    prefix = Path(prefix).expanduser().resolve()
    pvars = [Path(str(prefix)+suffix) for suffix in (".pvar", ".pvar.zst") if Path(str(prefix)+suffix).is_file()]
    if len(pvars) != 1: raise ValueError(f"Need exactly one .pvar or .pvar.zst for {prefix}")
    files = dict(pgen=str(prefix)+".pgen", pvar=str(pvars[0]), psam=str(prefix)+".psam")
    for path in files.values():
        if not Path(path).is_file(): raise ValueError(f"Missing PGEN input: {path}")
    for suffix in (".pgen.pgi", ".pgi"):
        if Path(str(prefix)+suffix).is_file(): files[suffix] = str(prefix)+suffix
    return files


def project_sample_ids(path):
    result, seen = [], set()
    for r in project_rows(path, {"IID"}):
        iid = r["IID"]; fid = r.get("FID", "0")
        if iid in PROJECT_NA or iid in seen: raise ValueError(f"PSAM requires unique nonmissing IID: {path}: {iid}")
        if fid in PROJECT_NA: raise ValueError(f"Invalid PSAM FID: {path}")
        seen.add(iid); result.append((fid, iid))
    if not result: raise ValueError(f"No PSAM samples: {path}")
    return result


def project_completed_outputs(basis_dir):
    """Verify the historical completion receipt without imposing new scoring QC."""
    from analysis_integrity import resolve_result_input
    root = Path(resolve_result_input(basis_dir)).expanduser().resolve()
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
        if project_sha(root/name) != digest.lower(): raise ValueError(f"Frozen PCA output changed: {name}")
        checked.add(name)
    if checked != expected: raise ValueError("Frozen PCA completed.sha256 lacks required outputs")
    return root, plan, expected


def project_validate_completed_coordinates(basis_dir):
    """Read-only validation for re-rendering; does not certify target projection."""
    root, plan, expected = project_completed_outputs(basis_dir)
    samples = project_sample_ids(root/"all_pruned.psam")
    old_ids, pcs = [], None
    for r in project_rows(root/"all.sscore", {"IID"}):
        if pcs is None:
            pcs = [key for key in r if re.fullmatch(r"PC\d+_AVG", key)]
            if len(pcs) < 2 or pcs != [f"PC{i}_AVG" for i in range(1,len(pcs)+1)]:
                raise ValueError("Completed scores require consecutive PC1_AVG..PCn_AVG columns")
        old_ids.append((r.get("FID","0"),r["IID"]))
        for p in pcs: project_numeric(r[p],p)
    if old_ids != samples:
        raise ValueError("Frozen all.sscore sample order differs from all_pruned PSAM")
    eigenvalues = [project_numeric(x,"eigenvalue") for x in (root/"ref.eigenval").read_text().split()]
    if len(eigenvalues) != len(pcs) or any(x <= 0 for x in eigenvalues):
        raise ValueError("Frozen eigenvalues disagree with completed score dimensions")
    return dict(status="verified_completed_coordinates", basis_dir=str(root),
                signature=plan["signature"], n_samples=len(samples), pcs=pcs,
                verified_outputs=sorted(expected), projection_eligibility="not_evaluated")


def project_basis(basis_dir, *, allow_duplicate_coordinates=False):
    root, plan, expected = project_completed_outputs(basis_dir)
    files = project_pgen_files(root/"all_pruned")
    positions, markers = {}, {}
    for index, r in enumerate(project_rows(files["pvar"], {"CHROM", "POS", "ID", "REF", "ALT"})):
        c, pos, ref, alt = project_chrom(r["CHROM"]), int(r["POS"]), r["REF"], r["ALT"]
        if c is None or pos < 1 or ref not in PROJECT_BASES or alt not in PROJECT_BASES or ref == alt:
            raise ValueError("Frozen all_pruned PVAR must contain autosomal biallelic A/C/G/T SNPs")
        if r["ID"] in PROJECT_NA or r["ID"] in markers or ((c,pos) in positions and not allow_duplicate_coordinates):
            raise ValueError("Frozen all_pruned requires unique variant IDs and coordinates")
        markers[r["ID"]] = dict(chrom=c, pos=pos, id=r["ID"], ref=ref, alt=alt, index=index, weights={})
        positions[c,pos] = r["ID"]
    pcs = None
    for r in project_rows(root/"ref.eigenvec.allele", {"CHROM", "ID", "REF", "ALT", "A1"}):
        names = [key for key in r if re.fullmatch(r"PC\d+", key)]
        if pcs is None:
            pcs = names
            if len(pcs) < 2 or pcs != [f"PC{i}" for i in range(1,len(pcs)+1)]: raise ValueError("Weights require consecutive PC1..PCn columns")
        marker = markers.get(r["ID"])
        if marker is None or (project_chrom(r["CHROM"]),r["REF"],r["ALT"]) != (marker["chrom"],marker["ref"],marker["alt"]):
            raise ValueError("Frozen weights disagree with original pruned PVAR identities")
        allele = r["A1"]
        if allele not in {marker["ref"],marker["alt"]} or allele in marker["weights"]:
            raise ValueError("Duplicate or invalid scored allele in frozen weights")
        marker["weights"][allele] = [project_numeric(r[p],p) for p in pcs]
    for r in project_rows(root/"ref.acount", {"CHROM", "ID", "REF", "ALT", "ALT_CTS", "OBS_CT"}):
        marker = markers.get(r["ID"])
        if marker is None: continue
        if "alt_count" in marker or (project_chrom(r["CHROM"]),r["REF"],r["ALT"]) != (marker["chrom"],marker["ref"],marker["alt"]):
            raise ValueError("Duplicate or discordant reference frequency identity")
        ac, an = project_numeric(r["ALT_CTS"],"ALT_CTS"), project_numeric(r["OBS_CT"],"OBS_CT")
        if not 0 < ac < an: raise ValueError("PCA variance standardization requires polymorphic reference frequencies")
        marker.update(alt_count=ac, obs_count=an)
    if not markers or any(not m["weights"] or "alt_count" not in m for m in markers.values()):
        raise ValueError("Every pruned marker needs exact allele weights and reference counts")
    eigenvalues = [project_numeric(x,"eigenvalue") for x in (root/"ref.eigenval").read_text().split()]
    if len(eigenvalues) != len(pcs) or any(x <= 0 for x in eigenvalues):
        raise ValueError("Frozen eigenvalues disagree with allele-weight dimensions")
    samples = project_sample_ids(files["psam"])
    old_ids = []
    for r in project_rows(root/"all.sscore", {"IID",*[p+"_AVG" for p in pcs]}):
        old_ids.append((r.get("FID","0"),r["IID"]))
        for p in pcs: project_numeric(r[p+"_AVG"],p+"_AVG")
    if old_ids != samples:
        raise ValueError("Frozen all.sscore sample order differs from all_pruned PSAM")
    inputs = {**files, **{name:str(root/name) for name in ["pca_plan.json", "PCA_COMPLETE", "completed.sha256", *sorted(expected)]}}
    return dict(root=str(root), plan=plan, files=files, markers=markers, pcs=pcs, samples=samples,
                inputs={key:project_state(path) for key,path in inputs.items()})


def project_validate_basis(basis_dir):
    """Public read-only basis check; no genotype decoding, writes, or target access."""
    basis = project_basis(basis_dir)
    return dict(status="verified_frozen_basis", basis_dir=basis["root"], build="GRCh37",
                signature=basis["plan"]["signature"], settings=basis["plan"].get("settings",{}),
                n_markers=len(basis["markers"]), n_samples=len(basis["samples"]),
                n_training_samples=basis["plan"].get("counts",{}).get("samples_basis"),pcs=basis["pcs"], inputs=basis["inputs"])


def project_quality_values(path, positions):
    result = {}
    if not path: return result
    for r in project_rows(path, {"chrom","pos","ref","alt","info"}):
        c, pos = project_chrom(r["chrom"]), int(r["pos"])
        if (c,pos) not in positions: continue
        key = c,pos,r["ref"],r["alt"]
        if key in result: raise ValueError(f"Duplicate INFO contrast: {key}")
        value = None if r["info"] in PROJECT_NA else project_numeric(r["info"],"INFO")
        if value is not None and not 0 <= value <= 1:
            raise ValueError(f"INFO must be in [0,1] for target contrast: {key}")
        result[key] = value
    return result


def project_intersect_targets(basis, targets, info_path, min_info, allow_unknown, allow_palindromic):
    positions = {(m["chrom"],m["pos"]):m for m in basis["markers"].values()}
    found = defaultdict(list); raw_counts = Counter()
    common_samples = None
    for number, target in enumerate(targets):
        ids = project_sample_ids(target["psam"])
        if common_samples is None: common_samples = ids
        elif common_samples != ids: raise ValueError("Target chromosome PSAM FID/IID sample order differs; never merge by row number across different people")
        for index,r in enumerate(project_rows(target["pvar"], {"CHROM","POS","ID","REF","ALT"})):
            raw_counts[str(number)] += 1
            c, pos = project_chrom(r["CHROM"]), int(r["POS"])
            if (c,pos) in positions:
                found[c,pos].append(dict(target_index=number, pvar_index=index, target_id=r["ID"], target_ref=r["REF"], target_alt=r["ALT"]))
        print(f"[projection plan] inspected PVAR {number+1}/{len(targets)}; {raw_counts[str(number)]:,} rows; no target genotypes decoded",flush=True)
    info = project_quality_values(info_path,positions)
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
            elif record["target_id"] in PROJECT_NA: record["status"] = "invalid_target_id"
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
        counts = Counter(r["ID"] for r in project_rows(targets[i]["pvar"], {"ID"}) if r["ID"] in ids)
        duplicates.update((i,vid) for vid,count in counts.items() if count != 1)
    for r in selected:
        if (r["target_index"],r["target_id"]) in duplicates: r["status"] = "duplicate_native_target_id"
    return audit, [r for r in selected if r["status"].startswith("matched_")], common_samples, dict(raw_counts)


def project_score_files(directory, records, basis, target=False):
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
    project_write_table(weight_path,weights,["ID","A1",*basis["pcs"]])
    project_write_table(count_path,counts,["CHROM","ID","REF","ALT","ALT_CTS","OBS_CT"])
    # PLINK treats the leading # as optional for imported frequency tables.
    project_atomic_text(count_path,"#"+count_path.read_text())
    project_atomic_text(id_path,"".join(vid+"\n" for vid in ids))
    return dict(weights=str(weight_path),frequencies=str(count_path),extract=str(id_path))


def project_score_command(args, files, prepared, output):
    return [str(args.plink2),"--threads",str(args.threads),"--memory",str(args.memory_mb),
            "--pgen",files["pgen"],"--pvar",files["pvar"],"--psam",files["psam"],
            "--extract",prepared["extract"],"--chr","1-22","--min-alleles","2","--max-alleles","2",
            "--snps-only","just-acgt","--read-freq",prepared["frequencies"],"--score",prepared["weights"],
            "1","2","header-read","variance-standardize","list-variants",
            "cols=fid,nallele,denom,dosagesum,scoreavgs,scoresums","--score-col-nums",f"3-{2+args.components}","--out",str(output)]


def project_prepare_plan(args):
    if not 0 <= args.min_info <= 1 or not 0 < args.min_overlap <= 1 or not 0 < args.min_call_rate <= 1:
        raise ValueError("Require INFO in [0,1], overlap/call-rate in (0,1]")
    if args.threads < 1 or args.memory_mb < 640: raise ValueError("Positive threads and at least 640 MiB PLINK memory required")
    if not args.target_info and not args.allow_unknown_info:
        raise ValueError("Supply --target-info (chrom,pos,ref,alt,info) or explicitly --allow-unknown-info; PGEN alone does not establish imputation QC")
    basis = project_basis(args.basis_dir); out = Path(args.out).expanduser().resolve()
    paths = sorted({Path(p).resolve() for p in glob.glob(os.path.expanduser(args.target_glob))})
    if not paths or any(p.suffix != ".pgen" for p in paths): raise ValueError("--target-glob must match one or more original .pgen files")
    targets = [project_pgen_files(str(path)[:-5]) for path in paths]
    coordinate_metadata = receipt_metadata_request(args, basis["plan"])
    inputs = dict(basis=basis["inputs"],targets=[{key:project_state(path) for key,path in t.items()} for t in targets],
                  info=project_state(args.target_info) if args.target_info else None, code=project_state(__file__),
                  coordinate_receipts_code=project_state(Path(__file__).with_name('pca.py')),
                  coordinate_metadata=coordinate_metadata)
    for path in [Path(basis["root"]),*paths,*([] if not args.target_info else [Path(args.target_info).resolve()])]:
        if out == path or out.is_relative_to(path) or path.is_relative_to(out): raise ValueError("Projection output must be separate from input/basis directories and files")
    args.components = len(basis["pcs"])
    settings = {key:getattr(args,key) for key in ("min_info","min_overlap","min_call_rate","allow_unknown_info","allow_palindromic","plink2","threads","memory_mb")}
    signature = hashlib.sha256(json.dumps(dict(version=PROJECT_VERSION,inputs=inputs,settings=settings),sort_keys=True).encode()).hexdigest()
    plan_path = out/"projection_plan.json"
    if plan_path.exists():
        prior = json.loads(plan_path.read_text())
        if prior.get("signature") != signature: raise ValueError("Projection inputs/settings/code changed; choose a new output directory")
        for path,digest in prior["prepared_sha256"].items():
            if not Path(path).is_file() or project_sha(path) != digest: raise ValueError(f"Prepared projection file changed: {path}")
        return prior
    if out.exists() and any(p.name != ".projection.lock" for p in out.iterdir()):
        raise ValueError("New projection plan requires an empty output directory")
    audit, selected, samples, raw_counts = project_intersect_targets(basis,targets,args.target_info,args.min_info,args.allow_unknown_info,args.allow_palindromic)
    eligible = sum(r["status"] != "palindromic_excluded" for r in audit)
    out.mkdir(parents=True,exist_ok=True)
    audit_path = out/"marker_audit.tsv.gz"
    fields = ["chrom","pos","id","ref","alt","status","target_index","pvar_index","target_id","target_ref","target_alt","target_coordinate_records","info","quality_status"]
    project_write_table(audit_path,audit,fields)
    if not selected or not eligible or len(selected)/eligible < args.min_overlap:
        project_save_json(out/"projection_manifest.json",dict(status="insufficient_marker_overlap",n_original_markers=len(audit),n_eligible=eligible,n_shared=len(selected),required_fraction=args.min_overlap,marker_audit=str(audit_path)))
        raise ValueError(f"Only {len(selected)}/{eligible} eligible markers shared; below --min-overlap {args.min_overlap}. Inspect marker_audit.tsv.gz; choose a new output for revised inputs")
    stages, generated = [], [audit_path]
    batches = [("reference",basis["files"],selected,False)]
    for i,target in enumerate(targets):
        subset = [r for r in selected if r["target_index"] == i]
        if subset: batches.append((f"target_{i+1:02d}",target,subset,True))
    for label,files,subset,is_target in batches:
        prepared = project_score_files(out/"prepared"/label,subset,basis,is_target)
        output = out/"scores"/label; output.parent.mkdir(exist_ok=True)
        stages.append(dict(label=label,target=is_target,files=files,prepared=prepared,output=str(output),n_markers=len(subset),
                           command=project_score_command(args,files,prepared,output)))
        generated.extend(Path(p) for p in prepared.values())
    script = out/"commands.sh"
    # The executable entry verifies receipts, markers, tool identity and merges
    # outputs; this script is an inspectable command list, not the resume engine.
    project_atomic_text(script,'#!/usr/bin/env bash\nset -euo pipefail\n# Inspect only: execute pca.py run for verification and final aggregation.\n'+
                "\n".join(shlex.join(s["command"]) for s in stages)+"\n")
    generated.append(script)
    plan = dict(schema=PROJECT_VERSION,status="planned",signature=signature,inputs=inputs,settings=settings,basis_dir=basis["root"],out=str(out),
        coordinate_metadata=coordinate_metadata,
        build="GRCh37",basis_signature=basis["plan"]["signature"],basis_settings=basis["plan"].get("settings",{}),pcs=basis["pcs"],
        n_reference_samples=len(basis["samples"]),n_target_samples=len(samples),n_training_samples=basis["plan"].get("counts",{}).get("samples_basis"),
        n_original_markers=len(audit),n_eligible_markers=eligible,n_shared_markers=len(selected),
        fraction_of_eligible=len(selected)/eligible,fraction_of_original=len(selected)/len(audit),marker_status_counts=dict(Counter(r["status"] for r in audit)),
        marker_quality_counts=dict(Counter(r["quality_status"] for r in selected)),pvar_rows_inspected=raw_counts,stages=stages,
        coordinate_kind="projected_average_same_common_panel",missing_rule="reference_frequency_mean_imputation; no target frequency estimation",
        allele_matching="exact GRCh37 coordinate and unordered two-base pair; explicit REF/ALT swaps, no strand complementation or rsID guess",
        mean_imputation_denominator="twice the common biallelic marker count, identical across cohorts; chromosome SUM and DENOM aggregated",
        prepared_sha256={str(p):project_sha(p) for p in generated},limitations=PROJECT_LIMITATIONS)
    project_save_json(plan_path,plan); project_save_json(out/"projection_manifest.json",dict(status="planned",signature=signature))
    return plan


def project_check_output_state(saved):
    for path,expected in saved.items():
        if not Path(path).is_file() or project_state(path) != expected: raise ValueError(f"Completed projection output changed: {path}")


def project_execute_stage(stage, plan, executable, tool):
    output = Path(stage["output"]); receipt = Path(str(output)+".receipt.json")
    identity = dict(signature=plan["signature"],label=stage["label"],tool=tool)
    if receipt.exists():
        prior = json.loads(receipt.read_text())
        if prior["identity"] != identity: raise ValueError("Completed projection tool/input identity changed; use a new output directory")
        project_check_output_state(prior["outputs"]); print(f"[projection resume] {stage['label']}",flush=True); return prior
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
    storage_path = Path(str(output)+".storage.txt"); project_atomic_text(storage_path,result.stdout+result.stderr)
    storage = "hardcalls_only" if "No dosages present" in result.stdout else "stored_dosages_present" if re.search(r"dosages? present",result.stdout,re.I) else "storage_status_not_recognized"
    paths = [Path(str(output)+suffix) for suffix in (".sscore",".sscore.vars",".log")]+[storage_path]
    receipt_data = dict(status="complete",identity=identity,storage_evidence=storage,outputs={str(p):project_state(p) for p in paths})
    project_save_json(receipt,receipt_data); print(f"[projection complete] {stage['label']}: {len(used):,} common markers",flush=True)
    return receipt_data


def project_aggregate_scores(paths, samples, pcs, n_markers, out_path, qc_path, min_call_rate):
    """Streaming aggregation: O(number of chromosomes * PCs) numeric memory."""
    expected_denom = 2*n_markers
    generators = [project_rows(path,{"FID","IID","ALLELE_CT","DENOM","NAMED_ALLELE_DOSAGE_SUM",*[p+"_SUM" for p in pcs]}) for path in paths]
    output, qc = Path(out_path),Path(qc_path)
    score_fields = ["#FID","IID","ALLELE_CT","DENOM","NAMED_ALLELE_DOSAGE_SUM",*[p+"_AVG" for p in pcs],*[p+"_SUM" for p in pcs]]
    n, low = 0,0
    with output.with_name(output.name+".partial").open("w") as dest, qc.with_name(qc.name+".partial").open("w") as quality:
        writer = csv.writer(dest,delimiter="\t",lineterminator="\n"); writer.writerow(score_fields)
        qwriter = csv.writer(quality,delimiter="\t",lineterminator="\n"); qwriter.writerow(["FID","IID","n_markers","observed_alleles","denominator","call_fraction","status"])
        for sample,parts in itertools.zip_longest(samples,itertools.zip_longest(*generators)):
            if sample is None or parts is None or any(r is None for r in parts): raise ValueError("Score row count differs from original PSAM")
            if any((r["FID"],r["IID"]) != sample for r in parts): raise ValueError("Score sample order differs from original PSAM")
            denom = sum(project_numeric(r["DENOM"],"DENOM") for r in parts)
            observed = sum(project_numeric(r["ALLELE_CT"],"ALLELE_CT") for r in parts)
            named = sum(project_numeric(r["NAMED_ALLELE_DOSAGE_SUM"],"NAMED_ALLELE_DOSAGE_SUM") for r in parts)
            if denom != expected_denom or not 0 <= observed <= denom: raise ValueError("Unexpected score denominator; cohort marker sets or mean-imputation rule differ")
            sums = [sum(project_numeric(r[p+"_SUM"],p+"_SUM") for r in parts) for p in pcs]
            writer.writerow([*sample,int(observed),int(denom),format(named,".12g"),*[format(v/denom,".12g") for v in sums],*[format(v,".12g") for v in sums]])
            status = "pass_call_fraction" if observed/denom >= min_call_rate else "low_call_fraction_projection_shrinkage_risk"
            qwriter.writerow([*sample,n_markers,int(observed),int(denom),format(observed/denom,".8g"),status])
            n += 1; low += int(observed/denom < min_call_rate)
    os.replace(output.with_name(output.name+".partial"),output); os.replace(qc.with_name(qc.name+".partial"),qc)
    return dict(n_samples=n,n_low_call_fraction=low,denominator=expected_denom)


def project_geometry_change(old_path, new_path, pcs, out):
    import numpy as np
    old = {r["IID"]:r for r in project_rows(old_path,{"IID",*[p+"_AVG" for p in pcs]})}
    pairs = [(old[r["IID"]],r) for r in project_rows(new_path,{"IID"}) if r["IID"] in old]
    result = []
    for p in pcs:
        x = np.asarray([float(a[p+"_AVG"]) for a,_ in pairs]); y = np.asarray([float(b[p+"_AVG"]) for _,b in pairs])
        centered = x-x.mean(); spread = float(centered@centered)
        slope = float(centered@(y-y.mean())/spread) if spread else None
        result.append(dict(pc=p,n=len(pairs),pearson_r=float(np.corrcoef(x,y)[0,1]) if spread and np.var(y)>0 else None,
                           slope_new_on_old=slope,intercept=float(y.mean()-slope*x.mean()) if slope is not None else None,
                           rms_difference=float(np.sqrt(np.mean((y-x)**2))),note="Descriptive geometry audit; no rescaling or correction applied"))
    project_write_table(out,result,["pc","n","pearson_r","slope_new_on_old","intercept","rms_difference","note"])
    return result


def project_run(args, *, prepared_plan=None):
    plan = project_prepare_plan(args) if prepared_plan is None else prepared_plan
    out = Path(plan["out"])
    executable = shutil.which(str(args.plink2))
    if not executable: raise ValueError(f"PLINK2 executable unavailable: {args.plink2}")
    version = subprocess.run([executable,"--version"],capture_output=True,text=True,check=True).stdout.strip()
    tool = dict(version=version,executable=project_state(executable))
    completion = out/"PROJECTION_COMPLETE.json"
    if completion.exists():
        saved = json.loads(completion.read_text())
        if saved["signature"] != plan["signature"] or saved["tool"] != tool: raise ValueError("Completed projection identity changed; use a new output directory")
        project_check_output_state(saved["outputs"])
        receipt_verify_projection_receipt(out/"reference.sscore", out/"reference.coordinate_receipt.json",
                                  [p+"_AVG" for p in plan["pcs"]], plan["basis_dir"])
        return json.loads((out/"projection_manifest.json").read_text())
    project_save_json(out/"projection_manifest.json",dict(status="running",signature=plan["signature"]))
    receipts = {s["label"]:project_execute_stage(s,plan,executable,tool) for s in plan["stages"]}
    refs = [s for s in plan["stages"] if not s["target"]]; targets = [s for s in plan["stages"] if s["target"]]
    summaries = {}
    for label,stages in (("reference",refs),("target",targets)):
        samples = project_sample_ids(stages[0]["files"]["psam"])
        summaries[label] = project_aggregate_scores([s["output"]+".sscore" for s in stages],samples,plan["pcs"],plan["n_shared_markers"],
            out/f"{label}.sscore",out/f"{label}_sample_qc.tsv",args.min_call_rate)
    geometry = project_geometry_change(Path(plan["basis_dir"])/"all.sscore",out/"reference.sscore",plan["pcs"],out/"reference_geometry_change.tsv")
    receipt_outputs = receipt_publish_projection_receipts(plan, receipts)
    coordinate_frame_id = json.loads((out/"coordinate_frame.json").read_text())["coordinate_frame_id"]
    manifest = dict(status="complete_frozen_basis_projection",signature=plan["signature"],basis_signature=plan["basis_signature"],
        coordinate_frame_id=coordinate_frame_id,
        coordinate_receipts={role:str(out/f"{role}.coordinate_receipt.json") for role in ("reference","target")},
        build="GRCh37",tool=tool,counts={k:plan[k] for k in ("n_original_markers","n_eligible_markers","n_shared_markers","n_reference_samples","n_target_samples","n_training_samples")},
        scores=summaries,coordinate_kind=plan["coordinate_kind"],marker_quality_counts=plan["marker_quality_counts"],
        dosage_storage={k:r["storage_evidence"] for k,r in receipts.items()},geometry_change=geometry,
        reference_coordinates=str(out/"reference.sscore"),target_coordinates=str(out/"target.sscore"),
        fit_cohort="frozen reference basis only; no UKB fit",imputation_pipeline_equivalence="not established",limitations=PROJECT_LIMITATIONS)
    project_save_json(out/"projection_manifest.json",manifest)
    outputs = [out/name for name in ("reference.sscore","target.sscore","reference_sample_qc.tsv","target_sample_qc.tsv","reference_geometry_change.tsv","projection_manifest.json")]
    outputs.extend(receipt_outputs)
    outputs.extend(Path(s["output"]+".receipt.json") for s in plan["stages"])
    states = {str(p):project_state(p) for p in outputs}
    for receipt in receipts.values(): states.update(receipt["outputs"])
    project_save_json(completion,dict(signature=plan["signature"],tool=tool,outputs=states))
    return manifest


def project_parser():
    p = argparse.ArgumentParser(description=_description_pca_project,allow_abbrev=False)
    p.add_argument("command",choices=("plan","run")); p.add_argument("--basis-dir",required=True,type=Path)
    p.add_argument("--target-glob",default=PROJECT_DEFAULT_GLOB); p.add_argument("--out",required=True,type=Path)
    p.add_argument("--reference-metadata",type=Path,help="Reference metadata; retain frozen historical provenance checks")
    p.add_argument("--target-metadata",type=Path,help="Optional explicit target metadata snapshot with sample_id")
    p.add_argument("--provenance-manifest",type=Path,help="Existing historical reference-metadata provenance receipt")
    p.add_argument("--target-info",type=Path,help="Exact target GRCh37 chrom,pos,ref,alt,info TSV; gzip/zstd accepted")
    p.add_argument("--min-info",type=float,default=.8); p.add_argument("--allow-unknown-info",action="store_true")
    p.add_argument("--allow-palindromic",action="store_true",help="Only when both inputs have independently verified common reference-strand labels")
    p.add_argument("--min-overlap",type=float,default=.9); p.add_argument("--min-call-rate",type=float,default=.9)
    p.add_argument("--plink2",default="plink2"); p.add_argument("--threads",type=int,default=4); p.add_argument("--memory-mb",type=int,default=8192)
    return p


def project_main(argv=None):
    args = project_parser().parse_args(argv); out = args.out.expanduser().resolve()
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
            plan = project_prepare_plan(args)
            try:
                result = plan if args.command == "plan" else project_run(args, prepared_plan=plan)
            except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
                if args.command == "run" and not (out/"PROJECTION_COMPLETE.json").exists():
                    project_save_json(out/"projection_manifest.json",dict(status="failed",signature=plan["signature"],reason=str(exc)))
                raise
        print(json.dumps({k:result[k] for k in ("status","signature")},indent=2)); return 0
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
        print(f"pca-project: {exc}",file=sys.stderr); return 2


def _entry_pca_project(): raise SystemExit(project_main())


# ========================================================================
# pca coordinate receipts
# ========================================================================
_description_pca_coordinate_receipts = 'Publish and verify common-panel PCA coordinate receipts.\n\nThis module reads metadata and existing scores only; it never fits PCA or\nchanges scores. Source-basis identity and scored-coordinate-frame identity\nremain separate, and reference metadata retains its existing provenance rule.'
"""Publish and verify common-panel PCA coordinate receipts.

This module reads metadata and existing scores only; it never fits PCA or
changes scores. Source-basis identity and scored-coordinate-frame identity
remain separate, and reference metadata retains its existing provenance rule.
"""
import io
from types import SimpleNamespace
from analysis_integrity import identity

RECEIPT_RECEIPT_SCHEMA = "allele3-pca-coordinate-receipt-1"
RECEIPT_FRAME_SCHEMA = "allele3-pca-coordinate-frame-1"


def receipt_digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def receipt_artifact(path):
    value = identity(path, full=True)
    if value.get("status") == "missing":
        raise ValueError(f"Missing coordinate evidence: {path}")
    return {key: value[key] for key in ("path", "size", "sha256")}


def receipt_verify_artifact(saved):
    current = receipt_artifact(saved["path"])
    if current != saved:
        raise ValueError("Coordinate evidence changed: " + saved["path"])
    return Path(saved["path"])


def receipt_metadata_request(args, basis_plan):
    """Bind inputs without promoting an unverified historical reference table.

    A target metadata file is a new explicit input snapshot.  A reference table
    remains subject to the pre-existing frozen-basis/provenance contract.
    Include this return value in the projection plan signature.
    """
    original = basis_plan.get("inputs", {}).get("metadata")
    explicit_reference = getattr(args, "reference_metadata", None)
    reference = explicit_reference or (original or {}).get("path")
    provenance = getattr(args, "provenance_manifest", None)
    if reference and Path(reference).is_file():
        history = frame_historical_identity(reference, provenance, original)
        verified = history["historical_metadata_verification"] == "verified"
        reference_binding = {
            "status": "verified_frozen_reference_metadata" if verified else "historically_unverified",
            "file": receipt_artifact(reference),
            # Only a prior verified digest can certify historical identity.
            "recorded": {"path": str(Path(reference).resolve()), "sha256": history["expected_sha256"]}
                        if verified else original,
            "historical_evidence": history,
        }
    elif explicit_reference:
        raise ValueError("Explicit --reference-metadata does not exist")
    else:
        reference_binding = {"status": "reference_metadata_not_available_here",
                             "file": None, "recorded": original}
    target = getattr(args, "target_metadata", None)
    target_snapshot = receipt_artifact(target) if target else None
    target_binding = ({"status": "explicit_new_target_metadata_snapshot", "file": target_snapshot,
                       "recorded": target_snapshot} if target else
                      {"status": "not_provided", "file": None, "recorded": None})
    return {"reference": reference_binding, "target": target_binding}


def receipt_by_id(path, required):
    result = {}
    for row in project_rows(path, required):
        if row["ID"] in result:
            raise ValueError("Duplicate prepared variant: " + row["ID"])
        result[row["ID"]] = row
    return result


def receipt_weights(path, pcs):
    result = {}
    for row in project_rows(path, {"ID", "A1", *pcs}):
        key = row["ID"], row["A1"]
        if key in result:
            raise ValueError("Duplicate scored allele")
        result[key] = [project_numeric(row[pc], pc) for pc in pcs]
    return result


def receipt_number_string(value):
    return format(float(value), ".17g")


def receipt_frame_spec(plan):
    """Reconstruct the actual common transform from existing prepared files.

    Target native IDs/orientation do not enter the frame identity.  Their
    equivalence to the reference transform is checked explicitly instead.
    """
    root, source_plan, _ = project_completed_outputs(plan["basis_dir"])
    if source_plan["signature"] != plan["basis_signature"]:
        raise ValueError("Projection source basis signature changed")
    if plan["build"] != "GRCh37" or plan["pcs"] != [f"PC{i+1}" for i in range(len(plan["pcs"]))]:
        raise ValueError("Unexpected coordinate build/PC order")
    for path, expected in plan["prepared_sha256"].items():
        if project_sha(path) != expected:
            raise ValueError("Prepared projection input changed: " + path)
    selected = [row for row in project_rows(Path(plan["out"]) / "marker_audit.tsv.gz")
                if row["status"].startswith("matched_")]
    selected.sort(key=lambda r: (int(r["chrom"]), int(r["pos"]), r["ref"], r["alt"], r["id"]))
    marker_ids = {r["id"] for r in selected}
    coordinates = {(r["chrom"], r["pos"]) for r in selected}
    if (not selected or len(marker_ids) != len(selected) or len(coordinates) != len(selected)
            or len(selected) != plan["n_shared_markers"]):
        raise ValueError("Common marker identities/count disagree")
    reference_stages = [s for s in plan["stages"] if not s["target"]]
    if len(reference_stages) != 1:
        raise ValueError("Expected one common-panel reference score stage")
    reference = reference_stages[0]
    ref_weights = receipt_weights(reference["prepared"]["weights"], plan["pcs"])
    ref_counts = receipt_by_id(reference["prepared"]["frequencies"], {"ID", "CHROM", "REF", "ALT", "ALT_CTS", "OBS_CT"})
    source_weights = receipt_weights(root / "ref.eigenvec.allele", plan["pcs"])
    source_counts = receipt_by_id(root / "ref.acount", {"ID", "CHROM", "REF", "ALT", "ALT_CTS", "OBS_CT"})
    if ref_weights != {key: value for key, value in source_weights.items() if key[0] in marker_ids}:
        raise ValueError("Prepared loadings differ from the frozen source axes")
    weights_by_id = {}
    for (vid, allele), value in ref_weights.items():
        weights_by_id.setdefault(vid, {})[allele] = value
    score_args = SimpleNamespace(**plan["settings"], components=len(plan["pcs"]))
    for stage in plan["stages"]:
        if stage["command"] != project_score_command(score_args, stage["files"], stage["prepared"], stage["output"]):
            raise ValueError("Executed score command does not implement the declared transform")
    ids = Path(reference["prepared"]["extract"]).read_text().splitlines()
    if set(ids) != marker_ids or len(ids) != len(selected) or set(ref_counts) != marker_ids:
        raise ValueError("Prepared reference panel differs from selected coordinates")
    weights, frequencies = [], []
    for marker in selected:
        vid, ref, alt = marker["id"], marker["ref"], marker["alt"]
        row = ref_counts[vid]
        if (row["CHROM"], row["REF"], row["ALT"]) != (marker["chrom"], ref, alt):
            raise ValueError("Prepared reference frequency identity differs")
        original_count = source_counts.get(vid)
        if original_count is None or any(row[key] != original_count[key] for key in ("CHROM", "REF", "ALT")):
            raise ValueError("Prepared reference contrast differs from frozen source")
        ac, an = project_numeric(row["ALT_CTS"], "AC"), project_numeric(row["OBS_CT"], "AN")
        if ac != project_numeric(original_count["ALT_CTS"], "AC") or an != project_numeric(original_count["OBS_CT"], "AN"):
            raise ValueError("Prepared reference counts differ from frozen source")
        if not 0 < ac < an:
            raise ValueError("Reference variance must be positive")
        alleles = sorted(weights_by_id.get(vid, {}))
        if not alleles or not set(alleles) <= {ref, alt}:
            raise ValueError("Missing/invalid common-panel allele weights")
        weights.extend([[vid, allele, *[receipt_number_string(v) for v in ref_weights[vid, allele]]]
                        for allele in alleles])
        frequencies.append([vid, ref, alt, receipt_number_string(ac), receipt_number_string(an)])
    if {v for v, _ in ref_weights} != marker_ids:
        raise ValueError("Unexpected reference weight variant")
    visited = set()
    for stage in [s for s in plan["stages"] if s["target"]]:
        match = re.fullmatch(r"target_(\d+)", stage["label"])
        if not match:
            raise ValueError("Unrecognized current projection target stage")
        index = int(match[1]) - 1
        subset = [r for r in selected if int(r["target_index"]) == index]
        native_ids = [r["target_id"] for r in subset]
        actual_ids = Path(stage["prepared"]["extract"]).read_text().splitlines()
        if (set(actual_ids) != set(native_ids) or len(actual_ids) != len(native_ids)
                or stage["n_markers"] != len(subset)):
            raise ValueError("Target prepared marker subset differs")
        tw = receipt_weights(stage["prepared"]["weights"], plan["pcs"])
        tc = receipt_by_id(stage["prepared"]["frequencies"], {"ID", "CHROM", "REF", "ALT", "ALT_CTS", "OBS_CT"})
        if set(tc) != set(native_ids):
            raise ValueError("Unexpected target frequency variants")
        expected_tw = {}
        for marker in subset:
            vid, native = marker["id"], marker["target_id"]
            if vid in visited:
                raise ValueError("Common marker scored in multiple target stages")
            visited.add(vid)
            if {marker["target_ref"], marker["target_alt"]} != {marker["ref"], marker["alt"]}:
                raise ValueError("Target allele contrast differs")
            reference_count = ref_counts[vid]
            ac, an = float(reference_count["ALT_CTS"]), float(reference_count["OBS_CT"])
            expected_ac = ac if marker["target_alt"] == marker["alt"] else an - ac
            row = tc[native]
            if ((row["CHROM"], row["REF"], row["ALT"]) !=
                    (marker["chrom"], marker["target_ref"], marker["target_alt"])
                    or float(row["ALT_CTS"]) != expected_ac or float(row["OBS_CT"]) != an):
                raise ValueError("Target does not use the same reference frequencies")
            expected_tw.update({(native, allele): value for allele, value in weights_by_id[vid].items()})
        if tw != expected_tw:
            raise ValueError("Target does not use the same allele loadings")
    if visited != marker_ids:
        raise ValueError("Target stages do not partition the common panel")
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(["chrom", "pos", "id", "ref", "alt"])
    writer.writerows([[r[k] for k in ("chrom", "pos", "id", "ref", "alt")] for r in selected])
    panel_text = buffer.getvalue()
    source_artifacts = {name: receipt_artifact(root / name) for name in
                        ("pca_plan.json", "PCA_COMPLETE", "completed.sha256",
                         "ref.eigenvec.allele", "ref.acount", "ref.eigenval")}
    transform = {
        "schema": RECEIPT_FRAME_SCHEMA, "build": plan["build"], "source_basis_id": plan["basis_signature"],
        "source_axes_sha256": source_artifacts["ref.eigenvec.allele"]["sha256"],
        "source_frequencies_sha256": source_artifacts["ref.acount"]["sha256"],
        "source_eigenvalues_sha256": source_artifacts["ref.eigenval"]["sha256"],
        "pcs": plan["pcs"], "n_markers": len(selected),
        "marker_sha256": hashlib.sha256(panel_text.encode()).hexdigest(),
        "used_loadings_sha256": receipt_digest_json(weights), "used_reference_frequencies_sha256": receipt_digest_json(frequencies),
        "scoring": {"centering": "reference_frequency", "variance_standardize": True,
                    "missing": "reference_mean_imputation", "coordinate_kind": "projected_average",
                    "denominator": 2 * len(selected), "chromosomes": "autosomes_1_22",
                    "aggregation": "sum_scores_and_denominators_before_division"},
    }
    return {"schema": RECEIPT_FRAME_SCHEMA, "coordinate_frame_id": "sha256:" + receipt_digest_json(transform),
            "transform": transform, "source_basis_dir": str(root), "source_artifacts": source_artifacts}, panel_text


def receipt_score_identity(path, pcs, expected_samples, n_markers):
    required = {"FID", "IID", "DENOM", *[pc + "_AVG" for pc in pcs]}
    samples = []
    for row in project_rows(path, required):
        samples.append((row["FID"], row["IID"]))
        if project_numeric(row["DENOM"], "DENOM") != 2 * n_markers:
            raise ValueError("Coordinate denominator differs from common panel")
        for pc in pcs:
            project_numeric(row[pc + "_AVG"], pc)
    if samples != expected_samples or len({sid for _, sid in samples}) != len(samples):
        raise ValueError("Coordinate samples differ from original PSAM")
    return {"file": receipt_artifact(path), "sample_ids_sha256": receipt_digest_json(samples), "n_samples": len(samples)}


def receipt_check_metadata(binding, sample_ids):
    if not binding.get("file"):
        return
    path = receipt_verify_artifact(binding["file"])
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "sample_id" not in (reader.fieldnames or []):
            raise ValueError("Coordinate metadata requires sample_id")
        ids = [r["sample_id"] for r in reader]
    if len(ids) != len(set(ids)) or not set(sample_ids) <= set(ids):
        raise ValueError("Coordinate metadata has duplicates or lacks scored samples")


def receipt_publish_projection_receipts(plan, stage_receipts):
    """Call after score aggregation, before publishing PROJECTION_COMPLETE.

    Return every newly generated path; add them to completion.outputs.  A
    consumer rejects these receipts until that final completion exists.
    """
    out = Path(plan["out"])
    frame, panel_text = receipt_frame_spec(plan)
    for stage in plan["stages"]:
        saved = stage_receipts[stage["label"]]
        if saved["status"] != "complete" or saved["identity"]["signature"] != plan["signature"]:
            raise ValueError("Scoring stage not verified for this projection")
        project_check_output_state(saved["outputs"])
        expected = Path(stage["prepared"]["extract"]).read_text().splitlines()
        used = Path(stage["output"] + ".sscore.vars").read_text().splitlines()
        if set(used) != set(expected) or len(used) != len(expected):
            raise ValueError("Actually scored stage markers differ")
    panel_path, frame_path = out / "common_panel.tsv", out / "coordinate_frame.json"
    project_atomic_text(panel_path, panel_text)
    frame["common_panel"] = receipt_artifact(panel_path)
    frame["projection_plan"] = receipt_artifact(out / "projection_plan.json")
    project_save_json(frame_path, frame)
    scores, samples = {}, {}
    for role in ("reference", "target"):
        stage = next(s for s in plan["stages"] if s["target"] == (role == "target"))
        samples[role] = project_sample_ids(stage["files"]["psam"])
        scores[role] = receipt_score_identity(out / f"{role}.sscore", plan["pcs"], samples[role], plan["n_shared_markers"])
    result = [panel_path, frame_path]
    bindings = plan["coordinate_metadata"]
    for role in ("reference", "target"):
        receipt_check_metadata(bindings[role], [sid for _, sid in samples[role]])
        proof = {
            "schema": RECEIPT_RECEIPT_SCHEMA, "status": "verified", "coordinate_kind": "projected_average",
            "role": role, "coordinates": scores[role]["file"], "coordinates_sha256": scores[role]["file"]["sha256"],
            "paired_scores": scores, "source_basis_id": plan["basis_signature"],
            "coordinate_frame_id": frame["coordinate_frame_id"], "basis_id": frame["coordinate_frame_id"],
            "marker_sha256": frame["transform"]["marker_sha256"], "frame": receipt_artifact(frame_path),
            "metadata": bindings[role], "projection_signature": plan["signature"],
            "completion_path": str(out / "PROJECTION_COMPLETE.json"),
        }
        path = out / f"{role}.coordinate_receipt.json"
        project_save_json(path, proof)
        result.append(path)
    return result


def receipt_verify_projection_receipt(coordinates, receipt_path, pcs, basis_dir=None):
    """Validate a producer receipt, its complete run, both scores and transform."""
    proof = json.loads(Path(receipt_path).read_text())
    if (proof.get("schema") != RECEIPT_RECEIPT_SCHEMA or proof.get("status") != "verified"
            or proof.get("coordinate_kind") != "projected_average" or proof.get("role") not in {"reference", "target"}):
        raise ValueError("Invalid projection coordinate receipt")
    if set(proof.get("paired_scores", {})) != {"reference", "target"}:
        raise ValueError("Coordinate receipt must bind both cohort score files")
    if Path(coordinates).resolve() != Path(proof["coordinates"]["path"]).resolve():
        raise ValueError("Coordinate receipt belongs to a different score file/role")
    completion = json.loads(Path(proof["completion_path"]).read_text())
    if completion["signature"] != proof["projection_signature"]:
        raise ValueError("Coordinate completion signature mismatch")
    project_check_output_state(completion["outputs"])
    if str(Path(receipt_path).resolve()) not in completion["outputs"]:
        raise ValueError("Coordinate receipt is not bound by the completed projection")
    frame_path = receipt_verify_artifact(proof["frame"])
    frame = json.loads(frame_path.read_text())
    plan_path = receipt_verify_artifact(frame["projection_plan"])
    plan = json.loads(plan_path.read_text())
    if plan["signature"] != proof["projection_signature"]:
        raise ValueError("Coordinate projection plan mismatch")
    rebuilt, panel_text = receipt_frame_spec(plan)
    if (rebuilt["transform"] != frame["transform"] or rebuilt["coordinate_frame_id"] != frame["coordinate_frame_id"]
            or rebuilt["source_artifacts"] != frame["source_artifacts"]
            or rebuilt["source_basis_dir"] != frame["source_basis_dir"]):
        raise ValueError("Recorded coordinate frame differs from actual shared transform")
    for saved in frame["source_artifacts"].values():
        receipt_verify_artifact(saved)
    panel_path = receipt_verify_artifact(frame["common_panel"])
    if panel_path.read_text() != panel_text:
        raise ValueError("Recorded common panel differs from actual scored markers")
    if (proof["coordinate_frame_id"] != frame["coordinate_frame_id"]
            or proof["basis_id"] != frame["coordinate_frame_id"]
            or proof["source_basis_id"] != frame["transform"]["source_basis_id"]
            or proof["marker_sha256"] != frame["transform"]["marker_sha256"]
            or frame["transform"]["pcs"][:len(pcs)] != [p.removesuffix("_AVG") for p in pcs]):
        raise ValueError("Coordinate frame aliases or PC order disagree")
    if basis_dir and Path(basis_dir).resolve() != Path(frame["source_basis_dir"]).resolve():
        raise ValueError("Explicit basis differs from coordinate receipt source")
    for role, score in proof["paired_scores"].items():
        path = receipt_verify_artifact(score["file"])
        stage = next(s for s in plan["stages"] if s["target"] == (role == "target"))
        samples = project_sample_ids(stage["files"]["psam"])
        if receipt_score_identity(path, plan["pcs"], samples, plan["n_shared_markers"]) != score:
            raise ValueError("Coordinate sample identity differs")
    if proof["coordinates"] != proof["paired_scores"][proof["role"]]["file"]:
        raise ValueError("Role-specific score identity differs")
    if proof["coordinates_sha256"] != proof["coordinates"]["sha256"]:
        raise ValueError("Role-specific score SHA alias differs")
    if proof["metadata"] != plan["coordinate_metadata"][proof["role"]]:
        raise ValueError("Coordinate metadata differs from planned binding")
    score_ids = [row["IID"] for row in project_rows(coordinates, {"IID"})]
    receipt_check_metadata(proof["metadata"], score_ids)
    audit = {"status": "verified", "basis_id": proof["coordinate_frame_id"],
             "coordinate_frame_id": proof["coordinate_frame_id"], "source_basis_id": proof["source_basis_id"],
             "source_basis_dir": frame["source_basis_dir"], "coordinates": identity(coordinates, full=True),
             "supplied_receipt": identity(receipt_path, full=True), "marker_sha256": proof["marker_sha256"],
             "metadata_recorded": proof["metadata"].get("recorded"),
             "metadata_binding_status": proof["metadata"]["status"], "projection_role": proof["role"]}
    if proof["role"] == "reference":
        # Original reference sample training membership remains meaningful.
        audit["basis_dir"] = frame["source_basis_dir"]
    return audit


def receipt_require_same_coordinate_frame(*contracts):
    """Use before comparing/concatenating independently verified coordinates."""
    if not contracts or any(c.get("status") != "verified" or not c.get("coordinate_frame_id") for c in contracts):
        raise ValueError("Comparison requires verified coordinate frame identities")
    if len({c["coordinate_frame_id"] for c in contracts}) != 1:
        raise ValueError("Different marker/weight/frequency/scale frames; re-score on one common panel")
    return contracts[0]["coordinate_frame_id"]


# ========================================================================
# pca contract
# ========================================================================
_description_pca_contract = 'Strict shared-coordinate contract; historical hashes are evidence, not guesses.'
"""Strict shared-coordinate contract; historical hashes are evidence, not guesses."""
frame_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(frame_paths_root / "0f"))
from analysis_integrity import sha256

FRAME_DEFAULT_PROVENANCE = configured_path('ANALYSIS_ROOT', 'allele3_migration_20261028/all_files_verified.json')


def frame_historical_identity(path, receipt=None, recorded=None):
    actual = identity(path, full=True)
    expected = (recorded or {}).get('sha256')
    provenance = None
    receipt = receipt or os.environ.get('ALLELE_PCA_PROVENANCE', FRAME_DEFAULT_PROVENANCE)
    if not expected and receipt and Path(receipt).is_file():
        evidence = json.loads(Path(receipt).read_text())
        if evidence.get('status') != 'verified':
            raise ValueError('Historical provenance receipt is not verified')
        matches = [r for r in evidence.get('records', [])
                   if str(Path(r.get('destination', '')).resolve()) == actual['path']
                   or (recorded and str(Path(evidence.get('source', '')) / r.get('source_relative', '')) == recorded.get('path'))]
        if len(matches) > 1: raise ValueError('Ambiguous historical provenance entry')
        if matches:
            expected = matches[0].get('sha256'); provenance = identity(receipt, full=True)
    # The original metadata adapter saved full output hashes alongside the
    # metadata. This receipt survives removal of D outputs / old audit folders.
    adapter = Path(path).resolve().parent / 'pca3d_metadata.json'
    if not expected and adapter.is_file():
        saved_outputs = json.loads(adapter.read_text()).get('outputs', {})
        matches = [(name, digest) for name, digest in saved_outputs.items()
                   if name == actual['path'] or (recorded and name == recorded.get('path'))
                   or Path(name).name == Path(path).name]
        if len(matches) != 1:
            raise ValueError('Ambiguous/missing original metadata-adapter output receipt')
        expected = matches[0][1]
        if not isinstance(expected, str) or not re.fullmatch(r'[0-9a-fA-F]{64}', expected):
            raise ValueError('Invalid original metadata-adapter digest')
        provenance = identity(adapter)
    if expected and actual['sha256'] != expected:
        raise ValueError(f'Historical content changed: {path}')
    return {'input': actual, 'historical_metadata_verification': 'verified' if expected else 'unverified',
            'expected_sha256': expected, 'provenance_receipt': provenance}


def frame_read_scores(path, count, *, basis_dir=None, receipt=None, require_verified=True):
    from analysis_integrity import resolve_result_input
    path = resolve_result_input(path)
    basis_dir = resolve_result_input(basis_dir)
    receipt = resolve_result_input(receipt)
    if count < 0: raise ValueError('PC count must be nonnegative')
    if count == 0:
        return {}, [], {'expected_pc_count': 0, 'actual_pc_count': 0, 'coordinate_kind': 'disabled_by_pcs_0',
                        'basis_id': None, 'status': 'explicitly_unadjusted'}
    if not path: raise ValueError(f'{count} PCs requested but no coordinates supplied')
    path = Path(path).resolve()
    opener = gzip.open if path.suffix == '.gz' else open
    scores, records, row_bases = {}, {}, set()
    with opener(path, 'rt') as handle:
        lines = (line for line in handle if line.strip() and not line.startswith('##'))
        try: names = [x.lstrip('#') for x in next(lines).split()]
        except StopIteration: raise ValueError('Empty PCA coordinates')
        if len(set(names)) != len(names) or 'IID' not in names:
            raise ValueError('Coordinates need unique headers and IID')
        if 'PC1_AVG' not in names and any(re.fullmatch(r'PC\d+_SUM', x) for x in names):
            raise ValueError('PC SUM scores are not accepted; supply verified projected AVG on a common marker panel')
        if any(re.fullmatch(r'PC\d+', x) for x in names):
            raise ValueError('Raw eigenvectors are not projected AVG coordinates')
        columns = [f'PC{i}_AVG' for i in range(1, count + 1)]
        missing = set(columns) - set(names)
        if missing: raise ValueError('Missing consecutive requested PCs: ' + ','.join(sorted(missing)))
        for line in lines:
            cells = line.split()
            if len(cells) != len(names): raise ValueError('Malformed PCA coordinate row')
            row = dict(zip(names, cells)); iid = row['IID']
            if iid in {'', '.', 'NA'} or iid in scores: raise ValueError('Duplicate/missing PCA IID: ' + iid)
            values = [float(row[k]) for k in columns]
            if not all(math.isfinite(v) for v in values): raise ValueError('Nonfinite PC coordinate: ' + iid)
            if row.get('coordinate_kind', 'projected_average') not in {'projected_average', 'projected_AVG'}:
                raise ValueError('Mixed coordinate kinds')
            if row.get('basis_id'): row_bases.add(row['basis_id'])
            scores[iid] = values; records[iid] = row
    if not scores: raise ValueError('No coordinate rows')
    if len(row_bases) > 1: raise ValueError('Mixed PCA basis IDs')
    audit = {'expected_pc_count': count, 'actual_pc_count': len(columns), 'coordinate_kind': 'projected_average',
             'coordinates': identity(path), 'basis_id': None, 'status': 'unverified', 'n_coordinates': len(scores)}
    # Automatic discovery is restricted to a completed receipt in the coordinate's
    # own directory, never a guess based on column names or another default basis.
    # Producer receipts may only be auto-discovered beside their exact score file.
    if not receipt:
        candidate = path.with_name(path.stem + ".coordinate_receipt.json")
        if candidate.is_file():
            receipt = candidate
    basis = Path(basis_dir).resolve() if basis_dir else path.parent
    plan_path = basis / 'pca_plan.json'
    if plan_path.is_file() and path == basis / 'all.sscore':
        verified = project_validate_completed_coordinates(basis)
        plan = json.loads(plan_path.read_text())
        audit.update(status='verified', basis_id=verified['signature'], basis_dir=str(basis),
                     source_basis_id=verified['signature'],
                     coordinate_frame_id='original-frozen-coordinates:'+verified['signature'],
                     basis_receipts={p: identity(basis / p, full=True) for p in ('pca_plan.json', 'completed.sha256', 'PCA_COMPLETE')},
                     metadata_recorded=plan.get('inputs', {}).get('metadata'))
    elif receipt:
        proof = json.loads(Path(receipt).read_text())
        if proof.get('schema') == 'allele3-pca-coordinate-receipt-1':
            audit.update(receipt_verify_projection_receipt(path, receipt, columns, basis_dir))
        elif proof.get('schema'):
            raise ValueError('Unsupported coordinate receipt schema')
        else:
            # Keep the existing explicitly supplied legacy receipt contract.
            audit['coordinates'] = identity(path, full=True)
            if (proof.get('status') != 'verified' or proof.get('coordinate_kind') != 'projected_average'
                    or proof.get('coordinates_sha256') != audit['coordinates']['sha256'] or not proof.get('basis_id')
                    or not proof.get('marker_sha256')):
                raise ValueError('Coordinate receipt does not verify AVG, basis, marker set and full score SHA')
            audit.update(status='verified', basis_id=proof['basis_id'], supplied_receipt=identity(receipt, full=True),
                         metadata_recorded={'sha256': proof.get('metadata_sha256')} if proof.get('metadata_sha256') else None)
    elif basis_dir:
        raise ValueError('Explicit basis does not bind these coordinates; supply their projection receipt')
    if row_bases and (audit['basis_id'] not in row_bases): raise ValueError('Row basis disagrees with verified receipt')
    if require_verified and audit['status'] != 'verified':
        raise ValueError('No verified shared-basis coordinate receipt; adjusted model unavailable')
    return scores, columns, audit


def frame_coordinate_join(samples, scores):
    ids = [s['sample_id'] for s in samples]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate metadata sample_id')
    return {'n_metadata': len(ids), 'n_coordinates': len(scores), 'n_joined': len(set(ids) & scores.keys()),
            'metadata_without_coordinates': sorted(set(ids) - scores.keys()),
            'coordinates_without_metadata': sorted(scores.keys() - set(ids))}


# ========================================================================
# pca basis
# ========================================================================
_description_pca_basis = 'Create a separate projection basis using unchanged axes and unique coordinates.'
#!/usr/bin/env python3
"""Create a separate projection basis using unchanged axes and unique coordinates."""


def basis_subset_table(source, destination, keep):
    """Keep exact original header/row bytes, including frequency orientation."""
    with source.open('rb') as src, destination.open('wb') as dest:
        header = src.readline()
        index = header.decode().lstrip('#').split().index('ID')
        dest.write(header)
        for line in src:
            if line.decode().split()[index] in keep:
                dest.write(line)


def basis_prepare(args):
    source, out = args.basis_dir.resolve(), args.out.resolve()
    if source == out or source in out.parents or out in source.parents:
        raise ValueError('New basis must be separate from the original basis')
    if args.threads < 1 or args.memory_mb < 640:
        raise ValueError('Positive threads and at least 640 MiB memory required')
    # Only this derivation step may read duplicate coordinates. IDs, alleles,
    # frequencies, loadings, receipts and sample order remain strictly checked.
    original = project_basis(source, allow_duplicate_coordinates=True)
    executable = shutil.which(str(args.plink2))
    if not executable:
        raise ValueError('PLINK2 executable unavailable')
    tool = dict(executable=project_state(executable), version=subprocess.check_output([executable,'--version'],text=True).strip())
    coordinate_counts = Counter((m['chrom'],m['pos']) for m in original['markers'].values())
    excluded = [m for m in original['markers'].values() if coordinate_counts[m['chrom'],m['pos']] > 1]
    keep = [m['id'] for m in original['markers'].values() if coordinate_counts[m['chrom'],m['pos']] == 1]
    if not keep:
        raise ValueError('No unique-coordinate reference markers remain')
    identity = dict(source=original['inputs'], tool=tool, policy='exclude_all_records_at_repeated_coordinates',
                    version=1, code_sha256=project_sha(__file__))
    signature = hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    out.mkdir(parents=True,exist_ok=True)
    with (out/'.basis.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        receipt = out/'derivation.json'
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved['signature'] != signature:
                raise ValueError('Basis derivation inputs/code/tool changed; use a new directory')
            if saved['status'] == 'complete':
                for name,digest in saved['outputs_sha256'].items():
                    if project_sha(out/name) != digest:
                        raise ValueError('Derived basis file changed: '+name)
                project_validate_basis(out)
                return saved
        elif any(p.name != '.basis.lock' for p in out.iterdir()):
            raise ValueError('New basis requires an empty output directory')
        record = dict(status='preparing',signature=signature,identity=identity,
                      original_basis=str(source),original_signature=original['plan']['signature'],
                      n_original_markers=len(original['markers']),n_retained_markers=len(keep),
                      n_excluded_markers=len(excluded),n_reference_samples=len(original['samples']),
                      axes_refitted=False,original_basis_modified=False)
        project_save_json(receipt,record)
        project_atomic_text(out/'retained_markers.ids',''.join(x+'\n' for x in keep))
        project_write_table(out/'excluded_markers.tsv',excluded,['chrom','pos','id','ref','alt'])
        for name in ('ref.acount','ref.eigenvec.allele'):
            basis_subset_table(source/name,out/name,set(keep))
        shutil.copy2(source/'ref.eigenval',out/'ref.eigenval')
        files=original['files']
        common=[executable,'--threads',str(args.threads),'--memory',str(args.memory_mb)]
        commands=[common+['--pgen',files['pgen'],'--pvar',files['pvar'],'--psam',files['psam'],
                  '--extract',str(out/'retained_markers.ids'),'--make-pgen','--out',str(out/'all_pruned')],
                  common+['--pfile',str(out/'all_pruned'),'--read-freq',str(out/'ref.acount'),
                  '--score',str(out/'ref.eigenvec.allele'),'2','5','header-read','variance-standardize','list-variants',
                  'cols=fid,nallele,denom,dosagesum,scoreavgs,scoresums','--score-col-nums',f"6-{5+len(original['pcs'])}",
                  '--out',str(out/'all')]]
        for label,command in zip(('extract','score'),commands):
            print('[projection basis] '+label,flush=True)
            with (out/(label+'.console.log')).open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        used=(out/'all.sscore.vars').read_text().splitlines()
        if set(used)!=set(keep) or len(used)!=len(keep):
            raise ValueError('PLINK did not score exactly the retained reference panel')
        for path,saved in original['inputs'].items():
            if project_state(saved['path'])!=saved:
                raise ValueError('Original basis changed during derivation')
        plan=dict(schema_version='derived-frozen-projection-basis-1',signature=signature,build='GRCh37',
                  inputs=original['plan'].get('inputs',{}),settings=original['plan'].get('settings',{}),
                  counts={**original['plan'].get('counts',{}),'retained_projection_markers':len(keep)},
                  coordinate_kind='unchanged_reference_axes_rescored_on_unique_coordinate_subset',
                  source_basis=str(source),source_signature=original['plan']['signature'],axes_refitted=False)
        project_save_json(out/'pca_plan.json',plan)
        expected=('all.sscore','ref.acount','ref.eigenvec.allele','ref.eigenval')
        project_atomic_text(out/'completed.sha256',''.join(f'{project_sha(out/name)}  {name}\n' for name in expected))
        project_atomic_text(out/'PCA_COMPLETE',signature+'\n')
        validation=project_validate_basis(out)
        geometry=project_geometry_change(source/'all.sscore',out/'all.sscore',original['pcs'],out/'reference_geometry_change.tsv')
        record.update(status='complete',commands=commands,validation=validation,geometry_change=geometry,
                      outputs_sha256={p.name:project_sha(p) for p in out.iterdir() if p.is_file() and p.name not in {'.basis.lock','derivation.json'}})
        project_save_json(receipt,record)
        return record


def basis_main(argv=None):
    p=argparse.ArgumentParser(description=_description_pca_basis)
    p.add_argument('--basis-dir',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--plink2',default='plink2');p.add_argument('--threads',type=int,default=4)
    p.add_argument('--memory-mb',type=int,default=8192)
    args=p.parse_args(argv)
    try:
        result=basis_prepare(args)
        print(json.dumps({k:result[k] for k in ('status','n_original_markers','n_retained_markers','n_excluded_markers','axes_refitted')},indent=2))
        return 0
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
        print('pca-basis: '+str(exc),file=sys.stderr)
        return 2


def _entry_pca_basis():raise SystemExit(basis_main())


# ========================================================================
# pca centroid
# ========================================================================
_description_pca_centroid = 'Shared dosage export and descriptive carrier centroids for PCA entry points.\n\nThis module preserves the legacy centroid interface. For joint nucleotide\noccupancy, explicit ploidy and modern ranks, use allele_space.py instead.'
"""Shared dosage export and descriptive carrier centroids for PCA entry points.

This module preserves the legacy centroid interface. For joint nucleotide
occupancy, explicit ploidy and modern ranks, use allele_space.py instead.
"""

from evolve import _coordinates, iter_pgen_variants

def centroid_op(path):return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)


def centroid_export_dosages(pgen,psam,selected,out,count_kind='hardcall'):
    """Stream original PVAR-indexed alleles, including explicit unavailable rows."""
    ids,variants=iter_pgen_variants(pgen,psam,selected,count_kind=count_kind)
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);n=0
    with gzip.open(out,'wt',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','variant_id','dosage','status','count_kind'])
        try:
            for r in variants:
                if r['dosages'] is None:
                    w.writerow(['.',r['variant_id'],'NA',r['status'],count_kind]);continue
                for iid,d in zip(ids,r['dosages']):
                    w.writerow([iid,r['variant_id'],'NA' if d is None else d,'missing_call' if d is None else 'ok',count_kind]);n+=1
        finally:variants.close()
    return dict(output=str(out),rows=n,count_kind=count_kind)


def centroid_run_centroids(a):
    if a.bin_years<=0:raise ValueError('Positive --bin-years required')
    coords,coordinate_kind=_coordinates(a.coordinates);meta={}
    with open(a.metadata) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            iid=r.get('sample_id') or r.get('IID')
            if iid in meta:raise ValueError(f'Duplicate metadata sample: {iid}')
            if not iid or r.get('include','1')!='1':continue
            cohort=r.get('cohort')
            if cohort not in {'ancient','modern'}:
                raise ValueError(f'{iid}: explicit ancient/modern cohort is required; age zero does not identify modern samples')
            try:age=float(r.get('age_bp') or r.get('age'))
            except (TypeError,ValueError):continue
            if not math.isfinite(age) or age<0:continue
            meta[iid]=(cohort,age)
    agg=defaultdict(lambda:[0,0.,0.,0.]);unavailable=defaultdict(set);measured=set();missing=Counter()
    current=None;seen_samples=set();finished=set()
    with centroid_op(a.dosages) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            iid=r.get('sample_id') or r.get('IID');vid=r.get('variant_id') or r.get('site_id')
            if not vid:raise ValueError('Missing dosage variant ID')
            if current!=vid:
                if vid in finished:raise ValueError('Dosage table must group rows by variant_id')
                if current is not None:finished.add(current)
                current=vid;seen_samples=set()
            if iid in seen_samples:raise ValueError('Duplicate dosage sample/variant')
            seen_samples.add(iid)
            if r.get('status','ok') not in {'ok','assayed_pass','measured'}:
                unavailable[vid].add(r.get('status','unavailable'));continue
            if iid not in coords or iid not in meta:missing[vid]+=1;continue
            try:d=float(r['dosage'])
            except (KeyError,ValueError):raise ValueError('Measured dosage must be numeric')
            if not math.isfinite(d) or not 0<=d<=2:raise ValueError('Invalid diploid dosage')
            cohort,age=meta[iid];k=-1 if cohort=='modern' else int(age//a.bin_years)
            p1,p2=coords[iid];z=agg[vid,cohort,k];z[0]+=1;z[1]+=d;z[2]+=d*p1;z[3]+=d*p2;measured.add(vid)
    rows=[]
    for (vid,cohort,k),(n,ac,p1,p2) in sorted(agg.items()):
        rows.append(dict(variant_id=vid,cohort=cohort,start_bp=0 if k<0 else k*a.bin_years,end_bp=0 if k<0 else (k+1)*a.bin_years,n=n,allele_count=ac,af=ac/(2*n),pc1=p1/ac if ac>0 else 'NA',pc2=p2/ac if ac>0 else 'NA',status='ok' if ac>0 else 'observed_zero'))
    for vid in sorted((set(unavailable)|set(missing))-measured):
        rows.append(dict(variant_id=vid,cohort='NA',start_bp='NA',end_bp='NA',n=0,allele_count='NA',af='NA',pc1='NA',pc2='NA',status='not_projectable:'+','.join(sorted(unavailable[vid] or {'missing_coordinates_or_age'}))))
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    fields=['variant_id','cohort','start_bp','end_bp','n','allele_count','af','pc1','pc2','status']
    with open(out/'allele_centroids.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
    valid=[r for r in rows if r['status']=='ok']
    if valid:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig=plt.figure(figsize=(9,7));ax=fig.add_subplot(projection='3d')
        for vid in sorted({r['variant_id'] for r in valid})[:12]:
            selected=sorted([r for r in valid if r['variant_id']==vid and r['cohort']=='ancient'],key=lambda r:r['start_bp'])
            if selected:
                ax.plot([r['pc1'] for r in selected],[r['pc2'] for r in selected],[(r['start_bp']+r['end_bp'])/2 for r in selected],marker='o',label=vid+' ancient')
            modern=[r for r in valid if r['variant_id']==vid and r['cohort']=='modern']
            if modern:
                ax.scatter([r['pc1'] for r in modern],[r['pc2'] for r in modern],[0.]*len(modern),marker='^',label=vid+' modern endpoint')
        ax.set(xlabel='PC1 carrier centroid',ylabel='PC2 carrier centroid',zlabel='Age BP');ax.legend(fontsize=6);fig.tight_layout()
        fig.savefig(out/'allele_centroid_trajectories.png',dpi=180)
        plt.close(fig)
    result=dict(status='complete',rows=len(rows),unavailable_variants=len((set(unavailable)|set(missing))-measured),coordinate_kind=coordinate_kind,bin_years=a.bin_years,rule='No measured individual dosage -> no PCA projection. Modern and ancient cohorts remain separate. Missing ages remain unavailable.')
    (out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result




# ========================================================================
# 0.pca3D
# ========================================================================
_description_0_pca3D = '0.pca3D'
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))

# The shared utility filename is sorted first; use a stable import name in Python.
import importlib.util as _common_import
import sys as _common_sys
from pathlib import Path as _CommonPath
if "allele_common" not in _common_sys.modules:
    _common_spec = _common_import.spec_from_file_location(
        "allele_common", _CommonPath(__file__).with_name("0.common.py"))
    _common_module = _common_import.module_from_spec(_common_spec)
    _common_sys.modules[_common_spec.name] = _common_module
    try:
        _common_spec.loader.exec_module(_common_module)
    except BaseException:
        _common_sys.modules.pop(_common_spec.name, None)
        raise


def main_centroid(argv=None):
    p=argparse.ArgumentParser(description=_description_0_pca3D);p.add_argument('--coordinates',required=True);p.add_argument('--metadata',required=True)
    p.add_argument('--dosages');p.add_argument('--pgen-prefix');p.add_argument('--selected-variants');p.add_argument('--count-kind',choices=['hardcall','dosage'],default='hardcall');p.add_argument('--out',required=True);p.add_argument('--bin-years',type=int,default=2000)
    a=p.parse_args(argv)
    if not a.dosages:
        if not a.pgen_prefix or not a.selected_variants:p.error('Supply --dosages or both --pgen-prefix and --selected-variants')
        a.dosages=str(Path(a.out)/'selected_variant_dosages.tsv.gz');centroid_export_dosages(a.pgen_prefix+'.pgen',a.pgen_prefix+'.psam',a.selected_variants,a.dosages,a.count_kind)
    print(json.dumps(centroid_run_centroids(a),indent=2))
from contextlib import nullcontext

from evolve import _atomic_text, _fingerprint, _pca_rows as _rows, _evolution_sha as _sha, prepare_pca_plan, render_pca_time

PREFIX = "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes"


def parser():
    p = argparse.ArgumentParser(description=_description_0_pca3D, allow_abbrev=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("command", choices=("run", "plan", "render", "atlas", "prepare-basis", "project-plan", "project"), nargs="?", default="run")
    p.add_argument("--out", "--dir-out", type=Path, default=Path(configured_path('ANALYSIS_ROOT', 'allele3/pca3d')))
    p.add_argument("--replace", choices=("TRUE", "FALSE"), default=os.environ.get("ALLELE_REPLACE", "TRUE"), help="Overwrite by default; FALSE reuses complete output")
    p.add_argument("--ancient-dir", "--dir-ancient", type=Path, default=Path(configured_path('DIR_GEN', 'ancient')))
    p.add_argument("--pgen-prefix", type=Path)
    p.add_argument("--metadata", type=Path, help="Canonical TSV; required for custom panels / render")
    p.add_argument("--hq-variants", type=Path, help="GRCh37 background variant TSV")
    p.add_argument("--supplement-xlsx", type=Path, help="Original ancient dataset metadata adapter")
    p.add_argument("--qpadm", type=Path)
    p.add_argument("--ancestry-model", choices=("3way", "4way"))
    p.add_argument("--dir-pgen-cache", type=Path)
    p.add_argument("--eigenvec", type=Path, help="Existing common-basis .sscore / .eigenvec for render")
    p.add_argument("--html", type=Path)
    p.add_argument("--source-html", type=Path, help="Existing complete allele3-data HTML for lossless atlas re-rendering")
    p.add_argument("--sample-features", type=Path, help="Optional per-person features joined by exact sample_id; requires --feature-schema")
    p.add_argument("--feature-schema", type=Path, help="TSV: feature_id,label,kind,units,provenance")
    p.add_argument("--atlas-theme", choices=("dark", "light"), default="dark")
    p.add_argument("--atlas-no-static", action="store_true", help="Offline interactive explorer and audits only; skip extra atlas static figures")
    p.add_argument("--selected-variants", type=Path, help="Optional indexed PVAR targets for allele overlays")
    p.add_argument("--target-variants", type=Path, help="Optional separate GRCh37 background exclusions")
    p.add_argument("--pca-basis", choices=("joint_balanced", "modern_projection", "ancient_balanced"), default="joint_balanced")
    p.add_argument("--pca-basis-label")
    p.add_argument("--pca-basis-min-age-bp", type=float)
    p.add_argument("--pca-basis-max-age-bp", type=float, help="Ancient training limit; ancient_balanced defaults to 10,000 BP")
    p.add_argument("--pca-maf", type=float, default=.05)
    p.add_argument("--pca-geno", type=float, default=.1)
    p.add_argument("--pca-components", type=int, default=10)
    p.add_argument("--pca-max-per-group", type=int, default=200)
    p.add_argument("--pca-time-bin-years", type=int, default=2000)
    p.add_argument("--pca-layer-years", type=int, default=1000, help="Display-only layer width; independent of PCA training time bins")
    p.add_argument("--pca-time-view", choices=("layered", "continuous"), default="layered")
    p.add_argument("--pca-color-by", choices=("region", "layer", "age", "cohort", "source"), default="region")
    p.add_argument("--pca-time-scale", choices=("linear", "log1p"), default="linear", help="Atlas display only; never changes original dates or genotype PCA")
    p.add_argument("--pca-display-max-age-bp", type=float, default=10000, help="Main atlas age limit; older samples remain in the full-range supplement")
    p.add_argument("--skip-atlas", action="store_true", help="Skip additional publication atlas after run/render")
    p.add_argument("--pca-target-flank-bp", type=int, default=1000000)
    p.add_argument("--pca-exclude-bed", type=Path)
    p.add_argument("--pca-threads", type=int, default=4)
    p.add_argument("--pca-memory-mb", type=int, default=16384)
    p.add_argument("--seed", type=int, default=20261021)
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--pca-max-variants", type=int, choices=range(1, 13), default=12)
    p.add_argument("--count-kind", choices=("hardcall", "dosage"), default="hardcall")
    p.add_argument("--title")
    p.add_argument("--simulated", action="store_true")
    p.add_argument("--basis-dir", type=Path, help="Completed frozen PCA directory for external projection")
    p.add_argument("--reuse-basis-dir", type=Path, help="Verify a completed basis before rendering its existing coordinates")
    p.add_argument("--target-glob", default=configured_path('DIR_GEN', 'ukb/37/imp/chr[0-9]*.pgen'))
    p.add_argument("--project-out", type=Path)
    p.add_argument("--reference-metadata", type=Path, help="External projection reference metadata")
    p.add_argument("--target-metadata", type=Path, help="Optional external target metadata snapshot")
    p.add_argument("--provenance-manifest", type=Path, help="Historical reference-metadata receipt")
    p.add_argument("--target-info", type=Path, help="Exact GRCh37 target variant quality table; see PCA guide")
    p.add_argument("--min-info", type=float, default=.8)
    p.add_argument("--min-overlap", type=float, default=.9)
    p.add_argument("--min-call-rate", type=float, default=.9)
    p.add_argument("--allow-unknown-info", action="store_true", help="Explicit exploratory projection without verified imputation INFO; audited")
    p.add_argument("--allow-palindromic", action="store_true", help="Explicitly retain A/T and C/G markers after external strand verification")
    return p


def metadata_path(a, prefix):
    if a.metadata:
        return a.metadata.expanduser().resolve()
    if a.command == "render":
        raise ValueError("render requires --metadata; no panel or ancient data are needed without overlays")
    if a.pgen_prefix and not a.supplement_xlsx:
        raise ValueError("Custom --pgen-prefix requires --metadata (or an explicit compatible --supplement-xlsx)")
    from evolve import prepare_metadata
    directory = a.out / "metadata"
    metadata = directory / "sample_metadata.tsv"
    receipt = directory / "pca3d_metadata.json"
    supplement = a.supplement_xlsx or a.ancient_dir / "41586_2026_10358_MOESM3_ESM.xlsx"
    psam = Path(str(prefix) + ".psam")
    inputs = {"psam": _fingerprint(psam), "supplement": _fingerprint(supplement),
              "qpadm": _fingerprint(a.qpadm) if a.qpadm else None, "ancestry_model": a.ancestry_model}
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved["inputs"] != inputs or any(not Path(path).is_file() or _sha(path) != digest
                for path, digest in saved["outputs"].items()):
            raise ValueError("Metadata inputs/output changed; use a new --out or supply verified --metadata")
        return metadata
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Existing metadata has no pca3d receipt; pass --metadata explicitly or use a new --out")
    result = prepare_metadata(psam, supplement, directory, qpadm_path=a.qpadm, ancestry_model=a.ancestry_model)
    outputs = {result[key]: _sha(result[key]) for key in ("metadata_path", "qc_path", "manifest_path")}
    _atomic_text(receipt, json.dumps(dict(inputs=inputs, outputs=outputs), indent=2) + "\n")
    return metadata


def overlay_options(a, prefix):
    if not a.selected_variants:
        return {}
    required = {"variant_id", "pvar_index", "chrom", "pos", "ref", "alt", "status"}
    selected, total, seen = [], 0, set()
    for row in _rows(a.selected_variants, required):
        total += 1
        if row["variant_id"] in seen:
            raise ValueError("Duplicate overlay variant_id: " + row["variant_id"])
        seen.add(row["variant_id"])
        if len(selected) < a.pca_max_variants:
            selected.append(row)
    if not selected:
        raise ValueError("--selected-variants is empty")
    folder = a.out / "pca"
    folder.mkdir(parents=True, exist_ok=True)
    bounded = folder / "overlay_targets.tsv"
    if bounded.resolve() == a.selected_variants.expanduser().resolve():
        raise ValueError("Overlay input must not be the generated overlay_targets.tsv; use a separate input file")
    with bounded.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(selected[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    _atomic_text(folder / "overlay_selection.json", json.dumps(dict(
        rule="first rows in explicit input order; bounded visualization only", total_available=total,
        displayed=len(selected), max_variants=a.pca_max_variants, source=_fingerprint(a.selected_variants)), indent=2) + "\n")
    return dict(pgen_path=Path(str(prefix) + ".pgen"), psam_path=Path(str(prefix) + ".psam"),
                selected_variants_path=bounded)


def workflow(a):
    if a.pca_layer_years < 1:
        raise ValueError('--pca-layer-years must be positive')
    if a.command == 'prepare-basis':
        if not a.basis_dir or not a.project_out:
            raise ValueError('prepare-basis requires original --basis-dir and separate --project-out')
        subprocess.run([sys.executable,str(Path(__file__).with_name('pca.py')),'basis','--basis-dir',str(a.basis_dir),'--out',str(a.project_out),'--plink2',a.plink2,
            '--threads',str(a.pca_threads),'--memory-mb',str(a.pca_memory_mb)],check=True)
        return dict(status='complete',operation='separate_unique_coordinate_projection_basis',output=str(a.project_out))
    if a.command in {"project-plan", "project"}:
        if not a.basis_dir:
            raise ValueError("External projection requires --basis-dir with a completed fixed PCA basis")
        destination = a.project_out or a.out / "ukb_projection"
        command = [sys.executable, str(Path(__file__).with_name('pca.py')),'project',"plan" if a.command == "project-plan" else "run", "--basis-dir", str(a.basis_dir),
            "--target-glob", a.target_glob, "--out", str(destination), "--plink2", a.plink2,
            "--threads", str(a.pca_threads), "--memory-mb", str(a.pca_memory_mb),
            "--min-info", str(a.min_info), "--min-overlap", str(a.min_overlap),
            "--min-call-rate", str(a.min_call_rate)]
        reference_metadata = a.reference_metadata or a.metadata
        if reference_metadata: command += ["--reference-metadata", str(reference_metadata)]
        if a.target_metadata: command += ["--target-metadata", str(a.target_metadata)]
        if a.provenance_manifest: command += ["--provenance-manifest", str(a.provenance_manifest)]
        if a.target_info: command += ["--target-info", str(a.target_info)]
        if a.allow_unknown_info: command += ["--allow-unknown-info"]
        if a.allow_palindromic: command += ["--allow-palindromic"]
        subprocess.run(command, check=True)
        return dict(status="planned" if a.command == "project-plan" else "complete",
                    operation="external_fixed_basis_projection", output=str(destination))
    if a.command == "atlas":
        command = [sys.executable, str(Path(__file__).with_name('pca_plot.py')),'atlas',"--out", str(a.out),
            "--replace", a.replace,
            "--time-scale", a.pca_time_scale, "--max-age-bp", str(a.pca_display_max_age_bp),
            "--layer-years", str(a.pca_layer_years), "--theme", a.atlas_theme,
            "--color-by", a.pca_color_by]
        if a.atlas_no_static:
            command.append("--no-static")
        if a.sample_features:
            command += ["--sample-features", str(a.sample_features), "--feature-schema", str(a.feature_schema)]
        if a.source_html:
            if a.eigenvec or a.metadata:
                raise ValueError("atlas uses either --source-html or --eigenvec plus --metadata")
            if a.pca_basis_label:
                raise ValueError("Source HTML retains its original basis label; --pca-basis-label is for shared coordinate input")
            command += ["--source-html", str(a.source_html)]
        elif a.eigenvec and a.metadata:
            command += ["--coordinates", str(a.eigenvec), "--metadata", str(a.metadata)]
            if a.pca_basis_label:
                command += ["--basis-label", a.pca_basis_label]
        else:
            raise ValueError("atlas needs --source-html or both --eigenvec and --metadata")
        if a.simulated:
            command.append("--simulated")
        subprocess.run(command, check=True)
        return dict(status="complete", operation="display_only_existing_coordinates", output=str(a.out))
    if a.source_html:
        raise ValueError("Use atlas for --source-html re-rendering")
    if a.command == "render" and not a.eigenvec:
        raise ValueError("render requires --eigenvec with shared PCA coordinates")
    if a.reuse_basis_dir:
        if a.command != "render" or a.eigenvec.resolve() != (a.reuse_basis_dir / "all.sscore").resolve():
            raise ValueError("Completed-basis reuse requires render with that basis's all.sscore")
        project_validate_completed_coordinates(a.reuse_basis_dir)
        reused_plan = json.loads((a.reuse_basis_dir / "pca_plan.json").read_text())
        saved = reused_plan["inputs"]["metadata"]
        if not a.metadata:
            raise ValueError("Completed-basis rendering requires its verified metadata")
        provenance = frame_historical_identity(a.metadata, recorded=saved)
        if provenance['historical_metadata_verification'] != 'verified':
            raise ValueError('Historical metadata lacks a trusted full SHA256 receipt; default verified render cannot proceed')
        _atomic_text(a.out / 'metadata_provenance.json', json.dumps(provenance, indent=2) + '\n')
        a.pca_basis_label = a.pca_basis_label or reused_plan["settings"]["basis"]
    if a.command != "render" and a.eigenvec:
        raise ValueError("Use render for existing --eigenvec coordinates")
    prefix = (a.pgen_prefix or a.ancient_dir / PREFIX).expanduser().resolve()
    metadata = metadata_path(a, prefix)
    plan = None
    if a.command != "render":
        if a.pgen_prefix and not a.hq_variants:
            raise ValueError("Custom panels require --hq-variants; background QC cannot be inferred from a PGEN")
        if a.dir_pgen_cache:
            from allele_common import cache_pgen_prefix
            prefix, cache_record = cache_pgen_prefix(prefix, a.dir_pgen_cache)
            _atomic_text(a.out / "pca_input_cache.json", json.dumps(cache_record, indent=2) + "\n")
        plan = prepare_pca_plan(prefix, metadata, a.hq_variants or a.ancient_dir / "HighQualityVariants.tsv",
            a.out / "pca", basis=a.pca_basis, maf=a.pca_maf, geno=a.pca_geno,
            n_components=a.pca_components, max_per_group=a.pca_max_per_group,
            time_bin_years=a.pca_time_bin_years, seed=a.seed, exclude_bed=a.pca_exclude_bed,
            target_variants=a.target_variants or a.selected_variants, target_flank_bp=a.pca_target_flank_bp,
            plink2=a.plink2, threads=a.pca_threads, memory_mb=a.pca_memory_mb,
            basis_min_age_bp=0 if a.pca_basis == "ancient_balanced" and a.pca_basis_min_age_bp is None else a.pca_basis_min_age_bp,
            basis_max_age_bp=10000 if a.pca_basis == "ancient_balanced" and a.pca_basis_max_age_bp is None else a.pca_basis_max_age_bp)
        if a.command == "plan":
            return dict(status="planned", manifest_path=plan["manifest_path"], script_path=plan["script_path"])
        subprocess.run(["bash", plan["script_path"]], check=True)
    source = plan["projected_path"] if plan else a.eigenvec
    label = a.pca_basis_label or (a.pca_basis if plan else "User-supplied shared PCA basis")
    overlay = overlay_options(a, prefix)
    result = None
    if a.replace == 'FALSE' and not overlay and not a.title:
        from pca_plot import display_reuse_display as reuse_display
        result = reuse_display(a.html or a.out / 'pca/pca_time.html',
            dict(coordinates=_fingerprint(source), metadata=_fingerprint(metadata)), label,
            layer_years=a.pca_layer_years, time_view=a.pca_time_view,
            color_by=a.pca_color_by, simulated=a.simulated)
    if result is None:
        result = render_pca_time(source, metadata, a.html or a.out / "pca" / "pca_time.html",
            basis_label=label, max_variants=a.pca_max_variants, simulated=a.simulated, title=a.title,
            count_kind=a.count_kind, layer_years=a.pca_layer_years, time_view=a.pca_time_view,
            color_by=a.pca_color_by, **overlay)
    if not a.skip_atlas:
        atlas_command = [sys.executable, str(Path(__file__).with_name('pca_plot.py')),'atlas',"--replace", a.replace,
            "--source-html", result["html_path"], "--out", str(a.out / "pca" / "atlas"),
            "--time-scale", a.pca_time_scale, "--max-age-bp", str(a.pca_display_max_age_bp),
            "--layer-years", str(a.pca_layer_years), "--theme", a.atlas_theme,
            "--color-by", a.pca_color_by]
        if a.atlas_no_static:
            atlas_command.append("--no-static")
        if a.sample_features:
            atlas_command += ["--sample-features", str(a.sample_features), "--feature-schema", str(a.feature_schema)]
        subprocess.run(atlas_command, check=True)
    if overlay:
        centroid = centroid_run_centroids
        from types import SimpleNamespace
        dosages=a.out/'pca'/'selected_variant_dosages.tsv.gz'
        centroid_export_dosages(overlay['pgen_path'],overlay['psam_path'],overlay['selected_variants_path'],dosages,a.count_kind)
        centroid(SimpleNamespace(coordinates=source,metadata=metadata,dosages=dosages,out=a.out/'pca'/'centroids',bin_years=a.pca_time_bin_years))
    return dict(status="complete", **{k: result[k] for k in ("html_path", "manifest_path", "samples", "variants")})


def main_pca(argv=None):
    a = parser().parse_args(argv)
    a.out = a.out.expanduser().resolve()
    a.ancient_dir = a.ancient_dir.expanduser().resolve()
    try:
        if bool(a.sample_features) != bool(a.feature_schema):
            raise ValueError("--sample-features and --feature-schema must be provided together")
        if a.sample_features and (a.skip_atlas or a.command not in {"run", "render", "atlas"}):
            raise ValueError("Sample features require run/render with atlas enabled, or the atlas command")
        context = nullcontext()
        if a.command == 'run':
            if a.eigenvec or a.reuse_basis_dir or a.source_html:
                raise ValueError('Use render/atlas for existing coordinates or HTML; run recalculates from raw genotypes')
            from analysis_integrity import output_run
            prefix = a.pgen_prefix or a.ancient_dir / PREFIX
            protected = [Path(str(prefix) + suffix) for suffix in ('.pgen', '.pvar', '.pvar.zst', '.psam')]
            protected += [a.metadata, a.sample_features, a.feature_schema, a.hq_variants or a.ancient_dir / 'HighQualityVariants.tsv',
                          a.supplement_xlsx or a.ancient_dir / '41586_2026_10358_MOESM3_ESM.xlsx',
                          a.qpadm, a.target_variants, a.selected_variants, a.pca_exclude_bed]
            context = output_run(a.out, a.replace == 'TRUE', protected)
        with context:
            a.out.mkdir(parents=True, exist_ok=True)
            with (a.out / ".pca3d.lock").open("a") as lock:
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    raise ValueError("Another pca3d job is already running in this output directory") from None
                result = workflow(a)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"pca3d: {exc}", file=sys.stderr)
        return 2

def _entry_0_pca3D():
    if len(sys.argv) > 1 and sys.argv[1] == 'centroid':
        raise SystemExit(main_centroid(sys.argv[2:]))
    raise SystemExit(main_pca())



if __name__ == '__main__':
    routes = {'project': _entry_pca_project, 'basis': _entry_pca_basis, 'core': _entry_0_pca3D}
    arguments = list(sys.argv[1:])
    action = arguments.pop(0) if arguments and arguments[0] in routes else 'core'
    if action is None:
        raise SystemExit('Choose a stage: ' + ', '.join(routes))
    sys.argv[1:] = arguments
    routes[action]()
