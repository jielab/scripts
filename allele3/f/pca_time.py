"""Fixed-basis PCA planning and offline PC1-PC2-time visualization.

Planning never executes PLINK or downloads inputs.  The generated workflow fits
one background-marker basis and projects every plotted sample through the same
score transform.  Time and target-allele dosage are overlays, not PCA features.
"""
from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shlex
import sys


META_COLUMNS = {"sample_id", "cohort", "age_bp", "include", "stratum"}
LONG_LD_SOURCE = "https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)"
# Explicit GRCh37 table, positions interpreted as 1-based inclusive intervals.
# These coarse historical masks are a sensitivity preset, not a completeness claim.
LONG_LD_GRCH37 = [
    (1, 48000000, 52000000), (2, 86000000, 100500000),
    (2, 134500000, 138000000), (2, 183000000, 190000000),
    (3, 47500000, 50000000), (3, 83500000, 87000000), (3, 89000000, 97500000),
    (5, 44500000, 50500000), (5, 98000000, 100500000),
    (5, 129000000, 132000000), (5, 135500000, 138500000),
    (6, 25000000, 35000000), (6, 57000000, 64000000), (6, 140000000, 142500000),
    (7, 55000000, 66000000), (8, 7000000, 13000000),
    (8, 43000000, 50000000), (8, 112000000, 115000000),
    (10, 37000000, 43000000), (11, 46000000, 57000000),
    (11, 87500000, 90500000), (12, 33000000, 40000000),
    (12, 109500000, 112000000), (20, 32000000, 34500000),
]


def _open(path):
    return (gzip.open(path, "rt", encoding="utf-8-sig", newline="")
            if str(path).endswith((".gz", ".bgz"))
            else open(path, "rt", encoding="utf-8-sig", newline=""))


def _rows(path, required):
    with _open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        names = reader.fieldnames or []
        if len(names) != len(set(names)) or not required.issubset(names):
            raise ValueError(f"{path}: requires unique exact columns {sorted(required)}")
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{line}: wrong number of fields")
            yield row


def _number(value, label):
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: finite number required") from exc
    if not math.isfinite(value):
        raise ValueError(f"{label}: finite number required")
    return value


def _metadata(path):
    result = {}
    for row in _rows(path, META_COLUMNS):
        sid = row["sample_id"]
        if not sid or sid == "." or any(x.isspace() for x in sid) or sid in result:
            raise ValueError(f"Duplicate or invalid metadata sample_id: {sid!r}")
        if row["include"] not in {"0", "1"} or row["cohort"] not in {"ancient", "modern"}:
            raise ValueError(f"{sid}: invalid include/cohort")
        age = _number(row["age_bp"], f"{sid} age_bp")
        if age < 0:
            raise ValueError(f"{sid}: age_bp must be nonnegative")
        if row["include"] == "1" and row["stratum"] in {"", "."}:
            raise ValueError(f"{sid}: included sample needs stratum")
        row["age"] = age
        result[sid] = row
    if not result or not any(r["include"] == "1" for r in result.values()):
        raise ValueError("No included metadata samples")
    return result


def _trajectory_module():
    try:
        from . import trajectory
    except ImportError:
        import trajectory
    return trajectory


def _fingerprint(path):
    path = Path(path).expanduser().resolve()
    stat = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        digest.update(handle.read(65536))
        if stat.st_size > 65536:
            handle.seek(max(65536, stat.st_size - 65536))
            digest.update(handle.read(65536))
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "edge_sha256": digest.hexdigest()}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def _atomic_text(path, contents):
    path = Path(path)
    partial = Path(str(path) + ".partial")
    if partial.exists():
        raise ValueError(f"Interrupted output exists: {partial}")
    try:
        with partial.open("w", encoding="utf-8", newline="") as handle:
            handle.write(contents)
        os.replace(partial, path)
    finally:
        if partial.exists():
            partial.unlink()


def _chrom(value):
    value = str(value)
    if value.lower().startswith("chr"):
        value = value[3:]
    if value.isdigit() and 1 <= int(value) <= 22:
        return int(value)
    return None


def _masks(exclude_bed, target_variants, flank, exclude_long_ld):
    intervals = defaultdict(list)
    if exclude_long_ld:
        for chrom, start, end in LONG_LD_GRCH37:
            intervals[chrom].append((start - 1, end))
    else:
        # Keep the MHC exclusion even in the reduced sensitivity preset.
        intervals[6].append((24999999, 35000000))
    if exclude_bed:
        with _open(exclude_bed) as handle:
            for line_no, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                values = line.split()
                if len(values) < 3:
                    raise ValueError(f"BED line {line_no}: need chrom,start0,end")
                chrom = _chrom(values[0])
                try:
                    start, end = int(values[1]), int(values[2])
                except ValueError as exc:
                    raise ValueError(f"BED line {line_no}: integer bounds required") from exc
                if start < 0 or end <= start:
                    raise ValueError(f"BED line {line_no}: invalid interval")
                if chrom is not None:
                    intervals[chrom].append((start, end))
    if target_variants:
        for row in _rows(target_variants, {"chrom", "pos"}):
            if row.get("build", "GRCh37") not in {"GRCh37", "hg19"}:
                raise ValueError("PCA target exclusions must already be harmonized to GRCh37")
            chrom = _chrom(row["chrom"])
            pos = int(row["pos"])
            if pos <= 0:
                raise ValueError("Target POS must be positive")
            if chrom is not None:
                intervals[chrom].append((max(0, pos - 1 - flank), pos + flank))
    merged = {}
    for chrom, values in intervals.items():
        compact = []
        for start, end in sorted(values):
            if compact and start <= compact[-1][1]:
                compact[-1] = (compact[-1][0], max(end, compact[-1][1]))
            else:
                compact.append((start, end))
        merged[chrom] = compact
    return merged


def _background_ids(hq_path, path, masks):
    counts = Counter()
    starts = {chrom: [pair[0] for pair in values] for chrom, values in masks.items()}
    partial = Path(str(path) + ".partial")
    if partial.exists():
        raise ValueError(f"Interrupted output exists: {partial}")
    try:
        with _open(hq_path) as source, partial.open("w", encoding="utf-8") as dest:
            reader = csv.DictReader(source, delimiter="\t")
            required = {"#CHROM", "POS", "ID", "REF", "ALT"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError("HQ table requires exact #CHROM,POS,ID,REF,ALT header")
            for line, row in enumerate(reader, 2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Malformed HQ row {line}")
                counts["hq_rows"] += 1
                chrom = _chrom(row["#CHROM"])
                if chrom is None or len(row["REF"]) != 1 or len(row["ALT"]) != 1 or row["REF"] not in "ACGT" or row["ALT"] not in "ACGT" or row["REF"] == row["ALT"]:
                    counts["non_autosomal_biallelic_snv"] += 1
                    continue
                pos = int(row["POS"])
                if pos <= 0:
                    raise ValueError(f"Invalid HQ POS at row {line}")
                j = bisect_right(starts.get(chrom, []), pos - 1) - 1
                if j >= 0 and pos - 1 < masks[chrom][j][1]:
                    counts["excluded_region"] += 1
                    continue
                vid = row["ID"]
                if not vid or vid == "." or any(c.isspace() for c in vid):
                    raise ValueError(f"Invalid HQ variant ID at row {line}")
                dest.write(vid + "\n")
                counts["background_ids"] += 1
        if counts["background_ids"] == 0:
            raise ValueError("No autosomal biallelic background SNVs remain in HQ table")
        os.replace(partial, path)
    finally:
        if partial.exists():
            partial.unlink()
    return dict(counts)


def pca_mode(n_basis):
    if not isinstance(n_basis,int) or n_basis<1: raise ValueError('Positive PCA basis size required')
    return 'standard' if n_basis<=5000 else 'approx'


def prepare_pca_plan(pgen_prefix, metadata_path, hq_variants_path, out_dir, *,
                     basis="joint_balanced", maf=0.05, geno=0.1, n_components=10,
                     max_per_group=200, time_bin_years=2000, seed=2026,
                     exclude_bed=None, target_variants=None, target_flank_bp=1000000,
                     plink2="plink2", exclude_long_ld=True, threads=4, memory_mb=8192):
    """Write keep files, streamed HQ IDs, a reviewable script and a manifest.

    ``joint_balanced`` caps each cohort/region/time group; it does not give every
    group equal sample size. ``modern_projection`` fits on all included modern
    samples. Both project all included people through one identical --score
    call. The script has safe complete-run reuse; interrupted stages are rerun.
    """
    if basis not in {"joint_balanced", "modern_projection"}:
        raise ValueError("basis must be joint_balanced or modern_projection")
    if not 0 < maf < .5 or not 0 <= geno < 1:
        raise ValueError("Require 0 < maf < 0.5 and 0 <= geno < 1")
    for label, value, minimum in [("n_components", n_components, 2), ("max_per_group", max_per_group, 1),
                                  ("time_bin_years", time_bin_years, 1), ("threads", threads, 1),
                                  ("memory_mb", memory_mb, 640), ("target_flank_bp", target_flank_bp, 0)]:
        if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
            raise ValueError(f"{label} must be an integer >= {minimum}")
    if not isinstance(seed, int) or not 0 <= seed <= 4294967295:
        raise ValueError("seed must be an unsigned 32-bit integer")
    prefix = Path(pgen_prefix).expanduser().resolve()
    meta_path = Path(metadata_path).expanduser().resolve()
    hq_path = Path(hq_variants_path).expanduser().resolve()
    files = {"metadata": meta_path, "high_quality": hq_path,
             **{ext: Path(str(prefix) + "." + ext) for ext in ("pgen", "pvar", "psam")}}
    if exclude_bed:
        files["exclude_bed"] = Path(exclude_bed).expanduser().resolve()
    if target_variants:
        files["target_variants"] = Path(target_variants).expanduser().resolve()
    inputs = {key: _fingerprint(path) for key, path in files.items()}
    metadata = _metadata(meta_path)
    psam_ids = _trajectory_module().read_psam_ids(files["psam"])
    if set(psam_ids) != set(metadata):
        raise ValueError("PCA metadata must cover exactly the complete original PSAM sample IDs")
    included = [sid for sid in psam_ids if metadata[sid]["include"] == "1"]
    groups = defaultdict(list)
    for sid in included:
        row = metadata[sid]
        epoch = "modern" if row["cohort"] == "modern" else str(int(row["age"] // time_bin_years))
        groups[(row["cohort"], row["stratum"], epoch)].append(sid)
    if basis == "modern_projection":
        basis_ids = [sid for sid in included if metadata[sid]["cohort"] == "modern"]
    else:
        chosen = set()
        for values in groups.values():
            chosen.update(sorted(values, key=lambda sid: hashlib.sha256(f"{seed}\t{sid}".encode()).digest())[:max_per_group])
        basis_ids = [sid for sid in included if sid in chosen]
    if len(basis_ids) < max(50, n_components + 1):
        raise ValueError("Fixed PCA basis needs at least max(50, n_components+1) included samples; do not fit tiny time-specific PCAs")
    mode = pca_mode(len(basis_ids))
    settings = dict(pca_mode=mode, basis=basis, maf=maf, geno=geno, n_components=n_components,
                    max_per_group=max_per_group, time_bin_years=time_bin_years, seed=seed,
                    target_flank_bp=target_flank_bp, exclude_long_ld=exclude_long_ld,
                    plink2=str(plink2), threads=threads, memory_mb=memory_mb,
                    ld_window_kb=200, ld_step=1, ld_r2=0.1)
    signature = hashlib.sha256(json.dumps({"schema": 2, "inputs": inputs, "settings": settings}, sort_keys=True).encode()).hexdigest()
    destination = Path(out_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "pca_plan.json"
    if manifest_path.exists():
        prior = json.loads(manifest_path.read_text())
        if prior.get("signature") != signature:
            raise ValueError("PCA plan inputs/settings changed; use a new output directory to preserve provenance")
        for key, digest in prior["plan_file_sha256"].items():
            path = Path(key)
            if not path.is_file() or _sha(path) != digest:
                raise ValueError(f"Existing PCA plan file is missing or modified: {path}")
        return prior
    # Existing work with no manifest cannot safely be assumed to belong here.
    if any(destination.iterdir()):
        raise ValueError("New PCA plan requires an empty output directory")
    masks = _masks(exclude_bed, target_variants, target_flank_bp, exclude_long_ld)
    all_keep, basis_keep = destination / "all.keep", destination / "basis.keep"
    marker_path, mask_path = destination / "background.ids", destination / "excluded_regions.GRCh37.bed"
    _atomic_text(all_keep, "#IID\n" + "\n".join(included) + "\n")
    _atomic_text(basis_keep, "#IID\n" + "\n".join(basis_ids) + "\n")
    _atomic_text(mask_path, "# GRCh37; 0-based half-open; historical LD/MHC plus target/custom masks\n" +
                 "".join(f"{chrom}\t{start}\t{end}\n" for chrom in sorted(masks) for start, end in masks[chrom]))
    counts = _background_ids(hq_path, marker_path, masks)
    if counts["background_ids"] < n_components + 1:
        raise ValueError("Fewer HQ background SNVs than required PCA dimensions")
    common, ld, pruned, ref, projected = [destination / x for x in ("basis_common", "ld", "all_pruned", "ref", "all")]
    shared = [str(plink2), "--threads", str(threads), "--memory", str(memory_mb)]
    commands = [
        shared + ["--pfile", str(prefix), "--keep", str(basis_keep), "--extract", str(marker_path),
                  "--chr", "1-22", "--snps-only", "just-acgt", "--min-alleles", "2", "--max-alleles", "2",
                  "--maf", str(maf), "--geno", str(geno), "--make-pgen", "--out", str(common)],
        shared + ["--pfile", str(common), "--indep-pairwise", "200kb", "1", "0.1", "--out", str(ld)],
        shared + ["--pfile", str(prefix), "--keep", str(all_keep), "--extract", str(ld) + ".prune.in",
                  "--make-pgen", "--out", str(pruned)],
        shared + ["--pfile", str(pruned), "--keep", str(basis_keep), "--freq", "counts", "--pca",
                  str(n_components), *(["approx"] if mode == "approx" else []), "allele-wts", "vcols=chrom,ref,alt", "--seed", str(seed), "--out", str(ref)],
        shared + ["--pfile", str(pruned), "--read-freq", str(ref) + ".acount", "--score",
                  str(ref) + ".eigenvec.allele", "2", "5", "header-read", "variance-standardize",
                  "cols=+scoresums", "--score-col-nums", f"6-{5+n_components}", "--out", str(projected)],
    ]
    script_path = destination / "commands.sh"
    projected_path = Path(str(projected) + ".sscore")
    checksums = destination / "completed.sha256"
    marker = destination / "PCA_COMPLETE"
    # Check original inputs before any reuse/execution, without rehashing huge PGEN files.
    verifier = """import hashlib,json,pathlib
p=json.loads(pathlib.Path(__import__('sys').argv[1]).read_text())
for spec in p['inputs'].values():
 f=pathlib.Path(spec['path']); s=f.stat(); h=hashlib.sha256()
 with f.open('rb') as r:
  h.update(r.read(65536))
  if s.st_size>65536:
   r.seek(max(65536,s.st_size-65536)); h.update(r.read(65536))
 if s.st_size!=spec['size'] or s.st_mtime_ns!=spec['mtime_ns'] or h.hexdigest()!=spec['edge_sha256']:
  raise SystemExit('PCA input changed; generate a new plan: '+str(f))
for name,want in p['plan_file_sha256'].items():
 h=hashlib.sha256()
 with open(name,'rb') as r:
  for b in iter(lambda:r.read(1048576),b''): h.update(b)
 if h.hexdigest()!=want: raise SystemExit('PCA plan modified: '+name)
"""
    quote = shlex.quote
    script = ["#!/usr/bin/env bash", "set -euo pipefail", "# Reviewable offline PCA workflow; no downloads.",
              "# Incomplete runs recompute stages; only a checksum-verified full result is reused.",
              f"python3 - {quote(str(manifest_path))} <<'ALLELE3_VERIFY_PCA_INPUTS'", verifier.rstrip(),
              "ALLELE3_VERIFY_PCA_INPUTS", f"if [[ -f {quote(str(marker))} ]]; then",
              f"  [[ $(cat {quote(str(marker))}) == {quote(signature)} ]] || exit 1",
              f"  sha256sum --check --status {quote(str(checksums))}", "  exit 0", "fi",
              shlex.join([str(plink2), "--version"])]
    for index, command in enumerate(commands, 1):
        script += [f"# Stage {index}: see pca_plan.json for scientific assumptions.", shlex.join(command)]
    completed = [projected_path, Path(str(ref) + ".acount"), Path(str(ref) + ".eigenvec.allele"), Path(str(ref) + ".eigenval")]
    script += ["sha256sum " + " ".join(quote(str(p)) for p in completed) + " > " + quote(str(checksums)),
               f"printf '%s\\n' {quote(signature)} > {quote(str(marker) + '.partial')}",
               f"mv {quote(str(marker) + '.partial')} {quote(str(marker))}"]
    _atomic_text(script_path, "\n".join(script) + "\n")
    script_path.chmod(0o755)
    basis_set = set(basis_ids)
    group_counts = [{"cohort": key[0], "stratum": key[1], "time_group": key[2],
                     "eligible": len(values), "basis": sum(sid in basis_set for sid in values)}
                    for key, values in sorted(groups.items())]
    manifest = dict(schema_version=1, signature=signature, build="GRCh37", inputs=inputs,
                    settings=settings, script_path=str(script_path), projected_path=str(projected_path),
                    manifest_path=str(manifest_path), commands=commands,
                    plan_file_sha256={str(p): _sha(p) for p in (all_keep, basis_keep, marker_path, mask_path, script_path)},
                    counts={**counts, "samples_eligible": len(included), "samples_basis": len(basis_ids)},
                    group_counts=group_counts, long_ld_source=LONG_LD_SOURCE,
                    coordinate_kind="all samples use PC*_AVG from the same --score transform; never mix with raw eigenvec units",
                    limitations=[
                        "Age is only an overlay; no separate per-epoch PCA is fitted.",
                        "joint_balanced caps cohort/region/time groups, rather than guaranteeing equal group weights.",
                        "Background common SNVs are not proven neutral; use target/MHC/long-LD exclusions as sensitivity analyses.",
                        "Projection can shrink out-of-basis samples; modern-only PCs may miss ancient variation.",
                        "MAF and genotype-missingness filters are evaluated on the fixed basis; LD pruning uses hardcall correlations.",
                        "PGEN imputation/batch differences may affect PCs. Inspect data-source, coverage and region overlays.",
                        "basis_common and all_pruned can occupy several GiB; incomplete stages are rerun, completed outputs are checksum-verified.",
                        "Huge inputs use size/mtime/edge hashes, not full-file cryptographic validation. Do not mutate original files during execution.",
                        "No target exclusion file supplied; provide GRCh37 selected variants before interpreting target overlays." if target_variants is None else "Targets and configured flanks were excluded from the background marker list.",
                    ])
    _atomic_text(manifest_path, json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    return manifest


def _coordinates(path):
    values, names, columns = {}, None, None
    with _open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("##"):
                continue
            cells = line.split()
            if names is None:
                names = [cell.lstrip("#") for cell in cells]
                if len(names) != len(set(names)) or "IID" not in names:
                    raise ValueError("PCA coordinates need a header with unique IID")
                present = [(a, b, kind) for a, b, kind in [
                    ("PC1_AVG", "PC2_AVG", "projected_average"),
                    ("PC1_SUM", "PC2_SUM", "projected_sum"),
                    ("PC1", "PC2", "eigenvector")]
                    if a in names and b in names]
                if not present:
                    raise ValueError("Coordinates require matching PC1/PC2, PC1_AVG/PC2_AVG, or PC1_SUM/PC2_SUM columns")
                columns = present[0]
                if any(kind == "eigenvector" for _, _, kind in present) and len(present) > 1:
                    raise ValueError("Do not mix eigenvector and projected score units in one coordinate table")
                continue
            if len(cells) != len(names):
                raise ValueError(f"PCA coordinate column count mismatch at line {line_no}")
            row = dict(zip(names, cells))
            sid = row["IID"]
            if sid in values or not sid or sid == ".":
                raise ValueError(f"Duplicate or invalid coordinate IID: {sid}")
            values[sid] = [_number(row[columns[0]], sid + " PC1"), _number(row[columns[1]], sid + " PC2")]
    if not values:
        raise ValueError("No PCA coordinates")
    return values, columns[2]


def render_pca_time(eigenvec_path, metadata_path, out_html, *, basis_label="fixed shared PCA",
                    pgen_path=None, psam_path=None, selected_variants_path=None,
                    max_variants=12, simulated=False, title=None, count_kind="hardcall", candidate_path=None):
    """Render one self-contained HTML.  All included sample IDs must have PCs.

    Optional PGEN overlays read only 1..12 preselected, verified original PVAR
    rows.  ALT dosage 0/1/2 never denotes three alleles. Unrepresented panel
    alleles and missing genotypes stay unavailable; no artificial third state.
    """
    if not isinstance(max_variants, int) or not 1 <= max_variants <= 12:
        raise ValueError("max_variants must be an integer from 1 to 12")
    metadata = _metadata(metadata_path)
    coords, coordinate_kind = _coordinates(eigenvec_path)
    extra = set(coords) - set(metadata)
    included = [sid for sid, row in metadata.items() if row["include"] == "1"]
    missing = set(included) - set(coords)
    if extra or missing:
        raise ValueError(f"Coordinate/metadata sample mismatch: {len(extra)} unknown, {len(missing)} included samples missing")
    requested = [pgen_path is not None, psam_path is not None, selected_variants_path is not None]
    if any(requested) and not all(requested):
        raise ValueError("PGEN overlays require pgen_path, psam_path and selected_variants_path together")
    overlays = []
    if all(requested):
        target_count = 0
        for _ in _rows(selected_variants_path, {"variant_id", "pvar_index", "chrom", "pos", "ref", "alt", "status"}):
            target_count += 1
            if target_count > max_variants:
                raise ValueError(f"Overlay file exceeds {max_variants} selected variants; supply a bounded target list")
        if target_count == 0:
            raise ValueError("Overlay target file is empty")
        sample_ids, iterator = _trajectory_module().iter_pgen_variants(pgen_path, psam_path, selected_variants_path, count_kind=count_kind)
        if set(sample_ids) != set(metadata):
            iterator.close()
            raise ValueError("Overlay PSAM must match the complete metadata sample set")
        indexes = {sid: i for i, sid in enumerate(sample_ids)}
        try:
            for row in iterator:
                if len(row["ref"]) != 1 or len(row["alt"]) != 1 or row["ref"] not in "ACGT" or row["alt"] not in "ACGT":
                    raise ValueError("PCA allele overlays currently accept biallelic SNV anchors only")
                dosage = None if row["dosages"] is None else [row["dosages"][indexes[sid]] for sid in included]
                if dosage is not None and any(x is not None and (not math.isfinite(x) or not 0 <= x <= 2) for x in dosage):
                    raise ValueError("ALT dosage must be finite in [0,2], or missing")
                overlays.append({k: row[k] for k in ("variant_id", "chrom", "pos", "ref", "alt", "status", "count_kind", "dosage_evidence")})
                overlays[-1]["values"] = dosage
        finally:
            iterator.close()
    points = []
    for sid in included:
        row = metadata[sid]
        points.append(dict(sample_id=sid, pc1=coords[sid][0], pc2=coords[sid][1], age_bp=row["age"],
                           cohort=row["cohort"], region=row.get("region", row["stratum"]),
                           source=row.get("data_source", row["cohort"]), data_types=row.get("data_types", "."),
                           coverage=row.get("target_coverage", "."), iqs=row.get("iqs", ".")))
    try:
        from plotly.offline import get_plotlyjs
        import plotly
    except ImportError as exc:
        raise RuntimeError("Self-contained PCA HTML needs plotly in the active environment") from exc
    payload = dict(points=points, variants=overlays, simulated=bool(simulated),
                   coordinate_kind=coordinate_kind, basis_label=str(basis_label))
    data_json = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    title = title or "PC1–PC2–时间：固定基底中的个体与等位基因"
    document = _HTML.replace("__TITLE__", html.escape(title)).replace("__PLOTLY__", get_plotlyjs()).replace("__PAYLOAD__", data_json)
    destination = Path(out_html).expanduser().resolve()
    input_paths = {Path(p).expanduser().resolve() for p in
                   (eigenvec_path, metadata_path, pgen_path, psam_path, selected_variants_path)
                   if p is not None}
    if destination in input_paths or Path(str(destination) + ".json") in input_paths:
        raise ValueError("HTML output must not overwrite inputs")
    destination.parent.mkdir(parents=True, exist_ok=True)
    _atomic_text(destination, document)
    manifest = dict(html_path=str(destination), samples=len(points), variants=len(overlays),
                    coordinate_kind=coordinate_kind, basis_label=str(basis_label), simulated=bool(simulated),
                    plotly_version=plotly.__version__, self_contained=True,
                    excluded_metadata_samples=len(metadata) - len(points),
                    third_allele_status="not_inferred_from_0_1_2; inspect_each_target_availability",
                    inputs={"coordinates": _fingerprint(eigenvec_path), "metadata": _fingerprint(metadata_path)},
                    warnings=["ALT copy number 0/1/2 is not three different alleles.",
                              "Coordinates use one fixed basis; positive BP points upward; REF is not assumed ancestral.",
                              "gnomAD aggregate counts are not individual coordinates."])
    from final_figures import static_pca
    manifest["static_figures"] = static_pca(payload, destination.parent, candidate_path)
    manifest_path = Path(str(destination) + ".json")
    _atomic_text(manifest_path, json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    manifest["manifest_path"] = str(manifest_path)
    return manifest


_HTML = r'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
:root{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;color:#152e45;background:#edf3f7;font-size:15px}body{margin:0}header{padding:24px 32px 16px;background:#102d45;color:#f8fbff}h1{font-size:26px;font-weight:650;margin:0 0 10px}header p{margin:7px 0;color:#cadce8;line-height:1.6}.badge{background:#f5b95f;color:#3b2506;padding:8px 14px;border-radius:6px;font-weight:800;margin:6px 0}.wrap{display:grid;grid-template-columns:minmax(450px,1fr) 310px;gap:16px;padding:18px 24px}.panel{background:white;border:1px solid #d5e0e9;border-radius:10px;box-shadow:0 3px 12px #102d4507}.chart{min-height:650px;height:76vh}aside{padding:20px}h2{font-size:17px;margin:0 0 15px}.field{margin:16px 0}.field label{display:block;font-weight:650;margin-bottom:7px}select,input[type=number]{box-sizing:border-box;padding:8px;border:1px solid #bdccd8;border-radius:6px;background:#f8fbfd;width:100%;color:#17344d}input[type=range]{width:100%;accent-color:#2478a2}.checks label{display:inline;font-weight:400;margin-right:14px}.note{font-size:13px;line-height:1.65;color:#496174;border-top:1px solid #e1e9ee;margin-top:18px;padding-top:14px}.status{padding:12px;background:#edf5f8;border-radius:6px;line-height:1.6}.warn{background:#fff3d9;color:#705015;padding:10px;border-radius:6px;line-height:1.5;font-size:13px}.numbers{display:grid;grid-template-columns:1fr 1fr;gap:12px}.sub{font-size:12px;color:#59758b}button{background:#146c94;color:white;border:0;padding:8px 12px;border-radius:6px;cursor:pointer}footer{padding:0 28px 25px;color:#456178;font-size:13px;line-height:1.6}@media(max-width:900px){.wrap{grid-template-columns:1fr;padding:12px}.chart{height:65vh;min-height:470px}header{padding:20px}aside{padding:18px}}
</style><script>__PLOTLY__</script></head><body>
<header><div id="simulated" class="badge" hidden>模拟示意数据 · 不是真实古代 DNA 结果</div><h1>__TITLE__</h1><p>每个点代表一个个体。横轴为共享 PC1、PC2；纵轴向上为更古老的 BP 年代。年代与目标等位基因只作叠加，不参与 PCA。</p><p id="basis"></p></header>
<main class="wrap"><section class="panel"><div id="plotError" class="warn" hidden></div><div id="plot" class="chart"></div></section><aside class="panel">
<h2>查看同一个坐标系</h2><div class="field"><label for="view">视图</label><select id="view"><option value="3d">3D 时间柱</option><option value="2d">2D 时间筛选截面</option></select></div><div class="field checks"><label><input id="ancient" type="checkbox" checked> 古代</label><label><input id="modern" type="checkbox" checked> 现代</label></div>
<div class="field"><label for="color">着色</label><select id="color"><option value="region">大区域</option><option value="cohort">古代 / 现代</option><option value="age">年代 BP</option></select></div>
<div class="field"><label for="region">区域</label><select id="region"><option value="all">全部区域</option></select></div>
<div class="field"><label for="source">数据来源</label><select id="source"><option value="all">全部来源</option></select></div>
<div class="field"><label>时间范围：<span id="rangeLabel"></span> BP</label><input id="young" type="range" min="0" step="100" value="0"><input id="old" type="range" min="0" step="100"><span class="sub">上滑条：年轻边界；下滑条：古老边界。仅改变可见点，坐标轴保持固定。</span></div>
<div class="field"><button id="reset">重置筛选</button></div><div id="counts" class="status"></div><div id="regionLegend" class="field"></div><div id="frequency" class="status" hidden></div><p id="alleleNote" class="warn"></p>
<div class="note">蓝 → 橙的 0/1/2 表示同一个 ALT 的拷贝数：REF/REF、REF/ALT、ALT/ALT。它们不是三个碱基。灰色表示缺失或当前目标未测定。第三碱基需要其独立目标记录和测定状态，不能从所选 ALT 的 0/1/2 推断。</div>
</aside></main><footer>REF 不等于祖先等位基因；最古老样本也可以携带 ALT。图中分离或频率差异可能来自迁移、抽样、imputation 或批次；不能仅凭图判断选择或 pleiotropy。gnomAD 汇总频率不被转换成虚构个体点。此 HTML 含全部绘图依赖，可离线打开。</footer>
<script id="allele3-data" type="application/json">__PAYLOAD__</script><script>
'use strict';
const DATA=JSON.parse(document.getElementById('allele3-data').textContent), pts=DATA.points;
const el=id=>document.getElementById(id), esc=v=>String(v).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
el('simulated').hidden=!DATA.simulated;
el('basis').textContent=DATA.basis_label+' · '+({projected_average:'全体统一投影分数（AVG）',projected_sum:'全体统一投影分数（SUM）',eigenvector:'共享 eigenvector 坐标'}[DATA.coordinate_kind]);
const maxBP=Math.max(100,...pts.map(p=>p.age_bp)), upper=Math.ceil(maxBP/100)*100;
for(const id of ['young','old'])el(id).max=upper;el('old').value=upper;
const regions=[...new Set(pts.map(p=>p.region))].sort(),sources=[...new Set(pts.map(p=>p.source))].sort();
for(const [id,values] of [['region',regions],['source',sources]])for(const v of values){let o=document.createElement('option');o.value=v;o.textContent=v;el(id).appendChild(o)}
DATA.variants.forEach((v,i)=>{let o=document.createElement('option');o.value='v'+i;o.textContent=v.variant_id+' · '+v.ref+'→'+v.alt;el('color').appendChild(o)});
const palette=['#2478a2','#e29135','#389b7c','#b25877','#7662a6','#798a99','#5c8d2c'];const regionColors=Object.fromEntries(regions.map((r,i)=>[r,i<palette.length?palette[i]:'hsl('+((i*137.508)%360)+',58%,45%)']));
function extent(values){let a=Math.min(...values),b=Math.max(...values),m=Math.max((b-a)*.06,Math.abs(a)*.01,.00001);return[a-m,b+m]}
const xRange=extent(pts.map(p=>p.pc1)),yRange=extent(pts.map(p=>p.pc2));
const coordLabel=DATA.coordinate_kind==='eigenvector'?'':' 投影分数';
const layout={margin:{l:0,r:0,b:20,t:20},paper_bgcolor:'#ffffff',font:{family:'system-ui',color:'#29475d'},legend:{orientation:'h',x:.04,y:.98},uirevision:'fixed-shared-pca',scene:{xaxis:{title:{text:'PC1'+coordLabel},range:xRange,autorange:false,backgroundcolor:'#f1f6f9',gridcolor:'#d9e5ed'},yaxis:{title:{text:'PC2'+coordLabel},range:yRange,autorange:false,backgroundcolor:'#f1f6f9',gridcolor:'#d9e5ed'},zaxis:{title:{text:'年代 BP ↑ 越高越古老'},range:[-upper*.025,upper*1.03],autorange:false,backgroundcolor:'#f9fbfc',gridcolor:'#d9e5ed'},aspectratio:{x:1,y:1,z:1.25},camera:{eye:{x:1.45,y:-1.65,z:1.1}}}};
const plotConfig={responsive:true,displaylogo:false,modeBarButtonsToRemove:['toImage']};
let renderedView=null,pendingDraw=null,drawing=false;
async function drawPlot(traces,displayLayout,view){
 pendingDraw={traces,displayLayout,view};if(drawing)return;drawing=true;
 try{
  while(pendingDraw){
   const job=pendingDraw;pendingDraw=null;
   // Plotly.react can retain incompatible subplot internals across 3D/2D.
   // Recreate the graph when its subplot family changes; retain fixed ranges.
   if(renderedView!==job.view){if(renderedView!==null)Plotly.purge('plot');await Plotly.newPlot('plot',job.traces,job.displayLayout,plotConfig)}
   else await Plotly.react('plot',job.traces,job.displayLayout,plotConfig);
   renderedView=job.view;el('plotError').hidden=true;
  }
 }catch(error){renderedView=null;el('plotError').hidden=false;el('plotError').textContent='绘图失败：'+String(error.message||error);console.error('PCA plot failed',error)}
 finally{drawing=false}
}
function update(){
 let lo=Number(el('young').value),hi=Number(el('old').value);if(lo>hi){[lo,hi]=[hi,lo]}
 el('rangeLabel').textContent=lo.toLocaleString()+' – '+hi.toLocaleString();
 const mode=el('color').value,variant=mode.startsWith('v')?DATA.variants[Number(mode.slice(1))]:null,is2d=el('view').value==='2d';
 const selected=pts.map((p,i)=>[p,i]).filter(([p])=>p.age_bp>=lo&&p.age_bp<=hi&&el(p.cohort).checked&&(el('region').value==='all'||p.region===el('region').value)&&(el('source').value==='all'||p.source===el('source').value));
 const traces=[];let present=0;
 for(const cohort of ['ancient','modern'])for(const unavailable of [false,true]){
  const subset=selected.filter(([p,i])=>p.cohort===cohort&&Boolean(variant&&(variant.values===null||variant.values[i]===null))===unavailable);if(!subset.length)continue;
  let numeric=Boolean(variant&&!unavailable)||mode==='age';const colors=subset.map(([p,i])=>unavailable?'#a5b3be':variant?variant.values[i]:mode==='age'?p.age_bp:mode==='cohort'?(cohort==='ancient'?'#2478a2':'#e29135'):regionColors[p.region]);
  if(variant&&!unavailable)present+=subset.length;
  const hover=subset.map(([p,i])=>{let s=esc(p.sample_id)+'<br>'+esc(p.cohort)+' · '+esc(p.region)+'<br>年代 '+p.age_bp.toLocaleString()+' BP<br>来源 '+esc(p.source)+'<br>类型 '+esc(p.data_types)+'<br>capture targets 覆盖 '+esc(p.coverage)+' · IQS '+esc(p.iqs);if(variant){const d=variant.values===null?null:variant.values[i];s+='<br>'+esc(variant.variant_id)+' '+esc(variant.ref)+'→'+esc(variant.alt)+'<br>ALT '+(variant.count_kind==='hardcall'?'拷贝数':'dosage')+': '+(d===null?'缺失 / 未测定':d);if(d!==null&&variant.count_kind==='hardcall')s+=' ('+[variant.ref+'/'+variant.ref,variant.ref+'/'+variant.alt,variant.alt+'/'+variant.alt][d]+')';s+='<br>0/1/2 不编码第三碱基；请查看独立目标测定状态'}return s});
  let marker={size:cohort==='modern'?(is2d?6:4):(is2d?4:3),opacity:.78,symbol:cohort==='modern'?'diamond':'circle',color:colors};if(numeric){marker.colorscale=variant?[[0,'#2452a2'],[.5,'#e7e7df'],[1,'#df7b25']]:'Viridis';marker.cmin=0;marker.cmax=variant?2:upper;marker.showscale=traces.every(t=>!t.marker.showscale);marker.colorbar={title:{text:variant?'ALT 拷贝数 / dosage':'年代 BP'},x:1.01,len:.6,tickvals:variant?[0,1,2]:undefined}}
  traces.push({type:is2d?'scattergl':'scatter3d',mode:'markers',name:(cohort==='ancient'?'古代':'现代')+(unavailable?' · 未知':''),x:subset.map(([p])=>p.pc1),y:subset.map(([p])=>p.pc2),...(is2d?{}:{z:subset.map(([p])=>p.age_bp)}),text:hover,hoverinfo:'text',marker});
 }
 const displayLayout=is2d?{...layout,scene:undefined,margin:{l:65,r:50,b:60,t:35},xaxis:{title:{text:'PC1'+coordLabel},range:xRange,autorange:false,gridcolor:'#e3ebf0'},yaxis:{title:{text:'PC2'+coordLabel},range:yRange,autorange:false,gridcolor:'#e3ebf0'},uirevision:'fixed-shared-pca-2d'}:layout;
 // Plotly may mutate layout objects. Keep the immutable axis template separate.
 drawPlot(traces,JSON.parse(JSON.stringify(displayLayout)),is2d?'2d':'3d');
 const anc=selected.filter(([p])=>p.cohort==='ancient').length;el('counts').textContent='显示 '+selected.length.toLocaleString()+' / '+pts.length.toLocaleString()+' 人：古代 '+anc.toLocaleString()+'，现代 '+(selected.length-anc).toLocaleString()+(variant?'；可用 ALT 值 '+present.toLocaleString():'');
 el('regionLegend').hidden=mode!=='region';el('regionLegend').innerHTML=mode==='region'?regions.map(r=>'<span style="display:inline-block;margin:4px 12px 4px 0"><span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:'+regionColors[r]+';margin-right:5px"></span>'+esc(r)+'</span>').join(''):'';
 el('frequency').hidden=!variant;if(variant){const vals=variant.values===null?[]:selected.map(([,i])=>variant.values[i]).filter(v=>v!==null);const n=vals.length,ac=vals.reduce((a,b)=>a+b,0),an=2*n;el('frequency').textContent='当前筛选个体中的 '+variant.alt+'：'+(n?'AF = '+(ac/an).toFixed(4)+'；'+(variant.count_kind==='dosage'?'expected AC':'AC')+' = '+Number(ac.toFixed(3))+'，AN = '+an:'AF / AC / AN 不可用')+'。这是所显示样本的描述性频率，不是总体频率估计或选择检验。'}
 el('alleleNote').textContent=variant?variant.variant_id+'：'+variant.ref+' / '+variant.alt+'；'+(variant.status==='assayed_pass'?'已通过面板位点筛选，灰色点为个体缺失。':'该目标未通过当前面板可测定条件：'+variant.status+'，不填 0。'):'选择目标位点后可叠加 ALT 拷贝数；横纵 PCA 坐标保持不变。第三碱基不能从 0/1/2 编码推断。';
}
for(const id of ['view','color','region','source','ancient','modern','young','old'])el(id).addEventListener('input',update);
el('reset').onclick=()=>{el('ancient').checked=el('modern').checked=true;el('region').value=el('source').value='all';el('young').value=0;el('old').value=upper;update()};
update();
</script></body></html>'''
