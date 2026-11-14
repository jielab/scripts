"""Strict, streaming input utilities for allele-oriented PHOLE summaries.

No allele is inferred from a variant identifier. Coordinates and allele columns
must refer to the same declared build; this module never performs liftover.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SCHEMA = "phole-signed-input-v2"
DATASETS = {"ukb": "prot", "ckb": "ckb", "decode": "decode", "scot": "scot"}
MISSING = {"", ".", "NA", "NAN", "NULL", "NONE"}
REQUIRED = {"CHR", "POS", "EA", "NEA", "BETA", "SE", "P"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def normal_chr(value):
    text = str(value).strip().upper().removeprefix("CHR")
    if text.endswith(".0"):
        text = text[:-2]
    return {"23": "X", "24": "Y", "25": "MT", "M": "MT"}.get(text, text)


def chrom_order(value):
    c = normal_chr(value)
    return (int(c), "") if c.isdigit() else ({"X": 23, "Y": 24, "MT": 25}.get(c, 999), c)


def normal_build(value):
    value = str(value).strip().lower().replace("grch", "").removeprefix("b")
    value = {"hg19": "37", "hg38": "38"}.get(value, value)
    if value not in {"37", "38"}:
        raise ValueError(f"Explicit BUILD=37 or 38 required, received {value!r}")
    return value


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    fd, name = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def file_identity(path):
    path = Path(path).resolve()
    st = path.stat()
    # Full GWAS can be hundreds of GB: do not scan it solely to make a cache key.
    # This stat fingerprint is explicitly distinguished from a content checksum.
    return {"path": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def summary_identity(assay):
    path=Path(assay["PATH"])
    possible=[Path(assay["INDEX"])] if assay.get("INDEX") else [Path(str(path)+".tbi"),Path(str(path)+".csi")]
    return {"summary":file_identity(path),"indices":[file_identity(p) for p in possible if p.is_file()]}


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix.lower() in {".gz", ".bgz", ".bgzip"} else path.open(encoding="utf-8-sig")


def read_header(path):
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            fields = line.rstrip("\r\n").lstrip("#").split()
            normalized = [x.upper() for x in fields]
            if not REQUIRED.issubset(normalized):
                raise ValueError(f"{path}: complete signed summary header requires {sorted(REQUIRED)}; found {fields}. MAGMA, sig.txt and selected-lead inputs are not accepted.")
            if len(normalized) != len(set(normalized)):
                raise ValueError(f"{path}: duplicate header columns")
            return normalized
    raise ValueError(f"{path}: empty input")


def safe_assay(value):
    value = str(value).strip()
    if value in {"", ".", ".."} or any(c in value for c in "/\\\x00\t\r\n"):
        raise ValueError(f"Unsafe assay identifier: {value!r}")
    return value


def number(value, optional=False):
    if str(value).strip().upper() in MISSING:
        if optional:
            return math.nan
        raise ValueError("missing_number")
    result = float(value)
    if not math.isfinite(result):
        if optional:
            return math.nan
        raise ValueError("nonfinite_number")
    return result


def parse_record(fields, header, default_n=None):
    if len(fields) != len(header):
        raise ValueError("column_count")
    row = dict(zip(header, fields))
    chrom = normal_chr(row["CHR"])
    raw_pos = number(row["POS"])
    if not raw_pos.is_integer() or raw_pos < 1 or chrom not in {*(str(x) for x in range(1, 23)), "X", "Y", "MT"}:
        raise ValueError("invalid_coordinate")
    ea, nea = row["EA"].upper(), row["NEA"].upper()
    if ea == nea:
        raise ValueError("identical_alleles")
    if not re.fullmatch("[ACGT]+", ea) or not re.fullmatch("[ACGT]+", nea):
        raise ValueError("non_sequence_allele")
    beta, se, p = number(row["BETA"]), number(row["SE"]), number(row["P"])
    if se <= 0 or not 0 <= p <= 1 or not math.isfinite(beta / se):
        raise ValueError("invalid_beta_se_p")
    eaf = number(row.get("EAF", "."), optional=True)
    raw_n = row.get("N", ".")
    if str(raw_n).strip().upper() in MISSING and default_n is not None:
        raw_n = default_n
    n = number(raw_n, optional=True)
    if math.isfinite(eaf) and not 0 <= eaf <= 1:
        raise ValueError("invalid_eaf")
    if math.isfinite(n) and n <= 0:
        raise ValueError("invalid_n")
    return {"chr": chrom, "pos": int(raw_pos), "ea": ea, "nea": nea, "beta": beta, "se": se,
            "p": p, "p_text": row["P"], "eaf": eaf, "n": n, "snp": row.get("SNP", "."), "z": beta / se}


def extract_positions(assay, positions, allow_scan=False, tabix="tabix"):
    """Return only the requested coordinates; query indexed files before scanning.

    positions is a set of (normalized chromosome, 1-based position). Empty and
    rejected sites remain absent. A scan is never silently substituted for an
    index that fails, or for an unavailable executable.
    """
    path = Path(assay["PATH"])
    header = read_header(path)
    stats = Counter()
    records = defaultdict(list)
    if not positions:
        return records, {"method": "no_requested_positions", "header": header}
    index_paths = [Path(assay["INDEX"])] if assay.get("INDEX") else [Path(str(path) + ".tbi"), Path(str(path) + ".csi")]
    index = next((p for p in index_paths if p.is_file()), None)
    exe = shutil.which(tabix)
    if index is not None and exe is not None:
        # Nonstandard index paths require htslib's explicit URI syntax. Do not
        # assume that a supplied index will be used by a plain tabix command.
        indexed_path = str(path) if index in [Path(str(path)+".tbi"), Path(str(path)+".csi")] else str(path)+"##idx##"+str(index)
        listed = subprocess.run([exe, "-l", indexed_path], text=True, capture_output=True, check=True).stdout.splitlines()
        contigs = defaultdict(list)
        for c in listed:
            contigs[normal_chr(c)].append(c)
        if any(len(v) > 1 for v in contigs.values()):
            raise ValueError(f"{path}: multiple contig labels normalize to the same chromosome")
        with tempfile.TemporaryDirectory(prefix="phole-query-") as temporary:
            region = Path(temporary)/"regions.tsv"
            stderr = Path(temporary)/"tabix.stderr"
            with region.open("w") as out:
                for chrom, pos in sorted(positions, key=lambda v:(chrom_order(v[0]),v[1])):
                    if chrom in contigs:
                        out.write(f"{contigs[chrom][0]}\t{pos}\t{pos}\n")
                    else:
                        stats["requested_positions_absent_contig"] += 1
            if region.stat().st_size:
                with stderr.open("w") as err:
                    proc = subprocess.Popen([exe, "-R", str(region), indexed_path], stdout=subprocess.PIPE, stderr=err, text=True)
                    try:
                        for line in proc.stdout:
                            _consume(line, header, positions, records, stats, assay)
                    finally:
                        proc.stdout.close()
                        rc = proc.wait()
                    if rc:
                        raise RuntimeError(f"tabix failed for {path}: {stderr.read_text()[-2000:]}")
        method = "tabix_exact_positions"
    else:
        if not allow_scan:
            reason = "no tabix executable" if exe is None else "no .tbi/.csi index"
            raise ValueError(f"{path}: {reason}; index the full summary or explicitly pass --allow-scan")
        with open_text(path) as stream:
            for line in stream:
                _consume(line, header, positions, records, stats, assay)
        method = "explicit_full_file_scan"
    return records, {"method": method, "header": header, "input_fingerprint": file_identity(path), "index": str(index) if index else None, **dict(stats)}


def _consume(line, header, positions, records, stats, assay):
    if not line.strip() or line.startswith("#"):
        return
    fields = line.split()
    if fields and fields[0].upper()==header[0] and len(fields)>1 and fields[1].upper()==header[1]:
        return
    stats["rows_read"] += 1
    # Screen coordinates before expensive allele/numeric parsing during scans.
    try:
        if len(fields)!=len(header):
            stats["malformed_rows"] += 1
            return
        chrom = normal_chr(fields[header.index("CHR")])
        pos = int(fields[header.index("POS")])
        if (chrom,pos) not in positions:
            return
        stats["requested_coordinate_rows"] += 1
        row = parse_record(fields,header,assay.get("N"))
    except (ValueError, IndexError) as exc:
        stats["rejected_"+str(exc)] += 1
        return
    records[chrom,pos].append(row)


def palindromic(ea, nea):
    return len(ea)==len(nea)==1 and {ea,nea} in ({"A","T"},{"C","G"})


def snp_tokens(value):
    return {x.strip() for x in str(value).split(",") if re.fullmatch(r"rs\d+",x.strip())}


def unique_record(rows):
    """Collapse literal repeats only, never silently pick one conflicting row."""
    if not rows:
        return None, "absent"
    first = rows[0]
    for other in rows[1:]:
        if (first["ea"],first["nea"])!=(other["ea"],other["nea"]):
            return None,"duplicate_allele_conflict"
        for k in ("beta","se","p","eaf","n"):
            a,b=first[k],other[k]
            if math.isnan(a) and math.isnan(b):
                continue
            if not math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-12):
                return None,"duplicate_statistic_conflict"
    return first,"unique" if len(rows)==1 else "identical_duplicates_collapsed"


def discovery_record(rows, candidate_snp):
    matching = [r for r in rows if r["snp"]==candidate_snp or snp_tokens(r["snp"]) & snp_tokens(candidate_snp)]
    if matching:
        return unique_record(matching)
    if len(rows)==1:
        if snp_tokens(candidate_snp) and snp_tokens(rows[0]["snp"]):
            return None,"discovery_variant_identifier_conflict"
        return rows[0],"unique_coordinate_id_differs"
    return unique_record(rows) if rows else (None,"absent")


def align_record(row, anchor, ancestry, *, palindrome_tolerance=0.10, palindrome_maf=0.42):
    """Return (oriented beta, oriented EAF, reason), or (None,None,reason)."""
    if row["chr"]!=normal_chr(anchor["CHR"]) or row["pos"]!=int(anchor["POS"]):
        return None,None,"coordinate_mismatch"
    a,b=row["ea"],row["nea"]
    x,y=anchor["EA"],anchor["NEA"]
    if palindromic(x,y):
        if {a,b}!={x,y}:
            return None,None,"allele_mismatch"
        f,g=row["eaf"],float(anchor.get("EAF",math.nan))
        if not math.isfinite(f) or not math.isfinite(g):
            return None,None,"palindrome_missing_frequency"
        if min(f,1-f)>palindrome_maf or min(g,1-g)>palindrome_maf:
            return None,None,"palindrome_near_half"
        source_ancestry=str(anchor.get("ANCESTRY","UNKNOWN")).upper()
        ancestry=str(ancestry).upper()
        if source_ancestry in {"UNKNOWN","",".","MIXED"} or ancestry != source_ancestry:
            return None,None,"palindrome_cross_or_unknown_ancestry"
        direct, flipped=abs(f-g),abs((1-f)-g)
        if min(direct,flipped)>palindrome_tolerance or abs(direct-flipped)<palindrome_tolerance:
            return None,None,"palindrome_frequency_ambiguous"
        sign=1 if direct<flipped else -1
        reason="palindrome_same_ancestry_frequency_"+("direct" if sign==1 else "flip")
    elif (a,b)==(x,y):
        sign,reason=1,"exact"
    elif (b,a)==(x,y):
        sign,reason=-1,"swap"
    elif len(a)==len(b)==len(x)==len(y)==1 and (a.translate(COMPLEMENT),b.translate(COMPLEMENT))==(x,y):
        sign,reason=1,"complement"
    elif len(a)==len(b)==len(x)==len(y)==1 and (b.translate(COMPLEMENT),a.translate(COMPLEMENT))==(x,y):
        sign,reason=-1,"complement_swap"
    else:
        return None,None,"allele_mismatch"
    eaf=row["eaf"] if sign==1 else 1-row["eaf"]
    return sign*row["beta"],eaf,reason


def aligned_at_anchor(rows, anchor, ancestry, **options):
    good=[];reasons=Counter()
    for row in rows:
        beta,eaf,reason=align_record(row,anchor,ancestry,**options)
        if beta is None:
            reasons[reason]+=1
        else:
            good.append({**row,"beta":beta,"eaf":eaf,"ea":anchor["EA"],"nea":anchor["NEA"],"alignment":reason})
    if not good:
        return None, next(iter(reasons),"variant_absent"),dict(reasons)
    row,reason=unique_record(good)
    return row, row["alignment"] if row is not None else reason,dict(reasons)


def validate_manifest(path):
    rows=read_tsv(path)
    required={"DATASET","ASSAY","PROT","PATH","BUILD","INPUT_SCOPE"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Manifest needs {sorted(required)}")
    seen=set();out=[]
    for raw in rows:
        row={k:str(v).strip() for k,v in raw.items() if k is not None}
        d=row["DATASET"].lower();d="ukb" if d=="prot" else d
        if d not in DATASETS:
            raise ValueError(f"Unknown dataset: {d}")
        row["DATASET"]=d;row["ASSAY"]=safe_assay(row["ASSAY"])
        if not row["PROT"] or row["PROT"]==".":
            raise ValueError(f"Missing PROT for {d}/{row['ASSAY']}")
        if (d,row["ASSAY"]) in seen:
            raise ValueError(f"Duplicate dataset/assay: {d}/{row['ASSAY']}")
        seen.add((d,row["ASSAY"]))
        if row["INPUT_SCOPE"].lower()!="full":
            raise ValueError(f"{d}/{row['ASSAY']}: INPUT_SCOPE must be full; selected significant summaries cannot supply negative/missing edges")
        row["BUILD"]=normal_build(row["BUILD"])
        p=Path(os.path.expandvars(os.path.expanduser(row["PATH"])))
        if not p.is_absolute():
            p=Path(path).resolve().parent/p
        p=p.resolve()
        if not p.is_file():
            raise ValueError(f"Summary not found: {p}")
        if p.name.endswith((".sig.txt",".genes.out",".genes.raw")) or re.search(r"(?:^|[._-])(thin|lead|sig|magma)(?:[._-]|$)",p.name,re.I):
            raise ValueError(f"Incomplete/MAGMA input is forbidden: {p}")
        row["PATH"]=str(p)
        if row.get("INDEX"):
            index=Path(row["INDEX"])
            row["INDEX"]=str((Path(path).resolve().parent/index if not index.is_absolute() else index).resolve())
        row.setdefault("ANCESTRY","UNKNOWN");row.setdefault("EFFECT_SCALE","unknown")
        row["ANCESTRY"]=row["ANCESTRY"].upper() or "UNKNOWN"
        row["SYNTHETIC"]=row.get("SYNTHETIC","FALSE").upper()
        if row["SYNTHETIC"] not in {"TRUE","FALSE"}:
            raise ValueError("SYNTHETIC must be TRUE or FALSE when supplied")
        out.append(row)
    return out
