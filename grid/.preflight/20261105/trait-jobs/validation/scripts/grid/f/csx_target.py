#!/usr/bin/env python3
"""Build the PRS-CSx target SNP list from chromosome PVAR/BIM metadata."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


@contextmanager
def metadata_text(path):
    if path.suffix != ".zst":
        with path.open() as stream:
            yield stream
        return
    # Read metadata only; no genotype decompression or conversion is needed.
    with subprocess.Popen(["zstdcat", str(path)], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True) as process:
        try:
            yield process.stdout
            error = process.stderr.read()
            if process.wait():
                raise ValueError(f"Cannot decompress {path}: {error.strip()}")
        finally:
            if process.poll() is None:
                process.terminate()
            process.stdout.close()
            process.stderr.close()


def variant_rows(path):
    with metadata_text(path) as stream:
        columns = None
        for number, line in enumerate(stream, 1):
            row = line.split()
            if not row or line.startswith("##"):
                continue
            if path.suffix == ".bim":
                if len(row) != 6:
                    raise ValueError(f"{path}:{number}: expected six BIM columns")
                yield row[0], row[3], row[1], row[5], row[4]
                continue
            if line.startswith("#"):
                names = [name.lstrip("#") for name in row]
                required = ("CHROM", "POS", "ID", "REF", "ALT")
                if not set(required) <= set(names) or len(set(names)) != len(names):
                    raise ValueError(f"{path}:{number}: invalid PVAR header")
                columns = [names.index(name) for name in required]
                continue
            if columns is None or len(row) <= max(columns):
                raise ValueError(f"{path}:{number}: missing header or incomplete PVAR row")
            yield tuple(row[i] for i in columns)


def build_target_bim(target_dir, chromosomes):
    chromosomes = sorted(map(int, chromosomes))
    if not chromosomes or len(set(chromosomes)) != len(chromosomes) or not set(chromosomes) <= set(range(1, 23)):
        raise ValueError("Use unique autosomes from 1 to 22")
    directory = Path(target_dir).resolve()
    prefix = directory / "ukb_array"
    bim, manifest = prefix.with_suffix(".bim"), directory / "ukb_array.bim.json"
    with (directory / "ukb_array.bim.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        # Keep one permanent list for all available autosomes. A run restricted
        # to one chromosome must not replace it with a partial SNP list.
        sources = []
        for chrom in range(1, 23):
            source_prefix = directory / f"chr{chrom}"
            path = next((Path(str(source_prefix) + suffix) for suffix in (".pvar", ".pvar.zst", ".bim")
                         if Path(str(source_prefix) + suffix).is_file()
                         and Path(str(source_prefix) + suffix).stat().st_size), None)
            if path is None:
                if chrom in chromosomes:
                    raise FileNotFoundError(f"Missing chr{chrom}.pvar, .pvar.zst or .bim under {directory}")
                continue
            sources.append(dict(chromosome=chrom, path=str(path.resolve()), sha256=sha256(path)))
        contract = dict(schema="grid.csx.target_bim.v1",
                        chromosomes=[source["chromosome"] for source in sources], sources=sources)
        try:
            saved = json.loads(manifest.read_text())
        except (OSError, ValueError):
            saved = {}
        if not isinstance(saved, dict):
            saved = {}
        owned = isinstance(saved.get("contract"), dict) and saved["contract"].get("schema") == contract["schema"]
        if (saved.get("contract") == contract and bim.is_file()
                and sha256(bim) == saved.get("bim_sha256")):
            print(f"[csx] Target BIM reused: {saved['snps']} SNPs; {bim}", file=sys.stderr)
            return prefix, saved
        counts, seen = [], set()
        fd, temporary = tempfile.mkstemp(prefix="ukb_array-", suffix=".bim.tmp", dir=directory)
        temporary = Path(temporary)
        try:
            with os.fdopen(fd, "w") as output:
                for source in sources:
                    chrom = source["chromosome"]
                    count = dict(chromosome=chrom, records=0, snps=0, excluded_non_snp=0, excluded_missing_id=0)
                    for code, position, snp, ref, alt in variant_rows(Path(source["path"])):
                        count["records"] += 1
                        code = code.removeprefix("chr").removeprefix("CHR")
                        if int(code) != chrom or int(position) < 1:
                            raise ValueError(f"{source['path']}: invalid chromosome/position for {snp}")
                        ref, alt = ref.upper(), alt.upper()
                        if ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"} or ref == alt:
                            count["excluded_non_snp"] += 1
                            continue
                        if snp in {".", "0"}:
                            count["excluded_missing_id"] += 1
                            continue
                        if snp in seen:
                            raise ValueError(f"Duplicate target SNP ID: {snp}; resolve ambiguous IDs before CSx")
                        seen.add(snp)
                        # PLINK BIM order: CHR, ID, cM, BP, A1(ALT), A2(REF).
                        # PRS-CSx consumes SNP/allele metadata; it does not use cM.
                        output.write(f"{chrom}\t{snp}\t0\t{int(position)}\t{alt}\t{ref}\n")
                        count["snps"] += 1
                    if not count["snps"]:
                        raise ValueError(f"{source['path']}: no usable biallelic A/C/G/T SNPs")
                    counts.append(count)
            for source in sources:
                if sha256(source["path"]) != source["sha256"]:
                    raise ValueError(f"Target metadata changed during conversion: {source['path']}")
            checksum = sha256(temporary)
            # Adopt an identical existing file, but do not overwrite a user's
            # different BIM unless our provenance identifies it as generated.
            if bim.exists() and bim.stat().st_size and not owned:
                if sha256(bim) != checksum:
                    raise ValueError(f"Existing BIM has no matching generated provenance: {bim}; "
                                     "use --bim-prefix to select it explicitly")
            else:
                os.replace(temporary, bim)
            saved = dict(contract=contract, snps=len(seen), chromosomes=counts, bim_sha256=checksum)
            marker = manifest.with_suffix(".json.tmp")
            marker.write_text(json.dumps(saved, indent=2) + "\n")
            os.replace(marker, manifest)
        finally:
            temporary.unlink(missing_ok=True)
    excluded = sum(c["records"] - c["snps"] for c in counts)
    print(f"[csx] Target BIM generated: {len(seen)} SNPs; excluded {excluded} unusable variants; {bim}",
          file=sys.stderr)
    return prefix, saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir-gen", required=True)
    parser.add_argument("--chrs", required=True, help="Space/comma separated autosomes")
    args = parser.parse_args()
    try:
        prefix, _ = build_target_bim(args.dir_gen, args.chrs.replace(",", " ").split())
    except (OSError, ValueError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")
    print(prefix)


if __name__ == "__main__":
    main()
