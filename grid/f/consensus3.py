"""Read one frozen GRCh38 consensus BED catalog, shared by every PRS channel.

This module consumes the supplied boundaries. It does not calculate LD, extend
the D grid, require convergence receipts, fill cut gaps, or change boundaries.
The Nov16 catalog provenance is documented in README/LD_AUDIT; BED coordinates
alone cannot establish genome build, so the input build is an explicit contract.
"""
from __future__ import annotations

import glob
from pathlib import Path
import re

import numpy as np
import pandas as pd


def consensus_files(source):
    """Resolve a directory, glob, file list, or repeated shared-map dictionary."""
    if isinstance(source, dict):
        sets=[consensus_files(v) for v in source.values()]
        if not sets or any(s != sets[0] for s in sets[1:]):
            raise ValueError("All score channels must use the same consensus BED files")
        return sets[0]
    if isinstance(source, (str, Path)):
        path=Path(source).expanduser()
        pattern=str(path/"chr*.bed") if path.is_dir() else str(path)
        files=[Path(p).resolve() for p in glob.glob(pattern)]
    else:
        files=[]
        for item in source:
            files.extend(consensus_files(item))
    files=sorted(set(files))
    if not files or any(not p.is_file() for p in files):
        raise FileNotFoundError(f"No complete consensus BED input: {source}")
    return files


def read_consensus(source, *, expected_chromosomes=None, build="GRCh38", require_sorted=True):
    suffix = str(build).removeprefix("GRCh").removeprefix("hg")
    if suffix not in {"37","38"}:
        raise ValueError("Expected GRCh37 or GRCh38 BED coordinates")
    chrom_col, start_col, end_col = (name+suffix for name in ("chrom","start","end"))
    files=consensus_files(source)
    rows=[]; file_records=[]
    for path in files:
        previous=None; n=0; observed=set()
        expected_file=re.fullmatch(r"chr([0-9]+)\.bed",path.name)
        with path.open() as stream:
            for line_no,line in enumerate(stream,1):
                if not line.strip() or line.lstrip().startswith(("#", "track ", "browser ")):
                    continue
                fields=line.split()
                if fields[0].lower() in {"chr", "chrom", "chromosome"}:
                    continue
                try:
                    chrom=int(re.sub(r"^chr", "", fields[0], flags=re.I))
                    start,end=int(fields[1]),int(fields[2])
                except (ValueError,IndexError) as exc:
                    raise ValueError(f"Malformed consensus BED {path}:{line_no}") from exc
                if chrom not in range(1,23) or not 0 <= start < end:
                    raise ValueError(f"Invalid 0-based half-open consensus BED {path}:{line_no}")
                if expected_file and chrom!=int(expected_file[1]):
                    raise ValueError(f"Chromosome in {path}:{line_no} disagrees with its filename")
                if require_sorted and previous and ((chrom,start)<(previous[0],previous[1]) or
                                 (chrom==previous[0] and start<previous[2])):
                    raise ValueError(f"Unsorted or overlapping consensus BED {path}:{line_no}")
                previous=(chrom,start,end);observed.add(chrom);n+=1
                rows.append(dict(**{chrom_col:chrom,start_col:start,end_col:end},
                    block_id=f"chr{chrom}:{start}-{end}",source_file=str(path),source_line=line_no))
        if not n:
            raise ValueError(f"Empty consensus BED: {path}")
        stat=path.stat()
        file_records.append(dict(path=str(path),size=stat.st_size,mtime_ns=stat.st_mtime_ns,
                                 chromosomes=sorted(observed),n_blocks=n))
    blocks=pd.DataFrame(rows).sort_values([chrom_col,start_col]).reset_index(drop=True)
    if blocks.block_id.duplicated().any():
        raise ValueError("Duplicated consensus intervals across BED files")
    summary=[]
    for chrom,g in blocks.groupby(chrom_col,sort=True):
        gaps=g[start_col].to_numpy(int)[1:]-g[end_col].to_numpy(int)[:-1]
        if np.any(gaps<0):
            raise ValueError("Overlapping consensus intervals across BED files")
        summary.append(dict(chromosome=int(chrom),n_blocks=len(g),
            n_internal_gaps=int((gaps>0).sum()),internal_gap_bp=int(gaps.sum()),
            first_start0=int(g[start_col].iloc[0]),last_end0=int(g[end_col].iloc[-1])))
    actual=set(blocks[chrom_col])
    if expected_chromosomes is not None and actual!=set(map(int,expected_chromosomes)):
        raise ValueError(f"Consensus chromosomes {sorted(actual)} do not match weights {sorted(set(expected_chromosomes))}")
    blocks.insert(0,"block_index",np.arange(len(blocks),dtype=int))
    provenance=dict(genome_build="GRCh"+suffix,build_source="declared_and_documented_input_provenance",
        coordinates="BED_0_based_half_open",n_blocks=len(blocks),chromosomes=sorted(actual),
        files=file_records,chromosome_summary=summary,shared_boundaries=True,
        gaps="preserved; unassigned SNPs are reported and excluded from regional features",
        ld_recalculated=False,larger_D_required=False,
        interpretation="Frozen supplied consensus catalog; no claim of unrestricted-D optimum or author-coordinate identity")
    return blocks,provenance


def unassigned_reasons(bp, blocks):
    suffix = "37" if "start37" in blocks else "38"
    pos=np.asarray(bp,dtype=int)-1
    return np.where(pos<blocks["start"+suffix].iloc[0],"before_catalog",
        np.where(pos>=blocks["end"+suffix].iloc[-1],"after_catalog","internal_cut_gap"))
