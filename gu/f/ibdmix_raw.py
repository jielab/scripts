"""Paper-time interpretation of native IBDmix calls; no segment exports."""
from dataclasses import dataclass
from pathlib import Path
import bisect
import json
import re
import numpy as np
import pandas as pd

REFS = ("Altai.2013", "Denisova.2013", "Denisova25", "Denisova3", "Chagyr", "Vindija", "Altai")
NEA = {"Altai.2013", "Altai", "Chagyr", "Vindija"}
AFR = {"ESN", "GWD", "LWK", "MSL", "YRI"}
COLUMNS = ["ID", "chrom", "start", "end", "length", "slod", "sites", "positive_lods", "negative_lods", "anc"]

@dataclass(frozen=True)
class RawUnit:
    root: Path
    name: str
    files: tuple
    config: dict

    @property
    def label(self):
        return ", ".join(p.name for p in self.files)


def input_units(folder):
    folder = Path(folder)
    root = folder.parent if folder.name == "raw" else folder
    if not (root / "raw").is_dir():
        raise ValueError(f"Missing native raw directory: {root / 'raw'}")
    found = {}
    for p in sorted((root / "raw").glob("*.txt.gz")):
        m = re.match(r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X))\.", p.name)
        if m: found.setdefault(m[1], []).append(p)
    result = {}
    for unit, files in found.items():
        meta = root / "jobs" / (unit + ".run.meta.tsv")
        if not meta.is_file(): raise ValueError(f"Missing raw calling provenance: {meta}")
        cfg = dict(line.split("\t", 1) for line in meta.read_text().splitlines() if "\t" in line)
        required = {"background_filter", "export_denisovan", "lod_cut", "len_cut", "refs", "genome_build"}
        if not required <= cfg.keys(): raise ValueError(f"Incomplete raw provenance: {meta}")
        spec = root / "metadata" / (unit + ".json")
        if spec.is_file(): cfg.update(json.loads(spec.read_text()).get("raw_interpretation", {}))
        observed = {reference(p,unit)[0] for p in files}
        if observed != set(cfg["refs"].split()): raise ValueError(f"Missing/unexpected references for {unit}: {observed}")
        result[unit[3:]] = RawUnit(root,unit,tuple(files),cfg)
    if not result: raise ValueError(f"No native chr*.txt.gz calls in {root / 'raw'}")
    return result


def reference(path, unit):
    value = path.name[len(unit)+1:].removesuffix(".txt.gz")
    for ref in REFS:
        if value == ref or value.startswith(ref + "."):
            return ref, value[len(ref):].lstrip(".") or "YGB"
    raise ValueError(f"Unknown reference in raw filename: {path}")


def subtract(lo,hi,mask):
    i=max(0,bisect.bisect_right(mask,(lo,float("inf")))-1)
    for left,right in mask[i:]:
        if left>=hi: break
        if right<=lo: continue
        if left>lo: yield lo,left
        lo=max(lo,right)
        if lo>=hi:return
    if lo<hi:yield lo,hi


def native_chunks(path, unit):
    for d in pd.read_csv(path,sep="\t",dtype={"ID":str,"chrom":str},keep_default_na=False,chunksize=200000):
        required={"ID","chrom","start","end","slod"}
        if not required<=set(d):raise ValueError(f"Invalid native header: {path}")
        d["chrom"]=d.chrom.str.removeprefix("chr").replace("23","X")
        d["start"]=pd.to_numeric(d.start,errors="raise")-1
        d["end"]=pd.to_numeric(d.end,errors="raise")-1
        d["slod"]=pd.to_numeric(d.slod,errors="raise")
        if (not d.chrom.eq(unit.name[3:]).all() or d.ID.eq("").any() or not d.start.ge(0).all()
                or not d.end.ge(d.start).all() or not np.isfinite(d.slod).all()):
            raise ValueError(f"Invalid native chromosome, coordinate or score: {path}")
        d=d[d.slod.ge(float(unit.config["lod_cut"])) & (d.end-d.start).ge(int(unit.config["len_cut"])) & d.end.gt(d.start)].copy()
        yield d


def frames(unit):
    cfg=unit.config
    for p in unit.files:
        ref,pop=unit.config.get("file_references",{}).get(str(p), (None,None))
        if ref is None: ref,pop=reference(p,unit.name)
        for d in native_chunks(p,unit):
            d["anc"]=ref
            core_end=int(cfg.get("core_end",0))
            if core_end and len(d):
                d["start"]=d.start.clip(lower=int(cfg.get("core_start",0)));d["end"]=d.end.clip(upper=core_end)
                d=d[d.end.gt(d.start)].copy()
            d["length"]=d.end-d.start
            for col in COLUMNS:
                if col not in d:d[col]=""
            yield d[COLUMNS]


def records(unit):
    for d in frames(unit):yield from d.to_dict("records")


def frame(unit):
    parts=list(frames(unit))
    return pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=COLUMNS)


def units_for_run(run):
    """Support both published flat cohorts and native locus calling runs."""
    import csv, gzip
    from ibdmix_layout import cohort_for, run_file
    run=Path(run);cohort=cohort_for(run)
    if cohort is not None and re.fullmatch(r"chr(?:[0-9]+|X)",run.name):
        return [input_units(cohort)[run.name[3:]]]
    meta=run/"run.meta.tsv"
    cfg=dict(line.split("\t",1) for line in meta.read_text().splitlines() if "\t" in line)
    units=[]
    for population in sorted((run/"samples").glob("*/populations.tsv")):
        native_unit=population.parent.name
        with population.open() as f:groups=list(csv.DictReader(f,delimiter="\t"))
        files=[];references={};chrom=None
        for ref in cfg["refs"].split():
            for group in groups:
                path=run_file(run,f"raw/{native_unit}/{ref}.{group['population']}.raw.txt.gz")
                if not path.is_file():raise ValueError(f"Missing native calls: {path}")
                files.append(path);references[str(path)]=(ref,group['population'])
                if chrom is None:
                    with gzip.open(path,"rt") as f:
                        row=next(csv.DictReader(f,delimiter="\t"),None)
                        if row:chrom=row["chrom"].removeprefix("chr").replace("23","X")
        descriptor=run/"metadata"/(native_unit+".json")
        if cohort is not None:
            from ibdmix_layout import unit_metadata_name
            descriptor=cohort/"metadata"/unit_metadata_name(native_unit,run.name)
        info=json.loads(descriptor.read_text()) if descriptor.is_file() else {}
        chrom=info.get("chromosome",chrom)
        if chrom is None:raise ValueError(f"Empty raw unit needs chromosome metadata: {run}/{native_unit}")
        options={**cfg,**info.get("raw_interpretation",{}),"file_references":references}
        units.append(RawUnit(cohort or run,"chr"+chrom,tuple(files),options))
    if not units:raise ValueError(f"No native sample units: {run}")
    return units
