#!/usr/bin/env python3
from __future__ import annotations
from collections import Counter, OrderedDict
from contextlib import ExitStack
import csv, gzip, io, json, re, shutil
from pathlib import Path
from common import DIR_LIB, fastq_records, paired_id, read_id, atomic_json
from sample_state import output_state
import reference_evidence

ROOTS={"human":9606,"animal":33208,"plant":33090,"fungi":4751,
       "bacteria":2,"archaea":2157,"virus":10239}
ORDER=("human","animal","plant","fungi","bacteria","archaea","virus")

class Taxonomy:
    def __init__(self,folder:Path):
        self.parent={}; self.rank={}; self.name={}
        with (folder/"nodes.dmp").open() as f:
            for line in f:
                p=[x.strip() for x in line.split("|")]
                self.parent[int(p[0])]=int(p[1]); self.rank[int(p[0])]=p[2]
        with (folder/"names.dmp").open() as f:
            for line in f:
                p=[x.strip() for x in line.split("|")]
                if len(p)>=4 and p[3]=="scientific name":
                    self.name[int(p[0])]=p[1]
        self.cache={}
    def lineage(self,tid):
        if tid in self.cache:return self.cache[tid]
        out=[]; seen=set()
        for _ in range(200):
            if tid in seen or tid not in self.parent: break
            seen.add(tid); out.append(tid)
            p=self.parent[tid]
            if p==tid: break
            tid=p
        self.cache[out[0] if out else tid]=tuple(out)
        return tuple(out)
    def group(self,tid):
        if tid<=0:return "unclassified"
        line=set(self.lineage(tid))
        if 9606 in line:return "human"
        if 33208 in line:return "animal"
        if 33090 in line:return "plant"
        if 4751 in line:return "fungi"
        if 2 in line:return "bacteria"
        if 2157 in line:return "archaea"
        if 10239 in line:return "virus"
        if 2759 in line:return "other_eukaryote"
        return "other"
    def species(self,tid):
        if tid<=0:return 0
        for x in self.lineage(tid):
            if self.rank.get(x)=="species":return x
        return 0

def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+","_",x).strip("._")[:80] or "unknown"

class Writers:
    def __init__(self,root,max_open=32):
        self.root=root; self.max_open=max_open; self.handles=OrderedDict(); self.paths=set()
    def write(self,key,rec):
        self.paths.add(key)
        path=self.root/key
        if key in self.handles:
            h=self.handles.pop(key)
        else:
            path.parent.mkdir(parents=True,exist_ok=True)
            h=gzip.open(path,"ab",compresslevel=1)
        self.handles[key]=h
        h.write(b"".join(rec))
        while len(self.handles)>self.max_open:
            _,old=self.handles.popitem(last=False); old.close()
    def close(self):
        for h in self.handles.values():h.close()
        self.handles.clear()

def mapping_lines(path):
    columns=None
    with path.open(errors="replace") as f:
        for line in f:
            if line.startswith("# TABLE_LAYOUT:"):
                columns=[x.strip() for x in line.split(":",1)[1].split("|")]
            if not line.strip() or line.startswith("#"):continue
            # MetaCache -separate-cols output uses | separators.
            parts=[x.strip() for x in line.rstrip("\n").split("|")]
            if columns is not None:
                if len(parts)!=len(columns) or not {"query_header","rank","taxid"}.issubset(columns):
                    raise RuntimeError(f"Unexpected MetaCache columns: {columns}")
                row=dict(zip(columns,parts))
                qid,rank,tax=row["query_header"].split()[0],row["rank"],row["taxid"]
            elif len(parts)==3:
                qid,rank,tax=parts[0].split()[0],parts[1],parts[2]
            else:
                raise RuntimeError(f"Unexpected MetaCache mapping: {line[:120]!r}")
            if tax=="--":
                tid=0
            elif re.fullmatch(r"-?\d+",tax):
                tid=int(tax)
                if tid<0: tid=0  # sequence-level synthetic ids are not NCBI taxids
            else:
                raise RuntimeError(f"Expected numeric MetaCache taxid, got {tax!r}")
            yield qid,rank,tid

def split(r1:Path,r2:Path|None,mapping:Path,out:Path,*,input_identity=None,before_publish=None,gather=None,evidence=None):
    tax=Taxonomy(DIR_LIB/"metacache"/"taxonomy")
    stage=out.with_name(out.name+".partial")
    if stage.exists():shutil.rmtree(stage)
    if out.exists():raise RuntimeError(f"Output already exists: {out}")
    stage.mkdir(parents=True)
    writers=Writers(stage)
    counts=Counter()
    source_counts=Counter()
    evidence_rows=iter(reference_evidence.grouped(evidence)) if evidence is not None else iter(())
    next_evidence=next(evidence_rows,None)
    stack=ExitStack()
    assignments=csv.DictWriter(stack.enter_context(gzip.open(stage/'read-assignments.tsv.gz','wt',compresslevel=1,newline='')),
        fieldnames=['fragment_index','fragment_id','read_id','mate','classification_unit','assigned_rank','assigned_taxid',
                    'species_taxid','species_name','group','species_status','reference_status','candidate_assembly_accessions',
                    'candidate_reference_sequences'],delimiter='\t')
    assignments.writeheader()
    sources=csv.DictWriter(stack.enter_context(gzip.open(stage/'read-reference-hits.tsv.gz','wt',compresslevel=1,newline='')),
        fieldnames=['fragment_index','fragment_id',*reference_evidence.FIELDS[2:]],delimiter='\t')
    sources.writeheader()
    a=fastq_records(r1); b=fastq_records(r2) if r2 else None
    n=0
    try:
        for qid,rank,tid in mapping_lines(mapping):
            rec1=next(a,None); rec2=next(b,None) if b else None
            if rec1 is None or (b and rec2 is None):
                raise RuntimeError("More MetaCache assignments than FASTQ records")
            rid=paired_id(rec1[0],rec2[0]) if rec2 else read_id(rec1[0])
            # MetaCache may preserve full header; compare first token after optional @.
            got=read_id(qid.encode())
            valid_ids={rid,read_id(rec1[0])}
            if got not in valid_ids:
                raise RuntimeError(f"MetaCache/FASTQ order-ID mismatch: {got!r} != {rid!r}")
            candidates=[]
            if next_evidence is not None:
                if next_evidence[0]<n+1:
                    raise RuntimeError('Reference evidence query IDs are out of order')
                if next_evidence[0]==n+1:
                    candidates=next_evidence[1]
                    next_evidence=next(evidence_rows,None)
            for row in candidates:
                if read_id(row['query_header'].encode()) not in valid_ids:
                    raise RuntimeError('Reference evidence/FASTQ read ID mismatch')
                sources.writerow(dict(fragment_index=n+1,fragment_id=rid,
                                      **{key:row[key] for key in reference_evidence.FIELDS[2:]}))
            group=tax.group(tid)
            species=tax.species(tid)
            status,accessions=reference_evidence.reference_status(candidates)
            if evidence is None:status='not_recorded'
            source_counts[status]+=1
            species_status='assigned_species' if rank=='species' and species else ('higher_rank' if tid else 'unclassified')
            for mate,record in ((1,rec1),(2,rec2)) if rec2 else ((0,rec1),):
                assignments.writerow(dict(fragment_index=n+1,fragment_id=rid,
                    read_id=record[0].decode().strip().split()[0].lstrip('@'),mate=mate,
                    classification_unit='paired_fragment' if rec2 else 'single_read',
                    assigned_rank=rank,assigned_taxid=tid,species_taxid=species,
                    species_name=tax.name.get(species,'') if species else '',group=group,
                    species_status=species_status,reference_status=status,
                    candidate_assembly_accessions=json.dumps(accessions,separators=(',',':')),
                    candidate_reference_sequences=len(candidates)))
            if rank=="species" and species:
                sname=safe_name(tax.name.get(species,str(species)))
                base=f"{group}/{species}_{sname}"
                dest="species"
                reason=""
            else:
                reason="unclassified" if tid==0 else "no_species_in_taxonomic_lineage"
                base=f"split.bad/{group}/{reason}"
                dest="split.bad"
            if rec2:
                writers.write(f"{base}/PE.R1.fastq.gz",rec1)
                writers.write(f"{base}/PE.R2.fastq.gz",rec2)
                reads=2
            else:
                writers.write(f"{base}/SE.fastq.gz",rec1); reads=1
            counts[(dest,group,species,rank,tid,reason)] += 1
            n+=1
        if next(a,None) is not None or (b and next(b,None) is not None):
            raise RuntimeError("Fewer MetaCache assignments than FASTQ records")
        if next_evidence is not None:
            raise RuntimeError('Reference evidence contains more queries than FASTQ')
    finally:
        writers.close()
        stack.close()
        if hasattr(evidence_rows,'close'):evidence_rows.close()

    rows=[]
    for (dest,group,species,rank,tid,reason),frags in sorted(counts.items()):
        rows.append(dict(destination=dest,group=group,species_taxid=species,
                         species_name=tax.name.get(species,"") if species else "",
                         assigned_rank=rank,assigned_taxid=tid,reason=reason,
                         fragments=frags,reads=frags*(2 if r2 else 1)))
    with (stage/"counts.tsv").open("w",newline="",encoding="utf-8") as f:
        fields=list(rows[0]) if rows else ["destination","group","species_taxid","species_name","assigned_rank","assigned_taxid","reason","fragments","reads"]
        w=csv.DictWriter(f,fields,delimiter="\t");w.writeheader();w.writerows(rows)
    reports={"counts.tsv","read-assignments.tsv.gz","read-reference-hits.tsv.gz"}
    if gather is not None:
        from sample_report import snapshot_references
        snapshot_references(gather,stage/"reference-matches.tsv",DIR_LIB)
        reports.add("reference-matches.tsv")
    if before_publish:
        before_publish()
    atomic_json(stage/"split.done",{
        "input_identity":input_identity,
        "outputs":{name:output_state(stage/name) for name in sorted(writers.paths|reports)},
        "fragments":n,"paired":bool(r2),"rows":rows,
        "reference_status_fragments":dict(source_counts),
        "reference_evidence":"All retained index feature hits; single_reported_assembly does not prove a unique biological origin. Paired mates share fragment evidence.",
        "classifier":"MetaCache merged competitive partition result"
    })
    stage.rename(out)
    print(f"PASS split: {n:,} fragments -> {out}")
    return out

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--r1",type=Path,required=True);ap.add_argument("--r2",type=Path)
    ap.add_argument("--mapping",type=Path,required=True);ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args(); split(a.r1,a.r2,a.mapping,a.out)
