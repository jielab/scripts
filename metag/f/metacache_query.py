#!/usr/bin/env python3
from __future__ import annotations
import json, os, shutil, subprocess
from pathlib import Path
from common import DIR_LIB, THREADS, atomic_json
from library import completed_library
from library import file_state
from metacache_build import complete_part
import reference_evidence

TOOLS = Path(os.environ.get("METAG_TOOLS", Path.home()/".cache/metag/tools"))
MC=TOOLS/"metacache"
CPU=MC/"metacache"
GPU=MC/"metacache_gpu"

def choose():
    if GPU.is_file() and shutil.which("nvidia-smi"):
        return GPU,"gpu"
    if CPU.is_file():
        return CPU,"cpu"
    raise RuntimeError("MetaCache not installed")

def database_parts():
    root=DIR_LIB/"metacache"/"db"
    saved=completed_library(root)
    parts=[Path(row["db"]) for row in saved.get("results", [])]
    if not saved.get("identity") or any(not complete_part(Path(row["db"]), row.get("identity", saved["identity"]))
                                           for row in saved.get("results", [])):
        raise RuntimeError("MetaCache partition differs from its completion receipt; run ./1.library.sh --metacache-only")
    if not parts:
        raise RuntimeError("No completed MetaCache database partitions")
    return parts


def completed_query(result,evidence,receipt):
    try:
        saved=json.loads(receipt.read_text())
        return (saved['mapping']==file_state(result) and saved['evidence']==file_state(evidence))
    except (OSError,ValueError,KeyError):
        return False

def query(r1:Path,r2:Path|None,name:str,out:Path):
    exe,kind=choose()
    tax=DIR_LIB/"metacache"/"taxonomy"
    parts=database_parts()
    out.mkdir(parents=True,exist_ok=True)
    result_files=[]
    evidence_files=[]
    backends=[]
    changed=False
    inventory=reference_evidence.load_inventory(DIR_LIB)
    for i,db in enumerate(parts,1):
        result=out/f"part{i:04d}.txt"
        evidence=out/f"part{i:04d}.references.tsv.gz"
        receipt=out/f"part{i:04d}.query.complete.json"
        result_files.append(result)
        evidence_files.append(evidence)
        if completed_query(result,evidence,receipt):
            backends.append("reused"); continue
        changed=True
        partial=out/f"part{i:04d}.native.partial.txt"
        partial.unlink(missing_ok=True)
        cmd=[str(exe),"query",str(db),str(r1)]
        if r2:
            cmd += [str(r2),"-pairfiles"]
        cmd += ["-threads",str(THREADS),"-batch-size","512","-refhits","-tophits","-queryids","-lowest","species",
                "-taxids-only","-separate-cols","-out",str(partial)]
        log=out/f"part{i:04d}.log"
        map_temp=result.with_name(result.name+'.partial')
        evidence_temp=evidence.with_name(evidence.name+'.partial')
        attempts=[(exe,kind)]+([(CPU,'cpu_fallback')] if kind=='gpu' else [])
        for attempt,(program,used) in enumerate(attempts):
            cmd[0]=str(program)
            partial.unlink(missing_ok=True)
            print(f'[MetaCache query {i}/{len(parts)}] {used}: {db.name}; log={log}',flush=True)
            with log.open('a' if attempt else 'w') as lf:
                if attempt:lf.write('\n### GPU FAILED/INCOMPLETE; CPU FALLBACK ###\n')
                p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT)
            try:
                if p.returncode:raise RuntimeError(f'exit status {p.returncode}')
                count=reference_evidence.extract(partial,map_temp,evidence_temp,inventory,paired=bool(r2))
                break
            except (RuntimeError,OSError,ValueError) as exc:
                if attempt+1==len(attempts):
                    raise RuntimeError(f'MetaCache query failed/incomplete for {db}: {exc}; see {log}') from exc
                print(f'GPU query failed/incomplete; retry with CPU: {exc}',flush=True)
        map_temp.replace(result)
        evidence_temp.replace(evidence)
        atomic_json(receipt,{'mapping':file_state(result),'evidence':file_state(evidence),'queries':count,'backend':used})
        partial.unlink()
        backends.append(used)

    merged=out/"merged.txt"
    if changed:merged.unlink(missing_ok=True)
    if not merged.is_file():
        partial=out/"merged.partial.txt"
        partial.unlink(missing_ok=True)
        log=out/"merge.log"
        if len(result_files)==1:
            shutil.copyfile(result_files[0],partial)
        else:
            cmd=[str(CPU),"merge",*[str(x) for x in result_files],
                 "-taxonomy",str(tax),"-threads",str(THREADS),"-lowest","species",
                 "-taxids-only","-separate-cols","-out",str(partial)]
            with log.open("w") as lf:
                p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT)
            if p.returncode:
                raise RuntimeError(f"MetaCache merge failed; see {log}")
        if not partial.is_file() or not partial.stat().st_size:
            raise RuntimeError(f"MetaCache produced no merged output; see {log}")
        partial.replace(merged)

    joined=out/'reference-evidence.sqlite'
    evidence_receipt=out/'reference-evidence.json'
    states={str(path):file_state(path) for path in evidence_files}
    try:
        previous=json.loads(evidence_receipt.read_text())
        reusable=previous.get('inputs')==states and previous.get('output')==file_state(joined)
    except (OSError,ValueError):
        reusable=False
    if not reusable:
        print('Collect per-read reference evidence across database partitions',flush=True)
        count=reference_evidence.build_index(evidence_files,joined)
        atomic_json(evidence_receipt,{'inputs':states,'output':file_state(joined),'reference_hit_rows':count})

    atomic_json(out/"query.json",{
        "name":name,"r1":str(r1),"r2":str(r2) if r2 else None,
        "parts":len(parts),"backends":backends,"merged":str(merged),"reference_evidence":str(joined),
        "evidence_kind":"retained index feature hits; not alignment-confirmed or guaranteed unique biological origins"
    })
    return merged

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--r1",type=Path,required=True); ap.add_argument("--r2",type=Path)
    ap.add_argument("--name",required=True); ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args(); print(query(a.r1,a.r2,a.name,a.out))
