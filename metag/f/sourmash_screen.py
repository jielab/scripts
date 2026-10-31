#!/usr/bin/env python3
from pathlib import Path
from common import DIR_LIB, SOURMASH_K, SOURMASH_SCALED, THREADS, run, require
from library import completed_library, file_state

def screen(r1:Path,r2:Path|None,name:str,out:Path):
    require("sourmash")
    db=DIR_LIB/"sourmash"/"references.zip"
    receipt = completed_library(db.parent)
    if not db.is_file() or receipt.get("output") != file_state(db):
        raise RuntimeError("Sourmash database differs from its completion receipt; run ./1.library.sh --sourmash-only")
    out.mkdir(parents=True,exist_ok=True)
    sig=out/f"{name}.sig.gz"
    gather=out/f"{name}.gather.csv"
    if not sig.is_file():
        partial=out/f"{name}.partial.sig.gz"
        partial.unlink(missing_ok=True)
        cmd=["sourmash","scripts","singlesketch",r1]
        if r2: cmd.append(r2)
        cmd += ["-p",f"k={SOURMASH_K},scaled={SOURMASH_SCALED},dna,abund","-o",partial,"--name",name]
        run(cmd)
        if not partial.is_file() or not partial.stat().st_size:
            raise RuntimeError("Sourmash did not produce a sample sketch")
        partial.replace(sig)
    if not gather.is_file():
        partial=out/f"{name}.partial.gather.csv"
        partial.unlink(missing_ok=True)
        run(["sourmash","scripts","fastgather",sig,db,"-o",partial,"--cores",str(THREADS)])
        if not partial.is_file():
            # branchwater 0.9.15 exits successfully without a CSV for no matches.
            partial.write_text("match_name,match_filename,match_md5,f_unique_to_query\n")
            print("No sourmash gather matches above the reporting threshold")
        partial.replace(gather)
    return gather

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--r1",type=Path,required=True); ap.add_argument("--r2",type=Path)
    ap.add_argument("--name",required=True); ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args(); screen(a.r1,a.r2,a.name,a.out)
