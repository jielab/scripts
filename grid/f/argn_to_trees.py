#!/usr/bin/env python3
"""Convert an ARG-Needle .argn object to a mutation-bearing tskit file."""
from __future__ import annotations
import argparse,gzip,importlib,json,sys
from pathlib import Path
import tskit

def is_ts(x): return hasattr(x,"dump") and hasattr(x,"num_trees") and hasattr(x,"samples")
def unwrap(x):
    if is_ts(x): return x
    if isinstance(x,(tuple,list)):
        for y in x:
            z=unwrap(y)
            if z is not None:return z
    if isinstance(x,dict):
        for y in x.values():
            z=unwrap(y)
            if z is not None:return z
    return None

def add_haps_mutations(ts, path):
    """Place HAPS alleles parsimoniously when arg_to_tskit omits mutations."""
    if not path or ts.num_sites:
        return ts, 0
    metadata = ts.metadata if isinstance(ts.metadata, dict) else {}
    offset = int(metadata.get("offset", 0))
    samples = ts.samples()
    tables = ts.dump_tables()
    tables.sites.metadata_schema = tskit.MetadataSchema.permissive_json()
    added = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as src:
        for line in src:
            z = line.split()
            if len(z) < 6:
                continue
            try:
                bp = int(float(z[2]))
            except ValueError:
                continue
            pos = bp - offset
            if not 0 <= pos < ts.sequence_length:
                continue
            genotypes = [int(x) if x in ("0", "1") else -1 for x in z[5:]]
            if len(genotypes) != len(samples):
                raise SystemExit(
                    f"HAPS sample count mismatch at {z[1]}: "
                    f"{len(genotypes)} haplotypes versus {len(samples)} ARG samples"
                )
            ancestral, mutations = ts.at(pos).map_mutations(
                genotypes, alleles=[z[3], z[4]]
            )
            site_id = tables.sites.add_row(
                position=pos,
                ancestral_state=ancestral,
                metadata={"bp": bp, "id": z[1], "allele0": z[3], "allele1": z[4]},
            )
            mutation_base = tables.mutations.num_rows
            for mutation in mutations:
                parent = (
                    tskit.NULL
                    if mutation.parent == tskit.NULL
                    else mutation_base + mutation.parent
                )
                tables.mutations.add_row(
                    site=site_id,
                    node=mutation.node,
                    derived_state=mutation.derived_state,
                    parent=parent,
                )
            added += 1
    tables.sort()
    return tables.tree_sequence(), added

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--argn",required=True); ap.add_argument("--haps",default=""); ap.add_argument("--out",required=True); ap.add_argument("--home",default=""); a=ap.parse_args()
    if a.home: sys.path.insert(0,a.home)
    import arg_needle_lib
    arg=arg_needle_lib.deserialize_arg(a.argn)
    # Newer builds may expose a method directly.
    for name in ("to_tskit","to_tree_sequence","as_tskit"):
        fn=getattr(arg,name,None)
        if callable(fn):
            try:
                ts=unwrap(fn())
                if ts is not None: break
            except Exception: pass
    else: ts=None
    modules=[arg_needle_lib]
    for mn in ("arg_needle_lib.convert","arg_needle_lib.utils"):
        try: modules.append(importlib.import_module(mn))
        except Exception: pass
    names=("arg_to_ts","arg_to_tree_sequence","arg_to_tskit","convert_arg_to_ts","convert_to_tskit","convert_arg")
    errors=[]
    if ts is None:
      for m in modules:
        for name in names:
          fn=getattr(m,name,None)
          if not callable(fn): continue
          attempts=[(arg,), (arg, None), (arg, []), (arg, None, None), (arg, [], [])]
          for aa in attempts:
            try:
              z=unwrap(fn(*aa))
              if z is not None: ts=z; break
            except Exception as e: errors.append(f"{m.__name__}.{name}{len(aa)}:{type(e).__name__}:{e}")
          if ts is not None: break
        if ts is not None: break
    if ts is None:
        avail=[]
        for m in modules: avail += [f"{m.__name__}.{x}" for x in dir(m) if "ts" in x.lower() or "tree" in x.lower() or "convert" in x.lower()]
        raise SystemExit("No compatible ARG-to-tskit converter found. Available="+",".join(avail[:50])+" Errors="+" | ".join(errors[:10]))
    if ts.num_samples<=0 or ts.num_trees<=0: raise SystemExit("Converted tree sequence is empty")
    ts, added_sites = add_haps_mutations(ts, a.haps)
    tables=ts.dump_tables()
    tables.time_units="generations"
    tables.provenances.add_row(record=json.dumps({
        "schema_version":"1.0.0",
        "software":{"name":"refGen.argn_to_trees","version":"1"},
        "parameters":{"source_method":"needle","node_time_units":"generations",
                      "ancestral_state_source":"tree_parsimony_not_external_AA",
                      "mutations":"HAPS alleles mapped by parsimony"},
    },separators=(",",":")))
    ts=tables.tree_sequence()
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    temp=Path(a.out+'.next'); ts.dump(temp); temp.replace(a.out)
    print(f"trees={ts.num_trees} nodes={ts.num_nodes} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations} added_sites={added_sites} sequence_length={ts.sequence_length}")
if __name__=="__main__": main()
