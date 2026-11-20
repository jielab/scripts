"""Independent chromosome scoring with deterministic, parent-only publication.

The harmonization/inference module stays unchanged, preserving CSx sampler
code receipts. Each child owns its files; no processes share writable arrays.
"""
from pathlib import Path
import gzip
import json
import shutil

import numpy as np
import pandas as pd

import genomics as g
from consensus3 import read_consensus
from parallel_runtime import ordered_map, worker_count


def _resources(threads, chromosomes, memory_mb):
    # PLINK does not use all its requested threads in every phase. Run several
    # chromosomes, giving each a small native pool and its own PLINK RAM budget.
    jobs=worker_count(max(1,int(threads)//4),len(chromosomes),
                      memory_per_worker_gib=memory_mb/1024+2.,reserve_gib=2.)
    inner=max(1,int(threads)//jobs)
    print(f"[scoring] chromosome processes={jobs}; PLINK threads/process={inner}; "
          f"PLINK memory/process={memory_mb} MiB",flush=True)
    return jobs,inner


def _region_chromosome(roster, weights, bed, genotype_dir, output, options):
    g.regional_scores(roster,weights,[bed],genotype_dir,output,**options)
    return str(output)


def _score_chromosome(roster,weights,genotype_dir,output,options):
    g.score_models(roster,weights,genotype_dir,output,**options)
    return str(output)


def score_models(roster_path,weights_paths,genotype_dir,output,*,common_variants=False,
                 plink="plink2",threads=4,memory_mb=4096,batch_size=32):
    catalog=g.index_weight_files(weights_paths)
    chromosomes=sorted(set().union(*(set(v) for v in catalog.values())))
    jobs,inner=_resources(threads,chromosomes,memory_mb)
    if common_variants:
        for chrom in chromosomes:
            if not any(len(w) for w in g.chromosome_weight_sets(catalog,chrom,True).values()):
                jobs=1;break
    if jobs==1:
        return g.score_models(roster_path,weights_paths,genotype_dir,output,common_variants=common_variants,
            plink=plink,threads=threads,memory_mb=memory_mb,batch_size=batch_size)
    roster=g.read_roster(roster_path);names=list(catalog)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    parts=output.parent/"tmp/parallel_scoring";parts.mkdir(parents=True,exist_ok=True)
    options=dict(common_variants=common_variants,plink=plink,threads=inner,memory_mb=memory_mb,batch_size=batch_size)
    arguments=((roster_path,{name:paths.get(chrom,[]) for name,paths in catalog.items()},
                genotype_dir,parts/f"chr{chrom}"/"scores.tsv.gz",options) for chrom in chromosomes)
    totals=np.lib.format.open_memmap(parts/"sums.npy",mode="w+",dtype="float64",shape=(len(roster),len(names)))
    totals[:]=0;audits=[]
    for filename in ordered_map(_score_chromosome,arguments,jobs,inner_threads=1):
        frame=pd.read_csv(filename,sep="\t",dtype={"id":str})
        if frame.id.tolist()!=roster.id.tolist():raise ValueError("Chromosome scoring roster mismatch")
        totals+=frame[names].to_numpy(float)
        audits.append(pd.read_csv(Path(filename).parent/"score_counts.tsv",sep="\t"))
    temporary=output.with_name(output.name+".partial")
    with (gzip.open(temporary,"wt") if output.suffix==".gz" else temporary.open("w")) as stream:
        for start in range(0,len(roster),50000):
            frame=pd.DataFrame(np.asarray(totals[start:start+50000]),columns=names)
            frame.insert(0,"id",roster.id.iloc[start:start+50000].to_numpy())
            frame.to_csv(stream,sep="\t",index=False,header=start==0)
    pd.concat(audits,ignore_index=True).to_csv(output.parent/"score_counts.tsv",sep="\t",index=False)
    temporary.replace(output);del totals;shutil.rmtree(parts)
    return output


def regional_scores(roster_path,weights_paths,block_paths,genotype_dir,output_dir,*,
                    plink="plink2",threads=4,memory_mb=4096,batch_size=32,block_build="GRCh38"):
    catalog=g.index_weight_files(weights_paths)
    chromosomes=sorted(set().union(*(set(v) for v in catalog.values())))
    if any(set(v)!=set(chromosomes) for v in catalog.values()):
        raise ValueError("Every score channel needs the same chromosomes")
    jobs,inner=_resources(threads,chromosomes,memory_mb)
    if jobs==1:
        return g.regional_scores(roster_path,weights_paths,block_paths,genotype_dir,output_dir,
            plink=plink,threads=threads,memory_mb=memory_mb,batch_size=batch_size,block_build=block_build)
    blocks,consensus=read_consensus(block_paths,expected_chromosomes=chromosomes,build=block_build)
    # The existing serial implementation permits an empty chromosome, but
    # rejects a whole channel with no covered SNPs. Keep that exact contract.
    for chrom in chromosomes:
        subset=blocks.loc[blocks.chrom38.eq(chrom)]
        for weights in g.chromosome_weight_sets(catalog,chrom).values():
            if not len(weights) or not np.any(g.assign_blocks(weights.BP38,subset)>=0):
                print("[scoring] empty chromosome/channel: retaining serial validation",flush=True)
                return g.regional_scores(roster_path,weights_paths,block_paths,genotype_dir,output_dir,
                    plink=plink,threads=threads,memory_mb=memory_mb,batch_size=batch_size,block_build=block_build)
    roster=g.read_roster(roster_path)
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    (out/"metadata.json").unlink(missing_ok=True)
    parts=out/"parallel_chromosomes";parts.mkdir(exist_ok=True)
    options=dict(plink=plink,threads=inner,memory_mb=memory_mb,batch_size=batch_size,block_build=block_build)
    arguments=[]
    for chrom in chromosomes:
        bed=parts/f"chr{chrom}.bed"
        blocks.loc[blocks.chrom38.eq(chrom),["chrom38","start38","end38"]].to_csv(bed,sep="\t",index=False,header=False)
        arguments.append((roster_path,{name:paths[chrom] for name,paths in catalog.items()},bed,
                          genotype_dir,parts/f"chr{chrom}",options))
    features=pd.concat([blocks.assign(map=name,population=name,name=name+":"+blocks.block_id)
                        for name in catalog],ignore_index=True)
    features.insert(0,"feature_index",np.arange(len(features)))
    features["n_snps"]=0
    index=dict(zip(features.name,features.feature_index))
    matrix=np.lib.format.open_memmap(out/"scores.npy",mode="w+",dtype="float32",shape=(len(roster),len(features)))
    totals=np.zeros((len(roster),len(catalog)),float);audits=[];checks=[]
    with gzip.open(out/"unassigned_variants.tsv.gz","wt") as missing:
        missing.write("map\tchromosome\tSNP\tTARGET_ID\tBP38\treason\n")
        for directory in ordered_map(_region_chromosome,arguments,jobs,inner_threads=1):
            part=Path(directory)
            ids=pd.read_csv(part/"ids.tsv",sep="\t",dtype=str).id
            if ids.tolist()!=roster.id.tolist():raise ValueError("Chromosome scoring roster mismatch")
            f=pd.read_csv(part/"features.tsv",sep="\t")
            columns=np.asarray([index[name] for name in f.name])
            features.loc[columns,"n_snps"]=f.n_snps.to_numpy()
            values=np.load(part/"scores.npy",mmap_mode="r")
            for start in range(0,len(roster),8192):matrix[start:start+8192,columns]=values[start:start+8192]
            del values
            total=pd.read_csv(part/"covered_global.tsv.gz",sep="\t",dtype={"id":str})
            if total.id.tolist()!=roster.id.tolist():raise ValueError("Covered scoring roster mismatch")
            totals+=total[list(catalog)].to_numpy(float)
            audits.append(pd.read_csv(part/"region_counts.tsv",sep="\t"))
            checks.append(pd.read_csv(part/"sum_check.tsv",sep="\t"))
            with gzip.open(part/"unassigned_variants.tsv.gz","rt") as stream:
                next(stream);shutil.copyfileobj(stream,missing)
    # Children validate sums independently. Publish metadata only after all pass.
    matrix.flush();del matrix
    audit=pd.concat(audits,ignore_index=True);check=pd.concat(checks,ignore_index=True)
    order={name:i for i,name in enumerate(catalog)}
    audit=audit.assign(_order=audit["map"].map(order)).sort_values(["_order","chromosome"]).drop(columns="_order")
    check=check.assign(_order=check["map"].map(order)).sort_values(["_order","chromosome"]).drop(columns="_order")
    for name in catalog:
        if not audit.loc[audit["map"].eq(name),"n_assigned"].sum():raise ValueError(f"No assigned SNPs for {name}")
    features["empty_block"]=features.n_snps.eq(0)
    features.to_csv(out/"features.tsv",sep="\t",index=False)
    blocks.to_csv(out/"catalog.tsv",sep="\t",index=False)
    roster[["id"]].to_csv(out/"ids.tsv",sep="\t",index=False)
    audit.to_csv(out/"region_counts.tsv",sep="\t",index=False)
    check.to_csv(out/"sum_check.tsv",sep="\t",index=False)
    with gzip.open(out/"covered_global.tsv.gz","wt") as stream:
        for start in range(0,len(roster),50000):
            frame=pd.DataFrame(totals[start:start+50000],columns=list(catalog))
            frame.insert(0,"id",roster.id.iloc[start:start+50000].to_numpy())
            frame.to_csv(stream,sep="\t",index=False,header=start==0)
    metadata=json.loads((parts/f"chr{chromosomes[0]}"/"metadata.json").read_text())
    metadata.update(shape=[len(roster),len(features)],consensus=consensus,n_common_blocks=len(blocks),
        n_empty_features=int(features.empty_block.sum()),n_unassigned_variants=int(audit.n_unassigned.sum()),
        chromosome_processes=jobs,plink_threads_per_process=inner)
    (out/"metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    shutil.rmtree(parts)
    return out
