"""Separate corrected group summaries; never rewrite signed split results."""
import argparse
from collections import Counter
from pathlib import Path
from result_output import managed_outputs, prepare as prepare_result, finish as finish_result, replace_bool
from common import atomic_json, load_config, sha256, code_identity
from qc_paths import configured_split_root
from storage_layout import metadata_dir
from split_summary import load_sample, table
from taxonomy import Taxonomy


@managed_outputs
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replace',type=replace_bool,default=True);p.add_argument('--config'); p.add_argument('--root'); p.add_argument('--sample',action='append',required=True)
    p.add_argument('--out',required=True); p.add_argument('--allow-local',action='store_true')
    a=p.parse_args(); cfg=load_config(a.config,a.root)
    out=Path(a.out).absolute()
    if not prepare_result(out,replace=a.replace,cfg=cfg):return
    tax=Taxonomy(cfg['taxonomy_dir']); corrected=[]; totals=[]; sources={}
    for sample in a.sample:
        if Path(sample).is_absolute() or '..' in Path(sample).parts: p.error('Unsafe sample ID')
        folder=configured_split_root(cfg)/sample
        summary,rows,_=load_sample(folder,sample)
        counts=Counter()
        for row in rows:
            info=tax.describe(int(row['assigned_taxid']))
            updated={**row,'sample_id':sample,'original_group':row['group'],'corrected_group':info['group']}
            corrected.append(updated); counts[(info['group'],row['status'])]+=int(row['fragments'])
        if sum(counts.values())!=int(summary['fragments']): raise RuntimeError('Count conservation failed')
        totals.extend(dict(sample_id=sample,group=g,status=s,fragments=n,fraction=n/summary['fragments']) for (g,s),n in sorted(counts.items()))
        sources[sample]=sha256(metadata_dir(folder)/'counts.tsv')
    table(out/'corrected-taxids.tsv',list(corrected[0]),corrected)
    table(out/'corrected-groups.tsv',['sample_id','group','status','fragments','fraction'],totals)
    atomic_json(out/'complete.json',dict(status='complete',code_identity=code_identity(),input_count_sha256=sources,
        taxonomy={n:sha256(Path(cfg['taxonomy_dir'])/n) for n in ('nodes.dmp','names.dmp')},
        policy='Display groups only; exact taxids, species counts and committed FASTQs preserved'))

    finish_result(out)

if __name__=='__main__': main()
