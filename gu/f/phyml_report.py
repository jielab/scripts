#!/usr/bin/env python3
"""Report predefined GWAS risk cores and cross-check their carriers with IBDmix."""
import argparse, hashlib, json, os, sqlite3
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import dual_lead as dl
from prepare_review import write_tsv, atomic_text, coverage, union
from phyml_gwas_input import read


def local_path(value):
    value=value or ''
    if os.name == 'nt' and str(value).startswith('/mnt/'):
        return Path(str(value)[5] + ':' + str(value)[6:])
    return Path(value)


def validate(copies,con,threshold=.8):
    """Same sample and actual tract; TRACE supports unknown ancestry, not lineage."""
    con.row_factory=sqlite3.Row
    runs={(r['dataset_id'],r['genome_build'],dl.chrom(r['chr']),r['method']):dict(r) for r in con.execute('SELECT * FROM method_runs')}
    for unit,run in runs.items():
        path=local_path(run.get('raw_file',''))
        if not path.is_file(): continue
        text=path.read_text(errors='replace').replace('\\t','\t')
        # A chromosome-level completion row must not turn an uncovered locus
        # into a negative when only a subset of that chromosome was analyzed.
        if any(line.startswith(('loci_file\t','loci\t')) for line in text.splitlines()):
            bed=path.parent/'request.loci.analysis.bed'
            run['_scope']=[]
            if bed.is_file():
                for line in bed.read_text().splitlines():
                    fields=line.split()
                    if len(fields)>=3 and dl.chrom(fields[0])==unit[2]:run['_scope'].append((int(fields[1]),int(fields[2])))
        if unit[3]=='ibdmix':
            meta=dict(line.split('\t',1) for line in text.splitlines() if '\t' in line)
            refs=meta.get('refs','').split()
            run['_tested_lineages']=set()
            if set(refs)&{'Altai','Vindija','Chagyr','Chagyrskaya'}:run['_tested_lineages'].add('Neanderthal')
            if any('denis' in ref.lower() for ref in refs) and (meta.get('background_filter')=='0' or meta.get('export_denisovan')=='1'):
                run['_tested_lineages'].add('Denisovan')
            chrom=unit[2]
            scope=('X_MALE' if 'male_haploid_nonpar' in text else 'X_PAR') if chrom=='X' else 'C'+chrom
            for roster in (path.parent/'samples'/scope/'ALL.txt',path.parent/'samples'/'ALL.txt'):
                if roster.is_file():
                    run['_tested_samples']=set(roster.read_text().split())
                    break
        if unit[3]=='trace':
            sm=path.parent/'samples/trace_sample_map.tsv'
            if sm.is_file():
                nodes=path.parent/'samples/tree_nodes.txt'
                chosen=set(nodes.read_text().split()) if nodes.is_file() else None
                run['_tested_copies']={str(r['sample'])+':'+str(r['haplotype']) for r in dl.rows(sm)
                    if chosen is None or str(r['tree_node_id']) in chosen}
    windows=defaultdict(list)
    for c in copies: windows[(c['dataset_id'],c['genome_build'],c['chr'])].append((c['candidate_start'],c['candidate_end']))
    segments=defaultdict(list)
    for unit,ww in windows.items():
        seen=set()
        for st,en in union(ww):
            for r in con.execute("SELECT * FROM segments WHERE dataset_id=? AND genome_build=? AND chr=? AND method IN ('ibdmix','trace') AND start<? AND end>?",(*unit,en,st)):
                r=dict(r);key=tuple(r.values())
                if key in seen: continue
                seen.add(key);segments[(*unit,r['sample_id'],r['method'])].append(r)
    result=[]
    for c in copies:
        unit=(c['dataset_id'],c['genome_build'],c['chr']);interval=(c['candidate_start'],c['candidate_end'])
        for method in ['ibdmix','trace']:
            run=runs.get((*unit,method),{});complete=run.get('status')=='complete'
            in_scope=('_scope' not in run or coverage(interval,run['_scope'])['union_fraction']>=1-1e-9)
            in_panel=(('_tested_copies' not in run or c['sample_id']+':'+str(c['haplotype']) in run['_tested_copies'])
                      and ('_tested_samples' not in run or c['sample_id'] in run['_tested_samples'])
                      and (method!='ibdmix' or c['lineage'] in run.get('_tested_lineages',set())))
            compared=complete and in_scope and in_panel
            eligible=compared and dl.truth(run.get('evidence_eligible'))
            pool=segments.get((*unit,c['sample_id'],method),[])
            pool=[s for s in pool if method=='trace' or s['source_class']==c['lineage']]
            samephase=[s for s in pool if str(s.get('haplotype'))==str(c['haplotype'])] if method=='trace' else []
            m=coverage(interval,[(s['start'],s['end']) for s in pool]);p=coverage(interval,[(s['start'],s['end']) for s in samephase])
            result.append(dict(c,method=method,method_complete=int(complete),comparison_available=int(compared),evidence_eligible=int(eligible),
                scope_status='in_scope' if in_scope else 'outside_or_unknown_run_scope',panel_status='included_or_panel_not_recorded' if in_panel else 'not_in_method_sample_map',
                availability_note=run.get('availability_note','no_completed_run'),coverage_threshold=threshold,
                overlap_fraction=m['best_single_fraction'] if complete or pool else None,
                union_fraction=m['union_fraction'] if complete or pool else None,
                overlap_pass=int(m['best_single_fraction']>=threshold) if compared else None,
                phase_overlap_pass=int(p['best_single_fraction']>=threshold) if compared and method=='trace' else None,
                phase_note='same_sample_only_unphased' if method=='ibdmix' else 'same_sample_and_stored_haplotype_index',
                lineage_note='same_lineage' if method=='ibdmix' else 'ghost_unknown_not_lineage_validation',
                callability_note='individual_callable_bases_not_available'))
    return result,runs


def attach_validation(rows,validation,runs,by_hap=False):
    grouped=defaultdict(list)
    for v in validation: grouped[v['candidate_id'] if by_hap else (v['locus_key'],v['lineage'])].append(v)
    for row in rows:
        rr=grouped[row.get('candidate_id') if by_hap else (row['locus_key'],row['lineage'])]
        for method in ['ibdmix','trace']:
            run=runs.get((row['dataset_id'],row['genome_build'],row['chr'],method),{})
            candidates=[v for v in rr if v['method']==method]
            vv=[v for v in candidates if v.get('comparison_available',1)];complete=run.get('status')=='complete'
            seen={v['sample_id'] for v in vv};passing={v['sample_id'] for v in vv if v['overlap_pass']==1}
            any_overlap={v['sample_id'] for v in vv if (v.get('overlap_fraction') or 0)>0}
            status='not_run' if not complete else 'no_candidate_to_test' if not candidates else 'not_evaluable_scope_or_panel' if not vv else 'exploratory' if not dl.truth(run.get('evidence_eligible')) else 'overlap_detected' if passing else 'partial_overlap' if any_overlap else 'not_detected'
            row.update({method+'_status':status,method+'_individuals':len(seen) if complete and vv else None,
                method+'_unassessed_copies':len(candidates)-len(vv) if complete else None,
                method+'_any_overlap_individuals':len(any_overlap) if complete and vv else None,
                method+'_supported_individuals':len(passing) if complete and vv else None,
                method+'_support_fraction':len(passing)/len(seen) if complete and seen else None,
                method+'_supported_copies':sum(v['phase_overlap_pass']==1 for v in vv) if method=='trace' and complete and vv else None,
                method+'_candidate_copies':len(vv) if complete and vv else None})



def gwas_rows(path):
    summaries=read(path); details=read(path.with_name('gwas_haplotypes.tsv')); copies=read(path.with_name('gwas_copies.tsv'))
    summaries=[s for s in summaries if s.get('locus_id')]
    result=[]; all_details=[]; all_copies=[]
    for s in summaries:
        s=dict(s)
        lineage=s.get('lineage') or 'Neanderthal'
        if lineage=='Denisova': lineage='Denisovan'
        for key in ('core_start','core_end','lead_pos','tree_pass','n_candidate_haplotypes','n_candidate_copies','n_sites','n_ld_sites'):
            if s.get(key) not in (None,''): s[key]=int(s[key])
        key=hashlib.sha256('|'.join(str(s.get(k,'')) for k in ('dataset_id','genome_build','locus_id','core_start','core_end')).encode()).hexdigest()[:20]
        s.update(locus_key=key,lineage=lineage,core_interval=f"{s['core_start']+1}–{s['core_end']}" if s.get('core_kb') else '未定义',
                 ld_rule='> 0.98 (EUR)',risk_haplotypes=f"{s['n_candidate_haplotypes']} / {s['n_candidate_copies']}" if s.get('n_candidate_haplotypes') not in (None,'') else None)
        result.append(s)
        for d in details:
            if d.get('locus_id')!=s['locus_id'] or d.get('lineage','Neanderthal')!=lineage: continue
            d=dict(d,lineage=lineage,locus_key=key,candidate_id=f"{key}|{lineage}|{d['hap_id']}")
            all_details.append(d)
            if d['role']=='risk':
                for c in copies:
                    if c.get('locus_id')!=s['locus_id'] or c.get('hap_id')!=d['hap_id'] or c.get('lineage','Neanderthal')!=lineage: continue
                    all_copies.append(dict(c,dataset_id=s['dataset_id'],genome_build=s['genome_build'],chr=s['chr'],
                        locus_key=key,lineage=lineage,candidate_id=d['candidate_id'],
                        candidate_start=s['core_start'],candidate_end=s['core_end']))
    return result,all_details,all_copies


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--normalize',required=True,type=Path);ap.add_argument('--database',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path);ap.add_argument('--loci',type=Path)
    args=ap.parse_args()
    summary=[];details=[];copies=[];files=[];trees=[]
    for path in sorted((args.normalize/'phyml').glob('**/final/gwas_loci.tsv')):
        ss,dd,cc=gwas_rows(path);summary+=ss;details+=dd;copies+=cc
        files.extend(path.parent.glob('gwas_*.tsv'));trees+=read(path.with_name('evidence_trees.tsv'))
    with sqlite3.connect(args.database.resolve().as_uri()+'?mode=ro',uri=True) as con:
        validation,runs=validate(copies,con)
    attach_validation(summary,validation,runs);attach_validation(details,validation,runs,True)
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    for name,rr in [('phyml_locus_report',summary),('phyml_haplotype_report',details),('phyml_copy_validation',validation),('phyml_lineage_trees',trees)]:
        write_tsv(out/(name+'.tsv'),rr)
    # This generated report no longer has two competing lead roles.
    (out/'phyml_lead_report.tsv').unlink(missing_ok=True)
    atomic_text(out/'phyml_report_manifest.json',json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
        workflow='gwas_lead_ld_core',n_loci=len({r['locus_key'] for r in summary}),n_lineage_tests=len(summary),n_haplotype_rows=len(details),ld_population='1KG EUR',ld_rule='r2 > 0.98',
        risk_definition='COJO refA and bJ; original lead retained',coverage_threshold=.8,
        validation_denominator='all recurrent risk-haplotype carriers, regardless of tree support; individuals counted once',
        tree_support_note='per-lineage predefined risk+archaic split, bootstrap >=70; not a high-confidence introgression classification',
        input_files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(set(files))]),indent=2))
    print(f"PhyML GWAS report ready: {len({r['locus_key'] for r in summary})} original leads, {len(summary)} lineage tests, {len(details)} recurrent haplotype rows",flush=True)


if __name__=='__main__':main()
