#!/usr/bin/env python3
"""Post-discovery, multi-lineage PhyML report and candidate-interval validation.

Consumes saved sequences, trees and normalized segment calls; never changes
discovery thresholds or uses a named SNP to select introgressed haplotypes.
"""
from __future__ import annotations
import argparse, hashlib, json, os, sqlite3
from datetime import datetime, timezone
from collections import defaultdict
from pathlib import Path
import numpy as np
import dual_lead as dl
from prepare_review import read_authoritative, write_tsv, atomic_text, coverage, union


def local_path(value):
    value=value or ''
    if os.name == 'nt' and str(value).startswith('/mnt/'):
        return Path(str(value)[5] + ':' + str(value)[6:])
    return Path(value)


def best_tag(data, best, target, variant_counts):
    """Rank SNPs against a fixed copy target using exact phased allele counts."""
    ranked=[]
    weights=np.array([len(h['_copies'] & target) for h in data['haps']])
    for v in dl.build_variants(data['_auth'], data, best):
        i=data['_position_index'][v['pos_1based']]
        called,positive=variant_counts[(i,v['tag_allele'])]
        alleles=data['_matrix'][:,i]
        target_called=int(weights[np.isin(alleles,[v['ref'],v['alt']])].sum())
        tp=int(weights[alleles==v['tag_allele']].sum())
        m=dl.binary_metrics(tp,positive-tp,target_called-tp,called-positive-target_called+tp,len(target),len(data['copies']))
        if dl.qualify(m,.95,2)=='eligible': ranked.append(dict(v,**m))
    ranked.sort(key=dl.rank_tag)
    b=ranked[0] if ranked else None
    ties=[r for r in ranked if b and all(round(r[k],12)==round(b[k],12) for k in ['r2','f1','call_rate'])]
    return b,len(ties)


def annotate_tags(row,data,best,target,counts,iv):
    b,nt=best_tag(data,best,target,counts)
    im=dl.marker_metrics(iv,target,data['copies'])
    row.update(input_pos_1based=iv['pos_1based'] if iv else None,
        input_ref=iv['ref'] if iv else None,input_alt=iv['alt'] if iv else None,
        input_allele_copy_counts=json.dumps({a:sum(iv['_alleles'].get(cp)==a for cp in target) for a in [iv['ref'],iv['alt']]}) if iv else None,
        input_tag_allele=im['tag_allele'],input_tag_r2=im['r2'],
        best_tag_snp=(b.get('snp_id') or b['variant_key']) if b else None,
        best_tag_pos_1based=b['pos_1based'] if b else None,
        best_tag_allele=b['tag_allele'] if b else None,best_tag_r2=b['r2'] if b else None,
        best_tag_ppv=b['ppv'] if b else None,best_tag_sensitivity=b['sensitivity'] if b else None,
        n_equivalent_best_tags=nt,
        best_tag_status=('strong_in_sample_proxy' if b['r2']>=.8 else 'weak_in_sample_proxy') if b else 'no_eligible_single_snp_tag',
        input_vs_best_snp_r2=dl.pair_r2(iv,b,data['copies'],im['tag_allele'],b['tag_allele']) if b else None)


def locus_rows(auth):
    ident=dict(locus_key=dl.coord_key(auth),dataset_id=auth['dataset_id'],genome_build=auth['genome_build'],
        locus_id=auth['locus_id'],chr=dl.chrom(auth['chr']),core_start=int(auth['core_start']),core_end=int(auth['core_end']),
        input_lead_snp=auth.get('_input_lead_override') or auth.get('name') or auth['locus_id'])
    data,state,source_files=dl.load_data(auth)
    auth['_report_source_files']=source_files
    if data is None:
        return [dict(ident,lineage='',call='not_evaluable',reason=state['status'],n_candidate_haplotypes=0,n_candidate_copies=0)],[],[],[]
    ident['n_tested_copies']=len(data['copies'])
    ident['n_tested_individuals']=len({cp.rsplit(':',1)[0] for cp in data['copies']})
    data['_auth']=auth
    data['_matrix']=np.array([list(h['seq']) for h in data['haps']])
    data['_position_index']={s['pos']:i for i,s in enumerate(data['sites'])}
    weights=np.array([len(h['_copies']) for h in data['haps']])
    counts={}
    for i in data['in_core']:
        s=data['sites'][i];valid=[s['ref'],s['alt']];aa=data['_matrix'][:,i]
        called=int(weights[np.isin(aa,valid)].sum())
        for a in valid: counts[i,a]=(called,int(weights[aa==a].sum()))
    # Use the original fourth-column label even when its VCF ID is CHR:POS:REF:ALT.
    auth=dict(auth,_input_lead_override=ident['input_lead_snp']);data['_auth']=auth
    _,iv,istatus=dl.input_variant(auth,data)
    ident['input_status']=istatus
    evidence={(e['hap_id'],e['diagnostic_lineage']):e for e in data['evidence']}
    trees={t['expected_lineage']:t for t in data['trees']}
    haps={h['hap_id']:h for h in data['haps']}
    comps={}
    for a in data['arch']:
        ix=[i for i in data['in_core'] if a['seq'][i] in dl.BASES]
        matrix=data['_matrix'][:,ix]
        ncalled=np.isin(matrix,list(dl.BASES)).sum(axis=1)
        matches=(matrix==np.array([a['seq'][i] for i in ix])).sum(axis=1)
        for hi,h in enumerate(data['haps']):
            nc=int(ncalled[hi]);nm=int(matches[hi])
            if not nc: continue
            e=evidence.get((h['hap_id'],a['lineage']),{})
            r=dict(hap_id=h['hap_id'],archaic=a['archaic'],lineage=a['lineage'],n_compared=nc,n_match=nm,
                prop_match=nm/nc,n_copies=len(h['_copies']),candidate_start=dl.integer(e.get('candidate_start')),
                candidate_end=dl.integer(e.get('candidate_end')),callable_fraction=nc/len(ix))
            k=(h['hap_id'],a['lineage'])
            if k not in comps or (r['prop_match'],r['n_compared'])>(comps[k]['prop_match'],comps[k]['n_compared']): comps[k]=r
    raw=[r for r in comps.values() if r['n_compared']>=10 and r['callable_fraction']>=.8]
    rawbest=max(raw,key=lambda r:(r['prop_match'],r['n_compared'],r['n_copies']),default=None)
    chosen={k for k,e in evidence.items() if dl.truth(e.get('diagnostic_candidate_pass'))}
    roles={k:'sequence_candidate' for k in chosen}
    if rawbest:
        k=(rawbest['hap_id'],rawbest['lineage']);chosen.add(k);roles[k]=roles.get(k,'')+';raw_similarity_best'
    # A common haplotype bearing the specified marker's archaic allele is a comparator.
    # It is never upgraded to a sequence candidate by this display-only rule.
    for lin in sorted({a['lineage'] for a in data['arch']}):
        eligible=[(k,e) for k,e in evidence.items() if k[1]==lin and dl.truth(e.get('named_anchor_matches_lineage')) and k in comps]
        if eligible:
            k,_=max(eligible,key=lambda ke:comps[ke[0]]['n_copies']);chosen.add(k);roles[k]=roles.get(k,'')+';input_marker_representative'
    details=[];copies=[];summaries=[];tree_rows=[]
    marker_map=Path(os.environ.get('PHYML_MARKER_MAP',Path(__file__).with_name('phyml_markers.tsv')))
    risk=next((r for r in dl.rows(marker_map) if
        r['genome_build']==ident['genome_build'] and r['marker_id']==ident['input_lead_snp']),{})
    risk_allele=risk.get('risk_allele')
    for k in sorted(chosen):
        if k not in comps: continue
        best=comps[k];e=evidence.get(k,{});t=trees.get(k[1],{});h=haps[k[0]];target=h['_copies']
        passed=dl.truth(e.get('diagnostic_candidate_pass'))
        in_tree=k[0] in t.get('candidate_tips_in_clade','').split(',')
        supported=passed and dl.truth(t.get('candidate_clade_pass')) and in_tree
        row=dict(ident,**best,role=roles[k].strip(';'),candidate_id=f"{ident['locus_key']}|{k[1]}|{k[0]}",
            diagnostic_candidate_pass=int(passed),call='supported_candidate' if supported else 'candidate_sequence_only' if passed else 'similarity_comparator_only',
            diagnostic_sites=dl.integer(e.get('n_lineage_diagnostic_sites')),diagnostic_matches=dl.integer(e.get('n_diagnostic_match')),
            diagnostic_match_prop=dl.num(e.get('diagnostic_match_prop')),ils_probability=dl.num(e.get('candidate_ils_probability')),
            tree_bootstrap=dl.num(t.get('candidate_clade_bootstrap')),tree_pass=int(supported),tree_purity=dl.num(t.get('candidate_purity')),
            tree_sensitivity=dl.num(t.get('candidate_sensitivity')),reason=t.get('tree_call_reason') if passed else e.get('diagnostic_candidate_reason'),
            risk_allele=risk_allele,trait=risk.get('trait'),n_individuals=len({cp.rsplit(':',1)[0] for cp in target}))
        if iv and risk_allele:
            rm=dl.count_variant(iv['_alleles'],target,target,risk_allele)
            row.update(n_risk_copies=rm['tp'],n_risk_callable=rm['n_callable_copies'],risk_fraction=rm['tp']/rm['n_callable_copies'] if rm['n_callable_copies'] else None)
        annotate_tags(row,data,best,target,counts,iv);details.append(row)
        if passed:
            for cp in sorted(target):
                sample,hap=cp.rsplit(':',1)
                copies.append(dict(candidate_id=row['candidate_id'],locus_key=ident['locus_key'],dataset_id=ident['dataset_id'],genome_build=ident['genome_build'],
                    locus_id=ident['locus_id'],chr=ident['chr'],lineage=k[1],hap_id=k[0],sample_id=sample,haplotype=hap,
                    candidate_start=best['candidate_start'],candidate_end=best['candidate_end']))
    for lin in sorted({a['lineage'] for a in data['arch']}):
        rr=[r for r in details if r['lineage']==lin and r['diagnostic_candidate_pass']]
        ee=[e for (hid,l),e in evidence.items() if l==lin]
        t=trees.get(lin,{})
        nd=max((dl.integer(e.get('n_lineage_diagnostic_sites'),0) for e in ee),default=0)
        call='supported_candidate' if any(r['tree_pass'] for r in rr) else 'candidate_sequence_only' if rr else 'no_candidate' if nd else 'not_evaluable'
        row=dict(ident,lineage=lin,call=call,diagnostic_sites=nd,n_candidate_haplotypes=len(rr),n_candidate_copies=sum(r['n_copies'] for r in rr),
            candidate_haplotypes=','.join(r['hap_id'] for r in rr),tree_bootstrap=dl.num(t.get('candidate_clade_bootstrap')),
            tree_pass=int(any(r['tree_pass'] for r in rr)),tree_purity=dl.num(t.get('candidate_purity')),tree_sensitivity=dl.num(t.get('candidate_sensitivity')),
            reason=t.get('tree_call_reason') if rr else 'no_candidate_passes_sequence_filters' if nd else 'no_lineage_diagnostic_sites',
            risk_allele=risk_allele,trait=risk.get('trait'),n_risk_copies=sum(r.get('n_risk_copies',0) for r in rr) if risk_allele and iv else None,
            n_risk_callable=sum(r.get('n_risk_callable',0) for r in rr) if risk_allele and iv else None,
            tree_file=Path(t.get('tree_file','')).name if t.get('tree_file') else None)
        if row['n_risk_callable']: row['risk_fraction']=row['n_risk_copies']/row['n_risk_callable']
        if rr:
            representative=max(rr,key=lambda r:(r['tree_pass'],r.get('diagnostic_match_prop') or 0,r['prop_match'],r['n_copies']))
            target=set().union(*(haps[r['hap_id']]['_copies'] for r in rr))
            annotate_tags(row,data,representative,target,counts,iv)
            row.update(best_haplotype=representative['hap_id'],candidate_start=min(r['candidate_start'] for r in rr),candidate_end=max(r['candidate_end'] for r in rr))
        summaries.append(row)
        if t:
            tree_rows.append(dict(ident,lineage=lin,**{k:v for k,v in t.items() if k not in ident and k not in ['tree_file','plot_file','stats_file','phy_file']}))
    summaries.sort(key=lambda r:({'supported_candidate':0,'candidate_sequence_only':1,'no_candidate':2,'not_evaluable':3}.get(r['call'],4),r['lineage']))
    return summaries,details,copies,tree_rows


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
            in_panel=('_tested_copies' not in run or c['sample_id']+':'+str(c['haplotype']) in run['_tested_copies'])
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


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--normalize',required=True,type=Path);ap.add_argument('--database',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path);ap.add_argument('--loci',type=Path)
    args=ap.parse_args()
    if list((args.normalize/'phyml').glob('**/final/evidence_loci.tsv')):
        auth,files=read_authoritative(args.normalize.resolve())
    else:
        auth,files=[],[]
    if args.loci: auth=dl.apply_loci_file(auth,args.loci)
    auth.sort(key=lambda a:(a['dataset_id'],a['genome_build'],dl.integer(dl.chrom(a['chr']),99),int(a['core_start'])))
    summary=[];details=[];copies=[];trees=[]
    for a in auth:
        ss,dd,cc,tt=locus_rows(a);summary+=ss;details+=dd;copies+=cc;trees+=tt
        files+=a.get('_report_source_files',[])
        print(f"PhyML report: {a['locus_id']}, {len(dd)} haplotype rows",flush=True)
    with sqlite3.connect(args.database.resolve().as_uri()+'?mode=ro',uri=True) as con: validation,runs=validate(copies,con)
    attach_validation(summary,validation,runs);attach_validation(details,validation,runs,True)
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    for name,rr in [('phyml_locus_report',summary),('phyml_haplotype_report',details),('phyml_copy_validation',validation),('phyml_lineage_trees',trees)]:
        write_tsv(out/(name+'.tsv'),rr)
    atomic_text(out/'phyml_report_manifest.json',json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),n_loci=len(auth),n_lineage_rows=len(summary),n_haplotype_rows=len(details),
        coverage_threshold=.8,copy_phase='IBDmix unphased; TRACE stored haplotype index',
        validation_denominator='candidate individuals/copies, not callable genomes',
        tag_scope='all stored biallelic core SNPs matching the representative archaic sequence; in-sample r2',
        input_files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(set(files))]),indent=2))
    print(f"PhyML report ready: {len(auth)} loci, {len(summary)} lineages, {len(details)} haplotypes",flush=True)


if __name__=='__main__':main()
