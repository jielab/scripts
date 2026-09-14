#!/usr/bin/env python3
"""GWAS-anchored LD-core phylogeny following Zeberg & Pääbo (2020).

EUR defines LD. All sampled populations supply recurrent modern haplotypes.
No archaic-similarity search, replacement index SNP, or candidate tip sampling.
"""
import argparse, hashlib, json, math, os, subprocess, sys, fcntl
from collections import Counter, defaultdict
from pathlib import Path
from phyml_gwas_input import read, write
from phyml_core import (BASES, SkipLocus, modern_data, query_rows, vcf_path, vcf_contig,
                        archaic_vcf, archaic_calls, called_base, called_haploid_base, read_sexes)
from phyml_tree_summary import parse_newick, bootstrap_values
from phyml_thresholds import ils_probability

REFS = ('Altai', 'Chagyr', 'Vindija')


def risk_from_effect(ref, alt, effect, beta):
    if effect not in (ref, alt) or ',' in alt or ref == alt:
        raise SkipLocus('effect allele does not match exact target variant', 'lead_allele_mismatch')
    return effect if beta > 0 else (alt if effect == ref else ref)


def exact_lead(row, records, n_samples, haploid=False):
    found=[]
    for f in records:
        if int(f[1]) != int(row['lead_pos']): continue
        ref, alt = f[3].upper(), f[4].upper()
        if row['ref']:
            if (ref,alt) != (row['ref'],row['alt']): continue
        elif row['index_snp'] not in f[2].split(';'): continue
        if ',' in alt or len(f) != 5+n_samples: continue
        found.append(f)
    if len(found) != 1:
        raise SkipLocus('exact original lead absent or ambiguous; no replacement SNP', 'lead_absent_or_ambiguous')
    f=found[0]; risk=risk_from_effect(f[3],f[4],row['effect_allele'],float(row['beta_j']))
    copies=[]
    for gt in f[5:]:
        alleles=(called_haploid_base(f[3],f[4],gt),'N') if haploid else called_base(f[3],f[4],gt,True)
        copies.extend(1 if x==risk else 0 if x in (f[3],f[4]) else -1 for x in alleles)
    return dict(pos=int(f[1]),ref=f[3],alt=f[4],vcf_id=f[2],risk_allele=risk,copies=copies)


def phased_ld(risk, site, indexes):
    n=sx=sy=sxy=0
    for i in indexes:
        a,b=risk[i],site['haps'][i]
        if a < 0 or b not in (site['ref'],site['alt']): continue
        y=int(b==site['alt']); n+=1; sx+=a; sy+=y; sxy+=a*y
    denom=sx*(n-sx)*sy*(n-sy)
    if n < 4 or denom==0: return None,n,''
    cov=n*sxy-sx*sy
    return min(1.,cov*cov/denom),n,site['alt'] if cov>0 else site['ref']


def define_core(sites, lead, indexes):
    ld=[]; high=[]
    for s in sites:
        r2,n,allele=phased_ld(lead['copies'],s,indexes)
        passed=r2 is not None and r2 > .98
        ld.append(dict(pos=s['pos'],id=s['vid'],ref=s['ref'],alt=s['alt'],ld_r2=r2,n_EUR_copies=n,
                       risk_linked_allele=allele,core_marker=int(passed)))
        if passed: high.append(s)
    if len(high)<2: raise SkipLocus('fewer than two SNPs with EUR r² > 0.98', 'insufficient_high_LD_markers',dict(ld=ld))
    start=min([lead['pos']]+[s['pos'] for s in high])-1
    end=max([lead['pos']+len(lead['ref'])-1]+[s['pos'] for s in high])
    return [s for s in sites if start<s['pos']<=end],ld,start,end


def risk_clade(newick, risk_tips, modern_tips):
    """Test the prespecified all-risk + three-Neanderthal split, either orientation."""
    root=parse_newick(newick); edges=[]
    def visit(node):
        tips={node.label} if not node.children else set().union(*(visit(c) for c in node.children))
        edges.append((node,tips)); return tips
    all_tips=visit(root)
    required=set(risk_tips)|set(REFS)
    if not risk_tips or not required<=all_tips or 'Ancestral' not in all_tips or not (set(modern_tips)-set(risk_tips)):
        return None
    matches=[]
    for node,tips in edges:
        if node is root: continue
        for side in (tips,all_tips-tips):
            if side==required:
                matches.append(dict(bootstrap=node.support, tips=','.join(sorted(side)), node=node.node_id))
    return max(matches,key=lambda x:-1 if x['bootstrap'] is None else x['bootstrap'],default=None)


def prepare_sequences(sites, calls, samples, lead, haploid):
    # Positions must be callable in all three Neanderthals. Missing ancestry
    # remains N (never substitute REF); constant columns are retained after
    # singleton removal, as required for HKY+I parameter estimation.
    sites=[s for s in sites if all(calls[r].get(s['pos'],'N') in BASES for r in REFS)]
    indexes=[2*i+h for i in range(len(samples)) for h in range(1 if haploid else 2)]
    grouped=defaultdict(list)
    for i in indexes:
        if lead['copies'][i] < 0: continue
        seq=''.join(s['haps'][i] for s in sites)
        grouped[seq].append(i)
    recurrent=[]
    for seq,ix in sorted(grouped.items(),key=lambda x:(-len(x[1]),x[0])):
        if len(ix)<2: continue
        flags={lead['copies'][i] for i in ix}
        role='risk' if flags=={1} else 'nonrisk' if flags=={0} else 'mixed'
        recurrent.append(dict(hap_id=f'H{len(recurrent)+1:05d}',seq=seq,indices=ix,role=role,n=len(ix)))
    return sites,recurrent,grouped


def run_verified_tree(command):
    result=subprocess.run(command)
    if result.returncode==0:
        result=subprocess.run(command+['--verify-only'])
    return result


def run_locus(a, row):
    out=a.out; final=out/'final'; loc=out/'loci'; final.mkdir(parents=True,exist_ok=True); loc.mkdir(parents=True,exist_ok=True)
    row=dict(row); lid=row['locus_id']; ch=row['chr']; pos=int(row['lead_pos'])
    # A rerun that becomes unassessable must not display earlier sequences.
    for name in ('sites.tsv','archaic.tsv','ancestral.tsv','haplotypes.tsv','ld.tsv'):
        (loc/name).unlink(missing_ok=True)
    haploid=ch=='X' and a.x_male_only
    ident=dict(locus_id=lid,index_snp=row['index_snp'],source_build=row['source_build'],source_pos=row['source_pos'],
        dataset_id=a.dataset,genome_build='GRCh37',chr=ch,lead_pos=pos,p_j=row['p_j'],beta_j=row['beta_j'],
        effect_allele=row['effect_allele'],ld_population='1KG EUR',ld_r2_threshold=.98,selection_method='gwas_lead_ld_core',
        analysis_start=int(row['search_start']),analysis_end=int(row['search_end']),core_start=pos-1,core_end=pos,
        lineage='Neanderthal',tree_pass=0,tree_bootstrap=None,status='not_evaluable',reason='',
        risk_allele='',n_ld_sites=None,n_sites=None,n_risk_copies=None,n_candidate_haplotypes=None,n_candidate_copies=None)
    locus=dict(chrom=ch,name=lid,start=int(row['search_start']),end=int(row['search_end']),core_start=pos-1,core_end=pos)
    haps=[]; detail=[]; copyrows=[]; allcopies=[]; tree=dict(locus_id=lid,tree_status='not_run',candidate_lineage='Neanderthal',
        expected_lineage='Neanderthal',candidate_clade_pass=0,candidate_clade_bootstrap=None,tree_call_reason='not_run',
        candidate_clade_tips='',candidate_tips_in_clade='',expected_archaic_tips_in_clade='',control_tips_in_clade='',candidate_context_tips_in_clade='',
        n_candidate_tips_in_clade=0,n_expected_archaic_tips_in_clade=0,candidate_clade_n_tips=0,
        candidate_clade_modern_tips=0,candidate_clade_archaic_tips=0,candidate_clade_specificity=None,
        tree_scope='gwas_risk_core',tree_newick='',tree_file='',stats_file='',plot_file='',tree_has_ancestral_outgroup=0)
    failed=False
    try:
        vcf=vcf_path(a.vcf_dir,ch); contig=vcf_contig(vcf,ch)
        panel=read(a.sample_file); meta={r['sample']:r for r in panel}
        samples, records=query_rows(vcf,f'{contig}:{pos}-{pos}')
        if any(s not in meta for s in samples): raise ValueError('sample panel must include every VCF sample and population')
        eur=[i for i,s in enumerate(samples) if meta[s].get('super_pop','')=='EUR' or meta[s].get('pop','') in {'CEU','GBR','FIN','IBS','TSI'}]
        if not eur: raise ValueError('No EUR samples in --sample-panel; LD cannot use all populations as a fallback')
        lead=exact_lead(row,records,len(samples),haploid)
        ident.update(risk_allele=lead['risk_allele'],target_variant=f"{ch}:{pos}:{lead['ref']}:{lead['alt']}",
                     n_EUR_individuals=len(eur),n_risk_copies=sum(x==1 for x in lead['copies']))
        indexes=[2*i+h for i in eur for h in range(1 if haploid else 2)]
        eur_called=[lead['copies'][i] for i in indexes if lead['copies'][i]>=0]
        ident['risk_frequency_EUR']=sum(eur_called)/len(eur_called) if eur_called else None
        _,sites,_=modern_data(vcf,contig,locus,2,read_sexes(a.sample_file),a.x_male_only,a.x_par_diploid,ancestral_any_base=True)
        # Multiple records at one coordinate cannot be matched to an archaic
        # base unambiguously; exclude them rather than merge unrelated alleles.
        counts=Counter(s['pos'] for s in sites); sites=[s for s in sites if counts[s['pos']]==1]
        core,ld,start,end=define_core(sites,lead,indexes)
        write(loc/'ld.tsv',ld)
        ident.update(core_start=start,core_end=end,core_kb=(end-start)/1000,n_ld_sites=sum(r['core_marker'] for r in ld),
                     n_search_sites=len(sites),anchor_pos=pos,selected_start=start,selected_end=end)
        high=[r for r in ld if r['core_marker']]
        # Search limits are explicit: the inferred span is conditional on this
        # finite window, and a marker near an edge requests a wider search.
        ident['search_edge_warning']=int(start-locus['start']<10000 or locus['end']-end<10000)
        calls={}
        for ref in REFS:
            av=archaic_vcf(a.archaic_root,ref,ch)
            calls[ref]=archaic_calls(av,vcf_contig(av,ch),dict(locus,start=start,end=end),core,allow_third_allele=True)
        markers={r['pos']:r['risk_linked_allele'] for r in high}
        for ref in REFS:
            called=[(p,b) for p,b in markers.items() if calls[ref].get(p,'N') in BASES]
            ident[ref+'_LD_matches']=sum(calls[ref][p]==b for p,b in called)
            ident[ref+'_LD_called']=len(called)
        sites,haps,grouped=prepare_sequences(core,calls,samples,lead,haploid)
        ident.update(n_sites=len(sites),n_compared=len(sites),n_ancestral_sites=sum(s['ancestral'] in BASES for s in sites),
                     n_tree_haplotypes=len(haps),n_singleton_copies=sum(len(ix) for ix in grouped.values() if len(ix)==1),
                     n_candidate_haplotypes=sum(h['role']=='risk' for h in haps),
                     n_candidate_copies=sum(h['n'] for h in haps if h['role']=='risk'),
                     n_nonrisk_haplotypes=sum(h['role']=='nonrisk' for h in haps))
        if len(sites)<2: raise SkipLocus('insufficient Neanderthal-callable core sites', 'insufficient_tree_sites')
        arch={r:''.join(calls[r][s['pos']] for s in sites) for r in REFS}
        ancestor=''.join(s['ancestral'] for s in sites)
        write(loc/'sites.tsv',[dict(chr=ch,pos=s['pos'],id=s['vid'],ref=s['ref'],alt=s['alt']) for s in sites])
        write(loc/'archaic.tsv',[dict(archaic=r,lineage='Neanderthal',seq=arch[r]) for r in REFS])
        write(loc/'ancestral.tsv',[dict(reference='Ancestral',n_callable=ident['n_ancestral_sites'],seq=ancestor)])
        for h in haps:
            h['copies']=';'.join(f'{samples[i//2]}:{i%2+1}' for i in h['indices'])
            comp=[]
            for ref,seq in arch.items():
                pairs=[(x,y) for x,y in zip(h['seq'],seq) if x in BASES and y in BASES]
                nm=sum(x==y for x,y in pairs); nc=len(pairs)
                comp.append((nm/nc if nc else 0,nc,nm,ref))
            prop,nc,nm,ref=max(comp)
            d=dict(ident,hap_id=h['hap_id'],role=h['role'],n_copies=h['n'],n_individuals=len({i//2 for i in h['indices']}),
                   n_compared=nc,n_match=nm,archaic=ref,prop_match=prop,candidate_start=start,candidate_end=end,
                   call='risk_haplotype' if h['role']=='risk' else 'nonrisk_control' if h['role']=='nonrisk' else 'risk_nonrisk_sequence_unresolved',tree_pass=0)
            d['superpopulation_copy_counts']=','.join(f'{p}:{n}' for p,n in sorted(Counter(meta[samples[i//2]].get('super_pop','UNKNOWN') for i in h['indices']).items()))
            detail.append(d)
            for i in h['indices']:
                cp=dict(locus_id=lid,hap_id=h['hap_id'],sample=samples[i//2],sample_id=samples[i//2],haplotype=i%2+1,
                        candidate_start=start,candidate_end=end,role=h['role'])
                allcopies.append(cp)
                if h['role']=='risk': copyrows.append(cp)
            h.update(locus_id=lid,genome_build='GRCh37',best_archaic=ref,best_lineage='Neanderthal',n_compared=nc,n_match=nm,prop_match=prop,direct_match_pass=0)
        write(loc/'haplotypes.tsv',[{k:v for k,v in h.items() if k!='indices'} for h in haps])
        if any(h['role']=='mixed' for h in haps):
            ident.update(n_candidate_haplotypes=None,n_candidate_copies=None)
            raise SkipLocus('lead alleles share identical callable core sequences', 'risk_nonrisk_sequence_unresolved')
        if not ident['n_candidate_haplotypes'] or not ident['n_nonrisk_haplotypes']:
            raise SkipLocus('both recurrent risk and nonrisk haplotypes required', 'insufficient_recurrent_haplotypes')
        if not ident['n_ancestral_sites']: raise SkipLocus('ancestral human bases absent; REF is not an outgroup', 'ancestral_sequence_unavailable')
        phy=loc/'haplotypes.phy'
        tree['phy_file']=str(phy)
        seqs=[(h['hap_id'],h['seq']) for h in haps]+list(arch.items())+[('Ancestral',ancestor)]
        phytext=f'{len(seqs)} {len(sites)}\n'+''.join(f'{label:<10} {seq}\n' for label,seq in seqs)
        if not phy.exists() or phy.read_text()!=phytext:
            from phyml_run import clean
            clean(phy)
            phy.write_text(phytext)
        write(loc/'haplotypes.phy.meta.tsv',[dict(phy_label=label,label=label,role=next((h['role'] for h in haps if h['hap_id']==label),'ancestral' if label=='Ancestral' else 'archaic')) for label,seq in seqs])
        ident.update(status='tree_not_requested',reason='sequence_prepared')
        if a.action=='run' and a.plot_phy=='TRUE':
            print(f"[GU PHYML] {lid}: EUR={len(eur)}; core={ch}:{start+1}-{end}; LD markers={ident['n_ld_sites']}; tree sites={len(sites)}; recurrent haplotypes={len(haps)}; bootstrap=100",flush=True)
            cmd=[sys.executable,str(Path(__file__).with_name('phyml_run.py')),'--phy',str(phy),'--scope','gwas_risk_core',
                 '--bootstrap','100','--timeout',str(a.timeout),'--cpus',str(a.cpus),'--mpi-fallback','serial']
            if os.environ.get('PHYML_REPLACE')=='TRUE':
                # Removing the completion receipt makes the runner recompute,
                # while its own lock continues to protect tree output files.
                Path(str(phy)+'.phyml.complete.json').unlink(missing_ok=True)
            proc=run_verified_tree(cmd)
            treepath=Path(str(phy)+'_phyml_tree.txt'); statpath=Path(str(phy)+'_phyml_stats.txt')
            if proc.returncode:
                failed=True
                logpath=Path(str(phy)+'.phyml.log')
                failure_log=logpath.read_text(errors='replace') if logpath.exists() else ''
                reason=('tree_timeout' if proc.returncode==124 else
                        'numerical_model_fit_failed' if 'Cannot work out eigen vectors' in failure_log else
                        'see_PhyML_run_log')
                ident.update(status='tree_failed',reason=reason)
                tree.update(tree_status='failed',tree_call_reason=reason)
            else:
                newick=treepath.read_text().strip(); risk={h['hap_id'] for h in haps if h['role']=='risk'}
                match=risk_clade(newick,risk,{h['hap_id'] for h in haps})
                bs=match['bootstrap'] if match else None; passed=bool(bs is not None and bs>=70)
                reason='risk_Neanderthal_split' if match else 'no_exclusive_risk_Neanderthal_split'
                ident.update(status='tree_supported' if passed else 'tree_not_supported',reason=reason,tree_pass=int(passed),tree_bootstrap=bs)
                tree.update(tree_status='complete',tree_newick=newick,tree_file=str(treepath),stats_file=str(statpath),
                    tree_has_ancestral_outgroup=1,candidate_clade_pass=int(passed),candidate_clade_bootstrap=bs,
                    candidate_clade_tips=match['tips'] if match else '',candidate_tips_in_clade=','.join(sorted(risk)) if match else '',
                    candidate_clade_rule='all_recurrent_risk_plus_three_Neanderthals_no_nonrisk_or_ancestor',tree_call_reason=reason,
                    n_bootstrap_nodes=len(bootstrap_values(newick)),
                    expected_archaic_tips_in_clade=','.join(REFS) if match else '',
                    n_candidate_tips_in_clade=len(risk) if match else 0,n_expected_archaic_tips_in_clade=3 if match else 0,
                    candidate_clade_n_tips=len(risk)+3 if match else 0,candidate_clade_modern_tips=len(risk) if match else 0,
                    candidate_clade_archaic_tips=3 if match else 0,candidate_clade_specificity=1 if match else None)
        tree['tree_call_reason']=ident['reason']
    except SkipLocus as e:
        ident.update(status=e.code,reason=str(e)); tree['tree_call_reason']=str(e)
        if e.code=='risk_nonrisk_sequence_unresolved': copyrows=[]
        if e.details and 'ld' in e.details:
            write(loc/'ld.tsv',e.details['ld'])
            ident.update(n_ld_sites=sum(r['core_marker'] for r in e.details['ld']),n_search_sites=len(e.details['ld']))
        print(f'[GU PHYML] SKIP {lid}: {e.code}: {e}',flush=True)
    ident['call']=ident['status']
    # A length-model sensitivity statistic, not a calibrated locus-specific P.
    if ident.get('core_kb'):
        ident['ils_probability']=ils_probability(ident['core_end']-ident['core_start'])
        ident['ils_model']='assumed_0.53cM/Mb_29yr_550k_split_50k_archaic_age;not_local_map;uncorrected'
    for d in detail:
        d.update(tree_bootstrap=ident['tree_bootstrap'],tree_pass=int(ident['tree_pass'] and d['role']=='risk'),
                 call=ident['status'] if d['role']=='risk' else 'nonrisk_control' if d['role']=='nonrisk' else 'risk_nonrisk_sequence_unresolved')
    write(final/'gwas_loci.tsv',[ident]); write(final/'gwas_haplotypes.tsv',detail)
    write(final/'gwas_copies.tsv',copyrows)
    write(final/'loci.tsv',[ident]); write(final/'trees.tsv',[tree]); write(final/'evidence_trees.tsv',[tree])
    completed_haps=[{k:v for k,v in h.items() if k!='indices'} for h in haps if 'copies' in h]
    write(final/'haplotypes.tsv',completed_haps, list(completed_haps[0]) if completed_haps else ['locus_id','hap_id','copies'])
    write(final/'haplotype_samples.tsv',allcopies,['locus_id','hap_id','sample','sample_id','haplotype','candidate_start','candidate_end','role'])
    write(final/'skipped_loci.tsv',[] if ident['status'].startswith('tree_') else [ident])
    write(final/'gwas_lead.tsv',[row])
    parameters=dict(workflow='gwas_lead_ld_core_v1',ld_population='1KG EUR',ld_rule='phased_r2 > 0.98',
        lead=row,tree_populations='all target 1KG samples',minimum_haplotype_copies=2,minimum_minor_allele_copies=2,
        ancestral_source='target VCF INFO/AA; unknown remains N; not verified as Ensembl release 100',
        references=list(REFS),heterozygous_archaic_policy='mask_as_N',bootstrap=100,model='HKY85+G4+I;estimated',
        tree_rule='all recurrent risk + three Neanderthals, no other modern tips or ancestor; reporting BS >=70',
        method_source='https://www.nature.com/articles/s41586-020-2818-3')
    (final/'gwas_parameters.json').write_text(json.dumps(parameters,indent=2)+'\n')
    if tree['tree_status']=='complete':
        plot=subprocess.run(['Rscript','--vanilla',str(Path(__file__).with_name('phyml_panel_b.R')),'--out',str(out)])
        if plot.returncode:
            print('[GU PHYML] WARNING: tree completed but plot export failed',flush=True)
            failed=True
    print(f"[GU PHYML] {lid}: {ident['status']}; {ident['reason']}",flush=True)
    return 1 if failed else 0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('lead-table','loci','vcf-dir','archaic-root','sample-file','out'): p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--dataset',default='1kg'); p.add_argument('--action',choices=['run','match','check'],default='run')
    p.add_argument('--plot-phy',choices=['TRUE','FALSE'],default='TRUE'); p.add_argument('--timeout',type=int,default=86400);p.add_argument('--cpus',type=int,default=4)
    p.add_argument('--x-male-only',action='store_true');p.add_argument('--x-par-diploid',action='store_true')
    a=p.parse_args()
    if a.dataset!='1kg': raise ValueError('This EUR LD workflow currently requires target 1kg; another target needs a separate fixed 1KG EUR LD reference')
    names=[line.split()[3] for line in a.loci.read_text().splitlines() if line.strip() and not line.startswith('#')]
    rows=[r for r in read(a.lead_table) if r['locus_id'] in names]
    if len(rows)!=1 or len(names)!=1: raise ValueError('Each worker requires exactly one original GWAS lead')
    if a.action=='check':
        print('[GU PHYML] original lead, EUR LD and recurrent-haplotype workflow configured'); return
    a.out.mkdir(parents=True,exist_ok=True)
    with (a.out/'.gwas.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        raise SystemExit(run_locus(a,rows[0]))

if __name__=='__main__': main()
