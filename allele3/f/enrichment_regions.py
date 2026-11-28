#!/usr/bin/env python3
"""Coordinate-resolved candidates from the main LOCO streaming pass.

Every finite O/E >= threshold is exported, including low-count candidates.
Conditional Poisson inference is exploratory: LD and fitted-background
uncertainty are not modeled. No case is counted without its expected rate.
"""
from __future__ import annotations
import bisect, gzip, heapq, html, itertools, json, math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from output_layout import data_path
import numpy as np
from research_common import (Intervals,atomic_json,bed_rows,bh,chrom_name,chrom_key,
    chromosomes,identity,number,write_tsv)

LENGTHS=dict(zip(['chr'+str(i) for i in range(1,23)]+['chrX','chrY','chrM'],[
    248956422,242193529,198295559,190214555,181538259,170805979,159345973,
    145138636,138394717,133797422,135086622,133275309,114364328,107043718,
    101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,
    156040895,57227415,16569]))
FIELDS=['region_id','region_kind','chrom','start0','end0','window_bp','annotation',
    'outcome','population','sensitivity','span_bp','region_bp','annotation_fraction',
    'observed','expected','oe','ci_low','ci_high','p','fdr','fdr_family','n_family_tests',
    'n_eligible_supported','n_unsupported','observed_total','support_fraction',
    'prior_only_sites','prior_only_fraction','passes_fold','ci_low_ge_fold',
    'support_status','mask_fraction','technical_status','stable_fraction','inference']


def ge(value,cutoff):
    """Inclusive boundary, allowing only accumulated binary64 rounding noise.

    The relative tolerance is far below the exported display precision; it
    does not relax a scientific fold or count-support threshold.
    """
    return value>=cutoff or math.isclose(value,cutoff,rel_tol=1e-12,abs_tol=0.)


def window_sizes(value):
    values=[int(x) for x in str(value).split(',') if x.strip()]
    if not values or min(values)<1 or len(values)!=len(set(values)):
        raise ValueError('--region-windows requires unique positive sizes')
    return sorted(values)


class Coverage:
    """Merged union with logarithmic-time covered-length queries."""
    def __init__(self,rows=()):
        self.by_chrom={}
        for (chrom,_),(starts,ends) in Intervals((c,s,e,'union') for c,s,e,*_ in rows).data.items():
            cumulative=[0]
            for s,e in zip(starts,ends):cumulative.append(cumulative[-1]+e-s)
            self.by_chrom[chrom]=(starts,ends,cumulative)
    def integral(self,chrom,x):
        if chrom not in self.by_chrom:return 0
        starts,ends,total=self.by_chrom[chrom];i=bisect.bisect_right(starts,x)-1
        return 0 if i<0 else total[i]+max(0,min(x,ends[i])-starts[i])
    def bp(self,chrom,start,end):return self.integral(chrom,end)-self.integral(chrom,start)
    def pieces(self,chrom,start,end):
        if chrom not in self.by_chrom:return []
        starts,ends,_=self.by_chrom[chrom];i=bisect.bisect_right(ends,start);answer=[]
        while i<len(starts) and starts[i]<end:
            answer.append((max(start,starts[i]),min(end,ends[i])));i+=1
        return answer


@dataclass(frozen=True,slots=True)
class Region:
    region_id:str
    region_kind:str
    chrom:str
    start0:int
    end0:int
    window_bp:int
    annotation:str
    segments:tuple=()


class LocusSweep:
    """Sorted-position sweep, with cost proportional to overlapping loci."""
    def __init__(self,regions):
        self.segments=defaultdict(list)
        for rid,r in regions.items():
            for start,end in r.segments:self.segments[r.chrom].append((start,end,rid))
        for c in self.segments:self.segments[c].sort()
        self.chrom=None;self.pos=-1;self.cursor=0;self.active=Counter();self.ends=[];self.last=()
    def lookup(self,chrom,pos):
        if chrom!=self.chrom:
            self.chrom=chrom;self.pos=-1;self.cursor=0;self.active=Counter();self.ends=[]
        if pos<self.pos:raise ValueError('Region enrichment requires coordinate-sorted sites')
        if pos==self.pos:return self.last
        records=self.segments.get(chrom,[])
        while self.cursor<len(records) and records[self.cursor][0]<=pos:
            start,end,rid=records[self.cursor];self.cursor+=1
            if end>pos:self.active[rid]+=1;heapq.heappush(self.ends,(end,rid))
        while self.ends and self.ends[0][0]<=pos:
            end,rid=heapq.heappop(self.ends);self.active[rid]-=1
            if not self.active[rid]:del self.active[rid]
        self.pos=pos;self.last=tuple(self.active);return self.last


def catalogue(a,reference_meta):
    if a.ref=='herv':path=Path(reference_meta['path'])
    elif a.ref=='archaic':path=Path(reference_meta['output'])
    else:path=Path(a.avi_dir)/'results/hotspots.bed'
    groups=defaultdict(list);allowed=set(chromosomes(a.chr))
    raw=bed_rows(path);families=bed_rows(path,family=True) if a.ref=='herv' else None
    for c,start,end,label in raw:
        c=chrom_name(c);family=next(families)[3] if families is not None else label
        if c not in allowed:continue
        if a.ref=='alphagenome' and not label.endswith('_w'+str(a.window_bp)):continue
        name=label.split('|',1)[1] if a.ref=='herv' and '|' in label else f'{c}:{start}-{end}'
        groups[(c,family,name)].append((start,end))
    answer={}
    for (c,annotation,name),segments in groups.items():
        merged=list(Intervals((c,s,e,'x') for s,e in segments).rows())
        pieces=tuple((s,e) for _,s,e,_ in merged);start=pieces[0][0];end=pieces[-1][1]
        rid=f'locus|{annotation}|{name}|{c}:{start}-{end}'
        answer[rid]=Region(rid,'locus',c,start,end,0,annotation,pieces)
    return answer,path


class RegionHighlighter:
    def __init__(self,a,out,reference_meta,landscape=None):
        self.a=a;self.out=Path(out);self.sizes=window_sizes(a.region_windows)
        if not math.isfinite(a.highlight_fold) or a.highlight_fold<1:raise ValueError('--highlight-fold must be finite and >=1')
        self.fold=a.highlight_fold;self.reference_meta=reference_meta
        self.regions,self.catalogue_path=catalogue(a,reference_meta)
        self.sweep=LocusSweep(self.regions)
        self.coverage=Coverage((r.chrom,s,e) for r in self.regions.values() for s,e in r.segments)
        self.mask=Coverage(itertools.chain.from_iterable(bed_rows(p) for p in a.exclude_bed))
        stable=Path(a.stable_herv_bed) if a.stable_herv_bed else Path(a.herv_scan_dir)/'results/stable_loci.bed'
        if a.stable_herv_bed and not stable.is_file():raise ValueError('Missing --stable-herv-bed: '+str(stable))
        self.stable_path=stable if a.ref=='herv' and stable.is_file() else None
        self.stable=Coverage(bed_rows(stable)) if self.stable_path else Coverage()
        self.count_modes=set();self.counts={}
        for rid in self.regions:
            for sensitivity in (['all','exclude_mask'] if a.exclude_bed else ['all']):
                self.counts[(rid,sensitivity)]=[0,0.,0,0,0,0]

    def blocks(self,r):
        if r.region_kind=='locus':return r.segments
        if r.region_kind=='reference_window':return self.coverage.pieces(r.chrom,r.start0,r.end0)
        return [(r.start0,r.end0)]

    def record(self,chrom,pos1,case,rate,sensitivity='all',prior_only=False,count_mode=None):
        chrom=chrom_name(chrom);pos=int(pos1)-1
        if pos<0 or (chrom in LENGTHS and pos>=LENGTHS[chrom]):raise ValueError('Invalid GRCh38 site coordinate')
        if rate is not None and (not math.isfinite(rate) or not 0<=rate<=1):raise ValueError('Expected site probability must lie in [0,1]')
        if count_mode is not None:self.count_modes.add(count_mode)
        hits=self.sweep.lookup(chrom,pos);identifiers=list(hits)
        for width in self.sizes:
            start=pos//width*width;end=min(start+width,LENGTHS.get(chrom,start+width))
            kinds=[('window','genome')]
            if hits:kinds.append(('reference_window',self.a.ref+'_any'))
            for kind,annotation in kinds:
                rid=f'{kind}|{annotation}|w{width}|{chrom}:{start}-{end}'
                if rid not in self.regions:self.regions[rid]=Region(rid,kind,chrom,start,end,width,annotation)
                identifiers.append(rid)
        for rid in identifiers:
            values=self.counts.setdefault((rid,sensitivity),[0,0.,0,0,0,0]);values[4]+=int(case)
            if rate is None:values[3]+=1
            else:values[0]+=int(case);values[1]+=rate;values[2]+=1;values[5]+=int(prior_only)

    def scientific_settings(self):
        a=self.a;effective=[]
        for mode in sorted(self.count_modes):
            if a.threshold_kind=='AF' or (a.threshold_kind=='auto' and mode=='AF_only'):rule=dict(kind='AF',min_af=a.min_af)
            else:rule=dict(kind='AC',min_ac=a.min_ac,min_an=a.min_an)
            if rule not in effective:effective.append(rule)
        masks=[identity(p) for p in a.exclude_bed]
        context=dict(mode='unadjusted') if a.no_context else dict(mode='trinucleotide',fasta=identity(a.fasta))
        ref_inputs=([identity(self.catalogue_path)] if a.ref=='herv' else
            self.reference_meta.get('files',[]) if a.ref=='archaic' else [identity(data_path(Path(a.avi_dir), 'qc')/'manifest.json')])
        definition=dict(ref=a.ref,inputs=ref_inputs,region_windows=self.sizes,
            avi_window_bp=a.window_bp if a.ref=='alphagenome' else None,
            archaic_min_carriers=a.archaic_min_carriers if a.ref=='archaic' else None,
            archaic_population=a.archaic_population if a.ref=='archaic' else None,
            archaic_reference=a.archaic_reference if a.ref=='archaic' else None,
            archaic_min_slod=a.archaic_min_slod if a.ref=='archaic' else None,
            archaic_min_bp=a.archaic_min_bp if a.ref=='archaic' else None)
        return dict(effective_thresholds=effective,count_modes=sorted(self.count_modes),mask_inputs=masks,
            context=context,ref_definition=definition,population=a.population,
            min_background_sites=a.min_background_sites,observation_rule='Count-bearing sources require positive filtered AC; AF-only sources require positive AF; complete declared threshold evidence')

    def finalize(self,category_rows=None):
        from scipy.stats import poisson,chi2
        rows=[]
        for (rid,sensitivity),v in self.counts.items():
            r=self.regions[rid];observed,expected,supported,unsupported,total,prior=v
            blocks=self.blocks(r);bp=sum(e-s for s,e in blocks);span=r.end0-r.start0
            masked=sum(self.mask.bp(r.chrom,s,e) for s,e in blocks)/bp if bp else None
            stable=sum(self.stable.bp(r.chrom,s,e) for s,e in blocks)/bp if bp and self.stable_path else None
            fraction=supported/(supported+unsupported) if supported+unsupported else None
            prior_fraction=prior/supported if supported else None
            valid=supported>0 and expected>0;oe=observed/expected if valid else None;warnings=[]
            if observed<10:warnings.append('low_observed')
            if not ge(expected,1.):warnings.append('low_expected')
            if fraction is None or fraction<.8:warnings.append('limited_background_support')
            if prior_fraction is not None and prior_fraction>.5:warnings.append('prior_dominated_background')
            family='|'.join([r.region_kind,str(r.window_bp),r.annotation,'multi_ALT',self.a.population,sensitivity])
            rows.append(dict(region_id=rid,region_kind=r.region_kind,chrom=r.chrom,start0=r.start0,end0=r.end0,
                window_bp=r.window_bp,annotation=r.annotation,outcome='multi_ALT',population=self.a.population,
                sensitivity=sensitivity,span_bp=span,region_bp=bp,annotation_fraction=self.coverage.bp(r.chrom,r.start0,r.end0)/span,
                observed=observed,expected=expected,oe=oe,ci_low=None,ci_high=None,p=None,fdr=None,
                fdr_family=family,n_family_tests=0,n_eligible_supported=supported,n_unsupported=unsupported,
                observed_total=total,support_fraction=fraction,prior_only_sites=prior,prior_only_fraction=prior_fraction,
                passes_fold=bool(valid and ge(observed,self.fold*expected)),ci_low_ge_fold=False,
                support_status='supported' if not warnings else ';'.join(warnings),mask_fraction=masked,
                technical_status=('not_assessed' if not self.a.exclude_bed else 'masked_analysis' if sensitivity=='exclude_mask' else 'mask_overlap' if masked and masked>0 else 'no_overlap_in_supplied_masks'),
                stable_fraction=stable,inference='conditional_Poisson_exploratory' if valid else 'no_supported_sites' if not supported else 'zero_expected'))
        valid=[i for i,r in enumerate(rows) if r['oe'] is not None]
        if valid:
            obs=np.array([rows[i]['observed'] for i in valid]);exp=np.array([rows[i]['expected'] for i in valid])
            p=poisson.sf(obs-1,exp);lo=np.zeros(len(obs));positive=obs>0
            lo[positive]=.5*chi2.ppf(.025,2*obs[positive])/exp[positive];hi=.5*chi2.ppf(.975,2*(obs+1))/exp
            for i,pv,lv,hv in zip(valid,p,lo,hi):rows[i].update(p=float(pv),ci_low=float(lv),ci_high=float(hv),ci_low_ge_fold=bool(lv>=self.fold))
        families=defaultdict(list)
        for i in valid:families[rows[i]['fdr_family']].append(i)
        for name,indices in families.items():
            for i,q in zip(indices,bh([rows[i]['p'] for i in indices])):rows[i]['fdr']=q
        for r in rows:r['n_family_tests']=len(families.get(r['fdr_family'],[]))
        rows.sort(key=lambda r:(chrom_key(r['chrom']),r['start0'],r['end0'],r['region_kind'],r['annotation'],r['sensitivity']))
        highlights=[r for r in rows if r['passes_fold']];supported=[r for r in highlights if r['support_status']=='supported']
        for name,items in [('regions.tsv.gz',rows),('highlights.tsv.gz',highlights),('highlights_supported.tsv.gz',supported)]:write_tsv(self.out/'results'/name,FIELDS,items)
        bed=[]
        for r in highlights:
            for start,end in self.blocks(self.regions[r['region_id']]):bed.append([r['chrom'],start,end,r['region_id'],r['oe'],r['region_kind'],r['sensitivity'],r['support_status']])
        write_tsv(self.out/'results/highlights.bed.gz',['#chrom','start0','end0','region_id','oe','region_kind','sensitivity','support_status'],bed)
        category_fields=['outcome','sensitivity','region','observed','expected','oe','ci_low','ci_high','p','fdr','n_eligible_supported','n_unsupported','n_blocks','inference']
        categories=[r for r in (category_rows or []) if number(r.get('oe')) is not None and ge(float(r['oe']),self.fold)]
        write_tsv(self.out/'results/category_highlights.tsv',category_fields,categories)
        self.plot(rows,highlights)
        counts=dict(tested_regions=len(valid),catalogued_regions=len(rows),highlight_regions=len(highlights),
            supported_highlight_regions=len(supported),herv_related_highlight_regions=sum(r['annotation_fraction']>0 for r in highlights) if self.a.ref=='herv' else 0,
            prior_dominated_highlights=sum(r['prior_only_fraction'] is not None and r['prior_only_fraction']>.5 for r in highlights),
            mask_overlap_highlights=sum(r['technical_status']=='mask_overlap' for r in highlights),category_highlights=len(categories))
        manifest=dict(schema='allele3-region-highlights-2',status='complete',build='GRCh38',ref=self.a.ref,target=self.a.target,
            population=self.a.population,threshold=self.fold,counts=counts,scientific_settings=self.scientific_settings(),
            selection='Every unrounded finite supported-site O/E >= threshold; no top-N or minimum-count filter in the all-highlights files.',
            support_label='O>=10, E>=1, background-supported fraction>=0.8, prior-only-site fraction<=0.5; labels do not remove raw highlights or certify significance.',
            inference='One-sided Poisson sf(O-1,E), Garwood 95% count interval divided by fixed E. Exploratory: within-region LD and fitted-background uncertainty are not modeled.',
            fdr='BH over all testable regions, including zero-case regions, within region_kind x scale x annotation category x outcome x population x sensitivity. Not a joint genome-wide FDR over every scale/reference.',
            denominator='The same LOCO-supported eligible observed sites contribute to O and E; unsupported cases are retained only in observed_total.',
            coordinate='GRCh38 BED0 half-open; BED exports actual annotation pieces for loci/reference_window. Physical bp is not cohort callability.',
            stable_annotation=str(self.stable_path) if self.stable_path else None,
            outputs=['results/regions.tsv.gz','results/highlights.tsv.gz','results/highlights_supported.tsv.gz','results/highlights.bed.gz',
                'results/category_highlights.tsv','figures/highlights_genome.png',
                'figures/highlights_ranked.png','reports/highlights.html'])
        self.report(highlights,manifest);atomic_json(data_path(self.out,'qc')/'region_highlights.json',manifest)
        return manifest

    def plot(self,rows,highlights):
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        folder=self.out/'figures';folder.mkdir(parents=True,exist_ok=True)
        groups=[('window',w) for w in self.sizes]+[('reference_window',w) for w in self.sizes]+[('locus',0)]
        present={(r['region_kind'],r['window_bp']) for r in rows};groups=[g for g in groups if g in present]
        if not groups:groups=[('window',self.sizes[0])]
        chroms=sorted({r['chrom'] for r in rows},key=chrom_key);offsets={};ticks=[];offset=0
        for c in chroms:
            offsets[c]=offset;length=LENGTHS.get(c,max(r['end0'] for r in rows if r['chrom']==c))
            ticks.append(offset+length/2);offset+=length+1000000
        fig,axes=plt.subplots(len(groups),1,figsize=(14,max(4,len(groups)*2.5)),sharex=True,squeeze=False)
        for ax,(kind,width) in zip(axes[:,0],groups):
            selected=[r for r in rows if (r['region_kind'],r['window_bp'])==(kind,width) and r['oe'] is not None]
            below=[r for r in selected if not r['passes_fold']]
            ax.scatter([offsets[r['chrom']]+(r['start0']+r['end0'])/2 for r in below],[r['oe'] for r in below],s=2,color='#b9c7cc',alpha=.35,rasterized=True)
            for sensitivity,marker in [('all','o'),('exclude_mask','x')]:
                for support,color in [('supported','#007d79'),('other','#d98724')]:
                    found=[r for r in selected if r['passes_fold'] and r['sensitivity']==sensitivity and (r['support_status']=='supported')==(support=='supported')]
                    if not found:continue
                    ax.scatter([offsets[r['chrom']]+(r['start0']+r['end0'])/2 for r in found],[r['oe'] for r in found],s=12,marker=marker,color=color,alpha=.65,
                        label=f'{sensitivity}: '+('supported counts' if support=='supported' else 'limited support'),rasterized=True)
            flagged=[r for r in selected if r['passes_fold'] and r['technical_status']=='mask_overlap']
            if flagged:ax.scatter([offsets[r['chrom']]+(r['start0']+r['end0'])/2 for r in flagged],[r['oe'] for r in flagged],s=22,facecolors='none',edgecolors='#b53340',linewidths=.7,label='Supplied-mask overlap',rasterized=True)
            ax.axhline(self.fold,color='#8a3844',ls='--',lw=1);ax.set_yscale('symlog',linthresh=1);ax.set_ylabel('O/E')
            ax.set_ylim(0,max([self.fold,*[r['oe'] for r in selected]])*1.4)
            ax.set_title(kind.replace('_',' ')+(f' — {width:,} bp' if width else '')+f' | {sum(r["passes_fold"] for r in selected):,} candidates',loc='left',fontsize=10);ax.grid(alpha=.15)
        legend={}
        for ax in axes[:,0]:
            handles,labels=ax.get_legend_handles_labels()
            for handle,label in zip(handles,labels):legend.setdefault(label,handle)
        if legend:axes[0,0].legend(list(legend.values()),list(legend),ncol=3,fontsize=7,frameon=False,loc='upper right')
        axes[-1,0].set_xticks(ticks,[c.removeprefix('chr') for c in chroms]);axes[-1,0].set_xlabel('GRCh38 chromosome; all candidate points retained')
        if chroms:axes[-1,0].set_xlim(0,offset-1000000)
        fig.suptitle(f'Every region with multi-ALT O/E ≥ {self.fold:g} — {self.a.target} / {self.a.ref}',fontsize=14)
        fig.tight_layout(rect=(0,0,1,.98))
        fig.savefig(folder/('highlights_genome.png'),dpi=160)
        plt.close(fig)
        ranked=sorted(highlights,key=lambda r:(r['support_status']!='supported',-r['observed'],-r['oe']))[:30]
        fig,ax=plt.subplots(figsize=(15,max(4,len(ranked)*.30+1.7)))
        if ranked:
            yy=np.arange(len(ranked));ax.scatter([r['oe'] for r in ranked],yy,c=['#007d79' if r['support_status']=='supported' else '#d98724' for r in ranked],s=22)
            for y,r in zip(yy,ranked):
                if r['ci_low'] is not None:ax.hlines(y,max(r['ci_low'],1e-8),r['ci_high'],color='#9aa9ae',lw=1,zorder=0)
            labels=[f'{r["chrom"]}:{r["start0"]+1:,}–{r["end0"]:,} | {r["region_kind"]}'+
                (f' {r["annotation"]}' if r['region_kind']=='locus' else '')+
                f' | {r["sensitivity"]} | O={r["observed"]}, E={r["expected"]:.2g}' for r in ranked]
            ax.set_yticks(yy,labels,fontsize=8);ax.invert_yaxis();ax.set_xscale('log')
        else:ax.text(.5,.5,f'No supported-site O/E reaches {self.fold:g}\nThis is a completed analysis, not missing input.',ha='center',va='center',transform=ax.transAxes)
        ax.axvline(self.fold,color='#8a3844',ls='--');ax.set_xlabel('O/E and conditional Poisson 95% interval (exploratory)')
        ax.set_title('Labelled overview: up to 30 candidates, ordered by support then observed count\nAll candidates remain in the genome plot and downloadable tables')
        fig.tight_layout()
        fig.savefig(folder/('highlights_ranked.png'),dpi=160)
        plt.close(fig)

    def report(self,highlights,manifest):
        columns=['region_kind','chrom','start0','end0','annotation','sensitivity','observed','expected','oe','ci_low','ci_high','fdr',
            'n_family_tests','support_status','technical_status','stable_fraction','prior_only_fraction','n_unsupported']
        records=[[r.get(k) for k in columns] for r in sorted(highlights,key=lambda r:-r['oe'])]
        payload=json.dumps(dict(columns=columns,rows=records),allow_nan=False,separators=(',',':')).replace('<','\\u003c')
        count=manifest['counts'];title=f'{self.a.target}: all regions with O/E ≥ {self.fold:g}'
        page='<!doctype html><html><meta charset="utf-8"><title>'+html.escape(title)+'</title><style>body{font:15px system-ui;color:#18394a;margin:28px auto;padding:0 22px;max-width:1450px}img{max-width:100%}table{border-collapse:collapse;font-size:12px;width:100%}th,td{padding:7px;border-bottom:1px solid #ddd;text-align:left}th{background:#edf4f4}input,button{padding:9px;margin:4px}.scroll{overflow:auto}small{color:#536572}</style><h1>'+html.escape(title)+'</h1>'
        page+=f'<p><b>{count["highlight_regions"]:,}</b> candidate rows; <b>{count["supported_highlight_regions"]:,}</b> with supported counts/background; {count["tested_regions"]:,} testable regions. Counts include separately labelled scales, kinds and mask sensitivities, not independent loci.</p>'
        page+='<p>The threshold selects effect-size candidates. Low counts and supplied-mask overlap remain visible. Conditional Poisson intervals/P values are exploratory; LD and fitted-background uncertainty are not modeled. FDR is calculated before candidate selection within each declared test family.</p>'
        page+='<p><a href="../results/highlights.tsv.gz">Every candidate TSV</a> · <a href="../results/highlights.bed.gz">Every candidate BED0</a> · <a href="../results/highlights_supported.tsv.gz">Supported-count subset</a> · <a href="../results/regions.tsv.gz">All regions, including non-candidates</a> · <a href="../results/category_highlights.tsv">Family/category candidates</a> · <a href="region_enrichment.html">Primary enrichment report</a></p>'
        page+='<h2>Genome-wide display — every candidate point</h2><img src="../figures/highlights_genome.png"><h2>Labelled overview</h2><img src="../figures/highlights_ranked.png">'
        page+='<h2>Search all candidates</h2><p>Coordinates below are zero-based, half-open. Genome windows count all eligible sites; reference windows count only annotation-overlapping sites; loci count their actual annotation pieces. Stable fraction is physical overlap with stable HERV annotations, not a population frequency.</p><input id="query" placeholder="Chromosome, annotation, support or mask status"><button id="previous">Previous</button><button id="next">Next</button><span id="status"></span><div class="scroll"><table><thead id="head"></thead><tbody id="body"></tbody></table></div>'
        page+='<script id="region-data" type="application/json">'+payload+'</script><script>const data=JSON.parse(document.getElementById("region-data").textContent);let page=0;const size=100;const query=document.getElementById("query");function cell(tag,text){let x=document.createElement(tag);x.textContent=text===null?"NA":typeof text==="number"&&!Number.isInteger(text)?text.toPrecision(6):String(text);return x}let h=document.createElement("tr");for(const c of data.columns)h.appendChild(cell("th",c));document.getElementById("head").appendChild(h);function render(){const q=query.value.toLowerCase();const found=data.rows.filter(r=>r.some(x=>String(x).toLowerCase().includes(q)));page=Math.max(0,Math.min(page,Math.max(0,Math.ceil(found.length/size)-1)));const body=document.getElementById("body");body.replaceChildren();for(const r of found.slice(page*size,(page+1)*size)){let tr=document.createElement("tr");for(const x of r)tr.appendChild(cell("td",x));body.appendChild(tr)}document.getElementById("status").textContent=`${found.length} matches; page ${page+1} / ${Math.max(1,Math.ceil(found.length/size))}`;}query.oninput=()=>{page=0;render()};document.getElementById("previous").onclick=()=>{page--;render()};document.getElementById("next").onclick=()=>{page++;render()};render();</script></html>'
        target=self.out/'reports/highlights.html';target.parent.mkdir(parents=True,exist_ok=True);target.write_text(page)
