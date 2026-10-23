"""Sidecar counters; original scanner alone defines eligible loci and ALT QC."""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import csv, gzip, math, json, shutil
try:
    from .common import canonical_chrom, af_bin, an_bin, open_text
    from .core import dump_json, fingerprint
    from .ldblocks import BlockSweep, load_block_sets, hotspot_consensus
    from .erv import ERVOverlapSweep, BoundarySweep, boundary_bin, parse_joined_rmsk
except ImportError:
    from common import canonical_chrom, af_bin, an_bin, open_text
    from core import dump_json, fingerprint
    from ldblocks import BlockSweep, load_block_sets, hotspot_consensus
    from erv import ERVOverlapSweep, BoundarySweep, boundary_bin, parse_joined_rmsk

SCHEMA=2
STRATA=['chrom','context','af_bin','an_bin','sensitivity','population_mode']
COUNTS=['threshold','n_sites','n_multi','n_three']
HEADERS={'global':STRATA+COUNTS,'ld':['ancestry','block_id','start','end']+STRATA+COUNTS,
         'erv':['label']+STRATA+COUNTS}

class ArchitectureAccumulator:
    def __init__(self,ld_block_sets=None,erv_elements=None,*,thresholds=(1,),out_prefix=None,
                 masks=None,flag_annotations=()):
        self.thresholds=tuple(thresholds)
        self.width=1+2*len(self.thresholds)
        self.ld={a:BlockSweep(b) for a,b in (ld_block_sets or {}).items()}
        self.erv_available=erv_elements is not None
        elements=erv_elements or []
        self.erv_sweep=ERVOverlapSweep(elements)
        self.boundaries={
            'erv_core':BoundarySweep(elements,('erv_core',)),
            'erv_related':BoundarySweep(elements,('erv_core','erv_related')),
            'ltr_all':BoundarySweep(elements,('erv_core','erv_related','ltr_other'))}
        self.global_strata=defaultdict(self.zeros)
        self.ld_strata=defaultdict(self.zeros)
        self.erv_strata=defaultdict(self.zeros)
        self.ld_current={}; self.ld_pending={a:defaultdict(self.zeros) for a in self.ld}
        self.out_prefix=out_prefix; self.ld_handle=None
        self.masks=masks or {}; self.flag_annotations=set(flag_annotations)
        self.sensitivities={'all'}
        self.classes=sorted({e.rep_class for e in elements})
        if out_prefix:
            Path(out_prefix).parent.mkdir(parents=True,exist_ok=True)
            self.ld_handle=gzip.open(str(out_prefix)+'.ld_strata.tsv.gz','wt',compresslevel=2,newline='')
            self.ld_writer=csv.writer(self.ld_handle,delimiter='\t'); self.ld_writer.writerow(HEADERS['ld'])
    def zeros(self): return [0]*self.width
    @staticmethod
    def stratum(context,af,an): return (str(context or 'NA'),af_bin(af),an_bin(an))
    def _add(self,table,key,values):
        z=table[key]
        for i,v in enumerate(values): z[i]+=v
    def _rows(self,data):
        for k,v in sorted(data.items()):
            for i,t in enumerate(self.thresholds): yield [*k,t,v[0],v[1+2*i],v[2+2*i]]
    def _flush_ld(self,anc):
        self.ld_writer.writerows(self._rows(self.ld_pending[anc])); self.ld_pending[anc].clear()
    def close(self):
        if self.ld_handle: self.ld_handle.close(); self.ld_handle=None
    def observe(self,chrom,pos1,n_alt,context,primary_alt_af,primary_alt_an,*,
                threshold_counts=None,population_mode='NA',region_flags=()):
        c=canonical_chrom(chrom); n=int(n_alt)
        counts=threshold_counts or {1:n}
        if any(t not in counts for t in self.thresholds): raise ValueError('Persistence counts must come from the original scanner')
        values=[1]+[v for t in self.thresholds for v in (int(counts[t]>=2),int(counts[t]>=3))]
        st=self.stratum(context,primary_alt_af,primary_alt_an)
        overlaps=self.erv_sweep.locate(c,pos1) if self.erv_available else ()
        sensitivities=['all']; flags=set(region_flags)
        exclusions={}
        for flag in ('lcr','segdup'):
            if flag in self.masks: exclusions[flag]=self.masks[flag].contains(c,pos1)
            elif flag in self.flag_annotations: exclusions[flag]=flag in flags
        if all(x in exclusions for x in ('lcr','segdup')) and not (exclusions['lcr'] or exclusions['segdup']):
            sensitivities.append('exclude_lcr_segdup')
        if self.erv_available:
            exclusions['simple_repeat']=any(e.rep_class in {'Simple_repeat','Low_complexity'} for e in overlaps)
            if not exclusions['simple_repeat']: sensitivities.append('exclude_simple_low_complexity')
        if 'low_mappability' in self.masks:
            exclusions['low_mappability']=self.masks['low_mappability'].contains(c,pos1)
            if not exclusions['low_mappability']: sensitivities.append('exclude_low_mappability')
        if len(exclusions)==4 and not any(exclusions.values()): sensitivities.append('exclude_combined')
        self.sensitivities.update(sensitivities)
        modes=['all']
        if population_mode!='NA' and n>=2: modes.append(population_mode)
        labels=[]
        if self.erv_available:
            ltr=[e for e in overlaps if e.rep_class=='LTR']
            labels=['all_sites','inside_any_ltr' if ltr else 'outside_ltr']
            if any(e.tier=='erv_core' for e in ltr): labels.append('inside_erv_core')
            if any(e.tier in {'erv_core','erv_related'} for e in ltr): labels.append('inside_erv_related')
            labels += ['family:'+fam for fam in sorted({e.family or 'unspecified' for e in ltr})]
            labels += ['negative_control:'+cl for cl in sorted({e.rep_class for e in overlaps}&{'LINE','SINE','DNA','Simple_repeat','Low_complexity'})]
            for tier,sweep in self.boundaries.items():
                for bt in ('5p','3p'):
                    near=sweep.nearest(c,pos1,bt)
                    if near:
                        label=f'{tier}:{bt}_boundary:{boundary_bin(near["distance"])}'
                        labels.append(label)
                        labels.append(label+(':upstream' if near['signed_transcription_distance']<0 else ':downstream'))
        locations=[]
        for anc,sweep in self.ld.items():
            b=sweep.locate(c,pos1)
            key=(anc,b.block_id,b.start,b.end) if b else (anc,'outside_partition',-1,-1)
            if self.ld_handle:
                if self.ld_current.get(anc)!=(c,key): self._flush_ld(anc); self.ld_current[anc]=(c,key)
                table=self.ld_pending[anc]
            else: table=self.ld_strata
            locations.append((table,key))
        for sens in sensitivities:
            for mode in modes:
                key=(c,*st,sens,mode)
                self._add(self.global_strata,key,values)
                for table,prefix in locations: self._add(table,(*prefix,*key),values)
                for label in labels: self._add(self.erv_strata,(label,*key),values)
    def write(self,out_prefix):
        Path(out_prefix).parent.mkdir(parents=True,exist_ok=True)
        for kind,data in [('global',self.global_strata),('ld',self.ld_strata),('erv',self.erv_strata)]:
            if kind=='ld' and self.ld_handle:
                for anc in self.ld: self._flush_ld(anc)
                self.close(); continue
            with gzip.open(str(out_prefix)+f'.{kind}_strata.tsv.gz','wt',compresslevel=2,newline='') as f:
                w=csv.writer(f,delimiter='\t'); w.writerow(HEADERS[kind]); w.writerows(self._rows(data))
        dump_json(str(out_prefix)+'.availability.json',dict(schema=SCHEMA,
            erv_status='measured' if self.erv_available else 'NA',repeat_classes=self.classes,
            ld_ancestries=sorted(self.ld),sensitivity_analyses=sorted(self.sensitivities),
            missing_sensitivity_annotations=sorted({'lcr','segdup','low_mappability'}-self.flag_annotations-set(self.masks)),
            denominator='Fixed original eligible loci at every AC tier; population profiles conditional on base multi-ALT population class',
            boundaries='Observed joined alignStart/alignEnd, BED edges, absolute distance <10000; signed distance follows repeat strand',
            body='Union of aligned joined blocks, excluding unaligned tails and nested insertion gaps'))

def _read_tsv(path):
    with open_text(path) as f: yield from csv.DictReader(f,delimiter='\t')

def annotation_identity(options):
    if not options.get('with_architecture'): return None
    files={}
    directory=options.get('ld_block_dir')
    if directory:
        load_block_sets(directory,options['architecture_build'])
        files.update({anc:fingerprint(Path(directory)/f'{anc}.{options["architecture_build"]}.bed') for anc in ('AFR','EAS','EUR','SAS')})
    for key in ('erv_joined','architecture_segdup_bed','architecture_lcr_bed','architecture_low_mappability_bed'):
        if options.get(key): files[key]=fingerprint(options[key])
    return dict(schema=SCHEMA,files=files)

def from_options(options,chrom,shard_dir,header):
    if not options.get('with_architecture'): return None
    blocks=load_block_sets(options['ld_block_dir'],options['architecture_build']) if options.get('ld_block_dir') else None
    elements=None
    if options.get('erv_joined'):
        source=options['erv_joined']
        if options.get('architecture_annotation_cache') and chrom:
            source=Path(options['architecture_annotation_cache'])/f'chr{chrom}.tsv.gz'
        elements=parse_joined_rmsk(source,chromosome=chrom)
    try: from .intervals import IntervalMask
    except ImportError: from intervals import IntervalMask
    masks={name:IntervalMask(options[key],[chrom] if chrom else None) for name,key in
        [('lcr','architecture_lcr_bed'),('segdup','architecture_segdup_bed'),('low_mappability','architecture_low_mappability_bed')]
        if options.get(key)}
    return ArchitectureAccumulator(blocks,elements,thresholds=options['thresholds'],
        out_prefix=Path(shard_dir)/options.get('architecture_prefix','architecture'),masks=masks,
        flag_annotations={'lcr','segdup'}&set(header['info']))

def loco_expected(global_path,target_path,kind,min_training_sites=1):
    """Streaming targets; report total observations even when LOCO training is missing."""
    globals_=defaultdict(lambda:[0,0,0]); by_chr=defaultdict(lambda:[0,0,0])
    dims=['context','af_bin','an_bin','sensitivity','population_mode','threshold']
    def stratum(r): return tuple(r.get(x,{'sensitivity':'all','population_mode':'all','threshold':'1'}.get(x,'NA')) for x in dims)
    for r in _read_tsv(global_path):
        st=stratum(r); vals=[int(r[x]) for x in ('n_sites','n_multi','n_three')]
        for table,key in ((globals_,st),(by_chr,(r['chrom'],*st))):
            for i,v in enumerate(vals): table[key][i]+=v
    ids=(['ancestry','block_id','chrom','start','end'] if kind=='ld' else ['label','chrom'])+['sensitivity','population_mode','threshold']
    grouped=defaultdict(lambda:[0,0,0,0.,0.,0.,0])
    for r in _read_tsv(target_path):
        gid=tuple(r.get(k,{'sensitivity':'all','population_mode':'all','threshold':'1'}.get(k,'NA')) for k in ids)
        z=grouped[gid]; n,m,t=[int(r[k]) for k in ('n_sites','n_multi','n_three')]
        z[0]+=n; z[1]+=m; z[2]+=t
        st=stratum(r); train=[a-b for a,b in zip(globals_[st],by_chr[(r['chrom'],*st)])]
        if train[0]<min_training_sites: z[6]+=n; continue
        pm=(train[1]+.5)/(train[0]+1); pt=(train[2]+.5)/(train[0]+1)
        z[3]+=n*pm; z[4]+=n*pt; z[5]+=n*pm*(1-pm)
    rows=[]
    for gid,(n,m,t,em,et,var,missing) in grouped.items():
        rows.append(dict(zip(ids,gid))|dict(n_sites=n,observed_multi=m,observed_three=t,
            expected_multi=em if not missing else 'NA',expected_three=et if not missing else 'NA',
            expected_multi_supported=em,n_sites_without_loco=missing,
            oe=m/em if em and not missing else 'NA',oe_three=t/et if et and not missing else 'NA',
            z_descriptive=(m-em)/math.sqrt(var) if var and not missing else 'NA',
            status='complete' if not missing else 'insufficient_loco_training'))
    return sorted(rows,key=lambda r:(r['oe']=='NA',-r['oe'] if r['oe']!='NA' else 0,str(r)))

def write_summary(rows,path,fields=None):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else ['status']
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)

def summarize_architecture(shards,out,options,min_training_sites=1):
    out=Path(out); prefix=options.get('architecture_prefix','architecture')
    for kind,header in HEADERS.items():
        with gzip.open(out/f'architecture.{kind}_strata.tsv.gz','wt',compresslevel=2,newline='') as f:
            w=csv.DictWriter(f,fieldnames=header,delimiter='\t'); w.writeheader()
            for shard in shards:
                source=Path(shard)/f'{prefix}.{kind}_strata.tsv.gz'
                if not source.exists(): raise ValueError(f'Architecture counter missing: {source}; old scans cannot be reused')
                w.writerows(_read_tsv(source))
    global_path=out/'architecture.global_strata.tsv.gz'
    for kind in ('ld','erv'):
        rows=loco_expected(global_path,out/f'architecture.{kind}_strata.tsv.gz',kind,min_training_sites)
        write_summary(rows,out/('ld_block_enrichment.tsv' if kind=='ld' else 'erv_enrichment.tsv'))
        if kind=='ld':
            ranked=[r for r in rows if r['block_id']!='outside_partition' and r['n_sites']>=options.get('architecture_min_sites',10000) and r['oe']!='NA']
            write_summary(ranked,out/'ld_block_ranked.tsv',list(rows[0]) if rows else None)
            consensus=hotspot_consensus(rows,options.get('architecture_min_sites',10000),options.get('architecture_top_fraction',.01))
            write_summary(consensus,out/'ld_hotspot_consensus.tsv',['chrom','start','end','overlap_bp','supporting_ancestries','n_ancestries','block_ids'])
    availability=[json.loads((Path(s)/f'{prefix}.availability.json').read_text()) for s in shards]
    dump_json(out/'architecture_manifest.json',dict(schema=SCHEMA,status='complete',shards=availability,
        inputs=annotation_identity(options),ld_primary='EUR',min_sites=options.get('architecture_min_sites',10000),
        missing_annotation='NA; unavailable annotation is never an observed zero',
        inference='Descriptive O/E and residuals, not independent-site GWAS P values'))
