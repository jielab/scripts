#!/usr/bin/env python3
"""Audit reported pQTL loci. No unreported edge is zero; no causal claims."""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, math, os, re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from summary_metadata import for_summary
DATASETS = {'ukb':'prot','ckb':'ckb','decode':'decode','scot':'scot'}
PLATFORMS = {'ukb':'Olink','ckb':'Olink / Chinese','decode':'SomaScan','scot':'LC-MS/MS'}
DEFAULT_SENTINELS = Path(os.environ.get('OS_SNP_LIST', str(Path(os.environ.get('DIR_UKB', str(Path(__file__).resolve().parents[3]/'phe/ukb')))/'common/snp.lst')))
PANEL = re.compile(r'^(Cardiometabolic|Inflammation|Neurology|Oncology)(?:_II)?$')
CODING = {'protein_coding','IG_C_gene','IG_D_gene','IG_J_gene','IG_V_gene','TR_C_gene','TR_D_gene','TR_J_gene','TR_V_gene'}
PRIMARY_CHROMS = {*(str(i) for i in range(1,23)),'X','Y','MT'}
def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path, fields, rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); tmp=Path(str(path)+'.tmp')
    with tmp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader()
        for r in rows: w.writerow({k:'.' if v is None else v for k,v in r.items()})
    tmp.replace(path)
def write_json(path,value):
    tmp=Path(str(path)+'.tmp'); tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)); tmp.replace(path)
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def chrom(v):
    v=re.sub(r'^chr','',str(v),flags=re.I).upper()
    return {'23':'X','24':'Y','25':'MT','M':'MT'}.get(v,v)
def chrom_key(v): return (0,int(v)) if v.isdigit() else (1,{'X':'01','Y':'02','MT':'03'}.get(v,v))
def pnum(v):
    try:
        p=Decimal(str(v))
        if not p.is_finite() or not 0<=p<=1: raise ValueError
        return p
    except (InvalidOperation,ValueError): raise ValueError(f'invalid P: {v}') from None
def significant(v,c): return pnum(v)<=Decimal(str(c))
def numeric(v):
    try:
        n=float(v); return n if math.isfinite(n) else None
    except (TypeError,ValueError): return None
def parse_signals(path):
    signals=[]; header=None; ended=False
    with Path(path).open(encoding='utf-8-sig') as f:
        for number,line in enumerate(f,1):
            line=line.strip()
            if line.lower().startswith('magma:'): ended=True; break
            if not line or line.lower()=='common:' or line.startswith('#'): continue
            fields=line.split('\t')
            if fields[0]=='SNP':
                header=fields
                if not {'SNP','CHR','POS','P','CLASS'}.issubset(header): raise ValueError(f'{path}:{number}: missing columns')
                continue
            if header is None or len(fields)!=len(header): raise ValueError(f'{path}:{number}: malformed common row')
            r=dict(zip(header,fields)); ch=chrom(r['CHR']); pos=int(r['POS']); kind=r['CLASS'].lower()
            if ch not in PRIMARY_CHROMS or pos<1 or kind not in ('cis','trans'):
                raise ValueError(f'{path}:{number}: invalid CHR/POS/CLASS')
            p=pnum(r['P']); flags=[]
            if p==0: flags.append('p_underflow_zero')
            if p>Decimal('5e-8'): flags.append('reported_nonsignificant')
            beta,se=numeric(r.get('BETA')),numeric(r.get('SE'))
            if beta is None or se is None or se<=0: flags.append('invalid_effect_or_se')
            m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([^:]+):([^:]+)',r['SNP'],flags=re.I)
            if m:
                if m[3].upper()==m[4].upper(): flags.append('identical_alleles_in_id')
                if chrom(m[1])!=ch or int(m[2])!=pos: flags.append('id_coordinate_differs')
            signals.append(dict(chr=ch,pos=pos,snp=r['SNP'],p=r['P'],beta=beta,se=se,gene=r.get('GENE','.'),
                                raw_kind=kind,kind=kind,flags=flags,ea=r.get('EA'),nea=r.get('NEA')))
    if header is None: raise ValueError(f'{path}: missing common header')
    meta=for_summary(path)
    if not ended and not meta: raise ValueError(f'{path}: truncated legacy summary')
    return signals,meta
def load_annotation(path):
    groups=defaultdict(list)
    for r in read_tsv(path):
        if r.get('CHR','.')=='.': continue
        t=dict(chr=chrom(r['CHR']),start=int(r['START']),end=int(r['END']),gene=r['gene_symbol'],
             gene_id=r.get('gene_id','.'),build=r.get('build','38'),match=r.get('match_method','provided'),
             status=r.get('annotation_status','single_gene'),biotype=r.get('biotype','protein_coding'),source=r.get('source',str(path)))
        for key in ('coordinate_status','coordinate_source','coordinate_window','supplied_CHR','supplied_POS','shared_anchor_targets'):
            if key in r: t[key]=r[key]
        if 'identity_verified' in r: t['identity_verified']=str(r['identity_verified']).lower() in ('1','true')
        for key in ('body_start','body_end','point_pos'):
            if r.get(key) not in (None,'','.'): t[key]=int(r[key])
        groups[r['prot']].append(t)
    return groups

def resolve_annotation(rows,reference,coordinate_window='body',map_build='38'):
    """Resolve authoritative mapping points without inventing gene identities.

    A supplied CHR/POS is a one-based point. A reference gene body is used only
    when that point lies within the named reference target. Explicit missing,
    invalid, alternate-contig, or conflicting coordinates never fall back to an
    unrelated primary assembly interval. Known complexes keep all components;
    a single supplied constituent cannot turn them into a single mother gene.
    Reference PAR X/Y representations are retained, including in point mode.
    """
    if coordinate_window not in ('body','point'): raise ValueError('coordinate_window must be body or point')
    if map_build not in ('37','38'): raise ValueError('map_build must be 37 or 38')
    if not rows: return {},[]
    has_coordinates='CHR' in rows[0] or 'POS' in rows[0]
    if has_coordinates and not {'CHR','POS'}.issubset(rows[0]): raise ValueError('Mapping coordinates require both CHR and POS')
    missing={'','.','NA','N/A','NAN','NONE'}
    shared=defaultdict(list)
    for r in rows:
        c=chrom(r.get('CHR','.')); raw=str(r.get('POS','.')).strip()
        if c not in missing and re.fullmatch(r'\d+',raw) and int(raw)>0:
            shared[(c,int(raw))].append(r['prot'])
    effective={}; audit=[]
    for r in rows:
        prot=r['prot']; refs=[dict(t) for t in reference.get(prot,[])]
        supplied_c=chrom(r.get('CHR','.')); supplied_p=str(r.get('POS','.')).strip()
        point=int(supplied_p) if re.fullmatch(r'\d+',supplied_p) and int(supplied_p)>0 else None
        matches=[]; chosen=[]; verified=False
        same_name=bool(refs) and len({t['gene'] for t in refs})==1
        ref_single=same_name and all(t['status']=='single_gene' for t in refs)
        # Two coordinate records for the same PAR target are one biological
        # gene, not two independent protein measurements or causal genes.
        par=ref_single and {t['chr'] for t in refs}=={'X','Y'}
        if ref_single and len({(t['chr'],t.get('gene_id',t['gene'])) for t in refs})>1 and not par:
            ref_single=False
        rec=dict(prot=prot,supplied_CHR=r.get('CHR','.'),supplied_POS=r.get('POS','.'),
                 supplied_build=map_build if has_coordinates else '.',coordinate_window=coordinate_window,
                 coordinate_status='unresolved',coordinate_source='none',identity_verified=False,
                 shared_anchor_targets=';'.join(p for p in shared.get((supplied_c,point),[]) if p!=prot) or '.',
                 matched_reference_genes='.',matched_reference_ids='.',reference_intervals=';',used_intervals='.',
                 anchor_matches_reference_midpoint=None,notes='.')
        rec['reference_intervals']=';'.join(f"{t['chr']}:{t['start']}-{t['end']}" for t in refs) or '.'
        if not has_coordinates:
            chosen=refs; verified=ref_single
            rec.update(coordinate_status='reference_only' if refs else 'unresolved',coordinate_source='reference_annotation',
                       notes='Mapping has no CHR/POS columns; retained legacy reference definition')
        elif supplied_c in missing and supplied_p.upper() in missing:
            rec.update(coordinate_status='supplied_missing',notes='Explicit missing coordinates suppress reference fallback')
        elif supplied_c in missing or point is None:
            rec.update(coordinate_status='supplied_invalid',notes='CHR and positive integer POS are both required')
        elif supplied_c not in PRIMARY_CHROMS:
            rec.update(coordinate_status='supplied_alternative_contig',
                       notes='Alternate-contig position cannot be compared with primary-assembly GWAS coordinates')
        elif refs:
            matches=[t for t in refs if str(t['build'])==map_build and t['chr']==supplied_c and t['start']<=point<=t['end']]
            if not matches or any(str(t['build'])!=map_build for t in refs):
                rec.update(coordinate_status='supplied_reference_conflict',
                           notes='Supplied point does not corroborate the named reference intervals/build; manual review required')
            else:
                chosen=refs; verified=ref_single
                rec.update(coordinate_status='matched_reference_par' if par else 'matched_reference_single' if ref_single else 'matched_reference_complex',
                           coordinate_source='supplied_point_corroborated_reference',
                           notes='All known complex/PAR intervals retained; supplied point matches one component' if len(refs)>1 else '.')
        else:
            # A valid point is useful for regional screening but does not certify
            # that a gene symbol, assay target, or coding sequence is correct.
            chosen=[dict(chr=supplied_c,start=point,end=point,gene=prot,gene_id='.',build=map_build,
                         match='supplied_point_only',status='point_unverified',biotype='unknown',source='protein_map_CHR_POS')]
            rec.update(coordinate_status='supplied_point_unverified',coordinate_source='supplied_point_only',
                       notes='Point window usable for screening; coding identity and gene body remain unverified')
        if matches:
            rec['matched_reference_genes']=';'.join(sorted({t['gene'] for t in matches}))
            rec['matched_reference_ids']=';'.join(sorted({t.get('gene_id','.') for t in matches}))
            rec['anchor_matches_reference_midpoint']=any(abs((t['start']+t['end'])/2-point)<=1 for t in matches)
        rec['identity_verified']=verified and all(t['biotype'] in CODING for t in chosen)
        used=[]
        for t in chosen:
            t=dict(t); body_start=t.get('body_start',t['start']); body_end=t.get('body_end',t['end'])
            anchor=point if has_coordinates and t['chr']==supplied_c and body_start<=point<=body_end else (body_start+body_end)//2
            t.update(body_start=body_start,body_end=body_end,point_pos=anchor,
                     start=anchor if coordinate_window=='point' else body_start,end=anchor if coordinate_window=='point' else body_end,
                     coordinate_window=coordinate_window,coordinate_status=rec['coordinate_status'],coordinate_source=rec['coordinate_source'],
                     supplied_CHR=rec['supplied_CHR'],supplied_POS=rec['supplied_POS'],identity_verified=rec['identity_verified'],
                     shared_anchor_targets=rec['shared_anchor_targets'])
            used.append(t)
        effective[prot]=used
        rec['used_intervals']=';'.join(f"{t['chr']}:{t['start']}-{t['end']}" for t in used) or '.'
        audit.append(rec)
    return effective,audit

def window_intervals(intervals,window):
    """The two definitions share the same target identity/build decisions."""
    return [{**t,'start':t.get('point_pos',(t['start']+t['end'])//2) if window=='point' else t.get('body_start',t['start']),
                 'end':t.get('point_pos',(t['start']+t['end'])//2) if window=='point' else t.get('body_end',t['end'])} for t in intervals]

def coordinates_used_rows(result):
    rows=[]
    for audit in result['coordinate_audit']:
        prot=audit['prot']; intervals=result['annotation'].get(prot,[])
        for t in intervals or [{}]:
            rows.append(dict(prot=prot,gene_symbol=t.get('gene','.'),gene_id=t.get('gene_id','.'),CHR=t.get('chr','.'),
                START=t.get('start','.'),END=t.get('end','.'),build=t.get('build',audit['supplied_build']),
                annotation_status=t.get('status','unresolved'),biotype=t.get('biotype','.'),source=t.get('source','.'),
                match_method=t.get('match','.'),coordinate_source=audit['coordinate_source'],coordinate_status=audit['coordinate_status'],
                coordinate_window=audit['coordinate_window'],supplied_CHR=audit['supplied_CHR'],supplied_POS=audit['supplied_POS'],
                identity_verified=int(audit['identity_verified']),shared_anchor_targets=audit['shared_anchor_targets'],
                body_start=t.get('body_start','.'),body_end=t.get('body_end','.'),point_pos=t.get('point_pos','.')))
    return rows
def infer_build(records,anchors,explicit=None):
    votes=Counter(); matched=defaultdict(set); conflicts=[]
    for alias,r in records.items():
        for s in r['signals']:
            ref=anchors.get(s['snp'])
            if not ref: continue
            hit=[b for b in ('37','38') if s['chr']==chrom(ref['CHR']) and s['pos']==int(ref['POS.'+b])]
            if len(hit)==1: votes[hit[0]]+=1; matched[hit[0]].add(s['snp'])
            elif not hit: conflicts.append(dict(assay=alias,snp=s['snp'],chr=s['chr'],pos=s['pos']))
    inferred=next((b for b in ('37','38') if len(matched[b])>=3 and not votes['37' if b=='38' else '38'] and not conflicts),None)
    if explicit and (votes['37' if explicit=='38' else '38'] or conflicts): raise ValueError('Explicit build conflicts with sentinels')
    return dict(build=explicit or inferred,source='explicit' if explicit else 'dataset_sentinels' if inferred else 'unknown',
                votes=dict(votes),unique_sentinels={b:len(matched[b]) for b in ('37','38')},conflicts=conflicts,
                scope='Dataset inference does not certify every assay')
def inside(s,intervals,flank):
    return any(s['chr']==t['chr'] and t['start']-flank<=s['pos']<=t['end']+flank for t in intervals)
def best(signals,alias=None):
    hits=[s for s in signals if s['kind']=='cis' and significant(s['p'],5e-8)]
    if not hits: return None
    s=min(hits,key=lambda x:(pnum(x['p']),x['pos'],x['snp']))
    return {k:s[k] for k in ('snp','chr','pos','p')}|{'assay':alias or ';'.join(s.get('assays',[]))}
def flag(code,label,level='info'): return dict(code=code,label=label,level=level)
def metadata_errors(path,signals,meta,alias):
    if not meta: return []
    errors=[]
    try:
        if int(meta.get('schema_version','0'))<2: errors.append('metadata_schema_invalid')
        if int(meta.get('common_rows','-1'))!=len(signals): errors.append('metadata_row_count_mismatch')
        if meta.get('sig_md5')!=hashlib.md5(Path(path).read_bytes()).hexdigest(): errors.append('metadata_checksum_mismatch')
        if meta.get('trait')!=alias: errors.append('metadata_assay_mismatch')
        if pnum(meta.get('p_threshold','.'))!=Decimal('5e-8'): errors.append('metadata_threshold_mismatch')
        if int(meta.get('cis_flank_bp','-1'))<0: errors.append('metadata_flank_invalid')
        if int(meta.get('cis_variant_rows','-1'))<0: errors.append('metadata_cis_coverage_invalid')
    except (ValueError,TypeError): errors.append('metadata_malformed')
    return errors

def complete_cis_search(meta,intervals,file_build,flank,errors):
    if not meta or errors: return False
    if any(t['status'] not in ('single_gene','complex_complete') for t in intervals): return False
    try:
        bounds=[]
        for s in meta.get('cis_intervals','').split(';'):
            m=re.fullmatch(r'([^:]+):(\d+)-(\d+)',s)
            if m: bounds.append((chrom(m[1]),int(m[2]),int(m[3])))
        covered=all(any(c==t['chr'] and start<=max(1,t['start']-flank) and end>=t['end']+flank
                    for c,start,end in bounds) for t in intervals)
        exact=set(bounds)=={(t['chr'],max(1,t['start']-flank),t['end']+flank) for t in intervals}
        no_hidden_lead=exact or int(meta.get('cis_significant_rows','-1'))==0
        return (covered and no_hidden_lead and bool(intervals) and meta.get('cis_annotation_available','').lower()=='true' and
            meta.get('trans_selection')=='all_position_clumped_loci' and meta.get('association_scope')=='full_gwas' and
            meta.get('p_source')=='raw_gwas' and meta.get('cis_window')=='gene_body_flank' and
            int(meta.get('cis_flank_bp','-1'))>=flank and int(meta.get('cis_variant_rows','0'))>0 and
            str(meta.get('cis_annotation_build',''))==file_build)
    except (ValueError,TypeError): return False
def collect(root,mapping,annotation,anchors,flank=500000,trans_rich=5,builds=None,coordinate_window='body',map_build='38'):
    rows=read_tsv(mapping)
    if not rows or not {'prot',*DATASETS}.issubset(rows[0]): raise ValueError('Mapping requires prot/ukb/ckb/decode/scot')
    if len({r['prot'] for r in rows})!=len(rows): raise ValueError('Duplicate prot in mapping')
    annotation,coordinate_audit=resolve_annotation(rows,annotation,coordinate_window,map_build)
    coordinate_lookup={r['prot']:r for r in coordinate_audit}
    cache={}; build_info={}; datasets={}
    file_audit=[]; coverage=[]; assay_rows=[]; signal_rows=[]; candidates=[]; repairs=[]
    all_aliases={d:{a for r in rows for a in r[d].split(';') if a and a!='.'} for d in DATASETS}
    for d,folder in DATASETS.items():
        records={}
        for path in sorted((root/folder/'mplot').glob('*.sig.txt')):
            alias=path.name.removesuffix('.sig.txt')
            try:
                signals,meta=parse_signals(path); error=''
            except (ValueError,OSError) as e: signals,meta,error=[],{},str(e)
            records[alias]=dict(signals=signals,meta=meta,error=error,path=path)
            file_audit.append(dict(dataset=d,assay=alias,file=str(path),sha256=digest(path),rows=len(signals),error=error or '.'))
        cache[d]=records; build_info[d]=infer_build(records,anchors,(builds or {}).get(d))
        raw=sum(s['raw_kind']=='cis' for r in records.values() for s in r['signals'])
        datasets[d]=dict(files=len(records),raw_cis_rows=raw,all_trans_labels=raw==0,build=build_info[d]['build'],build_source=build_info[d]['source'])
    resolved=[]
    for row in rows:
        row=dict(row)
        for d in DATASETS:
            names=[]
            for alias in [a.strip() for a in row[d].split(';') if a.strip() and a.strip()!='.']:
                if '/' in alias or '\\' in alias or alias in ('.','..'): raise ValueError(f'Unsafe assay: {alias}')
                hits=[]
                if d=='ukb' and alias not in cache[d]:
                    hits=sorted(a for a in cache[d] if a.startswith(alias+'.') and PANEL.fullmatch(a[len(alias)+1:]) and a not in all_aliases[d])
                if hits:
                    repairs.append(dict(prot=row['prot'],dataset=d,requested=alias,resolved=';'.join(hits),rule='exact_base_known_Olink_panel_suffix'))
                    names.extend(hits)
                else: names.append(alias)
            row[d]=';'.join(dict.fromkeys(names)) or '.'
        resolved.append(row)
    strict={d:5e-8/max(1,len({a for r in resolved for a in r[d].split(';') if a!='.'})) for d in DATASETS}
    # This certifies the supplied summaries, not coverage of an unmeasured or
    # missing protein. Missing assays remain separately visible in coverage.tsv.
    for d,records in cache.items():
        full=0
        for alias,r in records.items():
            m=r['meta']
            full+=int(not r['error'] and bool(m) and not metadata_errors(r['path'],r['signals'],m,alias)
                and m.get('association_scope')=='full_gwas' and m.get('p_source')=='raw_gwas'
                and m.get('beta_se_source')=='raw_gwas' and m.get('trans_selection')=='all_position_clumped_loci'
                and m.get('cis_window')=='gene_body_flank'
                and str(m.get('build','')).removeprefix('GRCh')==build_info[d]['build'])
        datasets[d].update(full_summary_files=full,complete_full_summaries=bool(records) and full==len(records),
                           completeness_scope='all supplied summaries; missing mapped assays are reported separately')
    proteins=[]
    for row in resolved:
        intervals=annotation.get(row['prot'],[]); single=bool(intervals) and all(t['status']=='single_gene' for t in intervals)
        coord=coordinate_lookup[row['prot']]; identity_verified=coord['identity_verified']
        p=dict(prot=row['prot'],target=intervals[0] if single else None,intervals=intervals,coordinate=coord,
               identity_verified=identity_verified,data={},flags=[],review_by_primary={d:[] for d in DATASETS})
        if not intervals: p['flags'].append(flag('no_target_annotation','目标基因坐标未解析','review'))
        elif any(t['status']=='point_unverified' for t in intervals):
            p['flags'].append(flag('target_identity_unverified','仅有提供的位置；编码基因身份和基因范围未核实','review'))
        elif not single: p['flags'].append(flag('composite_target','复合物/多成分目标不能指认为单一母基因','review'))
        if any(t['biotype'] not in CODING and t['biotype']!='unknown' for t in intervals):
            p['flags'].append(flag('noncoding_target','非编码/假基因映射需回查检测目标','review'))
        coordinate_labels={'supplied_missing':'提供的 CHR/POS 缺失；不以旧注释补回',
            'supplied_invalid':'提供的 CHR/POS 无效或不完整',
            'supplied_alternative_contig':'提供的是替代 contig 坐标，不能直接对应主组装 GWAS',
            'supplied_reference_conflict':'提供的位置与参考目标/版本冲突，暂不分类 cis/trans'}
        if coord['coordinate_status'] in coordinate_labels:
            p['flags'].append(flag(coord['coordinate_status'],coordinate_labels[coord['coordinate_status']],'review'))
        if coord['shared_anchor_targets']!='.':
            p['flags'].append(flag('shared_anchor_targets','同一提供坐标还有目标：'+coord['shared_anchor_targets']))
        if coord['anchor_matches_reference_midpoint'] is False:
            p['flags'].append(flag('supplied_point_not_reference_midpoint','提供点在参考基因内，但不是该参考范围的中点'))
        body_intervals=window_intervals(intervals,'body'); point_intervals=window_intervals(intervals,'point')
        for d in DATASETS:
            aliases=[a for a in row[d].split(';') if a!='.']; merged={}; assays=[]; missing=[]; invalid=[]; corrected=0
            for alias in aliases:
                r=cache[d].get(alias)
                a=dict(alias=alias,status='complete',cis=None,trans=0,raw_cis=0,best_cis=None,flags=[],
                       cis_search_complete=False,cis_search_status='unavailable',coordinate_sensitive_rows=0)
                if r is None:
                    a['status']='missing'; missing.append(alias); assays.append(a); continue
                if r['error']:
                    a['status']='invalid'; a['flags'].append(r['error']); invalid.append(alias); assays.append(a); continue
                meta=r['meta']; file_build=str(meta.get('grch',meta.get('build',''))).removeprefix('GRCh') or build_info[d]['build']
                meta_errors=metadata_errors(r['path'],r['signals'],meta,alias)
                a['flags'].extend(meta_errors)
                compatible=bool(intervals) and file_build is not None and all(t['build']==file_build for t in intervals)
                if meta and file_build!=build_info[d]['build'] and build_info[d]['build']:
                    compatible=False; a['flags'].append('file_dataset_build_conflict')
                classified=[]
                for orig in r['signals']:
                    s={**orig,'flags':list(orig['flags']),'assays':[alias]}
                    s['kind']=('cis' if inside(s,intervals,flank) else 'trans') if compatible else 'unknown'
                    s['kind_body']=('cis' if inside(s,body_intervals,flank) else 'trans') if compatible else 'unknown'
                    s['kind_point']=('cis' if inside(s,point_intervals,flank) else 'trans') if compatible else 'unknown'
                    s['coordinate_sensitive']=compatible and s['kind_body']!=s['kind_point']
                    if s['coordinate_sensitive']:
                        s['flags'].append('cis_window_definition_sensitive'); a['coordinate_sensitive_rows']+=1
                    if compatible and s['kind']!=s['raw_kind']:
                        corrected+=1; s['flags'].append('class_changed_by_coordinates')
                    if not significant(s['p'],5e-8):
                        for k in ('kind','kind_body','kind_point'): s[k]='unknown'
                    classified.append(s); key=(s['chr'],s['pos'],s['kind'],s['snp'])
                    if key not in merged: merged[key]=s.copy()
                    else:
                        merged[key]['assays']=sorted(set(merged[key]['assays']+[alias]))
                        merged[key]['flags']=sorted(set(merged[key]['flags']+s['flags']))
                        if pnum(s['p'])<pnum(merged[key]['p']):
                            for k in ('p','beta','se'): merged[key][k]=s[k]
                    signal_rows.append(dict(prot=row['prot'],dataset=d,assay=alias,CHR=s['chr'],POS=s['pos'],SNP=s['snp'],P=s['p'],BETA=s['beta'],SE=s['se'],
                           GENE=s['gene'],CLASS_ORIGINAL=s['raw_kind'],CLASS=s['kind'],CLASS_BODY=s['kind_body'],CLASS_POINT=s['kind_point'],
                           COORDINATE_SENSITIVE=int(s['coordinate_sensitive']),COORDINATE_STATUS=coord['coordinate_status'],
                           build=file_build,flags=';'.join(s['flags']) or '.'))
                    if significant(s['p'],5e-8):
                        candidates.append(dict(DATASET=d,ASSAY=alias,PROT=row['prot'],CHR=s['chr'],POS=s['pos'],SNP=s['snp'],P=s['p'],KIND=s['kind'],CLASS_STATUS='coordinate_checked' if compatible else 'unknown'))
                a['raw_cis']=int(any(s['raw_kind']=='cis' for s in classified)); a['best_cis']=best(classified,alias)
                a['trans']=sum(s['kind']=='trans' for s in classified); a['cis']=1 if a['best_cis'] else None
                complete_cis=complete_cis_search(meta,intervals,file_build,flank,meta_errors)
                a['cis_search_complete']=bool(compatible and complete_cis)
                a['cis_search_status']='complete_covered_search' if a['cis_search_complete'] else 'unresolved_coordinates_or_build' if not compatible else 'uncertified_summary'
                if a['cis'] is None and compatible:
                    if complete_cis: a['cis']=0
                    elif not meta and d in ('ukb','ckb'):
                        a['cis']=0; a['flags'].append('legacy_cis_window_not_fully_certified')
                    else: a['flags'].append('cis_search_censored_by_legacy_export')
                if a['coordinate_sensitive_rows']:
                    a['flags'].append('cis_window_definition_sensitive')
                if not compatible: a['flags'].append('target_or_build_unresolved')
                if any('identical_alleles_in_id' in s['flags'] for s in classified if s['kind']=='cis'): a['flags'].append('cis_variant_id_requires_allele_verification')
                if any('invalid_effect_or_se' in s['flags'] for s in classified if s['kind']=='cis'): a['flags'].append('cis_effect_invalid')
                assays.append(a)
            valid=[a for a in assays if a['status']=='complete']
            status='not_mapped' if not aliases else 'invalid' if invalid and not valid else 'missing' if not valid else 'partial' if missing or invalid else 'complete'
            signals=sorted(merged.values(),key=lambda s:(chrom_key(s['chr']),s['pos'],s['snp']))
            seen=any(a['cis']==1 for a in valid)
            cis=1 if seen else 0 if valid and all(a['cis']==0 for a in valid) and status=='complete' else None
            cs='observed' if cis==1 else 'not_reported' if cis==0 else 'unavailable' if not valid else 'censored' if intervals else 'unannotated'
            lead=best(signals)
            v=dict(aliases=aliases,status=status,cis=cis,cis_status=cs,raw_cis=int(any(a['raw_cis'] for a in valid)),
                raw_trans=sum(s['raw_kind']=='trans' for s in signals),trans=sum(s['kind']=='trans' for s in signals),
                best_cis=lead,strong=bool(lead and significant(lead['p'],strict[d])),assays=assays,signals=signals,
                corrected_rows=corrected,build=build_info[d]['build'],assay_heterogeneity=seen and any(a['cis']==0 for a in valid),
                assay_support_partial=seen and any(a['cis'] is None or a['status']!='complete' for a in assays),
                cis_search_complete=bool(aliases) and status=='complete' and all(a['cis_search_complete'] for a in assays),
                coordinate_sensitive_rows=sum(a['coordinate_sensitive_rows'] for a in valid))
            p['data'][d]=v
            if corrected: p['flags'].append(flag('coordinate_reclassification_'+d,f'{d}: {corrected} 条信号经坐标重分类'))
            if v['assay_heterogeneity'] or v['assay_support_partial']:
                p['flags'].append(flag('assay_discordance_'+d,f'{d}: 多测定 cis 支持不齐','review')); p['review_by_primary'][d].append('多个测定未全部检出 cis；逐测定核对')
            if status in ('partial','missing','invalid'): p['review_by_primary'][d].append('输入文件缺失或格式不完整')
            if any('cis_variant_id_requires_allele_verification' in a['flags'] for a in valid): p['review_by_primary'][d].append('cis 变异ID的两个等位基因相同，须回GWAS核对')
            if any('cis_effect_invalid' in a['flags'] for a in valid): p['review_by_primary'][d].append('cis 的效应值或标准误无效')
            if any(any(x.startswith('metadata_') for x in a['flags']) for a in valid): p['review_by_primary'][d].append('摘要与完成元数据不匹配，须重建')
            coverage.append(dict(prot=row['prot'],dataset=d,status=status,aliases=';'.join(aliases) or '.',files_read=';'.join(a['alias'] for a in valid) or '.',
                    files_missing=';'.join(missing) or '.',files_invalid=';'.join(invalid) or '.',cis_status=cs,
                    cis_search_complete=int(v['cis_search_complete']),coordinate_status=coord['coordinate_status'],
                    build=build_info[d]['build'],reclassified_rows=corrected))
            for a in assays:
                l=a['best_cis'] or {}; assay_rows.append(dict(prot=row['prot'],dataset=d,assay=a['alias'],status=a['status'],cis=a['cis'],raw_cis=a['raw_cis'],
                      trans=a['trans'],cis_P=l.get('p'),cis_SNP=l.get('snp'),cis_search_complete=int(a['cis_search_complete']),
                      cis_search_status=a['cis_search_status'],coordinate_sensitive_rows=a['coordinate_sensitive_rows'],flags=';'.join(a['flags']) or '.'))
        p['support_n']=sum(v['cis']==1 for v in p['data'].values()); p['strong_n']=sum(v['strong'] for v in p['data'].values())
        p['measured_n']=sum(bool(v['aliases']) for v in p['data'].values()); p['ms_support']=p['data']['scot']['cis']==1
        hard=not identity_verified or any(f['code'] in ('composite_target','noncoding_target','target_identity_unverified') for f in p['flags'])
        hard|=any(v['assay_heterogeneity'] or v['assay_support_partial'] for v in p['data'].values())
        hard|=any('cis_variant_id_requires_allele_verification' in a['flags'] for v in p['data'].values() for a in v['assays'])
        hard|=any('cis_effect_invalid' in a['flags'] or any(f.startswith('metadata_') for f in a['flags']) for v in p['data'].values() for a in v['assays'])
        p['priority']=p['strong_n']>=3 and single and not hard
        p['tier']='A' if p['priority'] else 'B' if p['support_n']>=2 else 'C' if p['support_n']==1 else 'D' if intervals else 'U'
        p['same_lead_pairs']=[]
        for x,y in itertools.combinations(DATASETS,2):
            a,b=p['data'][x]['best_cis'],p['data'][y]['best_cis']
            if a and b and a['snp'].startswith('rs') and a['snp']==b['snp'] and (a['chr'],a['pos'])==(b['chr'],b['pos']):
                p['same_lead_pairs'].append(f'{x}/{y}: {a["snp"]}')
        for d,v in p['data'].items():
            if v['status']=='not_mapped': continue
            if v['cis']!=1 and v['trans']>=trans_rich: p['review_by_primary'][d].append(f'未报告可核对的 cis，有 {v["trans"]} 个已导出 trans 位点')
            if v['cis']==1 and not v['strong']: p['review_by_primary'][d].append('cis 达5e-8，但未达蛋白组多重检验敏感性阈值')
            others=sum(x['cis']==1 for k,x in p['data'].items() if k!=d)
            if v['cis']==0 and others>=2: p['review_by_primary'][d].append(f'本摘要未报告 cis，另 {others} 队列有 cis 区域支持')
            if not single: p['review_by_primary'][d].append('目标基因未唯一解析，暂不进入可信单基因集合')
        proteins.append(p)
    for d in DATASETS:
        datasets[d].update(cis=sum(p['data'][d]['cis']==1 for p in proteins),strict_cis=sum(p['data'][d]['strong'] for p in proteins),
                           coverage=dict(Counter(p['data'][d]['status'] for p in proteins)))
    return dict(proteins=proteins,datasets=datasets,build_info=build_info,strict=strict,resolved=resolved,repairs=repairs,
                file_audit=file_audit,coverage=coverage,assay_rows=assay_rows,signal_rows=signal_rows,candidates=candidates,
                annotation=annotation,coordinate_audit=coordinate_audit,coordinate_window=coordinate_window,map_build=map_build)
def flatten(p):
    r=dict(prot=p['prot'],tier=p['tier'],priority=int(p['priority']),cis_cohorts=p['support_n'],strict_cis_cohorts=p['strong_n'],mapped_cohorts=p['measured_n'],MS_support=int(p['ms_support']),
           supplied_CHR=p['coordinate']['supplied_CHR'],supplied_POS=p['coordinate']['supplied_POS'],coordinate_status=p['coordinate']['coordinate_status'],
           coordinate_source=p['coordinate']['coordinate_source'],coordinate_window=p['coordinate']['coordinate_window'],identity_verified=int(p['identity_verified']),
           shared_anchor_targets=p['coordinate']['shared_anchor_targets'])
    for d,v in p['data'].items():
        r.update({f'cis-{d}':v['cis'],f'trans-{d}':v['trans'],f'raw_cis-{d}':v['raw_cis'],f'cis_P-{d}':(v['best_cis'] or {}).get('p'),
                  f'cis_SNP-{d}':(v['best_cis'] or {}).get('snp'),f'status-{d}':v['cis_status'],f'cis_search_complete-{d}':int(v['cis_search_complete']),
                  f'review-{d}':'; '.join(p['review_by_primary'][d]) or '.'})
    r['same_lead_pairs']='; '.join(p['same_lead_pairs']) or '.'; r['flags']='; '.join(f['label'] for f in p['flags']) or '.'
    return r
def trans_regions(proteins,strict,max_span=1000000,signal_rows=None):
    """Bounded-span coordinate correspondence, NEVER an LD/allelic replication test.

    Filter invalid assay receipts before choosing cohort leads. Per-assay rows
    permit a healthy duplicate to survive without inheriting an invalid probe's
    minimum P or unioned signal flags. With merged data alone, a mixed healthy /
    invalid signal is excluded because its individual contributions are lost.
    """
    raw=defaultdict(list)
    if signal_rows is not None:
        for row in signal_rows:
            if row['CLASS']=='trans': raw[row['prot']].append(row)
    output=[]
    for p in proteins:
        usable={d:{a['alias'] for a in v.get('assays',[]) if a.get('status')=='complete' and
                   not any(f.startswith('metadata_') or f in ('file_dataset_build_conflict','target_or_build_unresolved')
                           for f in a.get('flags',[]))} for d,v in p['data'].items()}
        eligible=defaultdict(list)
        if signal_rows is None:
            for d,v in p['data'].items():
                for s in v['signals']:
                    aliases=set(s.get('assays',[]))
                    if s['kind']=='trans' and aliases and aliases.issubset(usable[d]): eligible[d].append(s)
        else:
            merged=defaultdict(dict)
            for row in raw.get(p['prot'],[]):
                d=row['dataset']; alias=row['assay']
                if alias not in usable.get(d,set()): continue
                s=dict(chr=chrom(row['CHR']),pos=int(row['POS']),snp=row['SNP'],p=row['P'],kind='trans',
                       flags=[f for f in str(row.get('flags','')).split(';') if f and f!='.'],assays=[alias])
                key=(s['chr'],s['pos'],s['snp'])
                if key not in merged[d]: merged[d][key]=s
                else:
                    previous=merged[d][key]
                    previous['assays']=sorted(set(previous['assays']+[alias]))
                    previous['flags']=sorted(set(previous['flags']+s['flags']))
                    if pnum(s['p'])<pnum(previous['p']): previous['p']=s['p']
            for d,signals in merged.items():
                eligible[d]=sorted(signals.values(),key=lambda s:(chrom_key(s['chr']),s['pos'],s['snp']))
        hits=sorted([(d,s) for d in p['data'] for s in eligible[d]],
                    key=lambda ds:(chrom_key(ds[1]['chr']),ds[1]['pos']))
        groups=[]
        for d,s in hits:
            if not groups or s['chr']!=groups[-1][0][1]['chr'] or s['pos']-groups[-1][0][1]['pos']>max_span:
                groups.append([])
            groups[-1].append((d,s))
        p['trans_regions']=[]
        for group in groups:
            by=defaultdict(list)
            for d,s in group: by[d].append(s)
            if len(by)<2: continue
            c=group[0][1]['chr']; lo=min(s['pos'] for d,s in group); hi=max(s['pos'] for d,s in group)
            leads={d:min(v,key=lambda s:pnum(s['p'])) for d,v in by.items()}
            strong=sum(significant(s['p'],strict[d]) for d,s in leads.items())
            pairs=[]
            for a,b in itertools.combinations(by,2):
                ka={(s['snp'],s['chr'],s['pos']) for s in by[a] if re.fullmatch(r'rs\d+',s['snp'])}
                kb={(s['snp'],s['chr'],s['pos']) for s in by[b] if re.fullmatch(r'rs\d+',s['snp'])}
                if ka&kb: pairs.append(a+'/'+b+':'+','.join(sorted(t[0] for t in ka&kb)))
            qc=sorted({f for d,s in group for f in s['flags'] if f in ('identical_alleles_in_id','invalid_effect_or_se','reported_nonsignificant')})
            identity=p.get('identity_verified',False)
            record=dict(prot=p['prot'],CHR=c,START=lo,END=hi,N_COHORTS=len(by),STRICT_COHORTS=strong,
                 MS_SUPPORT=int('scot' in by),TARGET_CIS_TIER=p['tier'],TARGET_IDENTITY_VERIFIED=int(identity),
                 REGION_PRIORITY=int(strong>=3 and not qc and identity),
                 SAME_RSID_PAIRS='; '.join(pairs) or '.',FLAGS=';'.join(qc) or '.',
                 EVIDENCE='bounded-span regional correspondence; direction and colocalization untested')
            for d in DATASETS:
                s=leads.get(d,{})
                record.update({f'P-{d}':s.get('p'),f'SNP-{d}':s.get('snp'),f'POS-{d}':s.get('pos'),
                               f'ASSAYS-{d}':';'.join(sorted({a for h in by.get(d,[]) for a in h['assays']})) or '.'})
            output.append(record); p['trans_regions'].append(record)
    return sorted(output,key=lambda r:(-r['REGION_PRIORITY'],-r['STRICT_COHORTS'],-r['MS_SUPPORT'],r['prot'],chrom_key(r['CHR']),r['START']))
def outputs(out,result,args,annotation):
    out.mkdir(parents=True,exist_ok=True); annotation=result['annotation']; proteins=result['proteins']; tables=[flatten(p) for p in proteins]; fields=list(tables[0])
    for filename,selected in [('cross-check.tsv',tables),('priority.tsv',[r for r,p in zip(tables,proteins) if p['priority']]),
            ('priority.with-MS.tsv',[r for r,p in zip(tables,proteins) if p['priority'] and p['ms_support']]),
            ('review.tsv',[r for r,p in zip(tables,proteins) if any(p['review_by_primary'].values())])]: write_tsv(out/filename,fields,selected)
    for d in DATASETS:
        write_tsv(out/f'primary-{d}.tsv',fields,[r for r,p in zip(tables,proteins) if p['data'][d]['aliases']])
        write_tsv(out/f'review-{d}.tsv',fields,[r for r,p in zip(tables,proteins) if p['review_by_primary'][d]])
    regions=trans_regions(proteins,result['strict'],signal_rows=result['signal_rows'])
    region_fields=list(regions[0]) if regions else ['prot','CHR','START','END','N_COHORTS','STRICT_COHORTS','MS_SUPPORT','REGION_PRIORITY','EVIDENCE']
    write_tsv(out/'trans-region-support.tsv',region_fields,regions)
    write_tsv(out/'trans-region-priority.tsv',region_fields,[r for r in regions if r['REGION_PRIORITY']])
    for key,filename,fallback in [('coverage','coverage.tsv',['prot','dataset','status']),('assay_rows','assays.tsv',['prot','dataset','assay']),
            ('signal_rows','signals.tsv',['prot','dataset','assay']),('file_audit','input-files.tsv',['dataset','assay','file']),
            ('repairs','mapping-repairs.tsv',['prot','dataset','requested','resolved','rule']),
            ('candidates','candidate-loci.tsv',['DATASET','ASSAY','PROT','CHR','POS','SNP','P','KIND','CLASS_STATUS'])]:
        rs=result[key]; write_tsv(out/filename,list(rs[0]) if rs else fallback,rs)
    qc=[r for r in result['signal_rows'] if any(f in r['flags'].split(';') for f in
         ('identical_alleles_in_id','invalid_effect_or_se','reported_nonsignificant'))]
    write_tsv(out/'signal-qc.tsv',list(result['signal_rows'][0]) if result['signal_rows'] else ['prot','dataset','assay','flags'],qc)
    write_tsv(out/'prot.names.resolved.tsv',list(result['resolved'][0]),result['resolved'])
    coordinate_rows=coordinates_used_rows(result)
    write_tsv(out/'coordinates-used.tsv',list(coordinate_rows[0]),coordinate_rows)
    write_tsv(out/'coordinate-source-audit.tsv',list(result['coordinate_audit'][0]),result['coordinate_audit'])
    sensitivity=[r for r in result['signal_rows'] if r['COORDINATE_SENSITIVE']]
    write_tsv(out/'signal-window-sensitivity.tsv',list(result['signal_rows'][0]) if result['signal_rows'] else ['prot','dataset','assay','CLASS_BODY','CLASS_POINT'],sensitivity)
    comparisons=[]
    for p in proteins:
        for d,v in p['data'].items():
            if not v['aliases']: continue
            body=best([{**s,'kind':s['kind_body']} for s in v['signals']]); point=best([{**s,'kind':s['kind_point']} for s in v['signals']])
            comparisons.append(dict(prot=p['prot'],dataset=d,status=v['status'],coordinate_status=p['coordinate']['coordinate_status'],
                BODY_CIS_REPORTED=int(bool(body)),POINT_CIS_REPORTED=int(bool(point)),BEST_BODY_P=(body or {}).get('p'),
                BEST_POINT_P=(point or {}).get('p'),WINDOW_SENSITIVE_ROWS=v['coordinate_sensitive_rows'],
                meaning='0 means no cis lead in supplied summaries under this window; does not certify absence in full GWAS'))
    write_tsv(out/'coordinate-window-comparison.tsv',list(comparisons[0]) if comparisons else ['prot','dataset','BODY_CIS_REPORTED','POINT_CIS_REPORTED'],comparisons)
    export_builds={t['build'] for ts in annotation.values() for t in ts}|{result['map_build']}
    for d in DATASETS:
        for build in sorted(export_builds):
            path=out/'cis-bed'/f'{d}.GRCh{build}.bed'; path.parent.mkdir(exist_ok=True)
            with path.open('w') as f:
                f.write(f'# build={build}; coordinates=1-based-closed; columns=CHR START END ASSAY\n')
                for r in result['resolved']:
                    ts=annotation.get(r['prot'],[])
                    if not ts or any(t['status'] not in ('single_gene','complex_complete') for t in ts): continue
                    for a in r[d].split(';'):
                        if a!='.':
                            for t in ts:
                                if t['build']==build: f.write(f'{t["chr"]}\t{t["start"]}\t{t["end"]}\t{a}\n')
            Path(str(path)+'.grch').write_text(build+'\n')
    bins=defaultdict(set)
    for r in result['signal_rows']:
        if r['CLASS']=='trans': bins[(r['dataset'],r['CHR'],r['POS']//1000000)].add(r['prot'])
    hotspots=[dict(dataset=d,CHR=c,START=b*1000000,END=(b+1)*1000000-1,proteins=len(ps)) for (d,c,b),ps in bins.items()]
    write_tsv(out/'trans-hotspots.tsv',['dataset','CHR','START','END','proteins'],sorted(hotspots,key=lambda r:(-r['proteins'],r['dataset'],chrom_key(r['CHR']))))
    spans={}
    for r in result['signal_rows']: spans[r['CHR']]=max(spans.get(r['CHR'],0),r['POS'])
    axes=[]; offset=0
    for c in sorted(spans,key=chrom_key):
        span=int(spans[c]*1.02)+1; axes.append(dict(chr=c,span=span,offset=offset)); offset+=span
    warnings=['旧版 .sig.txt 每条染色体仅保留最强 trans 代表位点；trans 数不是完整独立 loci 数。',
        'deCODE/SCOT 旧 CLASS 全为 trans：坐标重注释可恢复阳性；未恢复者是 ?，须重建摘要。',
        '0 仅表示可读摘要未报告 cis；— 是未测，? 是不可判定。cis 缺失不证明无效应。',
        '提供的 CHR/POS 按声明的基因组版本解释；替代 contig、缺失或冲突不回填主组装坐标。',
        '默认沿用经提供点核对的基因体±窗口；另导出点±窗口敏感性比较，未报告不等于完整区域阴性。',
        '同一坐标可对应历史符号或复合物成分；坐标相同不能证明不同 assay 测到同一蛋白。',
        '基因组版本有数据集级 rsID 哨兵证据，但不能认证每个无哨兵文件。',
        'cis 区域重复出现不等同方向一致、共定位或因果作用；P=0表示数值下溢。',
        '严格阈值为 5e-8/映射测定数量，仅作保守敏感性分析。']
    summary=dict(datasets=result['datasets'],counts=dict(proteins=len(proteins),priority=sum(p['priority'] for p in proteins),
        priority_with_MS=sum(p['priority'] and p['ms_support'] for p in proteins),four_cohort_cis=sum(p['support_n']==4 for p in proteins),
        review=sum(any(p['review_by_primary'].values()) for p in proteins),trans_regions_2plus=len(regions),
        trans_regions_strict3plus=sum(r['REGION_PRIORITY'] for r in regions),
        coordinate_sensitive_signal_rows=len(sensitivity),coordinate_sensitive_proteins=len({r['prot'] for r in sensitivity}),
        coordinate_sensitive_cis_calls=sum(r['BODY_CIS_REPORTED']!=r['POINT_CIS_REPORTED'] for r in comparisons)),
        coordinates=dict(status_counts=dict(Counter(r['coordinate_status'] for r in result['coordinate_audit'])),
                         identity_verified=sum(r['identity_verified'] for r in result['coordinate_audit']),
                         shared_anchor_targets=sum(r['shared_anchor_targets']!='.' for r in result['coordinate_audit']),
                         window=result['coordinate_window'],map_declared_build=result['map_build']),
        build_audit=result['build_info'],mapping_repairs=result['repairs'],strong_thresholds=result['strict'])
    payload=dict(proteins=proteins,axes=axes,primary=args.primary,generated=datetime.now(timezone.utc).isoformat(timespec='seconds'),root=str(args.gwas_root),
                 summary=summary,warnings=warnings,config=dict(cis_flank=args.cis_flank,p_threshold=5e-8,strong_thresholds=result['strict'],
                 window='corroborated gene body +/- flank' if result['coordinate_window']=='body' else 'target point +/- flank',
                 coordinate_window=result['coordinate_window'],map_declared_build=result['map_build'],platforms=PLATFORMS,trans_region_max_span=1000000))
    write_json(out/'summary.json',summary); write_json(out/'cross-check.data.json',payload)
    template=Path(__file__).with_suffix('.html').read_text()
    body=json.dumps(payload,ensure_ascii=False,separators=(',',':'),allow_nan=False).replace('<','\\u003c')
    template=template.replace('__FINGERPRINT_CSS__',Path(__file__).with_name('fingerprint_view.css').read_text())
    template=template.replace('__FINGERPRINT_SCRIPT__',Path(__file__).with_name('fingerprint_view.js').read_text())
    (out/'cross-check.html').write_text(template.replace('__PAYLOAD__',body),encoding='utf-8')
    write_json(out/'provenance.json',dict(format='phole-cross-check-v3',mapping_sha256=digest(args.protein_map),annotation_sha256=digest(args.annotation),
              sentinels_sha256=digest(args.sentinels),cis_flank=args.cis_flank,trans_rich=args.trans_rich,generated=payload['generated'],
              coordinate_window=result['coordinate_window'],map_declared_build=result['map_build'],
              coordinates_used_sha256=digest(out/'coordinates-used.tsv'),coordinate_source_audit_sha256=digest(out/'coordinate-source-audit.tsv'),
              output_sha256={name:digest(out/name) for name in ('cross-check.data.json','signals.tsv','assays.tsv','input-files.tsv',
                    'coordinates-used.tsv','coordinate-source-audit.tsv','summary.json','prot.names.resolved.tsv',
                    'trans-region-priority.tsv','trans-region-support.tsv')},
              scope='reported locus screening; no allelic harmonization, causal mediation, or colocalization'))
    return summary
def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--gwas-root',type=Path,required=True); ap.add_argument('--protein-map',type=Path)
    ap.add_argument('--out',type=Path,required=True); ap.add_argument('--primary',choices=DATASETS,default='ukb')
    ap.add_argument('--annotation',type=Path); ap.add_argument('--sentinels',type=Path,default=DEFAULT_SENTINELS)
    ap.add_argument('--cis-flank',type=int,default=500000); ap.add_argument('--trans-rich',type=int,default=5); ap.add_argument('--builds',default='')
    ap.add_argument('--coordinate-window',choices=('body','point'),default='body',help='Corroborated gene body or supplied/reference point, before adding cis-flank')
    ap.add_argument('--map-build',choices=('37','38'),default='38',help='Declared build of optional CHR/POS columns; never inferred from a protein point')
    args=ap.parse_args(); args.gwas_root=args.gwas_root.expanduser().resolve()
    args.protein_map=(args.protein_map or args.gwas_root/'prot.names.tsv').expanduser().resolve()
    args.annotation=(args.annotation or args.gwas_root/'prot.annot.tsv').expanduser().resolve(); args.sentinels=args.sentinels.expanduser().resolve()
    args.out=args.out.expanduser().resolve()
    if args.cis_flank<0 or args.trans_rich<1: ap.error('cis-flank must be nonnegative; trans-rich positive')
    builds=dict(x.split('=',1) for x in args.builds.split(',') if x)
    if any(d not in DATASETS or b not in ('37','38') for d,b in builds.items()): ap.error('builds needs dataset=37|38')
    annotation=load_annotation(args.annotation); anchors={r['SNP']:r for r in read_tsv(args.sentinels)}
    result=collect(args.gwas_root,args.protein_map,annotation,anchors,args.cis_flank,args.trans_rich,builds,args.coordinate_window,args.map_build)
    print(json.dumps(outputs(args.out,result,args,annotation),ensure_ascii=False,indent=2)); print(f'Report: {args.out/"cross-check.html"}')
if __name__=='__main__': main()
