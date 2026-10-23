"""Exact GRCh37 allele matching with explicit unavailable/absent/swapped states."""
from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path
import csv, re, math, json, sqlite3, tempfile, os
try:
    from .common import canonical_chrom, open_text, chrom_sort_key
    from .core import dump_json, fingerprint, tsv_writer
except ImportError:
    from common import canonical_chrom, open_text, chrom_sort_key
    from core import dump_json, fingerprint, tsv_writer

VAR_RE=re.compile(r'^(?:chr)?(?P<chr>[0-9XYMT]+):(?P<pos>[0-9]+)[_:](?P<a1>[ACGT]+)[_:](?P<a2>[ACGT]+)$',re.I)
COORD_RE=re.compile(r'^(?:chr)?([0-9XYMT]+):([0-9]+)$',re.I)
COMPLEMENT=str.maketrans('ACGT','TGCA')

def _float(value,label,lo,hi,missing=False):
    if missing and value in {'.','NA','nan','-1'}: return None
    try: v=float(value)
    except ValueError as exc: raise ValueError(f'Invalid MFI {label}: {value}') from exc
    if not math.isfinite(v) or not lo<=v<=hi: raise ValueError(f'Invalid MFI {label}: {value}')
    return v

def parse_mfi_row(fields,chrom_hint=None):
    """7-col, 8-col (minor allele), or local 9-col split alternate ID.

    The local first real row verifies order; never silently skip malformed data.
    A 7-column MAF cannot identify ALT AF. Both minor/MAF orders are accepted
    only when the allele token is unambiguous.
    """
    n=len(fields)
    if n==9:
        coord,pair,rsid,pos,a1,a2,v1,v2,info=fields
        alt_id=coord+'_'+pair
    elif n==8: alt_id,rsid,pos,a1,a2,v1,v2,info=fields
    elif n==7:
        alt_id,rsid,pos,a1,a2,maf,info=fields; v1=maf; v2=None
    else: raise ValueError(f'Expected 7, 8 or 9 MFI columns, found {n}')
    a1=a1.upper(); a2=a2.upper()
    if not re.fullmatch('[ACGT]+',a1) or not re.fullmatch('[ACGT]+',a2) or a1==a2:
        raise ValueError('Invalid MFI alleles')
    try: pos=int(pos)
    except ValueError as exc: raise ValueError('Invalid MFI position') from exc
    if pos<1: raise ValueError('Invalid MFI position')
    minor=None
    if v2 is not None:
        if v1.upper() in {a1,a2}: minor=v1.upper(); maf=v2
        elif v2.upper() in {a1,a2}: minor=v2.upper(); maf=v1
        elif v2 in {'.','NA'}: maf=v1
        elif v1 in {'.','NA'}: maf=v2
        else: raise ValueError('MFI minor allele disagrees with REF/ALT')
    else: maf=v1
    maf=_float(maf,'MAF',0,.5,True); info=_float(info,'INFO',-1,1,True)
    match=VAR_RE.fullmatch(alt_id)
    coord=COORD_RE.fullmatch(alt_id)
    inferred=match.group('chr') if match else coord.group(1) if coord else None
    if chrom_hint is None and inferred is None: raise ValueError('MFI chromosome requires filename hint')
    chrom=canonical_chrom(chrom_hint or inferred)
    id_alleles=None
    if match:
        id_alleles=(match.group('a1').upper(),match.group('a2').upper())
        if canonical_chrom(match.group('chr'))!=chrom or int(match.group('pos'))!=pos or id_alleles!=(a1,a2):
            raise ValueError(f'MFI alternate_id disagrees with coordinates/alleles: {alt_id}')
    elif coord and (canonical_chrom(coord.group(1))!=chrom or int(coord.group(2))!=pos):
        raise ValueError('MFI alternate_id coordinate mismatch')
    alt_af=None if minor is None or maf is None else maf if minor==a2 else 1-maf
    return dict(chrom=chrom,pos=pos,ref=a1,alt=a2,maf=maf,info=info,minor_allele=minor,
        alt_af=alt_af,alternate_id=alt_id,rsid=rsid,id_alleles=id_alleles,columns=n)

def _records(path,chrom):
    with open_text(path) as f:
        for ln,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'): continue
            fields=line.split()
            if fields[0].lower() in {'alternate_id','alternate_ids'}: continue
            try: yield parse_mfi_row(fields,chrom)
            except ValueError as exc: raise ValueError(f'{path}:{ln}: {exc}') from exc

def find_mfi(directory,c):
    names=[f'ukb_mfi_chr{c}_v3.txt']
    if c.isdigit(): names.append(f'ukb_mfi_chr{int(c):02d}_v3.txt')
    for name in names:
        for suffix in ('','.gz'):
            p=Path(directory)/(name+suffix)
            if p.is_file(): return p
    return None

def target_rows(path):
    with open_text(path) as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'chrom37','pos37','ref37'}.issubset(reader.fieldnames or []):
            raise ValueError('MFI targets require explicit GRCh37 chrom37,pos37,ref37 columns')
        for row in reader:
            c=canonical_chrom(row['chrom37']); pos=int(row['pos37']); ref=row['ref37'].upper()
            if row.get('alt37') not in (None,'','.','NA'):
                yield c,pos,ref,row['alt37'].upper(),'anchor',row
            if row.get('extra_alt') not in (None,'','.','NA'):
                strand=row.get('strand37to38')
                if strand not in {'+','-'}: raise ValueError('extra_alt needs strand37to38 for GRCh37 matching')
                alt=row['extra_alt'].upper()
                if strand=='-': alt=alt.translate(COMPLEMENT)
                yield c,pos,ref,alt,'extra',row

def scan_mfi_dir(mfi_dir,target_tsv,out_tsv,info_thresholds=(.3,.8)):
    """Disk-backed U0 targets; preserves every contrast/role and every absence.

    Matching is pos+REF+ALT. A swapped row is flagged but never exact-present.
    Missing chromosome files mean unavailable, not absent. Output includes both
    MAF and ALT AF where a minor-allele column identifies frequency orientation.
    """
    out=Path(out_tsv); out.parent.mkdir(parents=True,exist_ok=True)
    marker=Path(str(out)+'.manifest.json')
    if marker.exists():
        old=json.loads(marker.read_text())
        if old.get('schema')!=2 or old.get('inputs',{}).get('targets')!=fingerprint(target_tsv) or old.get('mfi_directory')!=str(Path(mfi_dir).resolve()) or old.get('info_thresholds')!=list(info_thresholds):
            raise ValueError('UKB MFI inputs/settings changed; choose a new output')
        for c,audit in old['chromosomes'].items():
            path=find_mfi(mfi_dir,c)
            if audit['status']=='unavailable_file':
                if path is not None: raise ValueError('Previously unavailable MFI chromosome now exists; choose a new output')
            elif path is None or fingerprint(path)!=audit['source']: raise ValueError('MFI source changed; choose a new output')
        if not out.exists() or fingerprint(out)!=old.get('output'): raise ValueError('Completed MFI output changed')
        return old
    with tempfile.TemporaryDirectory(prefix='mfi-',dir=out.parent) as tmp:
        db=sqlite3.connect(str(Path(tmp)/'targets.sqlite'))
        db.execute('PRAGMA temp_store=FILE'); db.execute('PRAGMA cache_size=-32000')
        db.execute('CREATE TABLE target (chrom TEXT,pos INT,ref TEXT,alt TEXT,role TEXT,payload TEXT)')
        batch=[]; fields=None; target_count=0
        for c,p,r,a,role,row in target_rows(target_tsv):
            if fields is None: fields=list(row)
            if len(r)!=1 or len(a)!=1 or r==a or r not in 'ACGT' or a not in 'ACGT': raise ValueError('U0 target is not an exact SNV contrast')
            batch.append((c,p,r,a,role,json.dumps(row))); target_count+=1
            if len(batch)>=10000: db.executemany('INSERT INTO target VALUES (?,?,?,?,?,?)',batch); batch=[]
        if batch: db.executemany('INSERT INTO target VALUES (?,?,?,?,?,?)',batch)
        db.execute('CREATE INDEX target_key ON target(chrom,pos,ref,alt)')
        db.execute('CREATE TABLE matches (chrom TEXT,pos INT,ref TEXT,alt TEXT,status TEXT,payload TEXT,PRIMARY KEY(chrom,pos,ref,alt))')
        db.commit()
        qc=Counter(); audits={}; missing=[]
        for (c,) in db.execute('SELECT DISTINCT chrom FROM target ORDER BY chrom'):
            path=find_mfi(mfi_dir,c)
            if path is None: missing.append(c); audits[c]=dict(status='unavailable_file'); continue
            want={(p,r,a) for p,r,a in db.execute('SELECT DISTINCT pos,ref,alt FROM target WHERE chrom=?',(c,))}
            positions={k[0] for k in want}; found={}; first=None; n=0
            for x in _records(path,c):
                if first is None: first=x
                n+=1
                if x['pos'] not in positions: continue
                for key,status in [((x['pos'],x['ref'],x['alt']),'present'),((x['pos'],x['alt'],x['ref']),'swapped_alleles')]:
                    if key not in want: continue
                    if key in found and found[key][0]=='present':
                        if status=='present': raise ValueError(f'Duplicate exact MFI key: {c}:{key}')
                        continue
                    found[key]=(status,x)
            if first is None: raise ValueError(f'Empty MFI file: {path}')
            audits[c]=dict(status='scanned',source=fingerprint(path),first_real_row=first,records=n)
            db.executemany('INSERT INTO matches VALUES (?,?,?,?,?,?)',[(c,*key,status,json.dumps(x)) for key,(status,x) in found.items()]); db.commit()
            print(f'[UKB MFI] chr{c}: {n:,} rows, {len(found):,}/{len(want):,} target keys located',flush=True)
        added=['target_role','ukb_chrom','ukb_pos','ukb_ref','ukb_alt','ukb_status','ukb_allele_swap','ukb_maf','ukb_info','ukb_minor_allele','ukb_alt_af','ukb_alternate_id','ukb_rsid','gnomad_nfe_af','af_difference']+[f'ukb_info_ge_{str(t).replace(".","p")}' for t in info_thresholds]
        fields=(fields or ['chrom37','pos37','ref37','alt37'])+added
        temp=out.with_name(out.name+'.partial')
        sums=[0,0.,0.,0.,0.,0.]
        with temp.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); writer.writeheader()
            query='SELECT t.chrom,t.pos,t.ref,t.alt,t.role,t.payload,m.status,m.payload FROM target t LEFT JOIN matches m USING(chrom,pos,ref,alt) ORDER BY t.rowid'
            for c,p,r,a,role,payload,status,match in db.execute(query):
                row=json.loads(payload); x=json.loads(match) if match else {}
                status=status or ('unavailable_file' if c in missing else 'absent_from_MFI')
                qc[status]+=1
                nfe=row.get('extra_nfe_af') if role=='extra' else None
                try: nfe=float(nfe)
                except (ValueError,TypeError): nfe=None
                if nfe is not None and (not math.isfinite(nfe) or not 0<=nfe<=1): nfe=None
                alt_af=x.get('alt_af') if status=='present' else None
                diff=alt_af-nfe if alt_af is not None and nfe is not None else None
                row.update(dict(target_role=role,ukb_chrom=c,ukb_pos=p,ukb_ref=x.get('ref'),ukb_alt=x.get('alt'),ukb_status=status,
                    ukb_allele_swap=int(status=='swapped_alleles') if match else None,
                    ukb_maf=x.get('maf'),ukb_info=x.get('info'),ukb_minor_allele=x.get('minor_allele'),ukb_alt_af=alt_af,
                    ukb_alternate_id=x.get('alternate_id'),ukb_rsid=x.get('rsid'),gnomad_nfe_af=nfe,af_difference=diff))
                for t in info_thresholds: row[f'ukb_info_ge_{str(t).replace(".","p")}']=int(x['info']>=t) if x.get('info') is not None else None
                if diff is not None and (x.get('info') or 0)>=.8:
                    sums[0]+=1; sums[1]+=nfe; sums[2]+=alt_af; sums[3]+=nfe*nfe; sums[4]+=alt_af*alt_af; sums[5]+=nfe*alt_af
                writer.writerow({k:'NA' if v is None else v for k,v in row.items()})
        os.replace(temp,out); db.close()
    n,sx,sy,sxx,syy,sxy=sums; denom=(n*sxx-sx*sx)*(n*syy-sy*sy)
    result=dict(schema=2,target_rows=target_count,matched_rows=qc['present'],counts=dict(qc),chromosomes=audits,
        inputs=dict(targets=fingerprint(target_tsv)),mfi_directory=str(Path(mfi_dir).resolve()),info_thresholds=list(info_thresholds),output=fingerprint(out),concordance=dict(n=n,info_floor=.8,
        pearson_r=(n*sxy-sx*sy)/math.sqrt(denom) if n>1 and denom>0 else None),
        caveats=['Absence from MFI is not extinction.','Swapped alleles are not exact matches.','Seven-column MAF has unknown ALT direction.','Concordance conditional on exact orientation and INFO>=0.8.'])
    dump_json(str(out)+'.manifest.json',result); return result
