#!/usr/bin/env python3
"""Collapse COJO signals, lift complete regions and leads, emit GRCh37 BED + audit."""
import math
import argparse, csv, hashlib, json, os, re, shutil, subprocess, tempfile, urllib.request
from pathlib import Path
from collections import defaultdict
from comm import CHROM_LENGTHS, parse_bp


def rows(p):
    with open(p) as f: return list(csv.DictReader(f, delimiter='\t'))


def write(p, data, fields):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)


def lift(binary, chain, source, output, unmapped):
    subprocess.run([str(binary),'-minMatch=0.95','-multiple',str(source),str(chain),str(output),str(unmapped)],check=True)
    hits=defaultdict(list)
    for line in output.read_text().splitlines():
        a=line.split('\t');hits[a[3]].append(a)
    return hits


def bounded_target_loci(signals, output, args):
    """Reapply the identical size cap after liftOver; retain every mapped signal."""
    target = output / 'grch37_collapse'
    source = output / 'signals.GRCh37.tsv'
    write(source, signals, list(signals[0]))
    subprocess.run(['Rscript', str(Path(__file__).with_name('cojo_collapse.R')),
                    str(args.helper), str(source), str(target), args.size,
                    args.chr_column, args.pos_column, args.p_column], check=True)
    membership = rows(target / 'members.tsv')
    result = rows(target / 'collapsed.tsv')
    for r in result:
        r['start'], r['end'], r['lead_pos'] = (int(float(r[k])) for k in ('start', 'end', 'lead_pos'))
        r['end'] = min(r['end'], CHROM_LENGTHS['37'][r['chr']])
        r['source_loci'] = ','.join(sorted({m['GU_SOURCE_LOCUS'] for m in membership if m['locus_id'] == r['locus_id']}))
        assert 0 < r['end'] - r['start'] <= parse_bp(args.size)
        assert r['start'] < r['lead_pos'] <= r['end']
    for left, right in zip(result, result[1:]):
        assert left['chr'] != right['chr'] or left['end'] <= right['start']
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--build',choices=['37','38'],required=True);ap.add_argument('--size',default='1Mb')
    ap.add_argument('--chr-column',default='Chr');ap.add_argument('--pos-column',default='bp');ap.add_argument('--p-column',default='pJ')
    ap.add_argument('--helper',type=Path,default=Path('/mnt/d/scripts/0f/0phe.f.R'))
    ap.add_argument('--chain',type=Path,default=Path('/mnt/d/files/liftOver/hg38ToHg19.over.chain.gz'))
    # PATH may contain an older UCSC executable with different multi-map behavior.
    ap.add_argument('--liftover',default=os.environ.get('GU_LIFTOVER','/mnt/d/software/bin/liftOver'))
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    lines=[line.split() for line in a.input.read_text().splitlines() if line.strip()]
    required={'SNP',a.chr_column,a.pos_column,a.p_column}
    if not lines or not required.issubset(lines[0]):
        raise ValueError('COJO header missing required columns: '+','.join(sorted(required)))
    original=[]; conflicts=[]
    for i,fields in enumerate(lines[1:],1):
        signal=dict(zip(lines[0],fields));name=signal.get('SNP','')
        ch=re.sub(r'^chr','',signal.get(a.chr_column,''),flags=re.I).upper()
        ch='X' if ch=='23' else ch
        reason=''
        try:
            pos=float(signal.get(a.pos_column,''));pval=float(signal.get(a.p_column,''))
            if len(fields)!=len(lines[0]):reason='wrong_column_count'
            elif ch not in CHROM_LENGTHS[a.build]:reason='unsupported_chromosome'
            elif not math.isfinite(pos) or pos!=int(pos) or not 1<=pos<=CHROM_LENGTHS[a.build][ch]:reason='invalid_position'
            elif not math.isfinite(pval) or not 0<=pval<=1:reason='invalid_P'
            elif not name or name.upper() in ('.','NA','NAN'):reason='missing_SNP_ID'
            else:
                m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)',name,re.I)
                if ':' in name and not m:reason='malformed_coordinate_SNP_ID'
                elif m and (('X' if m[1]=='23' else m[1].upper())!=ch or int(m[2])!=int(pos)):
                    reason='SNP_ID_disagrees_with_CHR_POS'
                elif m and int(pos)+len(m[3])-1>CHROM_LENGTHS[a.build][ch]:reason='allele_outside_chromosome'
        except (ValueError,TypeError,OverflowError):
            reason='invalid_numeric_field'
        if reason:
            conflicts.append(dict(input_row=i,SNP=name,CHR=ch,POS=signal.get(a.pos_column,''),action='skip',reason=reason))
            print(f"[GU COJO] SKIP input_row={i} SNP={name} reason={reason}",flush=True)
            continue
        signal[a.chr_column]=ch
        signal['GU_ORIGINAL_SNP']=name
        signal['GU_ORIGINAL_INPUT_ROW']=i
        original.append(signal)
    source=a.output/'cojo.normalized.tsv'
    write(a.output/'input_coordinate.audit.tsv',conflicts,['input_row','SNP','CHR','POS','action','reason'])
    if not original:raise ValueError('No valid COJO SNPs remain after filtering; see input_coordinate.audit.tsv')
    write(source,original,list(original[0]))
    subprocess.run(['Rscript',str(Path(__file__).with_name('cojo_collapse.R')),str(a.helper),str(source),str(a.output),a.size,a.chr_column,a.pos_column,a.p_column],check=True)
    data=rows(a.output/'collapsed.tsv'); limits=CHROM_LENGTHS[a.build]
    for r in data:
        if r['chr'] not in limits:raise ValueError('Unsupported chromosome '+r['chr'])
        r['start']=int(float(r['start']));r['end']=min(int(float(r['end'])),limits[r['chr']]);r['lead_pos']=int(float(r['lead_pos']));r['lead_p']=float(r['lead_p']);r['lead_log_p']=float(r['lead_log_p'])
        if not r['start']<r['lead_pos']<=r['end']:raise ValueError('Lead outside chromosome: '+r['lead_snp'])
    # Collapsed BED is the source-build region, clipped to chromosome bounds.
    with (a.output/f'loci.GRCh{a.build}.source.bed').open('w') as f:
        for r in data:f.write(f"{r['chr']}\t{r['start']}\t{r['end']}\t{r['lead_snp']}\n")
    regions={}; leads={}; member_hits={}; liftover_info={}
    members=rows(a.output/'members.tsv'); width=parse_bp(a.size)
    if a.build=='38':
        if not a.chain.is_file():raise FileNotFoundError(f'GRCh38 -> GRCh37 chain missing: {a.chain}')
        binary=Path(a.liftover)
        if not binary.exists():
            binary.parent.mkdir(parents=True,exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=binary.parent,delete=False) as f:
                tmp=Path(f.name)
                try:
                    with urllib.request.urlopen('https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver',timeout=120) as src:shutil.copyfileobj(src,f)
                except Exception:tmp.unlink(missing_ok=True);raise
            tmp.chmod(0o755);tmp.replace(binary)
        liftover_info=dict(liftover_binary=str(binary.resolve()),liftover_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
        print('[GU COJO] liftOver='+liftover_info['liftover_binary']+' sha256='+liftover_info['liftover_sha256'],flush=True)
        with (a.output/'regions.hg38.bed').open('w') as f, (a.output/'leads.hg38.bed').open('w') as g:
            for r in data:
                m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)',r['lead_snp'],re.I)
                if m and (m[1].replace('23','X')!=r['chr'] or int(m[2])!=r['lead_pos']):raise ValueError('SNP ID disagrees with CHR/POS: '+r['lead_snp'])
                length=len(m[3]) if m else 1
                f.write(f"chr{r['chr']}\t{r['start']}\t{r['end']}\t{r['locus_id']}\t0\t+\n")
                g.write(f"chr{r['chr']}\t{r['lead_pos']-1}\t{r['lead_pos']-1+length}\t{r['locus_id']}\t0\t+\n")
        with (a.output/'members.hg38.bed').open('w') as f:
            for m in members:
                ch=str(m[a.chr_column]).removeprefix('chr');ch='X' if ch=='23' else ch
                pos=int(float(m[a.pos_column])); allele=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)',m['SNP'],re.I)
                length=len(allele[3]) if allele else 1
                if allele and (allele[1].removeprefix('chr') not in (ch, '23' if ch=='X' else ch) or int(allele[2])!=pos):raise ValueError('SNP ID disagrees with CHR/POS: '+m['SNP'])
                f.write(f"chr{ch}\t{pos-1}\t{pos-1+length}\tsignal_{m['input_row']}\t0\t+\n")
        member_hits=lift(binary,a.chain,a.output/'members.hg38.bed',a.output/'members.hg19.bed',a.output/'members.unmapped.bed')
        regions=lift(binary,a.chain,a.output/'regions.hg38.bed',a.output/'regions.hg19.bed',a.output/'regions.unmapped.bed')
        leads=lift(binary,a.chain,a.output/'leads.hg38.bed',a.output/'leads.hg19.bed',a.output/'leads.unmapped.bed')
    accepted=[]; audit=[]; target_signals=[]
    for r in data:
        status='accepted'; dest=dict(r);strand='+'; mapping_method='identity' if a.build=='37' else 'whole_region_liftover'
        if a.build=='38':
            rr=regions[r['locus_id']];ll=leads[r['locus_id']]
            if not rr and len(ll)==1:
                group=[m for m in members if m['locus_id']==r['locus_id']]
                hits=[member_hits['signal_'+m['input_row']] for m in group]
                if hits and all(len(h)==1 and h[0][0]==ll[0][0] and h[0][5]==ll[0][5] for h in hits):
                    ch=ll[0][0].removeprefix('chr')
                    if ch in CHROM_LENGTHS['37']:
                        lo=max(0,min(int(h[0][1]) for h in hits)-width//2)
                        hi=min(CHROM_LENGTHS['37'][ch],max(int(h[0][1]) for h in hits)-width//2+width)
                        rr=[[ll[0][0],str(lo),str(hi),r['locus_id'],'0',ll[0][5]]]
                        mapping_method='reconstructed_from_all_lifted_signals'
            if len(rr)!=1 or len(ll)!=1:status='unmapped_or_multiple_mapping'
            else:
                x,y=rr[0],ll[0];ch=x[0].removeprefix('chr');start,end=int(x[1]),int(x[2]);pos=int(y[1])+1;strand=y[5]
                if ch not in CHROM_LENGTHS['37'] or ch!=r['chr'] or x[0]!=y[0] or x[5]!=y[5] or not start<=int(y[1])<int(y[2])<=end:status='inconsistent_region_and_lead_mapping'
                elif not .5 <= (end-start)/(r['end']-r['start']) <= 2:status='large_length_change'
                else:
                    m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)',r['lead_snp'],re.I)
                    name=r['lead_snp']
                    if m:
                        if int(y[2])-int(y[1])!=len(m[3]):status='lead_allele_span_changed'
                        ref,alt=m[3].upper(),m[4].upper()
                        if strand=='-':ref=ref.translate(str.maketrans('ACGT','TGCA'))[::-1];alt=alt.translate(str.maketrans('ACGT','TGCA'))[::-1]
                        name=f'{ch}:{pos}:{ref}:{alt}'
                    elif not re.fullmatch(r'rs\d+',name):
                        # Unknown coordinate labels must not masquerade as GRCh37 markers.
                        name=f"cojo_{r['locus_id']}_GRCh37"
                    dest.update(chr=ch,start=start,end=end,lead_pos=pos,lead_snp=name)
        group_signals=[]
        if status=='accepted':
            for member in (m for m in members if m['locus_id']==r['locus_id']):
                signal=dict(member)
                signal['GU_SOURCE_LOCUS']=r['locus_id']
                signal['GU_INPUT_ROW']=member.get('GU_ORIGINAL_INPUT_ROW',member['input_row'])
                signal.pop('locus_id'); signal.pop('input_row')
                if a.build=='38':
                    hits=member_hits['signal_'+member['input_row']]
                    if len(hits)!=1 or hits[0][0].removeprefix('chr')!=dest['chr'] or hits[0][5]!=strand:
                        status='member_unmapped_or_inconsistent'; break
                    hit=hits[0]; pos=int(hit[1])+1
                    if not 1<=pos<=CHROM_LENGTHS['37'][dest['chr']]:
                        status='member_outside_chromosome'; break
                    name=member['SNP']
                    allele=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)',name,re.I)
                    if allele:
                        if int(hit[2])-int(hit[1])!=len(allele[3]):
                            status='member_allele_span_changed'; break
                        ref,alt=allele[3].upper(),allele[4].upper()
                        if strand=='-':
                            ref=ref.translate(str.maketrans('ACGT','TGCA'))[::-1]
                            alt=alt.translate(str.maketrans('ACGT','TGCA'))[::-1]
                        name=f"{dest['chr']}:{pos}:{ref}:{alt}"
                    elif not re.fullmatch(r'rs\d+',name):
                        name=f"cojo_signal_{member['input_row']}_GRCh37"
                    signal.update({a.chr_column:dest['chr'],a.pos_column:pos,'SNP':name})
                group_signals.append(signal)
        audit.append(dict(source_locus=r['locus_id'],source_lead=r['lead_snp'],source_chr=r['chr'],source_start=r['start'],source_end=r['end'],status=status,mapping_method=mapping_method,target_lead=dest['lead_snp'] if status=='accepted' else '',target_chr=dest['chr'] if status=='accepted' else '',target_start=dest['start'] if status=='accepted' else '',target_end=dest['end'] if status=='accepted' else '',strand=strand))
        if status!='accepted':print(f"[GU COJO] SKIP locus={r['locus_id']} lead={r['lead_snp']} reason={status}",flush=True)
        if status=='accepted':
            dest['source_loci']=r['locus_id'];accepted.append(dest);target_signals.extend(group_signals)
    write(a.output/'liftover.audit.tsv',audit,list(audit[0]))
    if not target_signals:raise ValueError('No loci survived coordinate conversion; inspect liftover.audit.tsv')
    merged=bounded_target_loci(target_signals,a.output,a)
    with (a.output/'loci.GRCh37.bed').open('w') as f:
        for r in merged:f.write(f"{r['chr']}\t{r['start']}\t{r['end']}\t{r['lead_snp']}\n")
    write(a.output/'loci.GRCh37.tsv',merged,list(merged[0]))
    manifest=dict(converter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),helper_sha256=hashlib.sha256(a.helper.read_bytes()).hexdigest(),source=str(a.input.resolve()),source_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),source_build=a.build,analysis_build='37',collapse_size=a.size,p_column=a.p_column,n_signals=sum(int(r['n_snps']) for r in data),n_source_loci=len(data),n_rejected=len(data)-len(accepted),n_analysis_loci=len(merged),n_reconstructed=sum(r['status']=='accepted' and r['mapping_method']=='reconstructed_from_all_lifted_signals' for r in audit),audit_directory=str(a.output.resolve()),chain=str(a.chain) if a.build=='38' else None,chain_sha256=hashlib.sha256(a.chain.read_bytes()).hexdigest() if a.build=='38' else None,n_input_signals=len(lines)-1,n_skipped_input_signals=len(conflicts),n_mapped_signals=len(target_signals),max_locus_bp=max(r['end']-r['start'] for r in merged),rule='strongest unassigned SNP first; centred windows clipped at assigned boundaries; width <= collapse_size in BOTH builds; regroup all accepted mapped signals in GRCh37; zero pJ ranked by joint Wald statistic; BED is 0-based half-open')
    manifest.update(liftover_info)
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest),flush=True)
    if manifest['n_rejected']:print('WARNING: rejected loci are excluded; inspect '+str(a.output/'liftover.audit.tsv'),flush=True)

if __name__=='__main__':main()
