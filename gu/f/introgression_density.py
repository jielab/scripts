#!/usr/bin/env python3
"""Cache per-individual genomic coverage, merging overlapping IBDmix references.

Fractions are lineage-specific union bp / diploid physical span. Unphased calls do
not resolve homozygous tract dosage. Male non-PAR X has a haploid denominator.
Missing chromosomes and untested female X rows are NA, never a negative call.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, re, sqlite3, tempfile, shutil, atexit
from pathlib import Path
import numpy as np
from comm import CHROM_LENGTHS, write_tsv_rows

SCHEMA = 10
SUMMARY_REFERENCES = ('Altai', 'Chagyr', 'Vindija', 'Denisova', 'Denisova25')


def whole_chromosome_run(text, chrom):
    """Recognize both native-VCF and cached modern_source run records."""
    records=[line.replace('\\t','\t').split('\t') for line in text.splitlines()]
    if any(row[0] in {'loci_file','loci'} for row in records):
        return False
    if any(row[0]=='modern_vcf' for row in records):
        return True
    source=[row[2:] for row in records if len(row)>=4 and row[:2]==['modern_source','chr'+str(chrom)]]
    scopes=[row[1] for row in source if row[0]=='loci']
    if scopes:
        return all(scope=='whole_chromosome' for scope in scopes)
    # Native VCF source manifests omit a loci row for whole chromosomes.
    return any(row[:2]==['format','vcf'] for row in source)

def neanderthal_summary(con, dataset, build, samples, lengths, targets):
    """Keep the Altai cache used by the map and Cell 2020 comparison."""
    return [dict(row, neanderthal_bp=row['archaic_bp']) for row in
            reference_summary(con, dataset, build, samples, lengths, targets, 'Altai')]

def reference_summary(con, dataset, build, samples, lengths, targets, reference):
    """Exact single-reference autosomal union; targets certify whole-chromosome runs.

    Untested individuals remain missing, including when provenance is unavailable.
    Do not infer a denominator from positive calls or extrapolate partial genomes.
    """
    burden = dict.fromkeys(samples, 0)
    autosomes=set(map(str,range(1,23)))
    sources = ('Chagyr', 'Chagyrskaya') if reference == 'Chagyr' else (reference,)
    placeholders = ','.join('?' for _ in sources)
    rows = con.execute(f"SELECT sample_id,chr,start,end FROM segments WHERE dataset_id=? AND genome_build=? AND method='ibdmix' AND source IN ({placeholders}) ORDER BY sample_id,chr,start,end", (dataset, build, *sources))
    for sample, chrom, left, right in merged_intervals(rows):
        if chrom in autosomes and sample in burden and sample in targets.get(chrom, set()):
            burden[sample] += max(0, min(lengths[chrom], right) - max(0, left))
    result = []
    for sample in samples:
        chroms = [ch for ch in map(str, range(1, 23)) if sample in targets.get(ch, set())]
        physical = sum(lengths[ch] for ch in chroms)
        denominator = 2 * physical
        result.append(dict(sample_id=sample, chromosomes=','.join(chroms), n_chromosomes=len(chroms),
                           tested_bp=denominator, haploid_bp=physical, archaic_bp=burden[sample] if denominator else '',
                           coverage_pct=100 * burden[sample] / denominator if denominator else ''))
    return result

def merged_intervals(rows):
    """Input sorted by sample, chromosome, start, end; yield union intervals."""
    key = None; start = end = 0
    for sample, chrom, left, right in rows:
        if right <= left: continue
        current = (sample, chrom)
        if current == key and left <= end:
            end = max(end, right)
        else:
            if key is not None: yield (*key, start, end)
            key, start, end = current, left, right
    if key is not None: yield (*key, start, end)

def add_interval(covered, row, offset, length, bin_bp, left, right):
    left=max(0,left);right=min(length,right)
    if right<=left:return
    for b in range(left//bin_bp,(right-1)//bin_bp+1):
        covered[row,offset+b]+=min(right,(b+1)*bin_bp)-max(left,b*bin_bp)

def prepare(database, output, sample_panel=None, bin_bp=5_000_000):
    stamp=database.stat()
    signature={'schema':SCHEMA,'database':str(database.resolve()),'mtime_ns':stamp.st_mtime_ns,'size':stamp.st_size,'bin_bp':bin_bp,
               'sample_panel':str(sample_panel) if sample_panel else '', 'panel_mtime_ns':sample_panel.stat().st_mtime_ns if sample_panel and sample_panel.exists() else 0}
    version=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    output.mkdir(parents=True,exist_ok=True)
    con=sqlite3.connect(f'file:{database.resolve()}?mode=ro',uri=True)
    units=con.execute("SELECT DISTINCT dataset_id,genome_build FROM method_runs WHERE method='ibdmix' AND status='complete'").fetchall()
    missing = any(not (output/re.sub(r'[^A-Za-z0-9_.-]','_',dataset+'.'+build)/version/'manifest.json').exists() for dataset,build in units)
    local_database = None
    if missing:
        # Random indexed reads on /mnt/* are much slower than one sequential
        # copy followed by local SQLite reads. Never mutate the source database.
        fd,name=tempfile.mkstemp(prefix='gu-density-',suffix='.sqlite');os.close(fd)
        local_database=Path(name);atexit.register(lambda:local_database.unlink(missing_ok=True))
        print('DENSITY staging SQLite for local interval aggregation',flush=True)
        shutil.copyfile(database,local_database)
        current=database.stat()
        if (current.st_mtime_ns,current.st_size)!=(stamp.st_mtime_ns,stamp.st_size):raise RuntimeError('database changed while preparing density; retry')
        con.close();con=sqlite3.connect(local_database)
    for dataset, build in units:
        lengths=CHROM_LENGTHS.get(build.replace('GRCh',''))
        if not lengths:continue
        slug=re.sub(r'[^A-Za-z0-9_.-]','_',dataset+'.'+build)
        root=output/slug; dest=root/version; pointer=root/'current.tsv'
        if (dest/'manifest.json').exists():
            write_tsv_rows(pointer,['directory'],[{'directory':version}]);print(f'DENSITY cached {dataset} {build}',flush=True);continue
        root.mkdir(parents=True,exist_ok=True)
        stage=Path(tempfile.mkdtemp(prefix='.building-',dir=root))
        metadata={r[0]:(r[1] or 'UNKNOWN',r[2] or 'UNKNOWN') for r in con.execute('SELECT sample_id,population,super_population FROM sample_populations WHERE dataset_id=?',(dataset,))}
        called_pairs=list(con.execute("SELECT DISTINCT sample_id,chr FROM segments WHERE dataset_id=? AND genome_build=? AND method='ibdmix'",(dataset,build)))
        called={r[0] for r in called_pairs}
        runs=con.execute("SELECT chr,evidence_eligible,raw_file,availability_note FROM method_runs WHERE dataset_id=? AND genome_build=? AND method='ibdmix' AND status='complete'",(dataset,build)).fetchall()
        tested={r[0] for r in runs};x_male=False;targets={};psam_cache={};summary_targets={};summary_refs=set();neanderthal_tested=set()
        lineage_targets={'Neanderthal':{},'Denisovan':{}}
        reference_targets={ref:{} for ref in SUMMARY_REFERENCES}
        # Use the database's provenance snapshot. Raw metadata may already have
        # changed during a rerun while this database still contains old calls.
        legacy_chromosomes=sorted({r[0] for r in runs if not str(r[3]).startswith('audited IBDmix ')})
        profiles=sorted({str(r[3]).split('profile=',1)[1] for r in runs if 'profile=' in str(r[3])})
        for chrom,_,raw,note in runs:
            p=Path(raw);p=p if p.is_absolute() else database.parent/p
            if not p.exists():continue
            text=p.read_text().replace("\\t", "\t")
            run_targets=set()
            if chrom=='X' and 'male_haploid_nonpar' in text:x_male=True
            for line in text.splitlines():
                fields=line.replace("\\t","\t").split("\t")
                if 'psam' not in fields:continue
                i=fields.index('psam')
                if i+1>=len(fields):continue
                psam=Path(fields[i+1].rsplit(':',2)[0])
                if not psam.exists():continue
                if psam not in psam_cache:
                    with psam.open() as f:
                        header=f.readline().split();ix=header.index('IID') if 'IID' in header else header.index('#IID')
                        psam_cache[psam]={r[ix] for r in (line.split() for line in f) if len(r)>ix}
                targets.setdefault(chrom,set()).update(psam_cache[psam])
                run_targets.update(psam_cache[psam])
            # Native VCF runs retain the actual caller sample roster here.
            roster=p.parent/'samples'/('male.txt' if chrom=='X' and 'male_haploid_nonpar' in text else 'ALL.txt')
            if not run_targets and roster.exists():
                run_targets={line.strip() for line in roster.read_text().splitlines() if line.strip()}
                targets.setdefault(chrom,set()).update(run_targets)
            unit=('X_MALE' if 'male_haploid_nonpar' in text else 'X_PAR') if chrom=='X' else 'C'+chrom
            actual_roster=p.parent/'samples'/unit/'ALL.txt'
            if actual_roster.exists():
                run_targets=set(actual_roster.read_text().split())
                targets[chrom]=run_targets
            refs=next((line.split('\t',1)[1].split() for line in text.splitlines() if line.startswith('refs\t')), [])
            # Locus-only calls cannot certify an entire chromosome denominator.
            whole=whole_chromosome_run(text,chrom)
            neand_refs=set(refs)&{'Altai','Chagyr','Chagyrskaya','Vindija'}
            if neand_refs:neanderthal_tested.add(chrom)
            # Background-only Denisova runs export no Denisovan ancestry calls.
            # Their absence in segments is unmeasured, not zero coverage.
            background_only=(any(line.split('\t')[:2]==['background_filter','1'] for line in text.splitlines())
                             and not any(line.split('\t')[:2]==['export_denisovan','1'] for line in text.splitlines()))
            if whole and run_targets and chrom in map(str,range(1,23)):
                for ref in refs:
                    ref = 'Chagyr' if ref == 'Chagyrskaya' else ref
                    if ref not in reference_targets or (background_only and ref.startswith('Denisova')):
                        continue
                    reference_targets[ref].setdefault(chrom,set()).update(run_targets)
            for lineage, has_refs in [('Neanderthal',bool(neand_refs)),('Denisovan',not background_only and any('denis' in ref.lower() for ref in refs))]:
                if has_refs and whole and run_targets:
                    lineage_targets[lineage].setdefault(chrom,set()).update(run_targets)
            if chrom in map(str,range(1,23)) and whole and 'Altai' in neand_refs and run_targets:
                summary_targets.setdefault(chrom,set()).update(run_targets)
                summary_refs.add('Altai')
        tested &= neanderthal_tested  # a Denisova-only run did not test Neanderthal ancestry
        # Use the actual caller target, not the larger annotation panel.
        samples=sorted(called|set().union(*targets.values()),key=lambda s:(metadata.get(s,('UNKNOWN','UNKNOWN'))[1],metadata.get(s,('UNKNOWN','UNKNOWN'))[0],s))
        sample_index={s:i for i,s in enumerate(samples)}
        for chrom in tested:
            if chrom not in targets:targets[chrom]={sample for sample,c in called_pairs if c==chrom}
        bins=[];offsets={}
        for chrom,length in lengths.items():
            offsets[chrom]=len(bins)
            for left in range(0,length,bin_bp):
                bins.append(dict(bin_index=len(bins),chr=chrom,start=left,end=min(length,left+bin_bp),tested=int(chrom in tested)))
        print(f'DENSITY building {dataset} {build}: {len(samples)} individuals x {len(bins)} bins',flush=True)
        widths=np.array([b['end']-b['start'] for b in bins])
        ploidy=np.array([1 if b['chr']=='X' and x_male else 2 for b in bins])
        unions={}
        for lineage, filename in [('Neanderthal','matrix.tsv.gz'),('Denisovan','denisovan_matrix.tsv.gz')]:
            covered=np.zeros((len(samples),len(bins)),dtype=np.int64)
            rows=con.execute("SELECT sample_id,chr,start,end FROM segments WHERE dataset_id=? AND genome_build=? AND method='ibdmix' AND source_class=? ORDER BY sample_id,chr,start,end",(dataset,build,lineage))
            unions[lineage]=0
            for sample,chrom,left,right in merged_intervals(rows):
                if chrom not in lengths:continue
                add_interval(covered,sample_index[sample],offsets[chrom],lengths[chrom],bin_bp,left,right);unions[lineage]+=1
            if np.any(covered>widths):raise ValueError('union coverage exceeds bin length')
            density=np.round(100*covered/(widths*ploidy),3)
            for j,b in enumerate(bins):
                eligible=lineage_targets[lineage].get(b['chr'],set())
                for i,sample in enumerate(samples):
                    if sample not in eligible:density[i,j]=np.nan
            with gzip.open(stage/filename,'wt') as f:np.savetxt(f,density,delimiter='\t',fmt='%.3f')
        write_tsv_rows(stage/'samples.tsv',['row_index','sample_id','population','super_population'],[
            dict(row_index=i,sample_id=s,population=metadata.get(s,('UNKNOWN','UNKNOWN'))[0],super_population=metadata.get(s,('UNKNOWN','UNKNOWN'))[1]) for i,s in enumerate(samples)])
        write_tsv_rows(stage/'bins.tsv',['bin_index','chr','start','end','tested'],bins)
        print(f'DENSITY summarizing Neanderthal autosomes: {len(summary_targets)}/22 chromosomes',flush=True)
        reference_rows=[]
        for ref in SUMMARY_REFERENCES:
            print(f'DENSITY summarizing {ref} autosomes: {len(reference_targets[ref])}/22 chromosomes',flush=True)
            reference_rows.extend(dict(row,reference=ref) for row in
                                  reference_summary(con,dataset,build,samples,lengths,reference_targets[ref],ref))
        write_tsv_rows(stage/'archaic_summary.tsv', ['reference','sample_id','chromosomes','n_chromosomes','tested_bp','haploid_bp','archaic_bp','coverage_pct'],reference_rows)
        write_tsv_rows(stage/'neanderthal_summary.tsv', ['sample_id','chromosomes','n_chromosomes','tested_bp','haploid_bp','neanderthal_bp','coverage_pct'],
                       [dict(row,neanderthal_bp=row['archaic_bp']) for row in reference_rows if row['reference']=='Altai'])
        filters=[]
        if con.execute("SELECT 1 FROM sqlite_master WHERE name='ibdmix_filter_runs'").fetchone():
            filters=[dict(chrom=c,status=s,daf_status=d) for c,s,d in con.execute('SELECT chr,status,daf_status FROM ibdmix_filter_runs WHERE dataset_id=? AND genome_build=?',(dataset,build))]
        (stage/'manifest.json').write_text(json.dumps(dict(signature,dataset=dataset,build=build,n_samples=len(samples),n_bins=len(bins),union_intervals=unions,x_male_only=x_male,
            legacy_chromosomes=legacy_chromosomes,profiles=profiles,diploid_autosome_bp=2*sum(lengths.get(str(ch),0) for ch in range(1,23)),
            lineages=['Neanderthal','Denisovan'],definition='Per-lineage union bp / diploid bin span (male non-PAR X: haploid); percent; untested entries are missing',
            summary_definition='Altai-only union bp / (2 x physical length of certified whole autosomes); percent; unphased calls cannot resolve homozygous dosage',summary_references=sorted(summary_refs),
            archaic_summary_references=list(SUMMARY_REFERENCES),
            archaic_summary_definition='Separate reference union bp on certified whole autosomes; no cross-reference summation; untested entries are missing',ibdmix_filters=filters),indent=2))
        os.replace(stage,dest)
        write_tsv_rows(pointer,['directory'],[{'directory':version}])
        print(f'DENSITY ready {dest}',flush=True)
    con.close()
    if local_database is not None:local_database.unlink(missing_ok=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--database',required=True,type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--sample-panel',type=Path);ap.add_argument('--bin-bp',type=int,default=5_000_000)
    args=ap.parse_args()
    if args.bin_bp<=0:ap.error('--bin-bp must be positive')
    prepare(args.database,args.output or args.database.parent/'normalize/density',args.sample_panel,args.bin_bp)
