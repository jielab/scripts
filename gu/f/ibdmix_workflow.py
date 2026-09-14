#!/usr/bin/env python3
"""IBDmix preparation and final filtering following the authors' Cell workflow.

The caller remains the unmodified upstream executable. BED files here are
0-based half-open; native IBDmix coordinates are converted exactly once.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import csv
import fcntl
import gzip
import hashlib
import http.client
import json
import os
from pathlib import Path
import re
import ssl
import subprocess
import tempfile
import time
import urllib.error
import urllib.request
import zlib
import numpy as np
from comm import CHROM_LENGTHS, write_tsv_rows

VERSION = '2026-09-12.2'
AFRICAN_CONTROLS = {'ESN', 'GWD', 'LWK', 'MSL', 'YRI'}
NEANDERTHALS = {'Altai', 'Chagyr', 'Chagyrskaya', 'Vindija'}

def opener(path):
    return gzip.open(path, 'rt') if str(path).endswith(('.gz', '.bgz')) else open(path)

def fingerprint(path):
    p = Path(path).resolve(); s = p.stat()
    return [str(p), s.st_size, s.st_mtime_ns]

def file_checksum(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with open(path,'rb') as handle:
        while block:=handle.read(1024*1024):h.update(block)
    return h.hexdigest()

def sha256_file(path):
    return file_checksum(path)

class IncompleteReferenceError(ValueError):
    """A response did not contain the complete requested reference."""


def download(url, path, attempts=5, chunk_bytes=4*1024*1024, expected_md5=None):
    """Validate and atomically cache references, retrying interrupted transfers.

    Fetch bounded ranges when supported so large AXT transfers do not depend
    on one long-lived connection. Resume only with an ETag/Last-Modified
    validator, a matching Content-Range and the original total size. A server
    ignoring Range starts a fresh file. A publisher checksum also pins content
    across mirror backends and allows partial files to survive interruptions.
    Attempts count consecutive failures.
    """
    if attempts < 1:
        raise ValueError('Reference download attempts must be positive')
    if expected_md5 is not None and not re.fullmatch(r'[0-9a-f]{32}', expected_md5):
        raise ValueError('Invalid publisher MD5 checksum')
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(str(path)+'.lock', 'w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        provenance = Path(str(path)+'.source.json')
        if path.exists() and path.stat().st_size:
            if expected_md5 and file_checksum(path, 'md5') != expected_md5:
                raise ValueError(f'Publisher checksum mismatch: {path}')
            if provenance.exists() and sha256_file(path)!=json.loads(provenance.read_text())['sha256']:
                raise ValueError(f'Reference checksum mismatch: {path}')
            if not provenance.exists() and path.suffix in {'.gz', '.bgz'}:
                with gzip.open(path, 'rb') as check:
                    while check.read(1024*1024): pass
            return path
        part = Path(str(path)+('.part' if expected_md5 else f'.part.{os.getpid()}'))
        meta_part = Path(str(provenance)+f'.part.{os.getpid()}')
        total = validator = None
        def publish():
            if path.suffix in {'.gz', '.bgz'}:
                try:
                    with gzip.open(part, 'rb') as check:
                        while check.read(1024*1024): pass
                except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
                    part.unlink(missing_ok=True)
                    raise IncompleteReferenceError(f'Invalid gzip: {exc}') from exc
            if expected_md5 and file_checksum(part, 'md5') != expected_md5:
                part.unlink(missing_ok=True)
                raise IncompleteReferenceError('Downloaded file does not match publisher MD5')
            meta = {'url':url,'sha256':sha256_file(part),'size_bytes':part.stat().st_size}
            if expected_md5:
                meta['publisher_md5'] = expected_md5
            meta_part.write_text(json.dumps(meta,indent=2))
            os.replace(part, path)
            os.replace(meta_part, provenance)
            return path
        try:
            # Interruption can occur after the last byte, before publication.
            # A complete retained part needs validation, not an invalid Range.
            if expected_md5 and part.exists() and file_checksum(part, 'md5') == expected_md5:
                return publish()
            attempt = 1
            while attempt <= attempts:
                offset = part.stat().st_size if part.exists() and (validator or expected_md5) else 0
                headers = {'Accept-Encoding': 'identity'}
                if chunk_bytes:
                    headers['Range'] = f'bytes={offset}-{offset+chunk_bytes-1}'
                elif offset:
                    headers['Range'] = f'bytes={offset}-'
                if offset and not expected_md5:
                    headers['If-Range'] = validator
                print(f'IBDMIX reference download: attempt={attempt}/{attempts} offset={offset} {url}', flush=True)
                try:
                    request = urllib.request.Request(url, headers=headers)
                    with urllib.request.urlopen(request, timeout=120) as response:
                        status = response.status
                        response_length = response.headers.get('Content-Length')
                        response_length = int(response_length) if response_length is not None else None
                        if status == 206:
                            match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                            if not match:
                                raise IncompleteReferenceError('Missing or invalid Content-Range')
                            start, end, size = map(int, match.groups())
                            if start != offset or not start <= end < size or (offset and total is not None and total != size):
                                raise IncompleteReferenceError('Content-Range does not match the requested reference')
                            if response_length is not None and response_length != end-start+1:
                                raise IncompleteReferenceError('Content-Length does not match Content-Range')
                            total = size
                        elif status == 200:
                            # Range may be unsupported, or If-Range may detect a
                            # new resource. Never append a complete HTTP 200 body.
                            offset = 0
                            total = response_length
                        else:
                            raise IncompleteReferenceError(f'Unexpected HTTP status: {status}')
                        etag = response.headers.get('ETag')
                        current_validator = etag if etag and not etag.startswith('W/') else response.headers.get('Last-Modified')
                        if offset and not expected_md5 and current_validator and current_validator != validator:
                            part.unlink(missing_ok=True)
                            total = validator = None
                            raise IncompleteReferenceError('Reference changed during resumed download')
                        validator = current_validator or (validator if offset else None)
                        with part.open('ab' if offset else 'wb') as out:
                            while True:
                                try:
                                    block = response.read(64*1024)
                                except http.client.IncompleteRead as exc:
                                    out.write(exc.partial)
                                    raise
                                if not block:
                                    break
                                out.write(block)
                    received = part.stat().st_size
                    if status == 206 and received == end+1 and received < total:
                        if not validator and not expected_md5:
                            # Without a validator the next range might belong
                            # to a different version; request a full body.
                            part.unlink(missing_ok=True)
                            chunk_bytes = total = None
                        attempt = 1
                        continue
                    if not received or (total is not None and received != total):
                        if total is not None and received > total:
                            part.unlink(missing_ok=True)
                            validator = None
                        raise IncompleteReferenceError(f'Incomplete body: received={received} expected={total} bytes')
                    return publish()
                except (urllib.error.URLError, http.client.HTTPException, TimeoutError, ConnectionError, IncompleteReferenceError) as exc:
                    if isinstance(exc, urllib.error.HTTPError):
                        if exc.code == 416 and expected_md5:
                            # An oversized/stale retained part must restart.
                            part.unlink(missing_ok=True)
                            total = validator = None
                        elif exc.code not in {408, 429, 500, 502, 503, 504}:
                            raise
                    if isinstance(getattr(exc, 'reason', None), ssl.SSLCertVerificationError):
                        raise
                    if attempt == attempts:
                        raise RuntimeError(f'Reference download failed after {attempts} consecutive attempts: {url}; {exc}. No incomplete file was published; rerun to retry.') from exc
                    delay = min(2**attempt, 30)
                    print(f'IBDMIX reference retry: {exc}; waiting {delay}s', flush=True)
                    time.sleep(delay)
                    attempt += 1
        finally:
            if not expected_md5:
                part.unlink(missing_ok=True)
            meta_part.unlink(missing_ok=True)
    return path

def samples(panel, actual, output, background_required=True):
    metadata = {}
    with open(panel) as handle:
        header = handle.readline().lower().split()
        def index(names):
            return next((header.index(n) for n in names if n in header), None)
        sample_col=index(['sample','sample_id','iid','#iid']); pop_col=index(['pop','population'])
        super_col=index(['super_pop','super_population','superpop'])
        if sample_col is None or pop_col is None:
            raise ValueError('IBDmix requires sample and pop columns; provide --sample-panel with population labels.')
        for line in handle:
            fields=line.split()
            if not fields: continue
            sample=fields[sample_col]; pop=fields[pop_col]
            if sample in metadata: raise ValueError(f'Duplicate metadata sample: {sample}')
            if not re.fullmatch(r'[A-Za-z0-9_-]+',pop) or pop.upper() in {'NA','NONE','UNKNOWN','ALL'}:
                raise ValueError(f'Missing or unsafe population label for {sample}: {pop}')
            metadata[sample]=(pop, fields[super_col] if super_col is not None else 'UNKNOWN')
    ids=Path(actual).read_text().split()
    if len(set(ids))!=len(ids) or not ids: raise ValueError('Target sample IDs must be nonempty and unique')
    groups={}
    for sample in ids:
        if sample not in metadata: raise ValueError(f'Target sample missing from metadata: {sample}')
        pop, superpop=metadata[sample]; groups.setdefault(pop,[]).append(sample)
    if background_required and not AFRICAN_CONTROLS.issubset(groups):
        raise ValueError('Cell 2020 filtering needs all five African control populations: ESN GWD LWK MSL YRI. Use a complete 1KG cohort, or an explicitly configured custom profile.')
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for pop, members in sorted(groups.items()):
        if len(members)<10: raise ValueError(f'Population {pop} has {len(members)} targets; IBDmix requires at least 10 individuals.')
        path=output/f'{pop}.txt'; path.write_text(''.join(s+'\n' for s in members))
        supers={metadata[s][1] for s in members}
        if len(supers)!=1: raise ValueError(f'Conflicting super-population labels for {pop}')
        rows.append(dict(population=pop,super_population=next(iter(supers)),n=len(members),sample_file=str(path.resolve())))
    write_tsv_rows(output/'populations.tsv',['population','super_population','n','sample_file'],rows)
    Path(output/'ALL.txt').write_text(''.join(s+'\n' for s in ids))
    return rows

def bed_intervals(path, chrom, length, columns=(0,1,2)):
    with opener(path) as handle:
        for line in handle:
            if line.startswith(('#','track','browser')) or not line.strip(): continue
            f=line.split()
            if f[columns[0]].removeprefix('chr')!=chrom: continue
            lo,hi=int(f[columns[1]]),int(f[columns[2]])
            if not 0<=lo<hi<=length: raise ValueError(f'Invalid BED interval in {path}: {line.strip()}')
            yield lo,hi

def write_boolean_bed(path, chrom, mask):
    # Emit runs without allocating one Python object per masked base.
    edges=np.flatnonzero(np.diff(mask.astype(np.int8),prepend=0,append=0))
    with open(path,'w') as out:
        for lo,hi in zip(edges[::2],edges[1::2]): out.write(f'{chrom}\t{lo}\t{hi}\n')

def cpg_sites(reference, variants, alignments, chrom):
    """Potential CpGs in hg19, 1KG alleles and the three ancestral alignments.

    Same union-of-C/G rule as upstream generate_cpg_mask.py, vectorized per
    chromosome. AXT target starts are 1-based inclusive.
    """
    sequence=np.frombuffer(reference.upper(),dtype='S1')
    can_c=sequence==b'C'; can_g=sequence==b'G'
    with opener(variants) as handle:
        for line in handle:
            ch,pos,ref,alt=line.split()
            if ch.removeprefix('chr')!=chrom: continue
            i=int(pos)-1
            if not 0<=i<len(sequence): raise ValueError('CpG variant outside chromosome')
            if ref not in {'N','M'} and sequence[i] not in (ref.encode(),b'N',b'M'):
                raise ValueError(f'CpG variant/reference mismatch at chr{chrom}:{pos}')
            if alt=='C':can_c[i]=True
            if alt=='G':can_g[i]=True
    for path in alignments:
        with opener(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):continue
                header=line.split(); human=handle.readline().strip().upper(); other=handle.readline().strip().upper()
                if len(header)!=9 or len(human)!=len(other):raise ValueError(f'Malformed AXT: {path}')
                if header[1].removeprefix('chr')!=chrom:raise ValueError(f'Wrong AXT chromosome: {path}')
                h=np.frombuffer(human.encode(),dtype='S1');o=np.frombuffer(other.encode(),dtype='S1')
                keep=h!=b'-'; mapped=o[keep]; start=int(header[2])-1;end=start+len(mapped)
                if end!=int(header[3]) or not 0<=start<end<=len(sequence):raise ValueError(f'Invalid AXT coordinates: {path}')
                can_c[start:end] |= mapped==b'C';can_g[start:end] |= mapped==b'G'
    pairs=can_c[:-1]&can_g[1:]
    result=np.zeros(len(sequence),dtype=bool);result[:-1]|=pairs;result[1:]|=pairs
    return result

def axt_directory():
    annot = Path(os.environ.get('GU_ANNOT_ROOT', '/mnt/i/annot'))
    return Path(os.environ.get('IBDMIX_AXT_DIR', str(annot/'axt/37')))


@contextmanager
def cached_mask_artifacts(root, signature, names):
    """Serialize builders and publish complete, checksum-verified cache entries."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
    destination = root/key
    with (root/f'{key}.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            meta = json.loads((destination/'cache.json').read_text())
            valid = meta['signature'] == signature and all(
                (destination/name).is_file() and
                sha256_file(destination/name) == meta['sha256'][name]
                for name in names)
        except (OSError, ValueError, KeyError, TypeError):
            valid = False
        if valid:
            yield destination, None
            return
        with tempfile.TemporaryDirectory(prefix=f'.{key}.', dir=root) as temporary:
            staging = Path(temporary)
            yield destination, staging
            checksums = {name:sha256_file(staging/name) for name in names}
            (staging/'cache.json').write_text(json.dumps(
                dict(signature=signature, sha256=checksums), indent=2)+'\n')
            destination.mkdir(exist_ok=True)
            for name in [*names, 'cache.json']:
                os.replace(staging/name, destination/name)


def cached_cpg_sites(fasta, variants, alignments, chrom, length, root):
    """Persist the cohort-independent CpG bitmap in compact packed-bit form."""
    signature = dict(schema=1, chrom=chrom, length=length,
                     code=sha256_file(__file__),
                     inputs=[fingerprint(p) for p in [fasta, str(fasta)+'.fai', variants, *alignments]])
    with cached_mask_artifacts(Path(root)/'cpg'/f'chr{chrom}', signature,
                               ['cpg.exclude.packed.npy']) as (cache, staging):
        if staging is None:
            packed = np.load(cache/'cpg.exclude.packed.npy', allow_pickle=False)
            if packed.dtype != np.uint8 or packed.shape != ((length+7)//8,):
                raise ValueError(f'Invalid cached CpG bitmap: {cache}')
            print(f'IBDMIX chr{chrom}: reusing permanent CpG cache {cache}', flush=True)
            return np.unpackbits(packed, count=length).astype(bool)
        with open(str(fasta)+'.fai') as handle:
            contigs = {line.split()[0]:int(line.split()[1]) for line in handle}
        contig = chrom if chrom in contigs else 'chr'+chrom
        if contigs.get(contig) != length:
            raise ValueError('Reference FASTA chromosome length does not match GRCh37')
        fa = subprocess.check_output(['samtools', 'faidx', str(fasta), contig])
        sequence = b''.join(fa.splitlines()[1:])
        del fa
        if len(sequence) != length:
            raise ValueError('Incomplete FASTA sequence')
        mask = cpg_sites(sequence, variants, alignments, chrom)
        np.save(staging/'cpg.exclude.packed.npy', np.packbits(mask), allow_pickle=False)
        print(f'IBDMIX chr{chrom}: built permanent CpG cache {cache}', flush=True)
    return mask


def link_cached_masks(cache, output, names):
    """Analysis directories contain disposable links, never the only BED copy."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.mask-links.', dir=output) as temporary:
        for name in names:
            link = Path(temporary)/name
            link.symlink_to((cache/name).resolve())
            os.replace(link, output/name)


def axt_reference(chrom, species, root=None):
    filename = f'chr{chrom}.hg19.{species}.synNet.axt.gz'
    url = f'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/vs{species[0].upper()+species[1:]}/syntenicNet/{filename}'
    # Explicit roots retain the flat-directory interface. The default cache
    # separates species, matching the user's 37/vsPanTro2 layout.
    directory = Path(root or axt_directory())
    if root is None and not os.environ.get('IBDMIX_AXT_DIR'):
        directory /= f'vs{species[0].upper()+species[1:]}'
    return url, directory/filename


def download_axt(chrom, species, root=None):
    url, path = axt_reference(chrom, species, root)
    checksums = path.parent/'md5sum.txt'
    if not checksums.is_file():
        checksums = download(url.rsplit('/', 1)[0]+'/md5sum.txt', path.parent/'checksums'/f'{species}.md5sum.txt')
    expected = {}
    for line in checksums.read_text().splitlines():
        digest, name = line.split()
        expected[Path(name.lstrip('*')).name] = digest
    if path.name not in expected:
        raise ValueError(f'AXT absent from publisher checksum list: {path.name}')
    return download(url, path, expected_md5=expected[path.name])


def prepare_masks(root, chrom, modern, fasta, upstream, output, axt_root=None):
    """Build hg19/GRCh37 masks, with an explicit non-blocking X extension.

    map35_50 mappability is already present in the published minimal masks;
    recomputing and intersecting that identical mask is unnecessary.
    """
    length=CHROM_LENGTHS['37'][chrom]
    root=Path(root); output=Path(output);output.mkdir(parents=True,exist_ok=True)
    assets=root/'sources'
    # This published accessibility BED contains autosomes only. Applying its
    # empty X subset would exclude all of X, independently of the AXT issue.
    strict = None if chrom == 'X' else download('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed',assets/'1kg.strict_mask.autosomes.bed')
    dups=download('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz',assets/'genomicSuperDups.txt.gz')
    minimal={}
    for ref,name in [('Altai','AltaiNea'),('Denisova','DenisovaPinky')]:
        filename=f'{name}.map35_50.MQ30.Cov.indels.TRF.bed.bgz'
        minimal[ref]=download('https://bioinf.eva.mpg.de/altai_minimal_filters/'+filename,assets/filename)
    alignments=[]
    species_list = ['panTro2','ponAbe2','rheMac2']
    missing_axt = [axt_reference(chrom, species, axt_root)[1] for species in species_list
                   if not axt_reference(chrom, species, axt_root)[1].is_file()
                   or not axt_reference(chrom, species, axt_root)[1].stat().st_size]
    cpg_enabled = chrom != 'X' or not missing_axt
    if cpg_enabled:
        for species in species_list:
            alignments.append(download_axt(chrom, species, axt_root))
        print(f'IBDMIX chr{chrom}: CpG filter enabled; verified AXT: '+', '.join(map(str,alignments)), flush=True)
    else:
        print('IBDMIX chrX: CpG filter SKIPPED (X AXT unavailable); continuing X analysis with archaic minimal, segmental-duplication and modern-indel masks.', flush=True)
    variants=Path(upstream)/'j.cell.2020.01.012-workflow/abridged_variants.gz'
    if cpg_enabled and not variants.exists():raise ValueError(f'Upstream CpG variant resource missing: {variants}')
    source_meta=Path(os.environ.get('GU_TARGET_TMP_DIR','/nonexistent'))/f'chr{chrom}'/'source.tsv'
    modern_signature=sha256_file(source_meta) if source_meta.exists() else fingerprint(modern)
    inputs = [dups,*minimal.values()]
    if strict is not None:inputs.append(strict)
    if cpg_enabled:inputs.extend([*alignments,variants,fasta])
    signature=dict(version=VERSION,chrom=chrom,modern=modern_signature,cpg_enabled=cpg_enabled,inputs=[fingerprint(p) for p in inputs])
    # Final masks also depend on modern indels; do not share them by chromosome alone.
    signature['cache_schema'] = 1
    signature['code'] = sha256_file(__file__)
    cache_root = Path(axt_root or axt_directory())/'ibdmix-cache'
    analysis_output = output
    names = [f'{ref}.exclude.bed' for ref in minimal]+['manifest.json']
    with cached_mask_artifacts(cache_root/'masks'/f'chr{chrom}', signature, names) as (cache, staging):
        if staging is not None:
            output = staging
            marker = output/'manifest.json'
            if cpg_enabled:
                exclude=cached_cpg_sites(fasta,variants,alignments,chrom,length,cache_root)
            else:
                exclude=np.zeros(length,dtype=bool)
            for lo,hi in bed_intervals(dups,chrom,length,(1,2,3)):exclude[lo:hi]=True
            with subprocess.Popen(['bcftools','query','-f','%CHROM\t%POS\t%REF\t%ALT\n',str(modern)],stdout=subprocess.PIPE,text=True) as proc:
                for line in proc.stdout:
                    ch,pos,ref,alt=line.split();pos=int(pos)
                    if ch.removeprefix('chr')!=chrom:raise ValueError('Modern VCF must contain only the requested chromosome')
                    if len(ref)>1:exclude[max(0,pos-6):min(length,pos+5+len(ref))]=True
                    if len(alt)>1:exclude[max(0,pos-6):min(length,pos+5)]=True
                if proc.wait():raise RuntimeError('bcftools failed while building indel mask')
            accessible=np.ones(length,dtype=bool) if strict is None else np.zeros(length,dtype=bool)
            if strict is not None:
                for lo,hi in bed_intervals(strict,chrom,length):accessible[lo:hi]=True
            accessible &= ~exclude
            totals={}
            for ref,path in minimal.items():
                allowed=np.zeros(length,dtype=bool)
                for lo,hi in bed_intervals(path,chrom,length):allowed[lo:hi]=True
                allowed &= accessible
                totals[ref]=int(allowed.sum())
                if not totals[ref]:raise ValueError(f'No callable bases after masking: {ref} chr{chrom}')
                write_boolean_bed(output/f'{ref}.exclude.bed','23' if chrom=='X' else chrom,~allowed)
            components=['archaic map35_50/MQ30/coverage/indel/TRF','segmental duplications','modern indels +/-5bp']
            if strict is not None:components.append('1KG strict accessibility')
            if cpg_enabled:components.append('CpG: hg19 + upstream modern variants + panTro2/ponAbe2/rheMac2')
            marker.write_text(json.dumps(dict(signature=signature,profile='cell2020_chrX_extension' if chrom=='X' else 'cell2020',callable_bp=totals,cpg_filter='applied' if cpg_enabled else 'skipped_missing_chrX_axt',strict_accessibility='applied' if strict is not None else 'unavailable_for_chrX',mask_semantics='excluded; BED0; native X contig=23',components=components),indent=2))
    link_cached_masks(cache, analysis_output, names)
    print(f'IBDMIX chr{chrom}: permanent masks {cache}', flush=True)


def union(intervals):
    result=[]
    for lo,hi in sorted(intervals):
        if result and lo<=result[-1][1]:result[-1]=(result[-1][0],max(hi,result[-1][1]))
        else:result.append((lo,hi))
    return result

def subtract(lo,hi,mask):
    import bisect
    i=max(0,bisect.bisect_right(mask,(lo,float('inf')))-1)
    for j in range(i,len(mask)):
        left,right=mask[j]
        if left>=hi:break
        if right<=lo:continue
        if left>lo:yield lo,left
        lo=max(lo,right)
        if lo>=hi:return
    if lo<hi:yield lo,hi

FIELDS=['ID','chrom','start','end','length','slod','sites','positive_lods','negative_lods','sample_set','super_pop','anc','locus_id','genome_build','parent_start','parent_end','score_scope']
def finalize(raw_dir, populations, refs, output, chrom, build, locus, core_start=0, core_end=0, min_bp=50000, lod=4, background=True):
    with open(populations) as handle:groups=list(csv.DictReader(handle,delimiter='\t'))
    calls=[];control=[]
    for ref in refs:
        for pop in groups:
            if background and ref=='Denisova' and pop['population'] not in AFRICAN_CONTROLS:continue
            path=Path(raw_dir)/f'{ref}.{pop["population"]}.raw.txt.gz'
            members=set(Path(pop['sample_file']).read_text().split())
            with opener(path) as handle:
                for row in csv.DictReader(handle,delimiter='\t'):
                    if row['ID'] not in members or row['chrom'].removeprefix('chr') not in {chrom,'23' if chrom=='X' else chrom}:
                        raise ValueError(f'Wrong population sample or chromosome in {path}')
                    start=int(row['start'])-1;end=int(row['end'])-1
                    score=float(row.get('LOD',row.get('slod','nan')))
                    if start<0 or end<=start:raise ValueError(f'Invalid native IBDmix interval: {path}')
                    if not np.isfinite(score):raise ValueError(f'Invalid LOD score in {path}')
                    if score<lod or end-start<min_bp:continue
                    if ref=='Denisova' and pop['population'] in AFRICAN_CONTROLS:control.append((start,end))
                    if ref in NEANDERTHALS or (ref=='Denisova' and not background):calls.append((ref,pop,row,start,end,score))
    mask=union(control) if background else []
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(output,'wt') as handle:
        writer=csv.DictWriter(handle,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
        for ref,pop,row,start,end,score in calls:
            for lo,hi in subtract(start,end,mask):
                if hi-lo<min_bp:continue  # authors reapply length filter after subtraction
                if core_end:
                    lo=max(lo,core_start);hi=min(hi,core_end)
                if hi<=lo:continue
                writer.writerow(dict(ID=row['ID'],chrom=chrom,start=lo,end=hi,length=hi-lo,slod=score,
                    sites=row.get('sites',''),positive_lods=row.get('positive_lods',''),negative_lods=row.get('negative_lods',''),
                    sample_set=pop['population'],super_pop=pop['super_population'],anc=ref,locus_id=locus,genome_build=build,
                    parent_start=start,parent_end=end,score_scope='parent_native_call'))
    Path(str(output)+'.afr_denisovan.bed').write_text(''.join(f'{chrom}\t{lo}\t{hi}\n' for lo,hi in mask))

def validate_genotypes(path, output):
    rows=ref_rows=0;last=None
    stream = gzip.open if str(path).endswith(('.gz', '.bgz')) else open
    with stream(path,'rb') as handle:
        header=handle.readline().strip().split(b'\t');expected=len(header)-4
        if expected<2 or len(set(header[4:]))!=expected:raise ValueError('Invalid genotype sample header')
        for line in handle:
            fields=line.split(b'\t',4)
            if len(fields)!=5:raise ValueError('Malformed genotype row')
            body=fields[4].strip();key=(fields[0],int(fields[1]))
            if body.count(b'\t')!=expected-1 or body.translate(None,b'0129\t'):
                raise ValueError(f'Unsupported genotype encoding at {key}; expected diploid 0/1/2/9 values')
            if last is not None and (key[0]!=last[0] or key[1]<=last[1]):raise ValueError('Genotype positions must be strictly increasing on one chromosome')
            last=key;rows+=1;ref_rows+=body.startswith(b'0\t')
    if not rows:raise ValueError('No informative genotype records')
    Path(output).write_text(json.dumps(dict(rows=rows,archaic_hom_ref_rows=ref_rows,modern_samples=expected-1),indent=2))

def main():
    parser=argparse.ArgumentParser();subs=parser.add_subparsers(dest='action',required=True)
    p=subs.add_parser('samples')
    for flag in ['panel','actual','output']:p.add_argument('--'+flag,required=True)
    p.add_argument('--no-background',action='store_true')
    p=subs.add_parser('mask')
    for flag in ['root','chrom','modern','fasta','upstream','output']:p.add_argument('--'+flag,required=True)
    p.add_argument('--axt-root')
    p=subs.add_parser('finalize')
    for flag in ['raw-dir','populations','output','chrom','build','locus']:p.add_argument('--'+flag,required=True)
    p.add_argument('--refs',nargs='+',required=True);p.add_argument('--core-start',type=int,default=0);p.add_argument('--core-end',type=int,default=0)
    p.add_argument('--min-bp',type=int,default=50000);p.add_argument('--lod',type=float,default=4);p.add_argument('--no-background',action='store_true')
    p=subs.add_parser('validate-mask')
    for flag in ['path','chrom','build']:p.add_argument('--'+flag,required=True)
    p=subs.add_parser('validate-genotypes');p.add_argument('--path',required=True);p.add_argument('--output',required=True)
    args=vars(parser.parse_args());action=args.pop('action')
    if action=='samples':args['background_required']=not args.pop('no_background');samples(**args)
    elif action=='mask':prepare_masks(**args)
    elif action=='validate-mask':
        length=CHROM_LENGTHS[args['build'].removeprefix('GRCh').removeprefix('b')][args['chrom']]
        previous=0
        for lo,hi in bed_intervals(args['path'],args['chrom'],length):
            if lo<previous:raise ValueError('Excluded mask must be sorted and merged')
            previous=hi
        # Reject wrong-contig rows rather than silently interpreting an include
        # mask or a mismatched chromosome as an empty excluded mask.
        with opener(args['path']) as handle:
            for line in handle:
                if line.strip() and not line.startswith(('#','track','browser')) and line.split()[0]!=args['chrom']:
                    raise ValueError('Use a chromosome-specific BED with numeric contig names (X for custom X)')
    elif action=='validate-genotypes':validate_genotypes(**args)
    else:args['background']=not args.pop('no_background');finalize(**args)

if __name__=='__main__':main()
