#!/usr/bin/env python3
"""Prepare original COJO leads, never distance-clump independent GWAS signals."""
import argparse, csv, hashlib, json, math, os, re, fcntl
from pathlib import Path
from comm import CHROM_LENGTHS
from cojo_loci import lift


def write(path, data, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for r in data for k in r)) or ['status']
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(data)
    tmp.replace(path)


def read(path):
    with path.open() as f: return list(csv.DictReader(f, delimiter='\t'))


def complement(a): return a.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def parse_cojo(path, build):
    lines = [x.split() for x in path.read_text().splitlines() if x.strip()]
    required = {'SNP', 'Chr', 'bp', 'refA', 'bJ', 'pJ'}
    if not lines or not required <= set(lines[0]):
        raise ValueError('COJO requires columns: ' + ', '.join(sorted(required)))
    accepted, audit, seen = [], [], set()
    for n, vals in enumerate(lines[1:], 1):
        r = dict(zip(lines[0], vals)); name = r.get('SNP', '')
        try:
            if len(vals) != len(lines[0]): raise ValueError('wrong_column_count')
            ch = re.sub('^chr', '', r['Chr'], flags=re.I).upper(); ch = 'X' if ch == '23' else ch
            pos = int(r['bp']); beta = float(r['bJ']); p = float(r['pJ']); effect = r['refA'].upper()
            if ch not in CHROM_LENGTHS[build] or not 1 <= pos <= CHROM_LENGTHS[build][ch]: raise ValueError('invalid_position')
            if not math.isfinite(beta) or beta == 0 or not math.isfinite(p) or not 0 <= p <= 1: raise ValueError('invalid_effect_or_P')
            if not re.fullmatch('[ACGT]+', effect): raise ValueError('invalid_effect_allele')
            m = re.fullmatch(r'(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)', name, re.I)
            ref, alt = ('', '')
            if m:
                mc = 'X' if m[1] == '23' else m[1].upper()
                if mc != ch or int(m[2]) != pos: raise ValueError('SNP_ID_disagrees_with_CHR_POS')
                ref, alt = m[3].upper(), m[4].upper()
                if ref == alt or effect not in (ref, alt): raise ValueError('effect_allele_disagrees_with_SNP_ID')
            elif not re.fullmatch('rs[0-9]+', name): raise ValueError('unsupported_SNP_ID')
            if name in seen: raise ValueError('duplicate_SNP_ID')
            seen.add(name)
            accepted.append(dict(input_row=n, locus_id=name, index_snp=name, source_build='GRCh'+build,
                source_chr=ch, source_pos=pos, chr=ch, lead_pos=pos, ref=ref, alt=alt,
                effect_allele=effect, beta_j=beta, p_j=p, strand='+'))
        except (ValueError, KeyError) as e:
            audit.append(dict(input_row=n, index_snp=name, status='skipped', reason=str(e)))
    return accepted, audit


def prepare(a):
    out = a.output; out.mkdir(parents=True, exist_ok=True)
    leads, audit = parse_cojo(a.input, a.build)
    provenance = dict(input=str(a.input.resolve()), input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),
        source_build=a.build, analysis_build='37', method='original_COJO_leads_no_distance_collapse', search_flank_bp=a.window)
    if a.build == '38' and leads:
        binary = Path(os.environ.get('GU_LIFTOVER', '/mnt/d/software/bin/liftOver'))
        print(f'[GU GWAS] liftOver={binary.resolve()} sha256={hashlib.sha256(binary.read_bytes()).hexdigest()}',flush=True)
        source = out/'leads.GRCh38.bed'
        source.write_text(''.join(f"chr{r['chr']}\t{r['lead_pos']-1}\t{r['lead_pos']-1+max(1,len(r['ref']))}\t{r['input_row']}\t0\t+\n" for r in leads))
        hits = lift(binary, a.chain, source, out/'leads.lifted.GRCh37.bed', out/'leads.unmapped.bed')
        provenance.update(liftover_binary=str(binary.resolve()), liftover_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                          chain=str(a.chain), chain_sha256=hashlib.sha256(a.chain.read_bytes()).hexdigest())
        mapped = []
        for r in leads:
            h = hits[str(r['input_row'])]
            reason = ''
            if len(h) != 1: reason = 'lead_unmapped_or_multiple_mapping'
            else:
                x = h[0]; ch = x[0].removeprefix('chr'); pos = int(x[1])+1
                if ch != r['chr'] or ch not in CHROM_LENGTHS['37'] or int(x[2])-int(x[1]) != max(1,len(r['ref'])):
                    reason = 'lead_mapping_changes_chromosome_or_allele_span'
                else:
                    r.update(chr=ch, lead_pos=pos, strand=x[5])
                    if x[5] == '-':
                        for key in ('ref', 'alt', 'effect_allele'): r[key] = complement(r[key])
            if reason: audit.append(dict(input_row=r['input_row'], index_snp=r['index_snp'], status='skipped', reason=reason))
            else: mapped.append(r)
        leads = mapped
    for r in leads:
        r.update(search_start=max(0,r['lead_pos']-1-a.window), search_end=min(CHROM_LENGTHS['37'][r['chr']],r['lead_pos']+a.window))
    write(out/'gwas_leads.GRCh37.tsv', leads)
    write(out/'input.audit.tsv', audit, ['input_row','index_snp','status','reason'])
    for r in audit: print(f"[GU GWAS] SKIP row={r['input_row']} SNP={r['index_snp']} reason={r['reason']}", flush=True)
    if not leads: raise ValueError('No valid original COJO leads remain')
    # These are extraction windows, not inferred core haplotypes or standalone inputs.
    (out/'search_windows.GRCh37.bed').write_text(''.join(f"{r['chr']}\t{r['search_start']}\t{r['search_end']}\t{r['locus_id']}\n" for r in leads))
    provenance.update(n_leads=len(leads), n_skipped=len(audit))
    (out/'manifest.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(f"[GU GWAS] {len(leads)} original leads retained; {len(audit)} skipped; LD=1KG EUR; cores will be defined by r² > 0.98", flush=True)
    print(f"[GU GWAS] lead audit: {out/'gwas_leads.GRCh37.tsv'}; extraction BED is NOT a standalone PhyML input", flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--build',choices=['37','38'],required=True); p.add_argument('--window',type=int,default=500000)
    p.add_argument('--chain',type=Path,default=Path('/mnt/d/files/liftOver/hg38ToHg19.over.chain.gz'))
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'.prepare.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        prepare(args)
