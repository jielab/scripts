#!/usr/bin/env python3
"""Stream complete native WGS genotypes into population-specific A/C/G/T counts.

This is an input bridge, not a variant caller. ``native_joint_wgs`` is an
explicit source contract: unsplit, jointly called WGS records, with complete
genotypes encoding every called nucleotide at that position. Other possible
bases have zero observed copies in THIS callset, not evidence of biological
absence. Repeated positions are discarded as a group; split ALT rows are never
combined or used to infer REF. Both biallelic and multiallelic source records
are retained. The sampling frame is the provided VCF's observed variant sites,
not all callable genome bases and not mutation opportunity per genomic base.
"""
from __future__ import annotations

# Joint Counts
# ======================================================================
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

from analysis_integrity import InputPath, identity, transactional

VERSION = 'joint-nucleotide-counts-1.0'
BASES = 'ACGT'
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
RANK = {chrom: index for index, chrom in enumerate(CHROMS)}
FIELDS = ['chrom', 'pos', 'ref', 'population', 'AN', 'A_AC', 'C_AC', 'G_AC', 'T_AC',
          'callable_all_four', 'call_rate', 'n_samples', 'n_called', 'n_missing',
          'n_haploid_called', 'n_diploid_called', 'build', 'context', 'count_semantics',
          'source_kind', 'sampling_frame', 'site_sampling_fraction', 'is_preview']
SITE_FIELDS = ['site_id', 'chrom', 'pos', 'ref', 'alts', 'candidate_alts', 'observed_alts_in_selected_genotypes', 'build', 'context',
               'sampling_frame', 'is_preview']


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith(('.gz', '.bgz', '.bgzf')) else open(path)


def canonical(value):
    value = value.removeprefix('chr').upper()
    return 'MT' if value == 'M' else value


def chromosome_selection(value):
    if value == 'all':
        return set(CHROMS)
    values = []
    for token in value.split(','):
        if token == '1-22':
            values.extend(str(i) for i in range(1, 23))
        else:
            values.append(canonical(token.strip()))
    if not values or set(values) - set(CHROMS) or len(values) != len(set(values)):
        raise ValueError('--chromosomes requires unique canonical names, 1-22, or all')
    return set(values)


def selected_by_hash(chrom, pos, fraction, seed):
    # Genotypes, ALT identities, annotation scores and population do not enter.
    digest = hashlib.blake2b(f'{seed}\t{chrom}\t{pos}'.encode(), digest_size=8).digest()
    return int.from_bytes(digest, 'big') < fraction * 2**64


def read_metadata(path):
    with open_text(path) as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith('##')), delimiter='\t')
        fields = reader.fieldnames or []
        normalized = {field: field.lstrip('#').strip().lower() for field in fields}
        if len(set(normalized.values())) != len(fields):
            raise ValueError('Duplicate metadata columns')
        sample = next((raw for raw, name in normalized.items() if name in {'sample_id', 'sample', 'iid'}), None)
        group = next((raw for wanted in ('population', 'super_pop', 'superpop', 'pop')
                      for raw, name in normalized.items() if name == wanted), None)
        include = next((raw for raw, name in normalized.items() if name == 'include'), None)
        if not sample or not group:
            raise ValueError('Metadata requires sample_id/sample/IID and population/super_pop/pop')
        result, seen = {}, set()
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('Malformed metadata row')
            name, population = row[sample].strip(), row[group].strip()
            if not name or name in seen:
                raise ValueError('Empty or duplicate sample in metadata: ' + name)
            seen.add(name)
            take = row[include].strip() if include else '1'
            if take not in {'0', '1'}:
                raise ValueError('Metadata include must be 0 or 1')
            if take == '0':
                continue
            if not population or population.upper() in {'NA', '.', 'NONE'}:
                raise ValueError('An included metadata sample lacks population: ' + name)
            result[name] = population
    if not result:
        raise ValueError('No included metadata samples')
    return result


def read_header(stream, build):
    raw, declared_builds = [], set()
    for line in stream:
        if line.startswith('##'):
            raw.append(line)
            if re.search(r'^##(?:source|.*Command|.*command).*\b(?:imput\w*|minimac\w*|beagle|glimpse\w*)\b', line, re.I):
                raise ValueError('VCF header declares imputation; native_joint_wgs is required')
            if re.search(r'bcftools.*norm', line, re.I) and re.search(r'(?:-m\s*-|--multiallelics[=\s]+-)', line):
                raise ValueError('VCF header declares multiallelic splitting; use unsplit native WGS')
            if line.lower().startswith(('##reference=', '##assembly=')):
                if re.search(r'(?i)(?:grch38|hg38)', line): declared_builds.add('GRCh38')
                if re.search(r'(?i)(?:grch37|hg19)', line): declared_builds.add('GRCh37')
            continue
        if line.startswith('#CHROM\t'):
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 10 or fields[8] != 'FORMAT':
                raise ValueError('Sites-only VCF is unsupported: FORMAT and sample genotypes are required')
            samples = fields[9:]
            if len(samples) != len(set(samples)) or any(not sample for sample in samples):
                raise ValueError('VCF sample IDs must be nonempty and unique')
            if declared_builds and declared_builds != {build}:
                raise ValueError(f'VCF assembly declaration {sorted(declared_builds)} conflicts with {build}')
            return samples, dict(header_sha256=hashlib.sha256(''.join(raw + [line]).encode()).hexdigest(),
                build_evidence='matching_vcf_header' if declared_builds else 'user_declaration_only')
        if line.strip():
            raise ValueError('Expected VCF meta-information and #CHROM header before variants')
    raise ValueError('VCF has no #CHROM header')


def loci(stream, n_samples, selected, audit):
    """Keep one record, plus a duplicate counter, for the current coordinate."""
    current, first, count, previous = None, None, 0, None
    for line in stream:
        if not line.strip():
            continue
        if line.startswith('#'):
            raise ValueError('Unexpected header inside the VCF variant stream')
        row = line.rstrip('\r\n').split('\t')
        audit['vcf_records_read'] += 1
        if len(row) != 9 + n_samples:
            raise ValueError('VCF row/sample cardinality mismatch')
        chrom = canonical(row[0])
        if chrom not in RANK:
            audit['noncanonical_records_excluded'] += 1
            continue
        if chrom not in selected:
            audit['out_of_scope_records_excluded'] += 1
            continue
        try:
            pos = int(row[1])
        except ValueError as exc:
            raise ValueError('Invalid VCF POS') from exc
        if pos < 1:
            raise ValueError('VCF POS must be positive')
        key = RANK[chrom], pos
        if previous is not None and key < previous:
            raise ValueError('VCF must be coordinate-sorted in canonical chromosome order 1..22,X,Y,MT')
        previous = key
        if re.search(r'(?:^|;)(?:IMPUTED|IMP)(?:=|;|$)', row[7], re.I):
            raise ValueError(f'Imputed record at {chrom}:{pos}; native joint WGS is required')
        if current is not None and key != current:
            yield first, count
            first, count = None, 0
        current = key
        count += 1
        if first is None:
            row[0] = chrom
            first = row
    if first is not None:
        yield first, count


def genotype(value, formats, n_alleles, min_gq, min_dp):
    tokens = value.split(':')
    if len(tokens) > len(formats):
        return None, 'malformed_format'
    values = dict(zip(formats, tokens))
    gt = values.get('GT', '.')
    alleles = re.split(r'[/|]', gt)
    if any(x in {'.', ''} for x in alleles):
        return None, 'missing_gt'
    if len(alleles) not in {1, 2} or any(not x.isdigit() for x in alleles):
        return None, 'invalid_gt'
    indices = [int(x) for x in alleles]
    if any(x >= n_alleles for x in indices):
        return None, 'invalid_gt'
    if values.get('FT', 'PASS') not in {'PASS', '.', ''}:
        return None, 'genotype_filter'
    for field, threshold in (('GQ', min_gq), ('DP', min_dp)):
        if threshold == 0:
            continue
        value = values.get(field, '.')
        try:
            numeric = float(value)
        except (ValueError, TypeError):
            return None, 'missing_' + field.lower()
        if not math.isfinite(numeric) or numeric < 0 or field == 'DP' and numeric != int(numeric):
            return None, 'invalid_' + field.lower()
        if numeric < threshold:
            return None, 'low_' + field.lower()
    return indices, 'called'


def reference_reader(path):
    # Reuse the project reader: existing .fai, bounded cache; compressed FASTA
    # additionally needs pysam/.gzi. No whole chromosome is loaded into memory.
    if 'allele_common' not in sys.modules:
        spec = importlib.util.spec_from_file_location('allele_common', Path(__file__).with_name('0.common.py'))
        module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module
        spec.loader.exec_module(module)
    return sys.modules['allele_common'].FastaContext(path)


def extra_inputs(raw):
    if raw.get('fasta'):
        for suffix in ('.fai', '.gzi'):
            path = Path(str(raw['fasta']) + suffix)
            if path.is_file():
                yield InputPath(path)


@transactional('qc/manifest.json', input_files=('vcf', 'metadata', 'fasta'), resolve_inputs=extra_inputs)
def run(args):
    if args.source_kind != 'native_joint_wgs' or args.sampling_frame != 'observed_variant_sites':
        raise ValueError('Requires explicit native_joint_wgs and observed_variant_sites contracts')
    if not 0 < args.site_fraction <= 1 or not 0 <= args.min_call_rate <= 1:
        raise ValueError('site-fraction must be (0,1]; min-call-rate must be [0,1]')
    if any(not math.isfinite(x) or x < 0 for x in (args.min_dp, args.min_gq)) or args.max_sites < 0:
        raise ValueError('Quality thresholds and preview max-sites must be nonnegative')
    if (args.min_gq == 0 or args.min_dp == 0) and not getattr(args, 'source_qc', '').strip():
        raise ValueError('Disabling GQ/DP thresholds requires --source-qc describing the upstream native WGS QC')
    selected = chromosome_selection(args.chromosomes)
    out = Path(args.out); (out/'results').mkdir(); (out/'qc').mkdir()
    counts_path, sites_path = out/'results/joint_counts.tsv.gz', out/'results/sites.tsv.gz'
    meta = read_metadata(args.metadata) if args.metadata else None
    audit = Counter(); populations = {}; header_info = {}; selected_samples = []
    with ExitStack() as stack:
        source = stack.enter_context(open_text(args.vcf))
        samples, header_info = read_header(source, args.build)
        slots = [(i, meta[sample] if meta else 'ALL', sample)
                 for i, sample in enumerate(samples) if meta is None or sample in meta]
        if not slots:
            raise ValueError('No selected samples match the VCF')
        selected_samples = [sample for _, _, sample in slots]
        populations = dict(sorted(Counter(pop for _, pop, _ in slots).items()))
        audit['vcf_samples'] = len(samples); audit['selected_samples'] = len(slots)
        audit['vcf_samples_not_selected'] = len(samples) - len(slots)
        audit['metadata_samples_not_in_vcf'] = len(set(meta or {}) - set(samples))
        fasta = stack.enter_context(reference_reader(args.fasta)) if args.fasta else None
        counts_out = stack.enter_context(gzip.open(counts_path, 'wt', newline=''))
        sites_out = stack.enter_context(gzip.open(sites_path, 'wt', newline=''))
        writer = csv.DictWriter(counts_out, FIELDS, delimiter='\t', lineterminator='\n'); writer.writeheader()
        site_writer = csv.DictWriter(sites_out, SITE_FIELDS, delimiter='\t', lineterminator='\n'); site_writer.writeheader()
        for row, duplicates in loci(source, len(samples), selected, audit):
            audit['source_positions_read'] += 1
            chrom, pos, ref = row[0], int(row[1]), row[3].upper()
            if duplicates != 1:
                audit['duplicate_positions_excluded'] += 1
                audit['duplicate_position_records_excluded'] += duplicates
                continue
            alleles = [ref] + row[4].upper().split(',')
            if any(base not in BASES or len(base) != 1 for base in alleles) or len(set(alleles)) != len(alleles):
                audit['non_snv_or_invalid_allele_positions_excluded'] += 1
                continue
            if row[6] != 'PASS' and not (args.allow_unfiltered and row[6] == '.'):
                audit['site_filter_positions_excluded'] += 1
                continue
            audit['eligible_native_snv_positions'] += 1
            if not selected_by_hash(chrom, pos, args.site_fraction, args.seed):
                audit['hash_unsampled_positions'] += 1
                continue
            if args.max_sites and audit['emitted_positions'] >= args.max_sites:
                audit['preview_stopped_before_end'] = 1
                break
            formats = row[8].split(':')
            if len(formats) != len(set(formats)) or any(not value for value in formats):
                raise ValueError(f'Duplicate/empty FORMAT labels at {chrom}:{pos}')
            context = fasta.context(chrom, pos) if fasta else 'NA'
            if fasta and context[1] != ref:
                raise ValueError(f'FASTA/VCF REF mismatch at {chrom}:{pos}: {context[1]} != {ref}')
            stats = {pop: Counter() for pop in populations}
            for index, pop, _ in slots:
                gt, reason = genotype(row[9+index], formats, len(alleles), args.min_gq, args.min_dp)
                audit['genotypes_' + reason] += 1
                if gt is None:
                    continue
                state = stats[pop]; state['n_called'] += 1
                state['n_haploid_called' if len(gt) == 1 else 'n_diploid_called'] += 1
                for index in gt:
                    state[alleles[index]] += 1
            for pop, n_samples in populations.items():
                state = stats[pop]; an = sum(state[base] for base in BASES)
                rate = state['n_called'] / n_samples
                callable_four = int(an > 0 and rate >= args.min_call_rate)
                writer.writerow(dict(chrom='chr' + chrom, pos=pos, ref=ref, population=pop, AN=an,
                    **{base + '_AC': state[base] for base in BASES}, callable_all_four=callable_four,
                    call_rate=rate, n_samples=n_samples, n_called=state['n_called'],
                    n_missing=n_samples-state['n_called'], n_haploid_called=state['n_haploid_called'],
                    n_diploid_called=state['n_diploid_called'], build=args.build, context=context,
                    count_semantics='joint_nucleotide_counts', source_kind=args.source_kind,
                    sampling_frame=args.sampling_frame, site_sampling_fraction=args.site_fraction,
                    is_preview=int(args.max_sites > 0)))
                audit['emitted_population_rows'] += 1
                audit['callable_population_rows' if callable_four else 'uncallable_population_rows'] += 1
            site_writer.writerow(dict(site_id=f'chr{chrom}:{pos}:{ref}', chrom='chr' + chrom,
                pos=pos, ref=ref, alts=','.join(alleles[1:]),
                candidate_alts=','.join(base for base in BASES if base != ref),
                observed_alts_in_selected_genotypes=','.join(base for base in BASES if base != ref
                    and sum(state[base] for state in stats.values()) > 0) or 'NA', build=args.build,
                context=context, sampling_frame=args.sampling_frame, is_preview=int(args.max_sites > 0)))
            audit['emitted_positions'] += 1
    with (out/'qc/audit.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n'); writer.writerow(['metric', 'count'])
        writer.writerows(sorted(audit.items()))
    status = 'complete_with_exclusions' if audit['callable_population_rows'] else 'not_available'
    return dict(status=status, version=VERSION, build=args.build, source_kind=args.source_kind,
        sampling_frame=args.sampling_frame, analysis_mode='preview_first_n_selected' if args.max_sites else 'declared_source_sampling_frame',
        inference_allowed=not bool(args.max_sites), site_fraction=args.site_fraction, seed=args.seed,
        output=str(counts_path), site_table=str(sites_path), counts=dict(audit), populations=populations,
        site_table_semantics='alts lists encoded source VCF ALT labels, not proof of observed copies in every population; observed_alts_in_selected_genotypes uses passing GT across selected samples; candidate_alts lists all three theoretical substitutions',
        selected_sample_order_sha256=hashlib.sha256('\n'.join(selected_samples).encode()).hexdigest(),
        population_assignment='explicit_metadata' if meta else 'ALL_not_ancestry_stratified',
        header=header_info, genotype_qc=dict(min_gq=args.min_gq, min_dp=args.min_dp,
            min_call_rate=args.min_call_rate, partial_genotypes='entire_genotype_missing',
            source_qc_declaration=getattr(args, 'source_qc', ''), source_qc_verified_by_this_tool=False,
            missing_quality='fails_when_threshold_positive', site_filter='PASS_or_dot' if args.allow_unfiltered else 'PASS'),
        chromosomes=[c for c in CHROMS if c in selected], count_semantics='joint_nucleotide_counts',
        zero_semantics='zero observed nucleotide copies among complete high-quality genotypes in the declared native joint callset; not absence from humanity',
        callable_all_four_semantics='Native unsplit joint-call WGS source declaration plus complete measured GT and population call-rate threshold; equal sensitivity to all four bases is not demonstrated',
        limitations=[
            'Only positions present as variant records in the provided source VCF are sampled. No genome-wide invariant/callable-base denominator is reconstructed.',
            'A candidate-only or prefiltered VCF remains a candidate-only sampling frame even if source-kind says native_joint_wgs; supply the full observed-variant callset for enrichment controls.',
            'All complete passing GT alleles are counted directly. INFO/AC, DS, GP and split-row REF reconstruction are never used.',
            'Low-call-rate and AN=0 rows remain in the table with callable_all_four=0 and must be excluded, not treated as negative observations.',
            'Repeated positions are excluded entirely, including duplicates caused by SNP/indel records sharing a position.',
            'The site hash depends only on chromosome, position and seed; max-sites is a nonrepresentative preview and is marked in every row.',
            'Call rate uses selected individuals as denominator. Sex chromosomes/MT use GT ploidy without inferring sex-specific eligibility; autosomes are the default.',
            'These counts cannot establish mutation age, independent population replication, uniform WGS sensitivity or extinction of an unobserved allele.'])


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf', required=True, help='Coordinate-sorted native joint WGS VCF[.gz/.bgz], with sample GT')
    p.add_argument('--source-kind', required=True, choices=('native_joint_wgs',))
    p.add_argument('--sampling-frame', required=True, choices=('observed_variant_sites',))
    p.add_argument('--build', required=True, choices=('GRCh37', 'GRCh38'))
    p.add_argument('--metadata', help='TSV with sample_id and population; optional include=0/1. Without metadata all VCF samples form ALL.')
    p.add_argument('--fasta', help='Optional matching indexed reference for trinucleotide context and REF checks')
    p.add_argument('--chromosomes', default='1-22', help='1-22 (default), comma-separated canonical names, or all')
    p.add_argument('--min-gq', type=float, default=20)
    p.add_argument('--min-dp', type=float, default=10)
    p.add_argument('--min-call-rate', type=float, default=.95)
    p.add_argument('--source-qc', default='', help='Upstream native WGS QC description; required if either GQ or DP threshold is set to zero')
    p.add_argument('--allow-unfiltered', action='store_true', help='Explicitly allow FILTER=. in addition to PASS')
    p.add_argument('--site-fraction', type=float, default=1., help='Outcome-independent chromosome/position hash sampling fraction')
    p.add_argument('--seed', type=int, default=20261110)
    p.add_argument('--max-sites', type=int, default=0, help='Positive values make a first-N preview, prohibited for formal enrichment')
    p.add_argument('--out', required=True)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    result = run(args)
    print(json.dumps({key: result[key] for key in ('status', 'output', 'site_table', 'counts')}, indent=2))
    return 0


# Site Regions
# ======================================================================
_description_site_regions = 'Attach same-build BED membership to an outcome-independent SNV site table.\n\nInput sites use 1-based POS; BED uses 0-based, half-open intervals. One output\nrow per site applies to every population in the joint-count table. Overlap is\nan annotation, not evidence that a variant originated in a virus or archaic\nhuman. No implicit coordinate conversion or missing-as-zero imputation.\n'

from bisect import bisect_right
import itertools


def regions_text_open(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def regions_chrom_key(chrom):
    c = re.sub("^chr", "", chrom, flags=re.I).upper()
    c = "MT" if c == "M" else c
    if c not in {*(str(i) for i in range(1, 23)), "X", "Y", "MT"}:
        return None
    return int(c) if c.isdigit() else {"X": 23, "Y": 24, "MT": 25}[c]


def regions_parse_beds(specifications):
    specs = {}
    for value in specifications:
        if "=" not in value:
            raise ValueError("--bed must be NAME=PATH")
        label, path = value.split("=", 1)
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", label) or label in specs:
            raise ValueError("BED labels must be distinct simple column names")
        if label in {"chrom", "pos", "ref", "build", "population", "sites"}:
            raise ValueError("Reserved BED label: " + label)
        specs[label] = Path(path).expanduser().resolve()
    return specs


def regions_bed_index(path):
    intervals, audit = defaultdict(list), Counter()
    with regions_text_open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("BED needs at least three tab-separated columns")
            c = regions_chrom_key(fields[0])
            if c is None:
                audit["noncanonical_contig_rows_ignored"] += 1
                continue
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError("BED intervals must satisfy 0 <= start < end")
            intervals[c].append((start, end)); audit["source_intervals"] += 1
    result = {}
    for c, rows in intervals.items():
        merged = []
        for start, end in sorted(rows):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        result[c] = ([s for s, e in merged], [e for s, e in merged])
        audit["union_bases"] += sum(e-s for s, e in merged)
        audit["union_intervals"] += len(merged)
    return result, audit


def regions_contains(index, chrom, pos):
    starts, ends = index.get(chrom, ([], []))
    i = bisect_right(starts, pos-1) - 1
    return int(i >= 0 and pos-1 < ends[i])


def regions_input_paths(raw):
    return [InputPath(p) for p in regions_parse_beds(raw.get("bed", [])).values()]


def regions_site_rows(path, build, keep):
    previous, last = None, None
    with regions_text_open(path) as f:
        reader = csv.DictReader((l for l in f if not l.startswith("##")), delimiter="\t")
        required = {"chrom", "pos", "ref", *keep}
        if required - set(reader.fieldnames or []):
            raise ValueError("Site table missing columns: " + str(sorted(required-set(reader.fieldnames or []))))
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed site TSV")
            c, pos, ref = regions_chrom_key(row["chrom"]), int(row["pos"]), row["ref"]
            if c is None or pos < 1 or ref not in "ACGT" or len(ref) != 1:
                raise ValueError("Expected canonical SNV site and positive 1-based POS")
            if row.get("build", build) != build:
                raise ValueError("Site build disagrees with --build")
            key = c, pos
            canon = str(c) if c <= 22 else {23:"X",24:"Y",25:"MT"}[c]
            current = dict(chrom="chr"+canon, pos=pos, ref=ref, build=build,
                           **{k:row[k] for k in keep})
            if previous is not None and key < previous:
                raise ValueError("Sites must be sorted by chromosome 1..22,X,Y,MT and POS")
            if previous == key:
                if current != last:
                    raise ValueError("Population rows disagree on site REF/build/retained context")
                continue
            previous, last = key, dict(current)
            yield c, current


@transactional("qc/manifest.json", input_files=("sites",), resolve_inputs=regions_input_paths)
def regions_run(args):
    if args.build != args.bed_build:
        raise ValueError("Sites and BED must use the same build; liftover is a separate allele-checked operation")
    specs = regions_parse_beds(args.bed)
    if not specs:
        raise ValueError("Supply at least one --bed NAME=PATH")
    keep = [v.strip() for v in args.keep_columns.split(",") if v.strip()]
    if len(keep) != len(set(keep)) or set(keep) & (set(specs) | {"chrom","pos","ref","build"}):
        raise ValueError("Retained columns must be distinct from coordinate/BED fields")
    before = {str(p): identity(p) for p in [Path(args.sites), *specs.values()]}
    indices, audits = {}, []
    for label, path in specs.items():
        indices[label], a = regions_bed_index(path)
        audits.append(dict(annotation=label, source=str(path), **a))
    out = Path(args.out); (out/"results").mkdir(parents=True, exist_ok=True)
    fields = ["chrom","pos","ref","build",*keep,*specs]
    totals = Counter()
    with gzip.open(out/"results/site_annotations.tsv.gz", "wt", newline="") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
        for c, row in regions_site_rows(args.sites, args.build, keep):
            for label, index in indices.items():
                row[label] = regions_contains(index, c, row["pos"])
                totals[label] += row[label]
            writer.writerow(row); totals["sites"] += 1
    if not totals["sites"]:
        raise ValueError("No site rows")
    if any(identity(p) != saved for p, saved in before.items()):
        raise ValueError("Input changed during annotation")
    manifest = dict(status="complete", settings=vars(args), inputs=before, bed_summary=audits,
                    totals=dict(totals), site_coordinates="1-based", bed_coordinates="0-based half-open",
                    interpretation="Exact position overlap only; each BED must be a complete declared annotation for its category",
                    zero_meaning="Site outside supplied BED union; no inference for omitted categories or unknown annotation coverage",
                    biological_limit="HERV/archaic overlap does not identify an SNV's historical origin")
    (out/"qc").mkdir(exist_ok=True)
    (out/"qc/manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def regions_parser():
    p = argparse.ArgumentParser(description=_description_site_regions)
    p.add_argument("--sites", required=True, help="Sorted joint counts or prespecified site TSV[.gz]")
    p.add_argument("--build", required=True, choices=["GRCh37","GRCh38"])
    p.add_argument("--bed-build", required=True, choices=["GRCh37","GRCh38"], help="Explicit common build of every BED")
    p.add_argument("--bed", action="append", default=[], help="NAME=PATH; repeat for independent overlapping annotations")
    p.add_argument("--keep-columns", default="", help="Optional shared per-site fields, e.g. context,mappability_bin")
    p.add_argument("--out", required=True)
    return p


def regions_main(argv=None):
    result=regions_run(regions_parser().parse_args(argv)); print(json.dumps(dict(status=result["status"], totals=result["totals"])))
    return 0


# Allele Balance
# ======================================================================
_description_allele_balance = 'Summarize actual nucleotide counts from genotype_validation.tsv.\n\nThis is conditional on the predefined target and genotype-encoded nucleotide\npanel. It neither infers REF counts from split population summaries nor dates\nmutations. Rarefaction samples chromosome copies, not unrelated individuals.\n'


BALANCE_BASES = "ACGT"
BALANCE_MISSING = {"", ".", "NA", "nan", "None"}


def balance_integer(value, label):
    x = float(value)
    if not math.isfinite(x) or x < 0 or x != int(x):
        raise ValueError(f"{label} must be a finite nonnegative integer")
    return int(x)


def balance_rarefaction(counts, n):
    """Exact finite-population P(>=3 observed bases), each with >=1 copy.

    Inclusion/exclusion uses integer binomial coefficients, so cancellation
    cannot create a negative tiny probability. No independent-ALT assumption.
    """
    values = [int(c) for c in counts if c > 0]
    N, k = sum(values), len(values)
    if n < 1 or n > N:
        return None, None
    choose = lambda x: math.comb(x, n) if x >= n else 0
    denominator = choose(N)
    expected_states = sum(denominator - choose(N-c) for c in values) / denominator
    if k < 3:
        return 0.0, expected_states
    # Each one-state outcome appears k-1 times in the pair sum.
    pair_sum = sum(choose(values[i]+values[j]) for i in range(k) for j in range(i+1, k))
    numerator = denominator - pair_sum + (k-2)*sum(choose(c) for c in values)
    if not 0 <= numerator <= denominator:
        raise AssertionError("Invalid multivariate rarefaction probability")
    return numerator/denominator, expected_states


def balance_summarize_row(row, minimum, target_an, min_call_rate):
    out = {k: row.get(k, "") for k in ("site_id", "chrom", "pos", "ref", "family", "group", "population")}
    out.update(status="unavailable", reason="", min_ac=minimum, rarefy_an=target_an,
               scope="predefined_targets_genotype_encoded_bases")
    if row.get("evaluable") != "1" or row.get("represented_status") != "all_expected_alts_encoded":
        out["reason"] = "target_not_fully_evaluable"
        return out
    if row.get("record_status") != "ok":
        raise ValueError("evaluable=1 conflicts with record_status")
    ref = row["ref"]
    if ref not in BALANCE_BASES or len(ref) != 1:
        raise ValueError("REF must be one DNA nucleotide")
    encoded = set(row["encoded_alts"].split(",")) - {"", "."}
    if not encoded <= set(BALANCE_BASES)-{ref}:
        raise ValueError("Invalid encoded ALT nucleotides")
    counts = {}
    for b in BALANCE_BASES:
        raw = row.get(b+"_AC", "NA")
        if raw in BALANCE_MISSING:
            if b in encoded | {ref}:
                raise ValueError("Missing count for an encoded nucleotide")
            continue
        counts[b] = balance_integer(raw, b+"_AC")
        if b not in encoded | {ref}:
            raise ValueError("Count supplied for an unencoded nucleotide")
    an = balance_integer(row["AN"], "AN")
    call_rate = float(row["call_rate"])
    if an < 1 or sum(counts.values()) != an or not math.isfinite(call_rate) or not 0 <= call_rate <= 1:
        raise ValueError("Counts, AN or call rate are inconsistent")
    if call_rate < min_call_rate:
        out["reason"] = "low_call_rate"
        return out
    order = sorted(counts, key=lambda b: (-counts[b], b))
    f = [counts[b]/an for b in order]
    n_supported = sum(c >= minimum for c in counts.values())
    p3, expected = balance_rarefaction(counts.values(), target_an)
    out.update(status="complete_descriptive", AN=an, call_rate=call_rate,
               n_observed_bases=sum(c > 0 for c in counts.values()),
               n_bases_ac_ge_min=n_supported, at_least_three_supported=int(n_supported >= 3),
               exactly_three_supported=int(n_supported == 3), four_supported=int(n_supported == 4),
               observed_ref_ac=counts[ref], observed_ref_af=counts[ref]/an,
               rank1_tied=int(len(f)>1 and f[0] == f[1]),
               rank3_tied=int(len(f)>3 and f[2] == f[3]),
               entropy_nats=-sum(x*math.log(x) for x in f if x),
               diversity_one_minus_sum_p2=1-sum(x*x for x in f),
               p_rarefied_ge3_each_ac_ge1=p3, expected_rarefied_states=expected,
               rarefaction_status="complete" if p3 is not None else "AN_below_target",
               rank_meaning="current_count_order_not_mutation_order")
    for i in range(4):
        out[f"rank{i+1}_base"] = order[i] if i<len(order) else None
        out[f"rank{i+1}_ac"] = counts[order[i]] if i<len(order) else None
        out[f"rank{i+1}_af"] = f[i] if i<len(order) else None
    return out


def balance_write_table(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({k: "NA" if r.get(k) is None else r.get(k, "NA") for k in keys} for r in rows)


@transactional("balance_manifest.json")
def balance_run(args):
    if args.min_ac < 1 or not 1 <= args.rarefy_an <= 5000 or not 0 <= args.min_call_rate <= 1:
        raise ValueError("Invalid min-ac, rarefy-an (1..5000) or min-call-rate")
    before = identity(args.genotypes, full=True)
    source = gzip.open(args.genotypes, "rt") if str(args.genotypes).endswith(".gz") else open(args.genotypes)
    found, results = {}, []
    with source as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"site_id", "chrom", "pos", "ref", "family", "group", "population", "AN", "call_rate",
                    "A_AC", "C_AC", "G_AC", "T_AC", "encoded_alts", "evaluable", "represented_status", "record_status"}
        if required - set(reader.fieldnames or []):
            raise ValueError("Input is not genotype_validation.tsv; missing: "+str(sorted(required-set(reader.fieldnames or []))))
        for row in reader:
            result = balance_summarize_row(row, args.min_ac, args.rarefy_an, args.min_call_rate)
            key = tuple(result[k] for k in ("chrom", "pos", "ref", "family", "group", "population"))
            if key in found:
                if result != found[key]:
                    raise ValueError("Conflicting duplicate site/population count")
                continue
            found[key] = result
            results.append(result)
    if not results:
        raise ValueError("Empty genotype table")
    groups = defaultdict(list)
    for row in results:
        groups[row["family"], row["group"], row["population"]].append(row)
    summary=[]
    for (fam, group, pop), rr in sorted(groups.items()):
        ok=[r for r in rr if r["status"]=="complete_descriptive"]
        equal=[r for r in ok if r["p_rarefied_ge3_each_ac_ge1"] is not None]
        summary.append(dict(family=fam, group=group, population=pop, n_targets=len(rr), n_evaluable=len(ok),
                            n_three_supported=sum(r["at_least_three_supported"] for r in ok),
                            n_rarefiable=len(equal), rarefy_an=args.rarefy_an,
                            mean_rarefied_ge3_probability=sum(r["p_rarefied_ge3_each_ac_ge1"] for r in equal)/len(equal) if equal else None,
                            interpretation="conditional_target_summary; no selection or genome-wide enrichment inference"))
    out=Path(args.out); out.mkdir(parents=True, exist_ok=True)
    balance_write_table(out/"allele_balance.tsv", results)
    balance_write_table(out/"balance_summary.tsv", summary)
    if before != identity(args.genotypes, full=True):
        raise ValueError("Genotype input changed during analysis")
    manifest=dict(status="complete", input=before, settings=vars(args),
                  estimand="Nucleotide-count diversity conditional on predefined targets and fully encoded callable genotypes",
                  rarefaction="exact multivariate hypergeometric over observed chromosome copies; >=1 per base, unrelatedness and LD not modeled",
                  caveat="Balanced frequencies or three states alone do not establish balancing selection; ranks do not imply ancestral or temporal order")
    (out/"balance_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def balance_main(argv=None):
    p=argparse.ArgumentParser(description=_description_allele_balance)
    p.add_argument("--genotypes", required=True, help="Existing genotypes/genotype_validation.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--min-ac", type=int, default=2)
    p.add_argument("--min-call-rate", type=float, default=.9)
    p.add_argument("--rarefy-an", type=int, default=100)
    balance_run(p.parse_args(argv))


def cli(argv=None):
    """Dispatch auxiliary commands owned by this module."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    routes = {'regions': regions_main, 'balance': balance_main}
    entry = routes[argv.pop(0)] if argv and argv[0] in routes else main
    return entry(argv)



if __name__ == '__main__':
    try:
        raise SystemExit(cli())
    except (ValueError, OSError, RuntimeError) as exc:
        raise SystemExit('joint-counts: ' + str(exc))

