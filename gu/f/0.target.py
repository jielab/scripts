#!/usr/bin/env python3
"""Resolve a direct PGEN cohort before scheduling; never read genotype data."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import sys

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('gu_ukb_input', Path(__file__).with_name('0.ukb.py'))
u = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = u
spec.loader.exec_module(u)


def psam_chunk(value):
	"""Slice sample rows first; PLINK applies sex filtering to this fixed slice."""
	parts = value.rsplit(',', 2)
	if len(parts) != 3 or not all(x.isdecimal() and int(x) > 0 for x in parts[1:]):
		raise ValueError('--keep-psam requires FILE,CHUNK_SIZE,CHUNK_INDEX with positive integers (index starts at 1)')
	path = Path(parts[0]).expanduser().resolve()
	size, index = map(int, parts[1:])
	digest = hashlib.sha256(path.read_bytes()).hexdigest()
	if os.environ.get('GU_CMD_WORKER') == '1' and os.environ.get('GU_KEEP_PSAM_SHA256', digest) != digest:
		raise ValueError('--keep-psam source changed after scheduling; regenerate the commands')
	samples = u.read_psam(path)
	start = (index - 1) * size
	if start >= len(samples):
		raise ValueError(f'--keep-psam chunk {index} is out of range: {len(samples)} samples, {(len(samples) + size - 1) // size} chunks')
	selected = samples[start:start + size]
	info = dict(source=str(path), source_sha256=digest, chunk_size=size, chunk_index=index,
		total_samples=len(samples), first_sample=start + 1, last_sample=start + len(selected),
		chunk_samples=len(selected), selection_order='PSAM rows, then optional --keep-males')
	key = hashlib.sha256(json.dumps(info, sort_keys=True).encode()).hexdigest()[:24]
	work = Path('/tmp/gu-psam-chunks') / key
	work.mkdir(parents=True, exist_ok=True)
	for name, content in [('chunk.keep', '#FID\tIID\n' + ''.join(s['fid'] + '\t' + s['sample'] + '\n' for s in selected)),
		('chunk.json', json.dumps(info, indent=2, sort_keys=True) + '\n')]:
		output = work / name
		if not output.exists() or output.read_text() != content:
			temporary = work / (name + f'.{os.getpid()}.tmp')
			temporary.write_text(content)
			temporary.replace(output)
	print(f'[GU samples] PSAM chunk {index}: rows {start + 1}–{start + len(selected)} of {len(samples)}; {len(selected)} samples before sex filtering', file=sys.stderr)
	return {'GU_KEEP': str(work / 'chunk.keep'), 'GU_KEEP_PSAM': f'{path},{size},{index}',
		'GU_KEEP_PSAM_INFO': str(work / 'chunk.json'), 'GU_KEEP_PSAM_SHA256': digest}


def context(prefix, keep=None, male=False, panel=None, dataset="ukb"):
	prefix = Path(prefix).expanduser().resolve()
	paths = [Path(str(prefix) + c + '.psam') for c in ['1', '22', 'X']]
	psam = next((p for p in paths if p.is_file()), None)
	if psam is None:
		psam = next(iter(sorted(prefix.parent.glob(prefix.name + '*.psam'))), None)
	if psam is None:
		raise ValueError('Direct PGEN input requires chromosome PSAM files')
	samples = u.read_psam(psam)
	if keep:
		samples = u.selected_samples(samples, u.read_ids(keep), '1')
	if male:
		samples = [s for s in samples if s['sex'] == 'male']
	if not samples:
		raise ValueError('No samples after --keep and PSAM --keep-males selection')
	samples = u.attach_panel(samples, panel)
	cohort = u.validate_cohort(samples, argparse.Namespace(sample_panel=panel)) if dataset == "ukb" else {"sample_metadata_sha256": u.cohort_hash(samples)}
	identity = {'male_only': male, 'cohort': cohort['sample_metadata_sha256'], 'keep_sha256': hashlib.sha256(Path(keep).read_bytes()).hexdigest() if keep else None}
	key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
	work = Path('/tmp/gu-cohorts') / key
	work.mkdir(parents=True, exist_ok=True)
	u.write_panel(work / 'samples.txt', samples)
	ids = '#IID\n' + ''.join(s['sample'] + '\n' for s in samples)
	if not (work / 'cohort.keep').exists() or (work / 'cohort.keep').read_text() != ids:
		(work / 'cohort.keep').write_text(ids)
	source = prefix.parent.name if prefix.parent.name in ('imp', 'hap', 'typ') else 'imp'
	if dataset != 'ukb':
		return {'GU_FILTERED_TARGET': '1', 'TARGET_INPUT': dataset + '-subset-' + key, 'GU_KEEP': str(Path(keep).resolve()) if keep else str(work / 'cohort.keep'), 'GU_SAMPLE_PANEL': str(work / 'samples.txt'), 'SAMPLE_PANEL_INPUT': str(work / 'samples.txt')}
	return {'GU_UKB_DIRECT': '1', 'GU_UKB_PGEN_PREFIX': str(prefix), 'GU_UKB_SOURCE': source,
		'TARGET_INPUT': 'ukb-cohort-' + key, 'GU_KEEP': str(Path(keep).resolve()) if keep else str(work / 'cohort.keep'),
		'GU_SAMPLE_PANEL': str(work / 'samples.txt'), 'SAMPLE_PANEL_INPUT': str(work / 'samples.txt')}


if __name__ == '__main__':
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument('--prefix')
	p.add_argument('--dataset', default='ukb')
	selection = p.add_mutually_exclusive_group()
	selection.add_argument('--keep', type=Path)
	selection.add_argument('--keep-psam')
	p.add_argument('--chunk-only', action='store_true')
	p.add_argument('--sample-panel', type=Path)
	p.add_argument('--keep-males', dest='male_only', action='store_true')
	a = p.parse_args()
	try:
		values = psam_chunk(a.keep_psam) if a.keep_psam else {}
		if a.chunk_only:
			if not a.keep_psam: raise ValueError('--chunk-only requires --keep-psam')
		else:
			if not a.prefix: raise ValueError('--prefix is required')
			values.update(context(a.prefix, values.get('GU_KEEP', a.keep), a.male_only, a.sample_panel, a.dataset))
		for key, value in values.items():
			print('export ' + key + '=' + shlex.quote(value))
	except (ValueError, OSError) as error:
		p.exit(2, f'ERROR: {error}\n')
