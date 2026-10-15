"""Resume without genotype export; reject changed scientific inputs."""
import contextlib
import gzip
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import ibdmix_resume as resume


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        # These are disposable paths of the same kind the real driver creates.
        cohort_root = Path('/tmp/gu-cohorts')
        cohort_root.mkdir(exist_ok=True)
        self.cohort_tmp = tempfile.TemporaryDirectory(dir=cohort_root)
        self.addCleanup(self.cohort_tmp.cleanup)
        self.selection = Path(self.cohort_tmp.name)
        self.run = self.root / 'run'
        self.run.mkdir()
        self.pgen_root = self.root / 'imp'
        self.pgen_root.mkdir()
        self.psam = self.pgen_root / 'chr1.psam'
        self.psam.write_text('#IID SEX\nA 1\nB 2\nC 1\n')
        self.pgen = self.pgen_root / 'chr1.pgen'
        self.pgen.write_bytes(b'genotypes')
        (self.pgen_root / 'chr1.pvar').write_text('variants')
        self.keep = self.selection / 'cohort.keep'
        self.keep.write_text('#IID\nA\nB\n')
        self.panel = self.selection / 'samples.txt'
        self.panel.write_text('sample pop super_pop sex\nA UKB ALL male\n')
        self.fasta = self.root / 'ref.fa'
        self.fasta.write_text('>1\nACGT\n')
        Path(str(self.fasta) + '.fai').write_text('1\t4\n')
        self.archaic = self.root / 'vcf'
        (self.archaic / 'Altai').mkdir(parents=True)
        source = self.archaic / 'Altai/chr1.vcf.gz'
        source.write_bytes(b'archaic')
        self.mask = self.root / 'mask'
        self.mask.mkdir()
        mask_input = self.mask / 'input.bed'
        mask_input.write_text('1\t0\t1\n')
        prepare, caller = self.root / 'prepare', self.root / 'caller'
        prepare.write_text('prepare binary')
        caller.write_text('caller binary')
        samples = [resume.ukb.read_psam(self.psam)[0]]
        self.prep = dict(inputs=resume.ukb.signature([
            self.pgen, self.psam, self.pgen_root / 'chr1.pvar', self.keep, self.panel,
            self.fasta, Path(str(self.fasta) + '.fai'), resume.F / 'ukb_prepare.py']),
            options=dict(existing_pgen=True, male_only=True, grch=37, ukb_source='imp', ukb_require_phase=False),
            pgen_root=str(self.pgen_root), sample_order_sha256=resume.ukb.sample_hash(samples),
            cohort_metadata_sha256=resume.ukb.cohort_hash(samples), scored_sample_metadata_sha256=resume.ukb.cohort_hash(samples))
        self.scoring = dict(contract=dict(schema='ukb-observed-sites-v1', chromosome='1', build=37,
            references=['Altai'], original_mask_root=str(self.mask), preparation=self.prep,
            postqc_sites_sha256='1' * 64, inputs=resume.ukb.signature([source, resume.F / 'ukb_prepare.py'])))
        self.scoring_path = self.run / 'inputs/ukb-scoring/chr1/target.json'
        self.scoring_path.parent.mkdir(parents=True)
        self.scoring_path.write_text(json.dumps(self.scoring))
        mask_manifest = self.run / 'qc/masks/C1.Altai.preparation.json'
        mask_manifest.parent.mkdir(parents=True)
        stat = mask_input.stat()
        mask_manifest.write_text(json.dumps(dict(signature=dict(code=resume.sha(resume.F / 'ibdmix.py'),
            refs=['Altai'], inputs=[[str(mask_input), stat.st_size, stat.st_mtime_ns]]))))
        fields = dict(pipeline_version='2026-10-12.shared-hts-v1.cuda-segment-v1', profile='multi_reference',
            refs='Altai', genome_build='b37', background_filter='0', custom_masks='', lod_cut='4', len_cut='50000',
            emit_lod_cut='4', minor_allele_count='1', archaic_error='0.01', modern_error_max='0.002',
            modern_error_proportion='2', genotype_binary_sha256=resume.sha(prepare), gpu_binary_sha256=resume.sha(caller),
            reference_fasta=f'{self.fasta}:{self.fasta.stat().st_size}:{int(self.fasta.stat().st_mtime)}')
        self.meta = ''.join(f'{k}\t{v}\n' for k, v in fields.items())
        self.meta += 'excluded_mask_sha256\tC1\tAltai\t' + '0' * 64 + '\n'
        self.meta += 'modern_source\tchr1\tpreparation_contract_sha256\t' + resume.digest(json.dumps(self.prep, sort_keys=True).encode()) + '\n'
        self.meta += 'sample_file\t/tmp/old/samples.txt:50:1\narchaic_vcf\t/tmp/old/Altai/chr1.vcf.gz:100:1\n'
        (self.run / 'run.meta.tsv').write_text(self.meta)
        (self.run / '.complete').touch()
        self.final = self.run / 'final/all_archaic_refs.lod4.len50000.segments.tsv.gz'
        self.final.parent.mkdir()
        self.final.write_bytes(gzip.compress(b'ID\tchr\tstart\tend\nA\t1\t1\t60000\n'))
        self.env = dict(GU_UKB_DIRECT='1', GU_UKB_SOURCE='imp', IBDMIX_OUT=str(self.run),
            IBDMIX_REFS='Altai', IBDMIX_PROFILE='multi_reference', IBDMIX_AFR_DENISOVAN_FILTER='0',
            IBDMIX_PREPARE_BIN=str(prepare), IBDMIX_GPU_BIN=str(caller), IBDMIX_REFERENCE_FASTA=str(self.fasta),
            UKB_REF_FASTA=str(self.fasta))
        self.opts = ['--chr', '1', '--grch', '37', '--target-dir', str(self.pgen_root / 'chr'),
            '--keep', str(self.keep), '--keep-males', '--sample-panel', str(self.panel),
            '--archaic-path', str(self.archaic), '--replace-ibdmix', 'FALSE']
        self.cmd = self.run / 'chr1.cmd'
        self.write_command()

    def write_command(self):
        self.cmd.write_text(''.join('export ' + k + '=' + shlex.quote(v) + '\n' for k, v in self.env.items())
                           + shlex.join(['exec', '/path/gu.sh', 'ibdmix-GPU', *self.opts]) + '\n')

    def check(self):
        with patch.dict(os.environ, {}, clear=True), contextlib.redirect_stdout(io.StringIO()):
            return resume.check(self.cmd)

    def test_recreated_temporary_selection_files_reuse_without_exports(self):
        self.assertEqual(self.check(), 0)
        for path in [self.panel, self.keep]:
            data = path.read_bytes()
            path.unlink()
            path.write_bytes(data)
        self.assertEqual(self.check(), 0)

    def test_changed_selection_is_rejected_even_with_same_file_size(self):
        self.keep.write_text('#IID\nC\nB\n')
        with self.assertRaisesRegex(ValueError, 'cohort|missing selected'):
            self.check()

    def test_changed_source_and_threshold_are_rejected(self):
        self.env['IBDMIX_LOD'] = '5'
        self.write_command()
        with self.assertRaisesRegex(ValueError, 'lod_cut'):
            self.check()
        del self.env['IBDMIX_LOD']
        self.write_command()
        self.pgen.write_bytes(b'new genotypes')
        with self.assertRaisesRegex(ValueError, 'input changed'):
            self.check()

    def test_changed_archaic_or_mask_inputs_are_rejected(self):
        source = self.archaic / 'Altai/chr1.vcf.gz'
        previous = source.read_bytes(), source.stat()
        source.write_bytes(b'new archaic')
        with self.assertRaisesRegex(ValueError, 'input changed'):
            self.check()
        source.write_bytes(previous[0])
        os.utime(source, ns=(previous[1].st_atime_ns, previous[1].st_mtime_ns))
        (self.mask / 'input.bed').write_text('1\t0\t2\n')
        with self.assertRaisesRegex(ValueError, 'input changed'):
            self.check()

    def test_corrupted_or_missing_final_is_not_skipped(self):
        self.assertEqual(self.check(), 0)
        self.final.write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'corrupted'):
            self.check()
        self.final.unlink()
        with self.assertRaisesRegex(ValueError, 'missing'):
            self.check()

    def test_partial_results_stay_pending_and_replace_bypasses_reuse(self):
        (self.run / '.complete').unlink()
        self.assertEqual(self.check(), 1)
        (self.run / '.complete').touch()
        self.assertEqual(self.check(), 0)
        self.opts[-1] = 'TRUE'
        self.write_command()
        self.assertEqual(self.check(), 1)

    def test_interrupted_final_without_completion_can_resume(self):
        (self.run / '.complete').unlink()
        self.final.write_bytes(b'interrupted gzip')
        self.assertEqual(self.check(), 1)

    def test_archive_survives_overwritten_working_scoring_receipt(self):
        paths = [p for p in self.run.rglob('*') if p.is_file() and p != self.cmd]
        manifest = {'files': [{'path': str(p.relative_to(self.run)), 'sha256': resume.sha(p)} for p in paths]}
        with tarfile.open(self.run / 'ibdmix.raw.tar.gz', 'w:gz') as bundle:
            for path in paths:
                bundle.add(path, arcname=str(path.relative_to(self.run)))
            data = json.dumps(manifest).encode()
            member = tarfile.TarInfo('GU-MANIFEST.json')
            member.size = len(data)
            bundle.addfile(member, io.BytesIO(data))
        self.scoring_path.write_text('{}')
        self.assertEqual(self.check(), 0)
        self.final.unlink()  # Compacted results are also reusable.
        self.assertEqual(self.check(), 0)

    def test_partial_provenance_accepts_only_exchange_identity_changes(self):
        (self.run / '.complete').unlink()
        self.assertEqual(self.check(), 1)
        old = self.run / 'run.meta.tsv'
        new = self.root / 'candidate.tsv'
        changed = json.loads(json.dumps(self.scoring))
        prep = changed['contract']['preparation']
        for item in prep['inputs']:
            if resume.disposable(item['path']):
                item['mtime_ns'] += 100000
        self.scoring_path.write_text(json.dumps(changed))
        old_hash = resume.digest(json.dumps(self.prep, sort_keys=True).encode())
        new_hash = resume.digest(json.dumps(prep, sort_keys=True).encode())
        current = self.meta.replace(old_hash, new_hash).replace('/tmp/old/', '/tmp/new/').replace(':100:1', ':104:1')
        new.write_text(current)
        self.assertEqual(resume.provenance(old, new, self.run), 0)
        new.write_text(current.replace('lod_cut\t4', 'lod_cut\t5'))
        self.assertEqual(resume.provenance(old, new, self.run), 1)
        changed['contract']['postqc_sites_sha256'] = '2' * 64
        self.scoring_path.write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError, 'scored sites changed'):
            resume.provenance(old, new, self.run)

    def scheduler(self):
        driver = (resume.F.parent / 'gu.sh').read_text()
        functions = driver[driver.index('\tgu_completed_analysis_cmd_output() {'):driver.index('\tgu_run_analysis_cmds_local() {')]
        script = 'set -euo pipefail\n' + functions + '\n'
        script += 'F=' + shlex.quote(str(resume.F)) + '\nMETHOD=ibdmix\nREPLACE_IBDMIX_INPUT=FALSE\n'
        script += 'bash(){ echo UNEXPECTED_WORKER_LAUNCH >&2; return 97; }\n'
        script += 'gu_run_one_analysis_cmd ' + shlex.quote(str(self.cmd)) + '\n'
        return subprocess.run(['/bin/bash', '-c', script], env={'PATH': os.defpath}, text=True, capture_output=True)

    def test_scheduler_skips_before_starting_any_worker(self):
        (self.run / 'chr1.err').write_text('stale error')
        result = self.scheduler()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('SKIP unit=chr1', result.stderr)
        self.assertNotIn('UNEXPECTED_WORKER_LAUNCH', result.stderr)
        self.assertFalse((self.run / 'chr1.err').exists())

    def test_scheduler_fails_changed_input_without_starting_worker(self):
        self.pgen.write_bytes(b'changed genotypes')
        result = self.scheduler()
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn('reason=resume_validation', result.stderr)
        self.assertNotIn('UNEXPECTED_WORKER_LAUNCH', result.stderr)


if __name__ == '__main__':
    unittest.main()
