"""Runtime tuning must preserve sample groups, masks and resume provenance."""
import ast
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import torch

import prsformer as p
from test_prsformer_cache import commands

ROOT = Path(__file__).resolve().parents[1]


class ResumeRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.current_source = hashlib.sha256(ast.dump(ast.parse(
            Path(p.__file__).read_text()), include_attributes=False).encode()).hexdigest()
        transition = json.loads((ROOT / 'f/prsformer_resume_compatibility.json').read_text())['transitions'][0]
        self.assertEqual(transition['to'], self.current_source)
        self.old_options = dict(epochs=30, patience=3, min_delta=1e-4, lr=5e-4,
            batch_size=1, accumulation_steps=64, amp='fp16', seed=42)
        self.old = dict(identity={'data_sha256': 'frozen data'}, architecture={'layers': 2},
            traits=['height'], preprocessing={'scale': 2.0}, sources={'model': 'frozen source'},
            training=self.old_options, adapter_logic_sha256=transition['from'])
        self.context = dict(training_options=self.old_options,
                            adapter_logic_sha256=transition['from'])
        self.current = copy.deepcopy(self.old)
        self.current['adapter_logic_sha256'] = self.current_source
        self.current['training'].update(batch_size=16, accumulation_steps=4)

    def validate(self, current):
        return p.validate_resume_signature(p.digest_json(self.old), self.context, current)

    def test_reviewed_upgrade_preserves_effective_batch_and_patience(self):
        self.assertEqual(self.validate(self.current), 3)
        self.current['training']['patience'] = 2
        self.assertEqual(self.validate(self.current), 3)

    def test_changes_to_scientific_inputs_and_schedule_still_fail(self):
        changes = [('identity', {'data_sha256': 'other'}), ('architecture', {'layers': 3}),
                   ('traits', ['ldl']), ('preprocessing', {'scale': 3.0}),
                   ('sources', {'model': 'other'})]
        for key, value in changes:
            with self.subTest(key=key):
                changed = copy.deepcopy(self.current)
                changed[key] = value
                with self.assertRaisesRegex(ValueError, 'Resume input'):
                    self.validate(changed)
        for key, value in [('epochs', 12), ('lr', .001), ('amp', 'bf16'), ('seed', 43), ('min_delta', .01)]:
            with self.subTest(key=key):
                changed = copy.deepcopy(self.current)
                changed['training'][key] = value
                with self.assertRaisesRegex(ValueError, 'Resume input'):
                    self.validate(changed)

    def test_effective_batch_and_unreviewed_source_changes_fail(self):
        changed = copy.deepcopy(self.current)
        changed['training']['batch_size'] = 4
        with self.assertRaisesRegex(ValueError, 'subjects per optimizer update'):
            self.validate(changed)
        changed = copy.deepcopy(self.current)
        changed['adapter_logic_sha256'] = 'unreviewed future edit'
        with self.assertRaisesRegex(ValueError, 'adapter source'):
            self.validate(changed)
        self.context['adapter_logic_sha256'] = 'unreviewed old code'
        with self.assertRaisesRegex(ValueError, 'adapter source'):
            self.validate(self.current)

    def test_same_source_resume_works_and_invalid_signature_fails(self):
        self.old['adapter_logic_sha256'] = self.current_source
        self.context['adapter_logic_sha256'] = self.current_source
        self.assertEqual(self.validate(self.current), 3)
        with self.assertRaisesRegex(ValueError, 'Resume input'):
            p.validate_resume_signature('corrupt', self.context, self.current)

    def test_wrapper_default_and_low_memory_override(self):
        for extra, batch, accumulation, workers in [([], 16, 4, 2),
                (['--batch-size', '1', '--accumulation-steps', '64', '--workers', '0'], 1, 64, 0)]:
            planned = commands('prsformer.sh', ['--mode', 'train', *extra])
            train = next(cmd for cmd in planned if len(cmd) > 2 and cmd[2] == 'train')
            for flag, value in [('--batch-size', batch), ('--accumulation-steps', accumulation), ('--workers', workers)]:
                self.assertEqual(train[train.index(flag) + 1], str(value))


class BatchEquivalenceTests(unittest.TestCase):
    def test_prefetch_batching_keeps_order_and_optimizer_groups_including_tail(self):
        torch.set_num_threads(2)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'genotypes.npy'
            rng = np.random.default_rng(17)
            matrix = rng.integers(-1, 3, size=(139, 5)).astype(np.float16)
            np.save(path, matrix)
            rows = np.arange(139)
            labels = rng.normal(size=(139, 2)).astype(np.float32)
            labels[::7, 1] = np.nan
            adjustment = torch.tensor(len(rows) / np.isfinite(labels).sum(axis=0), dtype=torch.float32)
            runs = []
            for batch_size, accumulation_steps, workers in [(1, 64, 0), (8, 8, 2), (16, 4, 2), (32, 2, 2)]:
                args = SimpleNamespace(genotypes=path, batch_size=batch_size,
                    workers=workers, seed=2030, device='cpu')
                torch.manual_seed(5)
                model = torch.nn.Linear(5, 2)
                optimizer = torch.optim.AdamW(model.parameters(), lr=.001)
                loader = p.make_loader(args, rows, torch, shuffle=True, seed=2031)
                seen, groups, group, losses = [], [], [], []
                optimizer.zero_grad(set_to_none=True)
                for batch, (x, indices) in enumerate(loader):
                    np.testing.assert_array_equal(x.numpy(), matrix[indices.numpy()].astype(np.float32))
                    ids = indices.numpy()
                    seen.extend(ids.tolist()); group.extend(ids.tolist())
                    start = batch // accumulation_steps * accumulation_steps * batch_size
                    effective_count = min(accumulation_steps * batch_size, len(rows) - start)
                    loss, mask = p.loss_matrix(model(x), ids, labels, np.zeros_like(labels),
                                               ['height', 'ldl'], torch, torch.device('cpu'))
                    ((loss * adjustment).sum() / effective_count / 2).backward()
                    losses.extend(loss.detach().numpy().tolist())
                    if (batch + 1) % accumulation_steps == 0 or batch + 1 == len(loader):
                        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.)
                        optimizer.step(); optimizer.zero_grad(set_to_none=True)
                        groups.append(group); group = []
                runs.append((seen, groups, losses, model.state_dict()))
            self.assertEqual(runs[0][0], runs[1][0])
            self.assertEqual(runs[0][1], runs[1][1])
            self.assertEqual([len(x) for x in runs[1][1]], [64, 64, 11])
            np.testing.assert_allclose(runs[0][2], runs[1][2], rtol=2e-5, atol=2e-6)
            for name in runs[0][3]:
                torch.testing.assert_close(runs[0][3][name], runs[1][3][name], rtol=2e-5, atol=2e-6)

    def test_training_heartbeat_remains_visible_after_midpoint(self):
        spec = importlib.util.spec_from_file_location('grid_console_test', ROOT.parent / '0f/console_run.py')
        module = importlib.util.module_from_spec(spec)
        with patch.object(sys, 'path', [str(ROOT.parent / '0f'), *sys.path]):
            spec.loader.exec_module(module)
        console = module.PrsformerConsole()
        for count in (1, 50, 51, 52, 60, 90):
            self.assertIsNotNone(console.feed(f'TRAIN epoch 10/30: {count}/100 samples; 20.00 samples/s'))
        self.assertIsNotNone(console.feed('RUNTIME: batch=8; accumulation=8; effective_batch=64'))


if __name__ == '__main__':
    unittest.main()
