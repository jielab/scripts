import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import blocknet as gpu
import block_inputs as blocks


def digest(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


class HashReuseTest(unittest.TestCase):
    def setUp(self):
        blocks._DIGESTS.clear()

    def test_metadata_scans_do_not_evict_genotype_hash(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            genotype = root/'genotypes.npy'
            genotype.write_bytes(b'genotype')
            with patch.object(blocks, 'sha256', wraps=blocks.sha256) as scan:
                expected = blocks.sha256_cached(genotype)
                for i in range(66):
                    p = root/str(i)
                    p.write_bytes(str(i).encode())
                    blocks.sha256_cached(p)
                self.assertEqual(blocks.sha256_cached(genotype), expected)
                self.assertEqual(sum(call.args[0] == genotype for call in scan.call_args_list), 1)
                # Same-size replacement with restored mtime must still invalidate.
                stamp = genotype.stat()
                replacement = root/'replacement'
                replacement.write_bytes(b'changed!')
                os.utime(replacement, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
                os.replace(replacement, genotype)
                self.assertNotEqual(blocks.sha256_cached(genotype), expected)
                self.assertEqual(sum(call.args[0] == genotype for call in scan.call_args_list), 2)

    def test_changed_during_hash_is_never_cached(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'matrix'
            p.write_bytes(b'original')
            def scan(path):
                result = hashlib.sha256(path.read_bytes()).hexdigest()
                path.write_bytes(b'changed')
                return result
            with patch.object(blocks, 'sha256', side_effect=scan):
                with self.assertRaisesRegex(RuntimeError, 'changed during'):
                    blocks.sha256_cached(p)
            self.assertNotIn(p, blocks._DIGESTS)

    def test_receipt_reuse_requires_exact_source_inputs_and_qc(self):
        payload = dict(source_sha256='current', genotype_sha256='verified',
                       weights={'AFR': 'weights'}, qc_config={'min_fraction': .95}, layout='layout')
        receipt = dict(signature=digest(payload), signature_payload=payload, status='complete')
        self.assertTrue(blocks.receipt_matches(receipt, payload, digest(payload)))
        for key, value in [('source_sha256', 'changed'), ('genotype_sha256', 'changed'),
                           ('weights', {'AFR': 'changed'}), ('qc_config', {'min_fraction': .99}),
                           ('layout', 'changed')]:
            changed = dict(payload, **{key: value})
            self.assertFalse(blocks.receipt_matches(receipt, changed, digest(changed)))


class RuntimeTest(unittest.TestCase):
    def test_scaler_uses_build_rows_only(self):
        x = np.arange(48, dtype=np.float32).reshape(8, 6, 1)
        before = gpu.feature_scaler(x, np.arange(4))
        x[4:] = 1e20
        after = gpu.feature_scaler(x, np.arange(4))
        for a, b in zip(before, after):
            np.testing.assert_array_equal(a, b)
        np.testing.assert_allclose(before[0], x[:4].mean(axis=0))
        np.testing.assert_allclose(before[1], x[:4].std(axis=0))

    def test_scaler_rejects_nonfinite_build_blocks(self):
        x = np.ones((4, 6, 1), dtype=np.float32)
        x[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            gpu.feature_scaler(x, np.arange(4))

    def test_transformer_forward_backward_cpu(self):
        torch = gpu.torch
        tokens = pd.DataFrame(dict(CHR=[1, 1, 2, 2],
            source_population=["AFR", "EAS", "EUR", "SAS"], midpoint=[1., 2., 3., 4.], length=1.))
        cfg = dict(d_model=8, heads=2, dropout=0., latent_tokens=4, layers=1)
        model = gpu.BlockTransformer(tokens, 1, 18, cfg)
        y, embedding = model(torch.randn(3, 4, 1), torch.randn(3, 18))
        self.assertEqual(tuple(y.shape), (3,))
        self.assertEqual(tuple(embedding.shape), (3, 8))
        y.square().mean().backward()
        self.assertTrue(all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters()))

    def test_ram_staging_preserves_values_and_limits_memory(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'matrix.npy'
            values = np.arange(96, dtype=np.float32).reshape(8, 12, 1)
            np.save(p, values)
            a = np.load(p, mmap_mode='r')
            with patch.object(gpu, 'available_memory', return_value=a.nbytes*8):
                staged = gpu.stage_block_matrix(a, {'block_ram_cache_gib': 1})
                self.assertIsNot(staged, a)
                np.testing.assert_array_equal(staged, values)
                self.assertFalse(staged.flags.writeable)
                self.assertIs(gpu.stage_block_matrix(a, {'block_ram_cache_gib': 0}), a)
            with patch.object(gpu, 'available_memory', return_value=a.nbytes*2):
                self.assertIs(gpu.stage_block_matrix(a, {'block_ram_cache_gib': 1}), a)



if __name__ == '__main__':
    unittest.main()
