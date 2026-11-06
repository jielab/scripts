import unittest
import numpy as np
import torch
from block_dosage import SupplementedDosage
from block_scoring import BatchReader, score_batches


class BlockScoringTests(unittest.TestCase):
    def test_bounded_reader_preserves_values_order_and_tail(self):
        original = np.arange(51, dtype=np.float16).reshape(17, 3)
        supplement = np.asfortranarray(np.arange(85, dtype=np.float32).reshape(17, 5))
        combined = SupplementedDosage(original, supplement)
        reader = BatchReader(combined, 4, buffer_bytes=160)
        expected = np.concatenate((original, supplement), axis=1)
        for start, end in [(0, 4), (4, 8), (8, 12), (16, 17), (0, 3), (12, 16)]:
            np.testing.assert_array_equal(reader.read(start, end), expected[start:end])
            self.assertLessEqual(reader.buffer.nbytes, 160)
            self.assertEqual(reader.host.dtype, np.float32)
        reader = BatchReader(combined, 4, buffer_bytes=0)
        np.testing.assert_array_equal(reader.read(2, 5), expected[2:5])
        self.assertIsNone(reader.buffer)

    def test_real_scoring_semantics_missing_ref_alt_duplicate_tokens_and_channels(self):
        values = np.array([[0, 1, 2, -1], [2, -1, .5, 1], [np.nan, .25, 1.5, 2]], dtype=np.float32)
        maps = [dict(channel=0, cache_col=np.array([0, 1, 2]), token_index=np.array([0, 0, 1]),
                     coef=np.array([.1, -.2, .3], np.float32), offset=np.array([0, .4, 0], np.float32)),
                dict(channel=1, cache_col=np.array([3, 0]), token_index=np.array([0, 1]),
                     coef=np.array([.5, .25], np.float32), offset=np.zeros(2, np.float32))]
        expected = np.array([[[.2, 0], [.6, 0]], [[.2, .5], [.15, .5]],
                             [[.35, 1], [.45, 0]]], dtype=np.float32)
        for layout in ['plain', 'F', 'C']:
            geno = values if layout=='plain' else SupplementedDosage(values[:, :2], np.array(values[:, 2:], order=layout))
            for device in ['cpu']+(['cuda'] if torch.cuda.is_available() else []):
                with self.subTest(layout=layout, device=device):
                    result = np.concatenate([v for _, _, v in score_batches(geno, maps, 2, 2, 2, device, buffer_bytes=64)])
                    np.testing.assert_allclose(result, expected, atol=1e-7, rtol=1e-6)

    def test_invalid_batch_and_out_of_range_read(self):
        with self.assertRaises(ValueError):
            BatchReader(np.ones((3, 2)), 0)
        reader = BatchReader(np.ones((3, 2)), 2)
        with self.assertRaises(ValueError):
            reader.read(0, 3)


if __name__ == '__main__':
    unittest.main()
