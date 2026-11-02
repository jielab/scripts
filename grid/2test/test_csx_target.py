import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from csx_target import build_target_bim


class TargetBimTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.target = self.root / "genotypes"
        self.target.mkdir()

    def pvar(self, chromosome=1, rows="1\t100\trs1\tA\tG\n"):
        path = self.target / f"chr{chromosome}.pvar"
        path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n" + rows)
        return path

    def test_pvar_conversion_uses_alt_ref_and_filters_unsupported_variants(self):
        self.pvar(rows="1\t100\trs1\tA\tG\n1\t101\tindel\tA\tAT\n"
                       "1\t102\tmulti\tA\tC,G\n1\t103\t.\tG\tA\n")
        self.pvar(2, "2\t200\trs2\tC\tT\n")
        prefix, meta = build_target_bim(self.target, [1, 2])
        self.assertEqual(prefix.with_suffix(".bim").read_text(), "1\trs1\t0\t100\tG\tA\n2\trs2\t0\t200\tT\tC\n")
        self.assertEqual(meta["snps"], 2)
        self.assertEqual(meta["chromosomes"][0]["excluded_non_snp"], 2)
        self.assertEqual(meta["chromosomes"][0]["excluded_missing_id"], 1)
        self.assertEqual(prefix, self.target / "ukb_array")
        self.assertTrue((self.target / "ukb_array.bim.json").is_file())

    def test_reordered_pvar_columns_and_chr_prefix(self):
        (self.target / "chr1.pvar").write_text("#CHROM\tID\tREF\tALT\tPOS\nchr1\trs1\tc\tt\t100\n")
        prefix, _ = build_target_bim(self.target, [1])
        self.assertEqual(prefix.with_suffix(".bim").read_text(), "1\trs1\t0\t100\tT\tC\n")

    def test_duplicate_ids_do_not_publish_partial_bim(self):
        self.pvar(rows="1\t100\trs1\tA\tG\n1\t200\trs1\tA\tC\n")
        with self.assertRaisesRegex(ValueError, "Duplicate target SNP ID"):
            build_target_bim(self.target, [1])
        self.assertFalse((self.target / "ukb_array.bim").exists())
        self.assertFalse(list(self.target.glob("*.tmp")))

    def test_missing_or_wrong_chromosome_is_rejected(self):
        self.pvar(rows="2\t100\trs1\tA\tG\n")
        with self.assertRaisesRegex(ValueError, "chromosome/position"):
            build_target_bim(self.target, [1])
        with self.assertRaises(FileNotFoundError):
            build_target_bim(self.target, [22])

    def test_permanent_reuse_content_change_and_corruption(self):
        path = self.pvar()
        prefix, _ = build_target_bim(self.target, [1])
        bim = prefix.with_suffix(".bim")
        before = bim.stat().st_mtime_ns
        self.assertEqual(build_target_bim(self.target, [1])[0], prefix)
        self.assertEqual(bim.stat().st_mtime_ns, before)
        bim.write_text("damaged\n")
        build_target_bim(self.target, [1])
        self.assertEqual(bim.read_text(), "1\trs1\t0\t100\tG\tA\n")
        # Same file size and timestamp must not disguise a changed allele.
        stamp = path.stat()
        path.write_text(path.read_text().replace("\tA\tG\n", "\tA\tT\n"))
        os.utime(path, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
        changed, _ = build_target_bim(self.target, [1])
        self.assertEqual(changed, prefix)
        self.assertEqual(changed.with_suffix(".bim").read_text(), "1\trs1\t0\t100\tT\tA\n")

    @unittest.skipUnless(shutil.which("zstd"), "zstd is unavailable")
    def test_zstd_pvar_and_legacy_bim(self):
        path = self.pvar()
        compressed = path.with_suffix(".pvar.zst")
        subprocess.run(["zstd", "-q", str(path), "-o", str(compressed)], check=True)
        path.unlink()
        (self.target / "chr2.bim").write_text("2\trs2\t0.5\t200\tT\tC\n")
        prefix, _ = build_target_bim(self.target, [1, 2])
        self.assertEqual(prefix.with_suffix(".bim").read_text(), "1\trs1\t0\t100\tG\tA\n2\trs2\t0\t200\tT\tC\n")

    def test_shell_automatic_and_explicit_prefix(self):
        self.pvar()
        script = '''
source "$1/f/0.common.sh"
grid_configure
grid_parse_args --dir-gen "$2" --output-root "$3" --chrs 1 "${@:4}"
grid_prepare_csx_bim
printf '%s\\n' "$GRID_CSX_BIM_PREFIX"
'''
        env = {k: v for k, v in os.environ.items() if not k.startswith("GRID_")}
        env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
        command = ["bash", "-c", script, "test", str(ROOT), str(self.target), str(self.root / "analysis")]
        completed = subprocess.run(command, env=env, check=True, capture_output=True, text=True)
        generated = Path(completed.stdout.strip())
        self.assertEqual(generated, self.target / "ukb_array")
        self.assertFalse((self.root / "analysis").exists())
        before = generated.with_suffix(".bim").stat().st_mtime_ns
        # An absent analysis/cache directory must not force BIM regeneration.
        subprocess.run(command, env=env, check=True, capture_output=True, text=True)
        self.assertEqual(generated.with_suffix(".bim").stat().st_mtime_ns, before)
        custom = self.root / "selected.bim"
        custom.write_text("1\trs_custom\t0\t123\tC\tA\n")
        completed = subprocess.run(command + ["--bim-prefix", str(custom.with_suffix(""))],
                                   env=env, check=True, capture_output=True, text=True)
        self.assertEqual(completed.stdout.strip(), str(custom.with_suffix("")))
        self.assertEqual(custom.read_text(), "1\trs_custom\t0\t123\tC\tA\n")

    def test_conversion_matches_plink_make_just_bim(self):
        plink = shutil.which("plink2")
        if plink is None:
            self.skipTest("plink2 unavailable")
        path = self.pvar(rows="1\t100\trs1\tA\tG\n1\t101\trs2\tT\tC\n1\t102\tindel\tA\tAT\n")
        prefix, _ = build_target_bim(self.target, [1])
        expected = self.root / "plink"
        subprocess.run([plink, "--pvar", str(path), "--make-just-bim", "--snps-only", "just-acgt",
                        "--max-alleles", "2", "--out", str(expected)], check=True, capture_output=True, text=True)
        self.assertEqual(prefix.with_suffix(".bim").read_text().split(), expected.with_suffix(".bim").read_text().split())

    def test_subset_runs_keep_all_available_autosomes(self):
        self.pvar()
        self.pvar(2, "2\t200\trs2\tC\tT\n")
        prefix, meta = build_target_bim(self.target, [1])
        self.assertEqual(meta["contract"]["chromosomes"], [1, 2])
        bim = prefix.with_suffix(".bim")
        before = bim.stat().st_mtime_ns
        _, meta = build_target_bim(self.target, [2])
        self.assertEqual(meta["snps"], 2)
        self.assertEqual(bim.stat().st_mtime_ns, before)

    def test_unmanaged_bim_is_preserved_and_failed_rebuild_keeps_previous_file(self):
        path = self.pvar()
        bim = self.target / "ukb_array.bim"
        original = "1\trs_custom\t0\t123\tC\tA\n"
        bim.write_text(original)
        with self.assertRaisesRegex(ValueError, "no matching generated provenance"):
            build_target_bim(self.target, [1])
        self.assertEqual(bim.read_text(), original)
        # An identical metadata-only file can acquire provenance after a crash.
        expected = "1\trs1\t0\t100\tG\tA\n"
        bim.write_text(expected)
        build_target_bim(self.target, [1])
        path.write_text("invalid input\n")
        with self.assertRaises(ValueError):
            build_target_bim(self.target, [1])
        self.assertEqual(bim.read_text(), expected)
        self.assertFalse(list(self.target.glob("*.tmp")))


if __name__ == "__main__":
    unittest.main()
