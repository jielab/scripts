"""Fresh and prepared HERV sources have different runtime dependencies."""
import contextlib
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
SPEC = importlib.util.spec_from_file_location("herv_source_audit_engine", ROOT / "f" / "1.herv.py")
HERV = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = HERV
SPEC.loader.exec_module(HERV)
CHROMS = ["chr" + str(i) for i in range(1, 23)]


class TestHervSourceAudit(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="herv-source-audit-")
        self.root = Path(self.temp.name)
        self.out = self.root / "result"
        self.out.mkdir()
        self.fasta = self.root / "current.fa"
        self.write_fasta(100)
        self.args = SimpleNamespace(
            out=self.out, effective_input="pgen", input_format="auto", sample="S1",
            pfile_dir=self.root / "unused-default-pgen", vcf_dir=self.root / "unused-default-vcf",
            fasta=self.fasta, repeatmasker=self.root / "unused-repeatmasker", jobs=1,
            chunk_bp=1000, halo=100, scratch=self.root / "scratch", library=self.root / "library.fa",
        )
        self.env = mock.patch.dict(os.environ, {
            "ALLELE_ACTIVE_HERV_BUILD": "38", "ALLELE_FASTA37": str(self.root / "unused-37.fa"),
            "ALLELE_FASTA": str(self.root / "unused-38.fa"),
        })
        self.env.start()

    def tearDown(self):
        self.env.stop()
        self.temp.cleanup()

    def write_fasta(self, length):
        # A valid existing FAI, created without indexing a large reference or
        # requiring a subprocess. Audit is allowed to read it but not edit it.
        offset, lines, index = 0, [], []
        for chrom in CHROMS:
            header = ">" + chrom + "\n"
            lines.extend([header, "A" * length + "\n"])
            offset += len(header)
            index.append(f"{chrom}\t{length}\t{offset}\t{length}\t{length + 1}\n")
            offset += length + 1
        self.fasta.write_text("".join(lines))
        Path(str(self.fasta) + ".fai").write_text("".join(index))

    def pvar(self, path, chrom="chr1", build="38"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"##reference=GRCh{build}\n#CHROM\tPOS\tID\tREF\tALT\n{chrom}\t10\tx\tA\tG\n")

    def prepared(self, kind="pgen", existing_first=False):
        self.args.effective_input = kind
        (self.out / "prepared.json").write_text(json.dumps(dict(
            sample="S1", input_source=kind, chromosomes=CHROMS, fasta=str(self.fasta),
        )))
        for chrom in CHROMS:
            folder = self.out / "variants" / chrom
            folder.mkdir(parents=True)
            path = self.root / "recorded-original" / (chrom + (".pvar" if kind == "pgen" else ".vcf.gz"))
            if existing_first and chrom == "chr1":
                if kind == "pgen":
                    self.pvar(path)
                else:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with gzip.open(path, "wt") as handle:
                        handle.write("##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t10\tx\tA\tG\n")
            size = path.stat().st_size if path.is_file() else 123
            mtime = path.stat().st_mtime_ns if path.is_file() else 1
            inputs = dict(length=100, sample="S1")
            if kind == "pgen":
                inputs["sources"] = {"pvar": dict(path=str(path), bytes=size, mtime_ns=mtime)}
            else:
                inputs.update(vcf=str(path), vcf_bytes=size, vcf_mtime_ns=mtime)
            (folder / "done.json").write_text(json.dumps(dict(inputs=inputs)))

    def audit(self):
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            HERV.audit_census_source(self.args)
        return stream.getvalue()

    def test_prepared_pgen_missing_originals_do_not_require_unused_filesets(self):
        self.prepared()
        before = {p: p.read_bytes() for p in self.out.rglob("*.json")}
        steps = HERV.commands(self.args)
        self.assertEqual(len(steps), 1)
        self.assertIn("run", steps[0])
        self.assertNotIn("prepare", steps[0])
        with mock.patch.object(HERV, "source_files", side_effect=AssertionError("unused CLI fileset opened")):
            text = self.audit()
        self.assertIn("deferred_missing", text)
        self.assertNotIn(str(self.args.pfile_dir), text)
        self.assertEqual(before, {p: p.read_bytes() for p in self.out.rglob("*.json")})

    def test_prepared_vcf_uses_receipt_path_not_unused_cli_directory(self):
        self.prepared("vcf", existing_first=True)
        text = self.audit()
        self.assertIn(str(self.root / "recorded-original" / "chr1.vcf.gz"), text)
        self.assertNotIn(str(self.args.vcf_dir), text)

    def test_prepared_pvar_metadata_does_not_need_original_pgen_or_psam(self):
        self.prepared(existing_first=True)
        text = self.audit()
        self.assertIn("metadata_supports_38", text)
        self.assertFalse((self.root / "recorded-original" / "chr1.pgen").exists())

    def test_changed_original_is_deferred_instead_of_applied_to_cached_model(self):
        self.prepared(existing_first=True)
        path = self.root / "recorded-original" / "chr1.pvar"
        # A replacement source is not evidence about the already prepared
        # model. Make the identity change definite even on coarse filesystems.
        self.pvar(path, build="37")
        with path.open("a") as handle:
            handle.write("\n")
        text = self.audit()
        self.assertIn("original source identity changed/unverified; deferred", text)

    def test_recorded_matching_original_build_conflict_still_blocks(self):
        self.prepared(existing_first=True)
        path = self.root / "recorded-original" / "chr1.pvar"
        self.pvar(path, build="37")
        receipt = self.out / "variants" / "chr1" / "done.json"
        value = json.loads(receipt.read_text())
        value["inputs"]["sources"]["pvar"].update(bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns)
        receipt.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "build conflicts"):
            self.audit()

    def test_current_fasta_lengths_still_bound_to_prepared_reference(self):
        self.prepared()
        self.write_fasta(101)
        with self.assertRaisesRegex(ValueError, "FASTA lengths differ"):
            self.audit()

    def test_prepared_canonical_lengths_conflicting_with_build_block(self):
        self.prepared()
        receipt = self.out / "variants" / "chr1" / "done.json"
        value = json.loads(receipt.read_text())
        value["inputs"]["length"] = 249250621  # canonical GRCh37 chr1
        receipt.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "Prepared reference chromosome lengths conflict"):
            self.audit()

    def test_fresh_selected_pgen_build_conflict_blocks_without_reading_vcf(self):
        self.args.pfile_dir.mkdir()
        for chrom in CHROMS:
            prefix = self.args.pfile_dir / chrom
            self.pvar(Path(str(prefix) + ".pvar"), chrom, "37" if chrom == "chr1" else "38")
            Path(str(prefix) + ".pgen").write_bytes(b"metadata audit only")
            Path(str(prefix) + ".psam").write_text("#IID\nS1\n")
        with self.assertRaisesRegex(ValueError, "build conflicts"):
            self.audit()
        self.assertFalse(self.args.vcf_dir.exists())

    def test_fresh_vcf_never_opens_unused_pgen(self):
        self.args.effective_input = "vcf"
        self.args.vcf_dir.mkdir()
        for chrom in CHROMS:
            with gzip.open(self.args.vcf_dir / (chrom + ".vcf.gz"), "wt") as handle:
                handle.write(f"##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\n{chrom}\t10\tx\tA\tG\n")
        with mock.patch.object(HERV, "source_files", side_effect=AssertionError("unused PGEN opened")):
            text = self.audit()
        self.assertIn("metadata_supports_38", text)

    def test_current_fasta_canonical_build_conflict_blocks_even_with_missing_originals(self):
        self.prepared()
        # Deliberately contradictory existing index metadata; the audit must
        # reject it before any attempt to use or rewrite the FASTA/index.
        Path(str(self.fasta) + ".fai").write_text("chr1\t249250621\t6\t100\t101\n")
        with self.assertRaisesRegex(ValueError, "FASTA canonical lengths conflict"):
            self.audit()


if __name__ == "__main__":
    unittest.main()
