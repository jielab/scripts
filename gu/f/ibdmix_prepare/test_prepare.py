#!/usr/bin/env python3
"""Compare the shared reader against the previous bcftools + upstream generate_gt path."""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import random
import shlex
import subprocess
import tempfile
import unittest
import sys

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[2]
PREPARE = os.environ.get("IBDMIX_PREPARE_BIN", "/mnt/d/software/gu/ibdmix-prepare/prepare")
UPSTREAM = os.environ.get("IBDMIX_RUNTIME", "/mnt/d/software/gu/ibdmix")
REFS = ["Altai", "Chagyr", "Denisova", "Denisova25", "Vindija"]
spec = importlib.util.spec_from_file_location("common", ROOT / "f/0.common.py")
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)


def run(cmd, **kwargs):
    return subprocess.run([str(x) for x in cmd], check=True, capture_output=True, **kwargs)


def legacy_stream(path, output, chrom="1", region="", build=None):
    view = ["bcftools", "view", "-Ou"]
    if region:
        view += ["-r", region]
    view += [str(path)]
    commands = [view, ["bcftools", "+setGT", "-Ou", "--", "-t", ".", "-n", "."],
                ["bcftools", "+fixploidy", "-Ou", "--", "-f", "2"],
                ["bcftools", "annotate", "-x", "INFO,^FORMAT/GT", "-Ov"]]
    if chrom == "X":
        fix = ["python3", str(ROOT / "f/0.common.py"), "fix-vcf-gt", "--chrom", "23", "--duplicate-haploid", "slash"]
        if build:
            fix += ["--x-nonpar-build", build]
        commands.append(fix)
    else:
        commands.append(["awk", 'BEGIN{OFS="\t"} /^#/{print;next}{sub(/^chr/,"",$1);print}'])
    with Path(output).open("wb") as handle:
        subprocess.run(["bash", "-o", "pipefail", "-c", " | ".join(shlex.join(c) for c in commands)],
                       stdout=handle, stderr=subprocess.PIPE, check=True)


def write_vcf(path, chrom, samples, rows):
    header = ("##fileformat=VCFv4.2\n"
              f"##contig=<ID={chrom},length=300000000>\n"
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
              '##FILTER=<ID=LowQual,Description="Low quality">\n'
              '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    content = header + "".join(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t{info}\tDP:GT\t" +
                             "\t".join("12:" + gt for gt in gts) + "\n"
                             for pos, ref, alt, gts, flt, info in rows)
    path = Path(path)
    with path.open("wb") as handle:
        subprocess.run(["bgzip", "-c"], input=content.encode(), stdout=handle, check=True)
    run(["bcftools", "index", "-t", path])
    return path


class PrepareTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gu-prepare-test-")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def fixture(self, chrom="1"):
        rng = random.Random(20261012)
        gts = ["0/0", "0|1", "1/1", "1|0", ".", "./.", "0/.", ".|1", "0", "1", "0/1/1", "1/1/."]
        variants = [("A", "G"), ("A", "."), ("AT", "A"), ("A", "C,G"), ("A", "<DEL>"), ("G", "T")]
        modern = [(i, *rng.choice(variants), [rng.choice(gts) for _ in range(12)], "PASS", ".")
                  for i in range(10, 401, 3)]
        # Include the boundaries of every X non-PAR component for both builds.
        if chrom == "X":
            for pos in [60000, 60001, 2699520, 2699521, 154931043, 154931044, 155260560, 155260561,
                        10000, 10001, 2781479, 2781480, 155701382, 155701383, 156030895, 156030896]:
                modern.append((pos, "A", "G", ["1"] * 12, "PASS", "."))
        modern.sort()
        m = write_vcf(self.root / "modern.vcf.gz", "chr" + chrom, [f"S{i}" for i in range(12)], modern)
        refs = {}
        for ref in REFS:
            rows = [(i, *rng.choice(variants), [rng.choice(gts)], "LowQual", ".") for i in range(1, 440)]
            if chrom == "X":
                rows += [(r[0], "A", "G", ["0/0" if r[0] % 2 else "1"], "PASS", ".") for r in modern if r[0] > 440]
            refs[ref] = write_vcf(self.root / f"{ref}.vcf.gz", "23" if chrom == "X" else chrom, [ref], sorted(rows))
        return m, refs

    def compare(self, chrom="1", build=None, region=""):
        modern, refs = self.fixture(chrom)
        out = self.root / "shared"
        command = [PREPARE, "--modern", modern, "--out", out, "--unit", "test", "--chrom", chrom, "--threads", "2"]
        mregion = "chr" + chrom + ":" + region if region else ""
        aregion = ("23" if chrom == "X" else chrom) + ":" + region if region else ""
        if region:
            command += ["--modern-region", mregion]
        if build:
            bed = self.root / "nonpar.bed"
            bed.write_text("".join(f"23\t{lo}\t{hi}\n" for lo, hi in common.x_nonpar_intervals(build)))
            command += ["--nonpar-bed", bed, "--build", build]
        for ref, path in refs.items():
            command += ["--archaic", ref, path, aregion]
        result = run(command)
        self.assertIn(b"modern_passes=1", result.stderr)
        legacy_stream(modern, self.root / "modern.txt", chrom, mregion, build)
        for ref, path in refs.items():
            legacy_stream(path, self.root / "archaic.txt", chrom, aregion, build)
            run([Path(UPSTREAM) / "build/src/generate_gt", "-m", self.root / "modern.txt", "-a", self.root / "archaic.txt", "-o", self.root / "expected.txt"])
            actual = out / f"{ref}.test.gt.txt.gz"
            self.assertEqual(gzip.decompress(actual.read_bytes()), (self.root / "expected.txt").read_bytes(), ref)
            qc = self.root / "qc.json"
            qargs = ["--x-nonpar-build", build] if build else []
            run(["python3", ROOT / "f/ibdmix.py", "validate-genotypes", "--path", actual, "--output", qc, *qargs])
            self.assertEqual(json.loads(qc.read_text()), json.loads((out / f"{ref}.test.qc.json").read_text()), ref)

    def test_autosomes(self):
        self.compare()

    def test_indexed_region(self):
        self.compare(region="90-300")

    def test_x_legacy(self):
        self.compare(chrom="X")

    def test_x_nonpar_37(self):
        self.compare(chrom="X", build="37")

    def test_x_nonpar_38(self):
        self.compare(chrom="X", build="38")

    def test_failure_does_not_commit(self):
        modern, refs = self.fixture()
        bad = self.root / "bad.vcf.gz"
        bad.write_bytes(b"not a VCF")
        out = self.root / "failed"
        result = subprocess.run([PREPARE, "--modern", str(modern), "--out", str(out), "--unit", "test", "--chrom", "1",
                                 "--archaic", "Altai", str(refs["Altai"]), "", "--archaic", "Chagyr", str(bad), ""], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(list(out.glob("*.gt.txt.gz")))
        self.assertFalse(list(out.glob("*.qc.json")))

    def pipeline(self, chrom):
        modern_dir = self.root / "modern"
        archaic_dir = self.root / "archaic"
        masks = self.root / "masks"
        modern_dir.mkdir()
        positions = [1000 + i * 25 for i in range(2100)]
        if chrom == "X":
            positions += [3000000 + i * 25 for i in range(2100)]
        rows = [(pos, "A", "G", ["1/1", "1/1"] + ["0/0"] * 10, "PASS", ".") for pos in positions]
        write_vcf(modern_dir / f"chr{chrom}.vcf.gz", chrom, [f"S{i}" for i in range(12)], rows)
        for ref in REFS:
            (archaic_dir / ref).mkdir(parents=True)
            (masks / ref).mkdir(parents=True)
            (masks / ref / f"chr{chrom}.bed").write_text(f"{chrom}\t20000\t20100\n")
            write_vcf(archaic_dir / ref / f"chr{chrom}.vcf.gz", chrom, [ref],
                      [(pos, "A", "G", ["1/1"], "PASS", ".") for pos in positions])
        panel = self.root / "samples.txt"
        panel.write_text("sample\tpop\tsuper_pop\tsex\n" + "".join(f"S{i}\tUKB\tEUR\tmale\n" for i in range(12)))
        env = {k:v for k,v in os.environ.items() if not k.startswith(("GU_", "IBDMIX_"))}
        env.update(PHE_F=str(ROOT.parent / "0f/phenotype.sh"), GU_TARGET_VCF_DIR=str(modern_dir),
                   dirarch=str(archaic_dir), sample_file=str(panel), dirsoft=UPSTREAM, chrs=chrom,
                   IBDMIX_PROFILE="custom", IBDMIX_MASK_DIR=str(masks), IBDMIX_AFR_DENISOVAN_FILTER="0",
                   refs=" ".join(REFS), len_cut="1", IBDMIX_PREPARE_BIN=PREPARE)
        outputs = []
        for backend in ["cpu", "cuda"]:
            output = self.root / backend
            current_env = dict(env, IBDMIX_BACKEND=backend, dirout=str(output))
            run(["bash", ROOT / "f/ibdmix.sh"], env=current_env)
            self.assertTrue((output / ".complete").is_file())
            outputs.append(output)
        unit = "X_MALE" if chrom == "X" else "C1"
        for ref in REFS:
            raw = Path("raw") / unit / f"{ref}.UKB.raw.txt.gz"
            cpu = gzip.decompress((outputs[0] / raw).read_bytes())
            gpu = gzip.decompress((outputs[1] / raw).read_bytes())
            self.assertEqual(cpu, gpu)
            self.assertGreater(len(cpu.splitlines()), 1, "Fixture must emit actual segments")
        final = Path("final/all_archaic_refs.lod4.len1.segments.tsv.gz")
        self.assertEqual(gzip.decompress((outputs[0]/final).read_bytes()), gzip.decompress((outputs[1]/final).read_bytes()))
        # Resume a partially completed chromosome: only the missing reference
        # should be prepared, and completed raw files must remain untouched.
        saved = outputs[0] / "raw" / unit / "Altai.UKB.raw.txt.gz"
        stamp = saved.stat().st_mtime_ns
        (outputs[0] / "raw" / unit / "Vindija.UKB.raw.txt.gz").unlink()
        (outputs[0] / ".complete").unlink()
        run(["bash", ROOT / "f/ibdmix.sh"], env=dict(env, IBDMIX_BACKEND="cpu", dirout=str(outputs[0])))
        self.assertEqual(stamp, saved.stat().st_mtime_ns)
        self.assertIn("references=1", (outputs[0]/"logs"/f"{unit}.prepare.log").read_text())

    def test_pipeline_cpu_gpu_and_resume(self):
        self.pipeline("1")

    def test_pipeline_x_cpu_gpu_and_resume(self):
        self.pipeline("X")


if __name__ == "__main__":
    unittest.main(verbosity=2)
