"""Exercise the real R plotter, reports and MAGMA colors with a shared protein map."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
class ProteinCis(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="protein-cis-")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.names = self.work / "shared protein names.tsv"
        self.names.write_text(
            "prot\tCHR\tPOS\tukb\tckb\tdecode\tscot\n"
            "CANON\t1\t2000000\tUKB\tCKB\tDECODE\tASSAY;ASSAY.1;NA\n"
            "OTHER\t2\t7000000\t.\tASSAY.1\t.\t.\n"
            "XGENE\tX\t2000000\t.\t.\t.\tXGENE\n"
            "ALT\tHSCHR22_1_CTG7\t2000000\t.\t.\t.\tALT\n"
            "UNKNOWN\t.\t.\t.\t.\t.\tUNKNOWN\n"
        )
        self.gwas = self.work / "input.tsv"
        self.gwas.write_text(
            "SNP\tCHR\tPOS\tEA\tNEA\tBETA\tSE\tP\tN\n"
            "rsleft\t1\t1500000\tA\tG\t0.1\t0.01\t1e-9\t1000\n"
            "rscenter\t1\t2000000\tA\tG\t0.1\t0.01\t1e-12\t1000\n"
            "rsright\t1\t2500000\tA\tG\t0.1\t0.01\t1e-10\t1000\n"
            "rsoutside\t1\t2500001\tA\tG\t0.1\t0.01\t1e-11\t1000\n"
            "rsbg\t1\t4000000\tA\tG\t0.1\t0.01\t0.05\t1000\n"
            "rstrans\t2\t7000000\tA\tG\t0.1\t0.01\t1e-8\t1000\n"
            "rsx\t23\t2000000\tA\tG\t0.1\t0.01\t1e-9\t1000\n"
        )
        self.genes = self.work / "input.genes.out"
        self.genes.write_text(
            "GENE CHR START STOP ZSTAT P\n"
            "NEAR 1 1950000 2050000 6 1e-8\n"
            "FAR 2 6990000 7010000 6 1e-8\n"
            "BACKGROUND 1 4000000 4001000 1 0.05\n"
        )
        self.plot = self.work / "capture.R"
        core = str(ROOT.parent / "0f/manhattan_plot.R")
        self.plot.write_text('source(' + repr(core) + ')\n' + '''
original_vectors <- .mh_plot_vectors
.mh_plot_vectors <- function(chr, pos, pvalue, annotate = NULL, mirror = NULL, ...) {
    dest <- Sys.getenv("CAPTURE_DIR")
    write.table(data.frame(CHR=chr, POS=pos, GROUP=as.character(annotate[[1]])),
        file.path(dest,"groups.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    settings <- annotate[-1]
    colors <- data.frame(GROUP=names(settings), COL=vapply(settings,function(x) x$col,character(1)))
    write.table(colors,file.path(dest,"colors.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    if (!is.null(mirror)) write.table(mirror,file.path(dest,"mirror.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    original_vectors(chr, pos, pvalue, annotate=annotate, mirror=mirror, ...)
}
''')

    def render(self, trait="ASSAY.1", cohort="scot", grch="38", success=True):
        args = ["Rscript", ROOT / "f/format.R", "mplot", self.gwas, self.work / "out.png",
                trait, "self", "magma", self.genes, self.plot, self.names, "", grch,
                self.work / "flags.tsv", "", "Protein", "[trait]", "Gene_chr,Gene_Start,Gene_end",
                "Gene,beta,p-value", "TRUE", self.work / "out.sig.txt", "", "8", "auto", "72"]
        env = dict(os.environ, GWAS_PROT_COHORT=cohort, CAPTURE_DIR=str(self.work),
                   GWAS_FULL="", GWAS_THIN="FALSE")
        result = subprocess.run(list(map(str, args)), env=env, text=True,
                                capture_output=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue((self.work / "out.png").read_bytes().startswith(b"\x89PNG"))
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def rows(self, file):
        with (self.work / file).open() as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def test_alias_window_red_points_and_magma_report(self):
        self.render()
        groups = self.rows("groups.tsv")
        self.assertEqual([r["GROUP"] for r in groups[:3]], ["ASSAY.1"] * 3)
        self.assertNotEqual(groups[3]["GROUP"], "ASSAY.1")
        self.assertEqual(next(r["COL"] for r in self.rows("colors.tsv") if r["GROUP"] == "ASSAY.1"), "#D62728")
        mirror = {r["GENE"]: r["CLASS"] for r in self.rows("mirror.tsv")}
        self.assertEqual(mirror, {"NEAR": "cis", "FAR": "trans", "BACKGROUND": "background"})
        report = (self.work / "out.sig.txt").read_text()
        self.assertTrue(any(line.startswith("rscenter\t") and line.endswith("\tcis") for line in report.splitlines()))
        self.assertTrue(any(line.startswith("NEAR\t") and line.endswith("\tcis") for line in report.splitlines()))

    def test_literal_na_alias_and_cohort_specific_match(self):
        for trait, cohort, chr_ in [("NA", "scot", "1"), ("ASSAY.1", "ckb", "2"),
                                   ("UKB", "prot", "1"), ("DECODE", "decode", "1"),
                                   ("CANON.Inflammation_II", "prot", "1"),
                                   ("XGENE", "scot", "23")]:
            with self.subTest(trait=trait, cohort=cohort):
                self.render(trait=trait, cohort=cohort)
                self.assertEqual({r["CHR"] for r in self.rows("groups.tsv") if r["GROUP"] == trait}, {chr_})

    def test_unmapped_and_alternative_contig_have_no_false_cis(self):
        for trait in ("ALT", "UNKNOWN"):
            with self.subTest(trait=trait):
                self.render(trait=trait)
                self.assertFalse(any(r["GROUP"] == trait for r in self.rows("groups.tsv")))
                self.assertFalse(any(r["CLASS"] == "cis" for r in self.rows("mirror.tsv")))

    def test_wrong_build_is_rejected(self):
        result = self.render(grch="37", success=False)
        self.assertIn("are GRCh38; GWAS is GRCh37", result.stderr)
        self.assertFalse((self.work / "out.png").exists())

    @unittest.skipUnless(shutil.which("bgzip") and shutil.which("tabix"), "htslib is required")
    def test_shell_worker_resume_tracks_table_changes_without_pheno_bed(self):
        project = self.work / "scot"
        gwas_dir = project / "common/ASSAY.1/gwas"
        gwas_dir.mkdir(parents=True)
        final = gwas_dir / "ASSAY.1.gz"
        with final.open("wb") as stream:
            subprocess.run(["bgzip", "-c", str(self.gwas)], stdout=stream, check=True)
        subprocess.run(["tabix", "-s", "2", "-b", "3", "-e", "3", "-S", "1", str(final)], check=True)
        # Tiny independent HM3 and annotation fixtures avoid workstation resources.
        hm3 = self.work / "hm3.txt"
        hm3.write_text("rscenter\n")
        hm3pos = self.work / "hm3.pos.tsv"
        hm3pos.write_text("SNP\tCHR\tPOS\nrscenter\t1\t2000000\n")
        bed = self.work / "genes.bed"
        bed.write_text("1\t1950000\t2050000\tCANON\n2\t6990000\t7010000\tOTHER\n")
        work = self.work / "coordinator"
        env = {k: v for k, v in os.environ.items() if not k.startswith(("GWAS_", "LSB_", "SCRIPT_"))}
        env.update(GWAS_WORK_ROOT=str(work), SCRIPT_CONSOLE_PASSTHROUGH="1",
                   SCRIPT_LOG_DIR=str(self.work / "console"))
        args = ["bash", ROOT / "format.sh", "mplot", "--submit", "local", "--label", "scot",
                "--dir-out", project, "--gwas", "ASSAY.1", "--grch", "38", "--prot-names", self.names,
                "--cis-bed", self.work / "does-not-exist.bed", "--mh-plot-bed", bed,
                "--hm3-file", hm3, "--hm3-pos", hm3pos, "--thin", "FALSE", "--local-io", "FALSE", "--write-sig", "TRUE",
                "--plot-width", "8", "--plot-res", "72", "--run-cmd", "TRUE", "--jobs", "1"]

        def run():
            r = subprocess.run(list(map(str, args)), env=env, capture_output=True, text=True, timeout=30)
            log = work / "logs/common/ASSAY.1/ASSAY.1.mplot.log"
            detail = log.read_text() if r.returncode and log.exists() else ""
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr + detail)
            return r

        run()
        png = project / "mplot/ASSAY.1.png"
        old = png.stat().st_mtime_ns
        state = work / "common/mplot/ASSAY.1.state"
        before = state.read_text()
        self.assertIn("prot_names_sha256\t", before)
        resumed = run()
        self.assertIn("SKIP completed Manhattan plot", resumed.stdout + resumed.stderr)
        self.assertEqual(png.stat().st_mtime_ns, old)
        self.names.write_text(self.names.read_text().replace("CANON\t1\t2000000", "CANON\t2\t7000000"))
        run()
        self.assertGreater(png.stat().st_mtime_ns, old)
        self.assertNotEqual(state.read_text(), before)
        report = (project / "mplot/ASSAY.1.sig.txt").read_text()
        self.assertTrue(any(line.startswith("rstrans\t") and line.endswith("\tcis") for line in report.splitlines()))


if __name__ == "__main__":
    unittest.main()
