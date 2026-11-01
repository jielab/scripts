"""Current result schema, actual-reference summaries, and denominator regression tests."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import csv
import json
import sqlite3
import sys

import pandas as pd
import pytest

GU = Path(__file__).resolve().parents[2]


def load(name, filename):
    spec = spec_from_file_location(name, GU / "f" / filename)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


results = load("gu_results_2026_tests", "results_normalize.py")
report = load("gu_report_2026_tests", "introgression_review.py")


def insert(con, table, row):
    columns = ",".join(row)
    con.execute(f"INSERT INTO {table} ({columns}) VALUES ({','.join('?' for _ in row)})", tuple(row.values()))


def complete_tables(con):
    results.build_catalog(con, "GRCh37", .5)
    results.build_burden(con)
    results.build_locus_support(con, 20)
    results.build_locus_evidence(con)
    results.build_database_summary(con)


def segment(sample="0001", method="ibdmix", source="Altai", start=100, end=250, chrom="1"):
    frame = pd.DataFrame([dict(sample_id=sample, source=source, chr=chrom, start=start, end=end,
                              haplotype="1" if method == "trace" else None,
                              score=8 if method == "ibdmix" else None,
                              posterior=.99 if method == "trace" else None)])
    return results.finish(frame, method, "GRCh37", Path("current.tsv"), dataset="1kg")


def test_current_schema_and_three_method_support():
    con = sqlite3.connect(":memory:")
    results.init_db(con)
    insert(con, "loci", dict(dataset_id="1kg", method="phyml", genome_build="GRCh37", locus_id="core",
           chr="1", start=100, end=300, input_start=100, input_end=300, analysis_start=50, analysis_end=350,
           source_class="Neanderthal", n_compared=20, n_match=20, prop_match=1, tree_status="complete",
           candidate_lineage="Neanderthal", tree_has_ancestral_outgroup=1, candidate_clade_pass=1,
           candidate_clade_bootstrap=95, candidate_clade_specificity=1, n_panel=2))
    for i, sample in enumerate(("0001", "0002")):
        insert(con, "sample_populations", dict(dataset_id="1kg", sample_id=sample, population="CEU", super_population="EUR"))
        insert(con, "phyml_haplotype_carriers", dict(dataset_id="1kg", genome_build="GRCh37", locus_id="core",
               sample_id=sample, haplotype="1", hap_id=f"H{i}", candidate_lineage="Neanderthal",
               tree_candidate_pass=int(i == 0), prop_match=1 if i == 0 else .3))
    for method, source in (("ibdmix", "Altai"), ("trace", "TRACE_ghost_or_unknown")):
        results.insert_df(con, "segments_raw", segment(method=method, source=source))
        insert(con, "method_runs", dict(dataset_id="1kg", genome_build="GRCh37", chr="1", method=method,
               status="complete", evidence_eligible=1, availability_note="current validated fixture"))
    complete_tables(con)
    assert con.execute("PRAGMA user_version").fetchone()[0] == 2026
    assert con.execute("SELECT sample_id,phyml,ibdmix,trace,n_methods FROM locus_sample_support ORDER BY sample_id").fetchall() == [
        ("0001", 1, 1, 1, 3), ("0002", 0, 0, 0, 0)]
    assert {r[0] for r in con.execute("SELECT DISTINCT method FROM locus_method_support")} == {"phyml", "ibdmix", "trace"}
    assert con.execute("SELECT independent_methods_available,independent_methods_supporting,ibdmix_max_lod,trace_max_posterior,n_panel FROM locus_evidence").fetchone() == (2, 2, 8.0, .99, 2)
    expected_tables = {"segments_raw", "loci", "coord_catalog", "segment_catalog", "sample_populations", "phyml_haplotype_carriers",
                       "method_runs", "database_summary", "segments", "carriers", "sample_burden",
                       "locus_sample_support", "locus_method_support", "locus_trajectory", "locus_evidence"}
    assert {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")} == expected_tables
    assert con.execute("PRAGMA quick_check").fetchone() == ("ok",)


def test_completed_zero_call_run_remains_visible():
    con = sqlite3.connect(":memory:")
    results.init_db(con)
    insert(con, "method_runs", dict(dataset_id="1kg", genome_build="GRCh37", chr="1", method="ibdmix", status="complete", evidence_eligible=1))
    complete_tables(con)
    summary = dict(con.execute("SELECT metric,value FROM database_summary"))
    assert summary["completed_method_units"] == 1
    assert summary["methods"] == 1
    assert con.execute("SELECT COUNT(*) FROM segments").fetchone() == (0,)


def test_reference_order_does_not_add_unrequested_references():
    assert results.summary_references(("Vindija", "Denisova25")) == ("Denisova25", "Vindija")
    assert results.summary_references(("Denisova.2013", "Altai.2013")) == ("Altai.2013", "Denisova.2013")
    assert results.summary_references(results.DEFAULT_ARCHAIC_REFERENCES) == ("Altai", "Chagyr", "Denisova", "Denisova25", "Vindija")
    with pytest.raises(ValueError, match="Unsupported archaic"):
        results.summary_references(("unrecognized",))


def test_exact_union_keeps_zero_and_untested_distinct():
    con = sqlite3.connect(":memory:")
    results.init_db(con)
    for start, end in ((100, 200), (150, 250)):
        results.insert_df(con, "segments_raw", segment(start=start, end=end))
    results.build_catalog(con, "GRCh37", .5)
    burden = results.reference_summary(con, "1kg", "GRCh37", ["0001", "0002", "0003"], {"1": 1000},
                                       {"1": {"0001", "0002"}}, "Altai")
    by_id = {r["sample_id"]: r for r in burden}
    assert by_id["0001"]["archaic_bp"] == 150
    assert by_id["0001"]["coordinate_coverage_pct"] == 15
    assert by_id["0001"]["coverage_pct"] == 7.5
    assert by_id["0002"]["archaic_bp"] == 0
    assert by_id["0002"]["tested_bp"] == 2000
    assert by_id["0003"]["archaic_bp"] == ""
    assert by_id["0003"]["tested_bp"] == 0


def test_x_nonpar_uses_one_physical_copy_and_clips_par():
    con = sqlite3.connect(":memory:")
    results.init_db(con)
    left = results.x_nonpar_intervals("GRCh37")[0][0]
    results.insert_df(con, "segments_raw", segment(start=left - 100, end=left + 100, chrom="X"))
    results.build_catalog(con, "GRCh37", .5)
    lengths = results.CHROM_LENGTHS["37"]
    burden = results.reference_summary(con, "1kg", "GRCh37", ["0001", "0002"], lengths, {"X": {"0001"}}, "Altai", x_male=True)
    assert burden[0]["archaic_bp"] == 100
    assert burden[0]["tested_bp"] == burden[0]["haploid_bp"]
    assert burden[0]["coverage_pct"] == burden[0]["coordinate_coverage_pct"]
    assert burden[1]["archaic_bp"] == ""


@pytest.mark.parametrize("references", [("Denisova25",), ("Altai.2013", "Denisova.2013")])
def test_density_exports_only_actual_tested_references(tmp_path, references):
    root = tmp_path / "analysis"
    run = root / "ibdmix" / "chr1"
    (run / "samples/C1").mkdir(parents=True)
    (run / "samples/C1/ALL.txt").write_text("0001\n0002\n")
    (run / ".complete").touch()
    (run / "run.meta.tsv").write_text(
        "pipeline_version\t2026-10-31.shared-hts-v1\ncoordinate_system\t0-based-half-open\n"
        "genome_build\tGRCh37\nrefs\t" + " ".join(references) + "\n"
        "background_filter\t0\nexport_denisovan\t1\nmodern_vcf\tcurrent/chr1.vcf.gz:10:10\n")
    db = root / "final/gu.sqlite"
    db.parent.mkdir()
    con = sqlite3.connect(db)
    results.init_db(con)
    runs = results.discover_method_runs(root, "GRCh37")
    assert len(runs) == 1
    con.executemany("INSERT INTO method_runs VALUES(?,?,?,?,?,?,?,?)", runs)
    for sample in ("0001", "0002"):
        insert(con, "sample_populations", dict(dataset_id="1kg", sample_id=sample, population="CEU", super_population="EUR"))
    for source in references:
        results.insert_df(con, "segments_raw", segment(source=source))
    complete_tables(con)
    con.close()
    output = root / "final/density"
    results.prepare(db, output, bin_bp=1_000_000_000)
    manifest_path, = output.glob("*/**/manifest.json")
    manifest = json.loads(manifest_path.read_text())
    assert manifest["archaic_summary_references"] == list(results.summary_references(references))
    with (manifest_path.parent / "archaic_summary.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["reference"] for r in rows} == set(references)
    assert len(rows) == 2 * len(references)
    assert all(r["archaic_bp"] == "0" for r in rows if r["sample_id"] == "0002")
    assert (manifest_path.parent / "archaic_summary_x.tsv").read_text().count("\n") == 1


def test_import_requires_current_coordinate_and_reference_records():
    good = dict(pipeline_version="2026-10-31.shared-hts-v1", coordinate_system="0-based-half-open", refs="Altai")
    results.require_ibdmix_metadata(good, Path("current"))
    for field in ("pipeline_version", "coordinate_system", "refs"):
        with pytest.raises(ValueError):
            results.require_ibdmix_metadata({k: v for k, v in good.items() if k != field}, Path("incomplete"))


def test_current_report_interval_and_tsv_helpers(tmp_path):
    covered = report.coverage((100, 300), [(50, 150), (140, 200), (250, 400)])
    assert covered["overlap_bp"] == 150
    assert covered["union_fraction"] == .75
    assert covered["n_covered_blocks"] == 2
    with pytest.raises(ValueError):
        report.union([(-1, 1)])
    path = tmp_path / "table.tsv"
    report.write_tsv(path, [{"sample_id": "0001", "value": None}])
    assert report.rows(path) == [{"sample_id": "0001", "value": ""}]
    assert report.chrom("chr23") == "X"
    assert report.truth("PASS")
