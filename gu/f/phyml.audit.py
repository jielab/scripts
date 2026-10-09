#!/usr/bin/env python3
"""Read-only audit of archived GWAS PhyML analyses; no tree fitting or extraction.

Invoked through phyml.py audit. Original risk-only calls stay unchanged.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter
from pathlib import Path
import subprocess
import tarfile


TABLES = ("final/gwas_loci.tsv", "final/evidence_trees.tsv", "loci/haplotypes.tsv",
	"loci/ld.tsv", "loci/site_qc.tsv", "loci/archaic.tsv")
TREE_PREFIX = "loci/haplotypes.phy"
PRODUCTS = (TREE_PREFIX + "_phyml_tree.txt", TREE_PREFIX + "_phyml_boot_trees.txt",
	TREE_PREFIX + ".phyml.log", "final/gwas_parameters.json", "GU-MANIFEST.json")


def archive_tables(path):
	"""Read each gzip stream once; validate the exact members used by the audit."""
	data = {}
	with tarfile.open(path, "r|gz") as archive:
		for member in archive:
			if member.isfile() and member.name in TABLES + PRODUCTS:
				if member.name in data:
					raise ValueError(f"duplicate archive member: {path}: {member.name}")
				data[member.name] = archive.extractfile(member).read()
	manifest = json.loads(data.get("GU-MANIFEST.json", b"{}"))
	checks = {r["path"]: r for r in manifest.get("files", [])}
	for name, blob in data.items():
		if checks and name != "GU-MANIFEST.json" and name not in checks:
			raise ValueError(f"unlisted archive member: {path}: {name}")
		if name in checks and hashlib.sha256(blob).hexdigest() != checks[name]["sha256"]:
			raise ValueError(f"archive checksum mismatch: {path}: {name}")
	tables = {name: list(csv.DictReader(io.StringIO(data.get(name, b"").decode()), delimiter="\t")) for name in TABLES}
	return tables, {k: v.decode() for k, v in data.items() if k not in TABLES}, bool(checks)


def affinity_edges(p, newick, haps, lineage):
	"""Exploratory pure-lineage edges, explicitly NOT genome-wide positive calls.

	Unlike the predefined whole-allele test, these can contain subsets of modern
	haplotypes and reference genomes. Scanning edges needs independent evidence.
	"""
	root = p.parse_newick(newick)
	edges = []
	def visit(node):
		tips = {node.label} if not node.children else set().union(*(visit(c) for c in node.children))
		edges.append((node, tips))
		return tips
	all_tips = visit(root)
	modern = {h["hap_id"] for h in haps}
	refs = set(p.LINEAGE_REFS[lineage])
	roles = {h["hap_id"]: h["role"] for h in haps}
	found = {}
	for node, tips in edges:
		if node is root or node.support is None or node.support < 70: continue
		for side in (tips, all_tips - tips):
			m = side & modern
			if not m or not side & refs or not modern - side or side - (modern | refs): continue
			key = ",".join(sorted(side))
			if key in found and found[key]["bootstrap"] >= node.support: continue
			found[key] = dict(bootstrap=node.support, tips=key, modern_tips=",".join(sorted(m)),
				archaic_tips=",".join(sorted(side & refs)), n_risk=sum(roles[h] == "risk" for h in m),
				n_nonrisk=sum(roles[h] == "nonrisk" for h in m),
				evidence="exploratory_affinity_only;not_introgression_call;edge_scan_not_multiplicity_calibrated")
	return list(found.values())


def audit(p, dataset, leads):
	summary, affinities, reference_qc, provenance = [], [], [], []
	seen = set()
	for index, archive in enumerate(sorted(dataset.glob("chr*/phyml.raw.tar.gz")), 1):
		t, raw, checksummed = archive_tables(archive)
		ss = t["final/gwas_loci.tsv"]
		if not ss: continue
		lid = ss[0]["locus_id"]
		if lid not in leads: continue
		if len(ss) != 2 or {r.get("lineage") for r in ss} != set(p.LINEAGE_REFS):
			raise ValueError(f"missing or duplicate lineage summaries: {archive}")
		if lid in seen: raise ValueError(f"duplicate analysis for {lid}")
		seen.add(lid)
		haps, ld, qc = (t[n] for n in ("loci/haplotypes.tsv", "loci/ld.tsv", "loci/site_qc.tsv"))
		arch = {r["archaic"]: r["seq"] for r in t["loci/archaic.tsv"]}
		for ref in p.REFS:
			for reason, count in Counter(r.get(ref + "_callability", "unknown") for r in qc).items():
				reference_qc.append(dict(locus_id=lid, reference=ref, reason=reason, n_sites=count))
		for source in ss:
			s = dict(source)
			lineage = s.get("lineage", "Neanderthal")
			tree = next((r for r in t["final/evidence_trees.tsv"] if r.get("expected_lineage") == lineage), {})
			complete = tree.get("tree_status") == "complete"
			nw = tree.get("tree_newick", "")
			if s["status"] in ("tree_supported", "tree_not_supported") and not complete:
				raise ValueError(f"saved decision without a complete tree: {archive} {lineage}")
			if complete:
				log = raw.get(TREE_PREFIX + ".phyml.log", "")
				boot = raw.get(TREE_PREFIX + "_phyml_boot_trees.txt", "").splitlines()
				if (raw.get(TREE_PREFIX + "_phyml_tree.txt", "").strip() != nw.strip()
					or len([x for x in boot if x.strip().endswith(";")]) != 100
					or "Time used" not in log or "Cannot work out eigen vectors" in log):
					raise ValueError(f"incomplete or inconsistent saved tree: {archive}")
				strict = p.risk_clade(nw, {h["hap_id"] for h in haps if h["role"] == "risk"},
					{h["hap_id"] for h in haps}, p.LINEAGE_REFS[lineage])
				strict_pass = bool(strict and strict["bootstrap"] is not None and strict["bootstrap"] >= 70)
				if strict_pass != (str(s["tree_pass"]) == "1"):
					raise ValueError(f"saved risk decision disagrees with tree: {lid} {lineage}")
				affinities.extend(dict(locus_id=lid, lineage=lineage, **edge) for edge in affinity_edges(p, nw, haps, lineage))
			s.update(p.allele_topology_audit(nw, haps, lineage, complete))
			s.update(p.locus_stage_fields(s))
			s.update(archive=str(archive), archive_members_checksummed=int(checksummed),
				ld_markers_r2_gt_08=sum(float(r.get("ld_r2") or 0) > 0.8 for r in ld),
				ld_markers_r2_gt_098=sum(r.get("core_marker") == "1" for r in ld),
				n_core_qc_sites=len(qc) if qc else None,
				all5_callable=sum(r.get("retained_alignment") == "1" for r in qc) if qc else None,
				nea3_callable=sum(all(r.get(ref + "_base") in p.BASES for ref in p.LINEAGE_REFS["Neanderthal"]) for r in qc) if qc else None,
				lead_nea3_callable=next((int(all(r.get(ref + "_base") in p.BASES for ref in p.LINEAGE_REFS["Neanderthal"])) for r in qc if r.get("lead") == "1"), None),
				unique_archaic_sequences=len(set(arch.values())) if arch else None,
				sequence_resolution="all_archaics_identical" if arch and len(set(arch.values())) == 1 else "see_site_qc")
			summary.append(s)
		provenance.append(dict(locus_id=lid, archive=str(archive), archive_sha256=p.phyml_run_digest(archive)))
		if index % 100 == 0: print(f"[GU AUDIT] {index} archives checked", flush=True)
	for lid in sorted(set(leads) - seen):
		for lineage in p.LINEAGE_REFS:
			summary.append(dict(leads[lid], lineage=lineage, status="no_saved_result", tree_run_status="not_run", tree_pass=None))
	return summary, affinities, reference_qc, provenance


def refresh_report(p, path, summary, output):
	"""Update diagnostic columns only, retaining original calls and validation."""
	data = p.read(path)
	key = lambda r: (r.get("dataset_id"), r["locus_id"], r["lineage"])
	audit_rows = {key(r): r for r in summary if r.get("dataset_id")}
	fields = list(p.allele_topology_audit("", [], "Neanderthal")) + list(p.locus_stage_fields({}))
	updated = 0
	for row in data:
		s = audit_rows.get(key(row))
		if s is None: continue
		if row.get("status", row.get("call")) != s["status"] or str(row["tree_pass"]) != str(s["tree_pass"]):
			raise ValueError("Report differs from archived results; run final first: " + str(key(row)))
		row.update({k: s.get(k) for k in fields})
		updated += 1
	if updated != len(audit_rows): raise ValueError("Report is missing audited locus/lineage rows")
	before = output / "report_before.tsv.gz"
	if not before.exists():
		with gzip.open(before, "wb") as handle: handle.write(path.read_bytes())
	p.write(path, data)
	return dict(path=str(path), updated_rows=updated, sha256=p.phyml_run_digest(path))


def main(p):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--dataset-dir", type=Path, required=True)
	parser.add_argument("--lead-table", type=Path, required=True)
	parser.add_argument("--output", type=Path, required=True)
	parser.add_argument("--refresh-report", type=Path, help="Explicitly refresh diagnostic fields in an existing report; preserve risk calls")
	a = parser.parse_args()
	leads = {r["locus_id"]: r for r in p.read(a.lead_table)}
	summary, affinities, qc, provenance = audit(p, a.dataset_dir, leads)
	a.output.mkdir(parents=True, exist_ok=True)
	stages = [dict(lineage=lineage, status=status, n_loci=n) for (lineage, status), n in sorted(Counter((r["lineage"], r["status"]) for r in summary).items())]
	supported = [r for r in summary if r.get("either_allele_tree_pass") == 1]
	input_audit = a.lead_table.parent / "input.audit.tsv"
	tables = dict(stages=stages, loci=summary, allele_tree_support=supported, exploratory_edges=affinities, reference_qc=qc,
		input_exclusions=p.read(input_audit) if input_audit.is_file() else [])
	for name, data in tables.items(): p.write(a.output / (name + ".tsv"), data)
	manifest = dict(lead_table=str(a.lead_table), lead_table_sha256=p.phyml_run_digest(a.lead_table), n_leads=len(leads),
		method_source="https://www.nature.com/articles/s41586-020-2818-3", code_sha256=p.phyml_run_digest(Path(__file__)),
		phyml_code_sha256=p.phyml_run_digest(Path(p.__file__)),
		note="No trees refitted. Whole-allele support and exploratory affinity are not confirmed introgression. LD 0.8 is a sensitivity count only; original threshold remains >0.98.",
		stages=stages, sources=provenance)
	if a.refresh_report: manifest["report_refresh"] = refresh_report(p, a.refresh_report, summary, a.output)
	(a.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
	# Use the repository's lossless workbook writer, with exact TSV attachments.
	numeric = {"bootstrap", "tree_pass", "tree_bootstrap", "nonrisk_tree_pass", "nonrisk_tree_bootstrap",
		"either_allele_tree_pass", "all5_callable", "nea3_callable", "lead_nea3_callable", "unique_archaic_sequences",
		"core_start", "core_end", "lead_pos", "core_kb"}
	spec = dict(tables=[], metadata=dict(method=manifest["note"]))
	for name in tables:
		path = a.output / (name + ".tsv")
		with path.open() as handle: fields = next(csv.reader(handle, delimiter="\t"))
		spec["tables"].append(dict(name=name, path=str(path.resolve()),
			types=["numeric" if k in numeric or k.startswith(("n_", "ld_markers_")) else "character" for k in fields]))
	job = a.output / "workbook.json"
	job.write_text(json.dumps(spec))
	writer = Path(__file__).resolve().parents[2] / "0f/results.py"
	subprocess.run([p.sys.executable, str(writer), "stream-workbook", str(job), str(a.output / "phyml_audit.xlsx")], check=True)
	print(json.dumps(dict(n_leads=len(leads), supported_allele_tests=len(supported), stages=stages), ensure_ascii=False), flush=True)
