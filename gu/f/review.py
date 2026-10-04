#!/usr/bin/env python3
"""GU review utilities. See --help for commands."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

if "gu_0_common" not in sys.modules:
	_spec = spec_from_file_location("gu_0_common", Path(__file__).with_name("0.common.py"))
	_common = module_from_spec(_spec)
	sys.modules[_spec.name] = _common
	try:
		_spec.loader.exec_module(_common)
	except BaseException:
		sys.modules.pop(_spec.name, None)
		raise
from gu_0_common import load_module


# 🚩 dual_lead
"""Region-first post-processing: original lead vs haplotype-derived tag.

Uses only saved, phased GU sequences and copy IDs. No GWAS input, risk allele,
lead-SNP LD window, environment changes, or rerunning of scientific callers.
The binary SNP--haplotype r^2 is calculated on chromosome copies, not people.
"""

import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from collections.abc import Mapping
from pathlib import Path

load_module("0.common.py")
from gu_0_common import resolve_tsv_path

BASES = set("ACGT")
dual_lead_VERSION = "2026-09-06.1"
dual_lead_MISSING = {"", ".", "NA", "nan", "None", "null"}


def num(v, default=None):
	try:
		x = float(v)
		return x if math.isfinite(x) else default
	except (ValueError, TypeError):
		return default


def dual_lead_integer(v, default=None):
	x = num(v)
	if x is None:
		return default
	if x != int(x):
		raise ValueError(f"Nonintegral value: {v!r}")
	return int(x)


def dual_lead_truth(v):
	return str(v).strip().upper() in {"1", "TRUE", "T", "YES", "PASS"}


def dual_lead_chrom(v):
	x = re.sub(r"^chr", "", str(v), flags=re.I)
	return "X" if x == "23" else x


def lineage(v):
	x = str(v).lower()
	if any(t in x for t in ("neander", "altai", "vindija", "chagyr")):
		return "Neanderthal"
	if "denis" in x:
		return "Denisovan"
	return str(v)


def rows(path):
	if path is None:
		return []
	path = resolve_tsv_path(Path(path))
	if not path.is_file():
		return []
	op = gzip.open if path.suffix == ".gz" else open
	with op(path, "rt", encoding="utf-8-sig", newline="") as f:
		out = list(csv.DictReader(f, delimiter="\t"))
	if any(None in r for r in out):
		raise ValueError(f"Malformed TSV fields: {path}")
	return out


class SiteAlleles(Mapping):
	"""Lazy shared sequence-backed map: avoids copies x SNPs dictionaries."""

	def __init__(self, copy_sequences, index, allowed):
		self.sequences = copy_sequences
		self.index = index
		self.allowed = set(allowed)

	def get(self, key, default=None):
		seq = self.sequences.get(key)
		if seq is None:
			return default
		base = seq[self.index]
		return base if base in self.allowed else default

	def __getitem__(self, key):
		value = self.get(key)
		if value is None:
			raise KeyError(key)
		return value

	def __iter__(self):
		return (cp for cp in self.sequences if self.get(cp) is not None)

	def __len__(self):
		return sum(self.get(cp) is not None for cp in self.sequences)

	def __bool__(self):
		return any(self.get(cp) is not None for cp in self.sequences)


def coord_key(r):
	st, en = dual_lead_integer(r.get("core_start")), dual_lead_integer(r.get("core_end"))
	if st is None or en is None or st < 0 or en <= st:
		raise ValueError("A valid core_start/core_end is required for a region-first locus key")
	return f"{r['dataset_id']}|{r['genome_build']}|{dual_lead_chrom(r['chr'])}:{st}-{en}"


def aliases(v):
	return {x.strip() for x in re.split("[;,]", str(v or "")) if x.strip() not in dual_lead_MISSING}


def copy_ids(h):
	tokens = [x.strip() for x in str(h.get("copies", "")).split(";") if x.strip()]
	if len(tokens) != len(set(tokens)) or any(":" not in t for t in tokens):
		raise ValueError(f"Duplicate or invalid sample:haplotype copy ID: {h.get('hap_id')}")
	if not tokens:
		raise ValueError(f"Copy IDs missing for haplotype {h.get('hap_id')}")
	return set(tokens)


def inside(st, en, pos1):
	return pos1 is not None and st <= pos1 - 1 < en


def binary_metrics(tp, fp, fn, tn, total_target=None, total=None):
	"""Target=haplotype membership; predicted=specified allele present.

	PPV/sensitivity here describe a tag for this computational target, NOT
	accuracy of introgression detection and NOT tree purity.
	"""
	a, b, c, d = tp, fp, fn, tn
	if min(a, b, c, d) < 0:
		raise ValueError("Negative contingency count")
	n = a + b + c + d
	denom = (a + b) * (c + d) * (a + c) * (b + d)
	signed = (a * d - b * c) / math.sqrt(denom) if denom else None
	return dict(
		n_callable_copies=n,
		n_target_callable=a + c,
		n_background_callable=b + d,
		tp=a,
		fp=b,
		fn=c,
		tn=d,
		r=signed,
		r2=signed * signed if signed is not None else None,
		sensitivity=a / (a + c) if a + c else None,
		specificity=d / (b + d) if b + d else None,
		ppv=a / (a + b) if a + b else None,
		f1=2 * a / (2 * a + b + c) if 2 * a + b + c else None,
		allele_frequency=(a + b) / n if n else None,
		haplotype_frequency=(a + c) / n if n else None,
		call_rate=n / total if total else None,
		target_call_rate=(a + c) / total_target if total_target else None,
	)


def count_variant(alleles, target, cohort, tag_allele):
	a = b = c = d = 0
	for cp in cohort:
		value = alleles.get(cp)
		if value is None:
			continue
		ish = cp in target
		ist = value == tag_allele
		if ish and ist:
			a += 1
		elif ist:
			b += 1
		elif ish:
			c += 1
		else:
			d += 1
	return binary_metrics(a, b, c, d, len(target & cohort), len(cohort))


def qualify(m, min_call_rate=0.95, min_target=2):
	if m["n_target_callable"] < min_target:
		return "too_few_callable_target_copies"
	if m["n_background_callable"] < 2:
		return "too_few_callable_background_copies"
	if m["call_rate"] is None or m["call_rate"] < min_call_rate:
		return "low_call_rate"
	if m["target_call_rate"] is None or m["target_call_rate"] < min_call_rate:
		return "low_target_call_rate"
	if m["r"] is None:
		return "monomorphic_allele_or_target"
	if m["r"] <= 0:
		return "allele_not_enriched_in_target"
	return "eligible"


def rank_tag(r):
	# Neither input lead position, input lead identity nor GWAS appears here.
	return (
		-round(r["r2"], 12),
		-round(r["f1"], 12),
		-round(r["call_rate"], 12),
		-int(r["inside_candidate"]),
		-int(r["diagnostic_for_target_lineage"]),
		r["pos_1based"],
		r["ref"],
		r["alt"],
		r["tag_allele"],
		r["variant_key"],
	)


def rank_evidence(r):
	return (
		-r["evidence_tier"],
		-(r.get("diagnostic_match_prop") or 0),
		-r.get("n_contiguous_diagnostic_match", 0),
		-r["prop_match"],
		-r["n_compared"],
		-r["n_copies"],
		r["hap_id"],
		r["lineage"],
		r["archaic"],
	)


def make_groups(copies, population_rows):
	meta = {}
	for p in population_rows:
		sid = p.get("sample_id", p.get("sample"))
		v = (p.get("population", p.get("pop")), p.get("super_population", p.get("super_pop")))
		if sid in meta and meta[sid] != v:
			raise ValueError(f"Conflicting population metadata for {sid}")
		meta[sid] = v
	groups = {"ALL": set(copies)}
	for cp in copies:
		pop, sup = meta.get(cp.rsplit(":", 1)[0], (None, None))
		for prefix, v in [("POP", pop), ("SUPER", sup)]:
			if v and v not in dual_lead_MISSING:
				groups.setdefault(prefix + ":" + v, set()).add(cp)
	return groups


def resolve_sequence_files(auth):
	f = Path(auth["_evidence_file"]).parent
	locus = f.parent / "loci"
	if not (locus / "sites.tsv").is_file():
		locus = locus / auth["locus_id"]
	# The region-wide catalogue is preferred; legacy tables may be LD-selected.
	region = resolve_tsv_path(f / "region_haplotypes.unfiltered.tsv.gz")
	hp = region if region.is_file() else f / "haplotypes.tsv"
	return f, locus, hp, hp == region


def load_data(auth):
	final, locus, hp, region_inventory = resolve_sequence_files(auth)
	paths = [hp, locus / "sites.tsv", locus / "archaic.tsv"]
	missing = [p.name for p in paths if not p.is_file()]
	if missing:
		return None, {"status": "sequence_files_missing", "detail": ",".join(missing)}, []
	sites = rows(paths[1])
	haps = rows(hp)
	arch = rows(paths[2])
	haps = [h for h in haps if not h.get("locus_id") or h["locus_id"] == auth["locus_id"]]
	if not sites or not haps or not arch:
		return None, {"status": "empty_sequence_catalogue"}, paths
	needed = {"chr", "pos", "id", "ref", "alt"}
	if not needed <= sites[0].keys():
		raise ValueError("sites.tsv lacks required chr/pos/id/ref/alt columns")
	nsite = len(sites)
	if any(len(str(h.get("seq", ""))) != nsite for h in haps + arch):
		raise ValueError("Sequence/site length mismatch; never truncate or infer SNP positions")
	if len({h["hap_id"] for h in haps}) != len(haps):
		raise ValueError("Duplicate hap_id in the same locus catalogue")
	st, en = dual_lead_integer(auth["core_start"]), dual_lead_integer(auth["core_end"])
	if region_inventory:
		if any(
			dual_lead_integer(h.get("core_start")) != st or dual_lead_integer(h.get("core_end")) != en for h in haps
		):
			raise ValueError("Unfiltered catalogue does not match the full input locus coordinates")
		scope = "region_wide_unfiltered_catalogue"
	else:
		lr = [r for r in rows(final / "loci.tsv") if r.get("locus_id") == auth["locus_id"]]
		if not lr:
			return (
				None,
				{
					"status": "cannot_verify_region_wide_scope",
					"detail": "Legacy haplotypes need loci.tsv core/selected bounds",
				},
				paths,
			)
		sel = lr[0]
		s0 = dual_lead_integer(sel.get("selected_start"), dual_lead_integer(sel.get("region_start")))
		e0 = dual_lead_integer(sel.get("selected_end"), dual_lead_integer(sel.get("region_end")))
		if (s0, e0) != (st, en):
			return (
				None,
				{
					"status": "legacy_lead_selected_region_not_used",
					"detail": "Run a region-wide PhyML inventory; an LD-pruned subset cannot produce a full-locus best tag",
				},
				paths + [final / "loci.tsv"],
			)
		paths.append(final / "loci.tsv")
		scope = "legacy_catalogue_full_core_bounds_verified"
	all_copies = set()
	for h in haps:
		h["_copies"] = copy_ids(h)
		n = dual_lead_integer(h.get("hap_n", h.get("n")), len(h["_copies"]))
		if n != len(h["_copies"]):
			raise ValueError(f"Copy count mismatch: {h['hap_id']}")
		if all_copies & h["_copies"]:
			raise ValueError("A chromosome copy belongs to multiple catalogue haplotypes")
		all_copies |= h["_copies"]
		h["seq"] = h["seq"].upper()
	for a in arch:
		a["seq"] = a["seq"].upper()
		a["lineage"] = lineage(a.get("lineage") or a.get("archaic"))
	in_core = []
	seen_sites = set()
	for i, s in enumerate(sites):
		s["pos"] = dual_lead_integer(s["pos"])
		if s["pos"] is None or s["pos"] < 1:
			raise ValueError("Invalid site position")
		sid = (dual_lead_chrom(s["chr"]), s["pos"], s["ref"], s["alt"])
		if sid in seen_sites:
			raise ValueError("Duplicate chromosome/position/REF/ALT site")
		seen_sites.add(sid)
		if dual_lead_chrom(s["chr"]) == dual_lead_chrom(auth["chr"]) and inside(st, en, s["pos"]):
			in_core.append(i)
	if not in_core:
		return None, {"status": "no_stored_sites_in_core"}, paths
	evpath = final / "evidence_haplotypes.tsv"
	ev = rows(evpath)
	ep = final / "evidence_sites.tsv"
	diag = rows(ep)
	tp = final / "evidence_trees.tsv"
	trees = rows(tp)
	for p in [evpath, ep, tp, final / "anchors.tsv", final / "anchor_copies.tsv"]:
		if p.is_file():
			paths.append(p)
	# Membership labels are validated against actual copy sets, not raw H labels.
	byid = {h["hap_id"]: h for h in haps}
	for e in ev:
		if e.get("locus_id") != auth["locus_id"] or e.get("hap_id") not in byid:
			continue
		e["_copies"] = copy_ids(e)
		if e["_copies"] != byid[e["hap_id"]]["_copies"]:
			raise ValueError("Evidence and region catalogue H labels map to different copies")
	return (
		dict(
			final=final,
			locus=locus,
			sites=sites,
			haps=haps,
			arch=arch,
			in_core=in_core,
			copies=all_copies,
			copy_sequences={cp: h["seq"] for h in haps for cp in h["_copies"]},
			evidence=ev,
			diagnostic=diag,
			trees=trees,
			scope=scope,
		),
		{"status": "computed"},
		paths,
	)


def choose_haplotype(auth, data, min_sites=10, min_callable=0.8):
	ev = {
		(r["hap_id"], r.get("diagnostic_lineage")): r
		for r in data["evidence"]
		if r.get("locus_id") == auth["locus_id"] and r.get("hap_id")
	}
	tree_members = {}
	for t in data["trees"]:
		if t.get("locus_id") == auth["locus_id"] and dual_lead_truth(t.get("candidate_clade_pass")):
			tree_members[t.get("expected_lineage")] = set(t.get("candidate_tips_in_clade", "").split(","))
	if dual_lead_truth(auth.get("candidate_tree_pass")):
		tree_members.setdefault(
			auth.get("evidence_lineage"), set(auth.get("candidate_tree_candidates_in_clade", "").split(","))
		)
	comparisons = []
	for h in data["haps"]:
		for a in data["arch"]:
			pos = [i for i in data["in_core"] if a["seq"][i] in BASES]
			called = [i for i in pos if h["seq"][i] in BASES]
			nc = len(called)
			nm = sum(h["seq"][i] == a["seq"][i] for i in called)
			if nc < min_sites or not pos or nc / len(pos) < min_callable:
				continue
			e = ev.get((h["hap_id"], a["lineage"]), {})
			diagnostic = dual_lead_truth(e.get("diagnostic_candidate_pass"))
			tree = diagnostic and h["hap_id"] in tree_members.get(a["lineage"], set())
			comparisons.append(
				dict(
					hap_id=h["hap_id"],
					archaic=a["archaic"],
					lineage=a["lineage"],
					n_compared=nc,
					n_match=nm,
					prop_match=nm / nc,
					n_copies=len(h["_copies"]),
					archaic_callable_fraction=nc / len(pos),
					evidence_tier=3 if tree else 2 if diagnostic else 1,
					diagnostic_match_prop=num(e.get("diagnostic_match_prop")),
					n_contiguous_diagnostic_match=dual_lead_integer(e.get("n_contiguous_diagnostic_match"), 0),
					candidate_start=dual_lead_integer(e.get("candidate_start")),
					candidate_end=dual_lead_integer(e.get("candidate_end")),
				)
			)
	if not comparisons:
		return None, None, ev, tree_members
	best = sorted(comparisons, key=rank_evidence)[0]
	raw = sorted(
		comparisons, key=lambda r: (-r["prop_match"], -r["n_compared"], -r["n_copies"], r["hap_id"], r["archaic"])
	)[0]
	return best, raw, ev, tree_members


def build_variants(auth, data, best):
	h = next(x for x in data["haps"] if x["hap_id"] == best["hap_id"])
	a = next(x for x in data["arch"] if x["archaic"] == best["archaic"])
	passed_diag = {
		(dual_lead_integer(x.get("pos")), x.get("lineage"), x.get("lineage_consensus"))
		for x in data["diagnostic"]
		if dual_lead_truth(x.get("diagnostic_site_pass")) and x.get("locus_id") == auth["locus_id"]
	}
	pool = []
	for i in data["in_core"]:
		s = data["sites"][i]
		ref, alt = s["ref"].upper(), s["alt"].upper()
		# An A/C/G/T alignment cannot reconstruct indels or multiallelic alleles.
		if ref not in BASES or alt not in BASES or ref == alt:
			continue
		allele = h["seq"][i]
		if allele not in {ref, alt} or a["seq"][i] != allele:
			continue
		amap = SiteAlleles(data["copy_sequences"], i, {ref, alt})
		cs, ce = best["candidate_start"], best["candidate_end"]
		pool.append(
			dict(
				chr=dual_lead_chrom(s["chr"]),
				pos_1based=s["pos"],
				snp_id=s["id"] if s["id"] not in dual_lead_MISSING else None,
				ref=ref,
				alt=alt,
				tag_allele=allele,
				other_allele=alt if allele == ref else ref,
				variant_key=f"{auth['genome_build']}:{dual_lead_chrom(s['chr'])}:{s['pos']}:{ref}:{alt}",
				archaic_ref=best["archaic"],
				archaic_matching_allele=allele,
				inside_candidate=int(cs is not None and ce is not None and inside(cs, ce, s["pos"])),
				diagnostic_for_target_lineage=int((s["pos"], best["lineage"], allele) in passed_diag),
				_alleles=amap,
			)
		)
	return pool


def input_variant(auth, data):
	requested = auth.get("_input_lead_override", auth.get("named_anchor_id", auth.get("name", auth["locus_id"])))
	if not requested or requested in dual_lead_MISSING:
		return requested, None, "input_not_provided"
	hits = [s for s in data["sites"] if requested in aliases(s.get("id"))]
	# Exact CHR:POS:REF:ALT is also accepted; no nearest-SNP substitution.
	if not hits:
		m = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([^:]+):([^:]+)", requested, re.I)
		if m:
			hits = [
				s
				for s in data["sites"]
				if (dual_lead_chrom(s["chr"]), s["pos"], s["ref"], s["alt"])
				== (dual_lead_chrom(m[1]), int(m[2]), m[3], m[4])
			]
	if len(hits) > 1:
		return requested, None, "ambiguous_input_variant_id"
	if hits:
		s = hits[0]
		i = data["sites"].index(s)
		valid = {s["ref"], s["alt"]}
		if len(valid) == 2 and valid <= BASES:
			amap = SiteAlleles(data["copy_sequences"], i, valid)
			return (
				requested,
				dict(
					chr=dual_lead_chrom(s["chr"]),
					pos_1based=s["pos"],
					snp_id=requested,
					ref=s["ref"],
					alt=s["alt"],
					variant_key=f"{auth['genome_build']}:{dual_lead_chrom(s['chr'])}:{s['pos']}:{s['ref']}:{s['alt']}",
					_alleles=amap,
				),
				"observed_in_phased_alignment",
			)
	anchors = [
		r
		for r in rows(data["final"] / "anchors.tsv")
		if r.get("locus_id") == auth["locus_id"]
		and requested in aliases(r.get("anchor_id")) | aliases(r.get("requested_anchor_id"))
	]
	if len(anchors) > 1:
		return requested, None, "ambiguous_input_anchor_record"
	if anchors and dual_lead_truth(anchors[0].get("exact_anchor_found")):
		s = anchors[0]
		pos = dual_lead_integer(s.get("pos"))
		ref = s.get("ref", "")
		alt = s.get("alt", "")
		rec = dict(
			chr=dual_lead_chrom(s.get("chr", auth["chr"])),
			pos_1based=pos,
			snp_id=requested,
			ref=ref,
			alt=alt,
			variant_key=f"{auth['genome_build']}:{dual_lead_chrom(s.get('chr', auth['chr']))}:{pos}:{ref}:{alt}",
			_alleles={},
		)
		if "," in alt or not ref or not alt:
			return requested, rec, "multiallelic_input_not_scored"
		for c in rows(data["final"] / "anchor_copies.tsv"):
			if c.get("locus_id") != auth["locus_id"] or not (
				aliases(c.get("anchor_id")) & (aliases(s.get("anchor_id")) | {requested})
			):
				continue
			if dual_lead_integer(c.get("anchor_pos")) != pos or not dual_lead_truth(c.get("assigned_to_hap_id")):
				continue
			cp = str(c.get("sample", c.get("sample_id"))) + ":" + str(c.get("haplotype"))
			val = c.get("allele")
			if cp in data["copies"] and val in {ref, alt}:
				old = rec["_alleles"].get(cp)
				if old is not None and old != val:
					raise ValueError("Conflicting phased input-marker allele")
				rec["_alleles"][cp] = val
		return requested, rec, "observed_in_anchor_copies" if rec["_alleles"] else "input_found_phase_unavailable"
	return requested, None, "input_not_found_no_substitution"


def clean(r):
	return {k: v for k, v in r.items() if not k.startswith("_")}


def empty_metrics():
	return {k: None for k in binary_metrics(0, 0, 0, 0)}


def marker_metrics(record, target, cohort, allele=None):
	if record is None or not record.get("_alleles"):
		return dict(tag_allele=allele, **empty_metrics())
	options = [allele] if allele else [record["ref"], record["alt"]]
	scored = [dict(tag_allele=v, **count_variant(record["_alleles"], target, cohort, v)) for v in options]
	# Freeze input orientation from ALL and re-use it in population comparisons.
	return max(scored, key=lambda x: (x["r"] if x["r"] is not None else -2, x["f1"] or 0, x["tag_allele"]))


def pair_r2(a, b, cohort, aa, ba):
	if a is None or b is None or not aa or not ba:
		return None
	common = cohort & set(a["_alleles"]) & set(b["_alleles"])
	first = {cp for cp in common if a["_alleles"][cp] == aa}
	return count_variant(b["_alleles"], first, common, ba)["r2"]


def analyse_locus(auth, population_rows, options):
	lk = coord_key(auth)
	ident = dict(
		locus_key=lk,
		dataset_id=auth["dataset_id"],
		genome_build=auth["genome_build"],
		locus_id=auth["locus_id"],
		chr=dual_lead_chrom(auth["chr"]),
		core_start=dual_lead_integer(auth["core_start"]),
		core_end=dual_lead_integer(auth["core_end"]),
	)
	requested = auth.get("_input_lead_override", auth.get("named_anchor_id", auth.get("name", auth["locus_id"])))
	summary = dict(
		ident,
		input_lead_snp=requested,
		best_haplotype_id=None,
		best_haplotype_tier=None,
		best_haplotype_lineage=None,
		best_haplotype_archaic=None,
		best_tag_snp=None,
		best_tag_pos_1based=None,
		best_tag_allele=None,
		best_tag_r2=None,
		best_tag_ppv=None,
		best_tag_sensitivity=None,
		best_tag_quality=None,
		input_tag_r2=None,
		input_tag_allele=None,
		input_pos_1based=None,
		input_vs_best_snp_r2=None,
		same_input_and_best=None,
		n_equivalent_best_tags=0,
		family_tag_snp=None,
		family_tag_r2=None,
		tag_status="not_computed",
		tag_detail="",
	)
	result = dict(
		summary=summary, comparisons=[], ranked=[], population=[], by_population=[], haplotypes=[], gwas=[], files=[]
	)
	data, state, files = load_data(auth)
	result["files"] = files
	if data is None:
		summary.update(tag_status=state["status"], tag_detail=state.get("detail", ""))
		return result
	best, raw, ev, members = choose_haplotype(auth, data, options.min_match_sites, options.min_match_callable)
	if best is None:
		summary.update(
			tag_status="no_haplotype_meets_sequence_information_minimum",
			tag_detail="No SNP or lead-SNP substitution is made",
		)
		return result
	targets = {}
	selected = next(h for h in data["haps"] if h["hap_id"] == best["hap_id"])
	targets["best_haplotype"] = selected["_copies"]
	fam = set()
	for h in data["haps"]:
		e = ev.get((h["hap_id"], best["lineage"]), {})
		if best["evidence_tier"] == 3:
			ok = dual_lead_truth(e.get("diagnostic_candidate_pass")) and h["hap_id"] in members.get(
				best["lineage"], set()
			)
		else:
			ok = dual_lead_truth(e.get("diagnostic_candidate_pass"))
		if ok:
			fam |= h["_copies"]
	if fam:
		targets["candidate_family"] = fam
	groups = make_groups(data["copies"], population_rows)
	pool = build_variants(auth, data, best)
	req, iv, input_status = input_variant(auth, data)
	summary.update(
		best_haplotype_id=best["hap_id"],
		best_haplotype_tier={3: "tree_supported_candidate", 2: "sequence_candidate", 1: "raw_similarity_only"}[
			best["evidence_tier"]
		],
		best_haplotype_lineage=best["lineage"],
		best_haplotype_archaic=best["archaic"],
		best_haplotype_copies=best["n_copies"],
		best_haplotype_prop_match=best["prop_match"],
		n_stored_region_sites=len(data["in_core"]),
		n_archaic_matching_snp_candidates=len(pool),
		n_tested_haplotype_copies=len(data["copies"]),
		sequence_scope=data["scope"],
		raw_best_haplotype_id=raw["hap_id"],
		raw_best_archaic=raw["archaic"],
		raw_best_prop_match=raw["prop_match"],
		input_status=input_status,
		input_pos_1based=iv["pos_1based"] if iv else None,
		raw_match_is_not_introgression=int(best["evidence_tier"] == 1),
	)
	result["haplotypes"] = [
		dict(ident, role="evidence_prioritized_best", **best),
		dict(ident, role="raw_similarity_best", **raw),
	]
	for target_name, target in targets.items():
		ranked = []
		for v in pool:
			m = count_variant(v["_alleles"], target, data["copies"], v["tag_allele"])
			if qualify(m, options.min_tag_call_rate, options.min_tag_target_copies) == "eligible":
				ranked.append(dict(v, **m))
		ranked.sort(key=rank_tag)
		chosen = ranked[0] if ranked else None
		ties = []
		if chosen:
			ties = [
				r
				for r in ranked
				if (round(r["r2"], 12), round(r["f1"], 12), round(r["call_rate"], 12))
				== (round(chosen["r2"], 12), round(chosen["f1"], 12), round(chosen["call_rate"], 12))
			]
		im = marker_metrics(iv, target, data["copies"])
		sm = marker_metrics(chosen, target, data["copies"], chosen["tag_allele"] if chosen else None)
		inp = dict(
			ident,
			target_definition=target_name,
			role="input_lead",
			requested_snp=req,
			status=input_status,
			**(clean(iv) if iv else {}),
			**im,
		)
		if iv:
			inp["inside_locus"] = int(inside(ident["core_start"], ident["core_end"], iv["pos_1based"]))
		bst = {
			**ident,
			"target_definition": target_name,
			"role": "best_tag",
			"requested_snp": None,
			"status": (
				"strong_in_sample_proxy"
				if chosen and chosen["r2"] >= options.strong_tag_r2
				else "weak_in_sample_proxy"
				if chosen
				else "no_eligible_single_snp_tag"
			),
			**(clean(chosen) if chosen else {}),
			**sm,
			"n_equivalent_best_tags": len(ties),
		}
		result["comparisons"] += [inp, bst]
		for rank, r in enumerate(ranked[: options.top_tags], 1):
			result["ranked"].append(
				dict(
					ident,
					target_definition=target_name,
					rank=rank,
					**clean(r),
					tied_on_statistical_metrics=int(r in ties),
				)
			)
		# Fixed ALL-selected alleles/tags: these rows quantify portability.
		for g, cohort in sorted(groups.items()):
			for role, v, allele in [
				("input_lead", iv, im["tag_allele"]),
				("best_tag", chosen, chosen["tag_allele"] if chosen else None),
			]:
				met = marker_metrics(v, target, cohort, allele)
				result["population"].append(
					dict(
						ident,
						target_definition=target_name,
						population=g,
						role=role,
						snp_id=v.get("snp_id") if v else req if role == "input_lead" else None,
						variant_key=v.get("variant_key") if v else None,
						**met,
						metric_status=qualify(met, options.min_tag_call_rate, options.min_tag_target_copies)
						if v and v.get("_alleles")
						else "not_evaluable",
						selection_population="ALL",
						tag_reselected=0,
					)
				)
		# Optional per-superpopulation reselection: same global haplotype target,
		# not a new ancestry-specific haplotype definition.
		for g, cohort in sorted(groups.items()):
			if (
				g == "ALL"
				or (g.startswith("POP:") and options.tag_population_level != "all")
				or options.tag_population_level == "none"
			):
				continue
			candidates = []
			for v in pool:
				met = count_variant(v["_alleles"], target, cohort, v["tag_allele"])
				if qualify(met, options.min_tag_call_rate, options.min_tag_target_copies) == "eligible":
					candidates.append(dict(v, **met))
			candidates.sort(key=rank_tag)
			bb = candidates[0] if candidates else None
			result["by_population"].append(
				dict(
					ident,
					target_definition=target_name,
					population=g,
					status="computed" if bb else "not_evaluable_or_no_single_snp_tag",
					**(clean(bb) if bb else {}),
					tag_reselected=1,
				)
			)
		for role, v, met in [("input_lead", iv, im), ("best_tag", chosen, sm)]:
			if v and v.get("pos_1based"):
				result["gwas"].append(
					dict(
						ident,
						target_definition=target_name,
						role=role,
						snp_id=v.get("snp_id"),
						variant_key=v.get("variant_key"),
						pos_1based=v["pos_1based"],
						ref=v["ref"],
						alt=v["alt"],
						tag_allele=met["tag_allele"],
						tag_haplotype_r2=met["r2"],
						gwas_effect_allele=None,
						gwas_other_allele=None,
						gwas_beta=None,
						gwas_se=None,
						gwas_p=None,
						gwas_status="not_queried_user_followup",
						allele_note="tag_allele_is_not_assumed_to_be_risk_or_effect_allele",
					)
				)
		if target_name == "best_haplotype":
			summary.update(
				best_tag_snp=(chosen["snp_id"] or chosen["variant_key"]) if chosen else None,
				best_tag_pos_1based=chosen["pos_1based"] if chosen else None,
				best_tag_allele=sm["tag_allele"],
				best_tag_r2=sm["r2"],
				best_tag_ppv=sm["ppv"],
				best_tag_sensitivity=sm["sensitivity"],
				best_tag_quality=bst["status"],
				input_tag_r2=im["r2"],
				input_tag_allele=im["tag_allele"],
				input_vs_best_snp_r2=pair_r2(iv, chosen, data["copies"], im["tag_allele"], sm["tag_allele"]),
				same_input_and_best=int(iv["variant_key"] == chosen["variant_key"]) if iv and chosen else None,
				n_equivalent_best_tags=len(ties),
				tag_status="computed",
				tag_detail="Full stored locus SNP scan; no GWAS or input-lead constraint. r2 is SNP--haplotype, not SNP--SNP. In-sample, not externally validated.",
			)
		else:
			summary.update(
				family_tag_snp=(chosen["snp_id"] or chosen["variant_key"]) if chosen else None, family_tag_r2=sm["r2"]
			)
	return result


def analyse_all(authoritative, population_rows, options):
	out = {
		k: []
		for k in ["summary", "comparisons", "ranked", "population", "by_population", "haplotypes", "gwas", "files"]
	}
	for auth in authoritative:
		try:
			r = analyse_locus(
				auth,
				[p for p in population_rows if p.get("dataset_id", auth["dataset_id"]) == auth["dataset_id"]],
				options,
			)
		except (OSError, ValueError, KeyError, IndexError, csv.Error) as exc:
			# Do not take down the other loci. Invalid inputs must not produce tags.
			r = {k: [] for k in out}
			r["summary"] = dict(
				locus_key=coord_key(auth),
				dataset_id=auth["dataset_id"],
				genome_build=auth["genome_build"],
				locus_id=auth["locus_id"],
				chr=dual_lead_chrom(auth["chr"]),
				core_start=dual_lead_integer(auth["core_start"]),
				core_end=dual_lead_integer(auth["core_end"]),
				input_lead_snp=auth.get("_input_lead_override", auth.get("named_anchor_id", auth["locus_id"])),
				tag_status="invalid_sequence_inputs",
				tag_detail=str(exc),
			)
		out["summary"].append(r["summary"])
		for k in out:
			if k != "summary":
				out[k] += r[k]
	out["files"] = sorted(set(out["files"]))
	return out


def apply_loci_file(authoritative, path, coordinate_format="bed", allow_missing=False):
	"""Exact coordinate matching only. Column 4 is output metadata, never a filter.

	Existing sequence results must already cover the requested region. An input
	file is not allowed to silently redefine or widen a previously computed locus.
	"""
	requests = {}
	with Path(path).open(encoding="utf-8-sig") as f:
		for ln, text in enumerate(f, 1):
			if not text.strip() or text.lstrip().startswith(("#", "track", "browser")):
				continue
			fields = text.split()
			if len(fields) < 3:
				raise ValueError(f"Invalid loci line {ln}")
			try:
				st, en = int(fields[1]), int(fields[2])
			except ValueError:
				if ln == 1 and fields[0].lower() in {"chr", "chrom", "chromosome"}:
					continue
				raise ValueError(f"Invalid coordinates on loci line {ln}")
			if coordinate_format == "1based":
				st -= 1
			if st < 0 or en <= st:
				raise ValueError(f"Invalid interval on loci line {ln}")
			ck = (dual_lead_chrom(fields[0]), st, en)
			name = fields[3] if len(fields) > 3 else ""
			if ck in requests and requests[ck] != name:
				raise ValueError("Two different input lead labels for the same coordinates")
			requests[ck] = name
	result = []
	matched = set()
	for r in authoritative:
		ck = (dual_lead_chrom(r["chr"]), dual_lead_integer(r["core_start"]), dual_lead_integer(r["core_end"]))
		if ck in requests:
			r = dict(r)
			r["_input_lead_override"] = requests[ck]
			result.append(r)
			matched.add(ck)
	missing = set(requests) - matched
	if missing and not allow_missing:
		raise ValueError("Requested regions lack an exact saved PhyML coordinate match: " + str(sorted(missing)))
	if not result and not allow_missing:
		raise ValueError("No loci in input file")
	return result


# 🚩 prepare_review
"""Read-only GU evidence review. Standard-library Python; never edits GU inputs.

The legacy support table is labelled ANY-LOCUS OVERLAP throughout. Optional
per-copy comparisons use the actual candidate interval in evidence_haplotypes.tsv.
IBDmix is unphased: a match in the same individual does not establish phase.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from typing import Iterable

csv.field_size_limit(100_000_000)
prepare_review_VERSION = "2026-09-05.2"
sys.path.insert(0, str(Path(__file__).resolve().parent))

prepare_review_MISSING = {"", "NA", "NaN", "nan", "None", "null"}


def number(v, default=None):
	if v is None or str(v).strip() in prepare_review_MISSING:
		return default
	try:
		x = float(v)
		return x if math.isfinite(x) else default
	except (ValueError, TypeError):
		return default


def prepare_review_integer(v, default=None):
	x = number(v)
	if x is None:
		return default
	if x != int(x):
		raise ValueError(f"Expected integer, got {v!r}")
	return int(x)


def prepare_review_truth(v):
	return str(v).strip().lower() in {"1", "true", "t", "yes", "pass"}


def prepare_review_chrom(v):
	c = str(v).removeprefix("chr").removeprefix("CHR")
	return "X" if c == "23" else c


def key(r):
	return (r["dataset_id"], r["genome_build"], r["locus_id"])


def read_tsv(path: Path):
	opener = gzip.open if path.suffix == ".gz" else open
	with opener(path, "rt", encoding="utf-8-sig", newline="") as f:
		reader = csv.DictReader(f, delimiter="\t")
		if not reader.fieldnames:
			raise ValueError(f"Missing header: {path}")
		for line, r in enumerate(reader, 2):
			if None in r:
				raise ValueError(f"Unexpected extra TSV fields at {path}:{line}")
			yield r


def atomic_text(path: Path, text: str):
	path.parent.mkdir(parents=True, exist_ok=True)
	fd, name = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
	try:
		with os.fdopen(fd, "w", encoding="utf-8", newline="") as f:
			f.write(text)
		os.replace(name, path)
	finally:
		if os.path.exists(name):
			os.unlink(name)


def write_tsv(path: Path, rows: list[dict], fields=None):
	import io

	fields = fields or list(dict.fromkeys(k for r in rows for k in r))
	if not fields:
		fields = ["status"]
	s = io.StringIO(newline="")
	w = csv.DictWriter(s, fieldnames=fields, delimiter="\t", extrasaction="raise", lineterminator="\n")
	w.writeheader()
	for row in rows:
		w.writerow({k: ("" if v is None else v) for k, v in row.items()})
	atomic_text(path, s.getvalue())


def union(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
	out = []
	for st, en in sorted(set(intervals)):
		if st < 0 or en <= st:
			raise ValueError(f"Invalid half-open interval: {st}, {en}")
		if out and st <= out[-1][1]:
			out[-1] = (out[-1][0], max(en, out[-1][1]))
		else:
			out.append((st, en))
	return out


def overlap(a, b):
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def contains_pos1(interval, pos):
	"""VCF position is 1-based; GU candidate interval is 0-based half-open."""
	return pos is not None and interval[0] <= pos - 1 < interval[1]


def coverage(candidate, intervals):
	if candidate[0] < 0 or candidate[1] <= candidate[0]:
		raise ValueError(f"Invalid candidate interval {candidate}")
	ints = list(intervals)
	clipped = [(max(candidate[0], s), min(candidate[1], e)) for s, e in ints if overlap(candidate, (s, e))]
	reduced = union(clipped)
	bp = sum(e - s for s, e in reduced)
	size = candidate[1] - candidate[0]
	return dict(
		overlap_bp=bp,
		union_fraction=bp / size,
		best_single_fraction=max((overlap(candidate, x) / size for x in ints), default=0),
		n_covered_blocks=len(reduced),
		longest_block_fraction=max(((e - s) / size for s, e in reduced), default=0),
	)


def wilson(n, total, z=1.959963984540054):
	if total == 0:
		return None, None
	if not 0 <= n <= total:
		raise ValueError(f"Invalid binomial counts {n}/{total}")
	p = n / total
	den = 1 + z * z / total
	mid = (p + z * z / (2 * total)) / den
	rad = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / den
	return max(0.0, mid - rad), min(1.0, mid + rad)


def contingency(n, p, b, both):
	vals = (both, p - both, b - both, n - p - b + both)
	if min(vals) < 0:
		raise ValueError(f"Inconsistent support counts N={n}, P={p}, B={b}, both={both}")
	lo, hi = wilson(both, p)
	return dict(
		both=vals[0],
		phyml_only=vals[1],
		ibdmix_only=vals[2],
		neither=vals[3],
		ibdmix_given_phyml=both / p if p else None,
		conditional_ci_low=lo,
		conditional_ci_high=hi,
		phyml_given_ibdmix=both / b if b else None,
		jaccard=both / (p + b - both) if p + b - both else None,
	)


def read_authoritative(root):
	rows = []
	files = []
	for path in sorted((root / "phyml").glob("**/final/evidence_loci.tsv")):
		parts = path.relative_to(root).parts
		# normalize/phyml/<dataset>/<scope>/final/evidence_loci.tsv
		if len(parts) < 5:
			raise ValueError(f"Cannot resolve dataset from {path}")
		for r in read_tsv(path):
			r["dataset_id"] = r.get("dataset_id") or parts[1]
			r["_evidence_file"] = str(path)
			if not r.get("genome_build"):
				raise ValueError(f"Genome build missing in {path}; refusing to infer it")
			rows.append(r)
		files.append(path)
	seen = set()
	for r in rows:
		region_key = coord_key(r)
		if region_key in seen:
			raise ValueError(f"Duplicate authoritative region key {region_key}")
		seen.add(region_key)
	if not rows:
		raise ValueError("No phyml/<dataset>/<scope>/final/evidence_loci.tsv found")
	return rows, files


def anchor_assessment(r):
	if not prepare_review_truth(r.get("named_anchor_found")):
		return "not_evaluable_anchor_missing"
	if r.get("introgression_call") == "not_evaluable":
		return "not_evaluable_phyml"
	cs, ce, pos = (prepare_review_integer(r.get(x)) for x in ["candidate_start", "candidate_end", "named_anchor_pos"])
	if cs is None or ce is None:
		return "no_supported_candidate_for_anchor"
	if not prepare_review_truth(r.get("candidate_tree_pass")):
		return "tree_not_supported_for_anchor"
	if not contains_pos1((cs, ce), pos):
		return "outside_candidate_envelope_allele_origin_not_established"
	status = r.get("named_anchor_support_status", "")
	if status in {"named_anchor_not_phase_assigned", "lineage_anchor_uninformative_or_representation_mismatch"}:
		return "not_evaluable_anchor_phase_or_representation"
	if number(r.get("named_anchor_candidate_match")) == 0:
		return "candidate_does_not_match_lineage_anchor"
	if prepare_review_truth(r.get("named_anchor_candidate_match")):
		return "allele_linked_candidate_inside_envelope_not_final_origin_proof"
	return "not_evaluable_anchor_allele"


def prepare_summary(authoritative, legacy_rows, support_rows):
	legacy = {key(r): r for r in legacy_rows}
	idx = {(key(r), r["method"], r["population"]): r for r in support_rows}
	overview = []
	label_counts = __import__("collections").Counter(key(r) for r in authoritative)
	for r in authoritative:
		k = key(r)
		ambiguous = label_counts[k] > 1
		old = {} if ambiguous else legacy.get(k, {})
		ibd = {} if ambiguous else idx.get((k, "ibdmix", "ALL"), {})
		n = prepare_review_integer(ibd.get("n_tested"), prepare_review_integer(old.get("n_tested"), 0))
		p = prepare_review_integer(ibd.get("n_phyml_carriers"), prepare_review_integer(old.get("n_phyml_carriers"), 0))
		available = prepare_review_truth(ibd.get("method_available"))
		eligible = prepare_review_truth(ibd.get("evidence_eligible"))
		b = prepare_review_integer(ibd.get("n_carriers"), 0) if available else None
		a = prepare_review_integer(ibd.get("n_phyml_supported"), 0) if available else None
		c = contingency(n, p, b, a) if available else {x: None for x in contingency(0, 0, 0, 0)}
		call = r.get("introgression_call", "not_evaluable")
		if call == "not_evaluable":
			state = "phyml_not_evaluable"
		elif prepare_review_truth(r.get("candidate_tree_pass")):
			state = "tree_plus_broad_ibdmix_overlap" if eligible and a else "tree_candidate_needs_tract_review"
		elif available and eligible and b:
			state = "discordant_phyml_ibdmix"
		elif available and eligible and b == 0:
			state = "no_call_under_current_filters_not_proven_absent"
		else:
			state = "phyml_no_candidate_ibdmix_not_evaluable"
		cs, ce = prepare_review_integer(r.get("candidate_start")), prepare_review_integer(r.get("candidate_end"))
		pos = prepare_review_integer(r.get("named_anchor_pos"))
		row = dict(
			legacy_join_status="ambiguous_legacy_label_not_joined" if ambiguous else "unique_legacy_label",
			dataset_id=k[0],
			genome_build=k[1],
			locus_id=k[2],
			chr=prepare_review_chrom(r["chr"]),
			core_start=prepare_review_integer(r.get("core_start")),
			core_end=prepare_review_integer(r.get("core_end")),
			candidate_start=cs,
			candidate_end=ce,
			candidate_envelope_bp=(ce - cs if cs is not None and ce is not None else None),
			candidate_span_definition="min_start_to_max_end_across_candidate_types_not_shared_contiguous_tract",
			evidence_lineage=r.get("evidence_lineage") or None,
			ibdmix_comparison_lineage=ibd.get("source_class"),
			phyml_call=call,
			review_state=state,
			diagnostic_sites=prepare_review_integer(r.get("n_lineage_diagnostic_sites")),
			sequence_candidate_types=prepare_review_integer(r.get("n_candidate_haplotypes")),
			sequence_candidate_copies=prepare_review_integer(r.get("n_candidate_copies")),
			tree_bootstrap=number(r.get("candidate_tree_bootstrap")),
			candidate_purity=number(r.get("candidate_tree_candidate_purity")),
			candidate_sensitivity=number(r.get("candidate_tree_candidate_sensitivity")),
			tree_pass=int(prepare_review_truth(r.get("candidate_tree_pass"))),
			n_tested=None if ambiguous else n,
			phyml_carriers=None if ambiguous else p,
			ibdmix_carriers=b,
			ibdmix_available=int(available),
			ibdmix_eligible=int(eligible),
			ibdmix_availability_note=ibd.get("availability_note"),
			anchor_pos_1based=pos,
			anchor_consensus_allele_index=prepare_review_integer(r.get("named_anchor_lineage_consensus_allele_index")),
			anchor_inside_candidate_envelope=(
				int(contains_pos1((cs, ce), pos)) if cs is not None and ce is not None and pos is not None else None
			),
			anchor_assessment=anchor_assessment(r),
			upstream_anchor_call=r.get("anchor_linked_introgression_call"),
			locus_wide_ibdmix_max_lod=number(old.get("ibdmix_max_lod")),
			raw_best_lineage_descriptive_only=old.get("best_lineage"),
			raw_sequence_identity_descriptive_only=number(old.get("prop_match")),
			upstream_reason=r.get("evidence_reason"),
			support_definition="legacy_same_individual_same_lineage_any_locus_overlap",
			**c,
		)
		overview.append(row)

	def order(r):
		return (
			r["dataset_id"],
			r["genome_build"],
			99 if r["chr"] == "X" else prepare_review_integer(r["chr"], 98),
			r["core_start"] or 0,
		)

	overview.sort(key=order)
	population = []
	for r in support_rows:
		if r["method"] != "ibdmix" or label_counts[key(r)] != 1:
			continue
		n, p, b, a = (
			prepare_review_integer(r.get(x), 0)
			for x in ["n_tested", "n_phyml_carriers", "n_carriers", "n_phyml_supported"]
		)
		available = prepare_review_truth(r.get("method_available"))
		c = contingency(n, p, b, a) if available else {x: None for x in contingency(0, 0, 0, 0)}
		population.append(
			dict(
				dataset_id=r["dataset_id"],
				genome_build=r["genome_build"],
				locus_id=r["locus_id"],
				population=r["population"],
				source_class=r["source_class"],
				method_available=int(available),
				evidence_eligible=int(prepare_review_truth(r.get("evidence_eligible"))),
				n_tested=n,
				phyml_carriers=p,
				ibdmix_carriers=b if available else None,
				**c,
			)
		)
	return overview, population


def candidate_copies(authoritative):
	rows = []
	files = []
	seen = set()
	for r in authoritative:
		path = Path(r["_evidence_file"]).with_name("evidence_haplotypes.tsv")
		if not path.exists():
			continue
		files.append(path)
		tips = set(filter(None, r.get("candidate_tree_candidates_in_clade", "").split(",")))
		for h in read_tsv(path):
			if h.get("locus_id") != r["locus_id"] or h.get("diagnostic_lineage") != r.get("evidence_lineage"):
				continue
			if not prepare_review_truth(h.get("diagnostic_candidate_pass")):
				continue
			st, en = prepare_review_integer(h.get("candidate_start")), prepare_review_integer(h.get("candidate_end"))
			if st is None or en is None or en <= st or st < 0:
				raise ValueError(f"Invalid per-haplotype candidate interval at {path}: {h.get('hap_id')}")
			tokens = [x.strip() for x in h.get("copies", "").split(";") if x.strip()]
			if not tokens:
				raise ValueError(f"Candidate with no copy IDs: {path} {h.get('hap_id')}")
			for token in tokens:
				if ":" not in token:
					raise ValueError(f"Invalid sample:haplotype token {token!r} in {path}")
				sample, hap = token.rsplit(":", 1)
				ident = (coord_key(r), sample, hap, h["hap_id"], st, en)
				if ident in seen:
					continue
				seen.add(ident)
				rows.append(
					dict(
						locus_key=coord_key(r),
						dataset_id=r["dataset_id"],
						genome_build=r["genome_build"],
						locus_id=r["locus_id"],
						chr=prepare_review_chrom(r["chr"]),
						sample_id=sample,
						haplotype=hap,
						hap_id=h["hap_id"],
						candidate_start=st,
						candidate_end=en,
						source_class=r["evidence_lineage"],
						candidate_group=(
							"tree_supported_copy"
							if prepare_review_truth(r.get("candidate_tree_pass")) and h["hap_id"] in tips
							else "sequence_candidate_copy"
						),
						input_lead_snp=r.get("_input_lead_override", r.get("named_anchor_id", r["locus_id"])),
						upstream_anchor_pos_1based=prepare_review_integer(r.get("named_anchor_pos")),
						anchor_pos_1based=prepare_review_integer(
							r.get("_display_input_pos1", r.get("named_anchor_pos"))
						),
					)
				)
	return rows, files


def selected_segments(root, database, candidates, offset):
	# One pass through a compressed file. With SQLite, use only indexed locus windows.
	windows = defaultdict(list)
	for r in candidates:
		windows[(r["dataset_id"], r["genome_build"], r["chr"])].append((r["candidate_start"], r["candidate_end"]))
	windows = {k: union(v) for k, v in windows.items()}
	targets = {(r["dataset_id"], r["genome_build"], r["chr"], r["sample_id"]) for r in candidates}
	kept = {}
	nread = 0

	def accept(r):
		nonlocal nread
		nread += 1
		if r.get("method") != "ibdmix":
			return
		unit = (r.get("dataset_id"), r.get("genome_build"), prepare_review_chrom(r.get("chr")))
		if (*unit, r.get("sample_id")) not in targets or unit not in windows:
			return
		st, en = prepare_review_integer(r.get("start")), prepare_review_integer(r.get("end"))
		if st is None or en is None or st + offset < 0 or en <= st:
			raise ValueError("Invalid IBDmix interval; refusing silent row removal")
		st += offset
		en += offset
		if not any(overlap((st, en), q) for q in windows[unit]):
			return
		out = {
			x: r.get(x)
			for x in [
				"dataset_id",
				"genome_build",
				"sample_id",
				"method",
				"source",
				"source_class",
				"locus_id",
				"raw_file",
			]
		}
		out.update(chr=unit[2], start=st, end=en, score=number(r.get("score")))
		ident = (*unit, out["sample_id"], out["source"], st, en, out["score"])
		kept[ident] = out

	if database:
		uri = Path(database).resolve().as_uri() + "?mode=ro"
		with sqlite3.connect(uri, uri=True) as con:
			con.row_factory = sqlite3.Row
			cols = {r[1] for r in con.execute("PRAGMA table_info(segments)")}
			need = {"dataset_id", "genome_build", "chr", "start", "end", "sample_id", "method", "source_class"}
			if not need <= cols:
				raise ValueError("Database segments schema is incompatible")
			for (d, b, c), ints in windows.items():
				for st, en in ints:
					for row in con.execute(
						"SELECT * FROM segments WHERE dataset_id=? AND genome_build=? AND chr=? AND method='ibdmix' AND end>? AND start<?",
						(d, b, c, st - offset, en - offset),
					):
						accept(dict(row))
		source = str(Path(database).resolve())
	else:
		path = root / "summary" / "segments.tsv.gz"
		if not path.exists():
			return None, {"status": "not_computed_segments_missing"}
		for r in read_tsv(path):
			accept(r)
		source = str(path)
	return list(kept.values()), dict(
		status="computed", source=source, rows_read=nread, selected_calls=len(kept), ibdmix_coordinate_offset=offset
	)


def exact_support(candidates, segments, runs, populations, threshold):
	bysample = defaultdict(list)
	for s in segments:
		bysample[(s["dataset_id"], s["genome_build"], s["chr"], s["sample_id"])].append(s)
	run = {}
	for r in runs:
		if r["method"] == "ibdmix":
			k = (r["dataset_id"], r["genome_build"], prepare_review_chrom(r["chr"]))
			prev = run.get(k)
			if prev and prepare_review_truth(prev.get("evidence_eligible")) != prepare_review_truth(
				r.get("evidence_eligible")
			):
				raise ValueError(f"Conflicting IBDmix eligibility: {k}")
			run[k] = r
	pop = {(r["dataset_id"], r["sample_id"]): (r.get("population"), r.get("super_population")) for r in populations}
	out = []
	for c in candidates:
		unit = (c["dataset_id"], c["genome_build"], c["chr"])
		pool = bysample.get((*unit, c["sample_id"]), [])
		interval = (c["candidate_start"], c["candidate_end"])
		same = [s for s in pool if s["source_class"] == c["source_class"] and overlap(interval, (s["start"], s["end"]))]
		other = [
			s for s in pool if s["source_class"] != c["source_class"] and overlap(interval, (s["start"], s["end"]))
		]
		metric = coverage(interval, [(s["start"], s["end"]) for s in same])
		info = run.get(unit)
		completed = bool(info and info.get("status") == "complete")
		eligible = bool(completed and prepare_review_truth(info.get("evidence_eligible")))
		if same:
			status = (
				"observed_same_individual_overlap"
				if eligible
				else "observed_overlap_eligibility_unconfirmed_or_exploratory"
			)
		else:
			status = (
				"not_detected_individual_callability_unknown" if completed else "not_evaluable_no_completion_metadata"
			)
		known = int(bool(info))
		if not same and not completed:
			metric = {x: None for x in metric}
		pp, sp = pop.get((c["dataset_id"], c["sample_id"]), (None, None))
		r = dict(
			c,
			population=pp,
			super_population=sp,
			**metric,
			n_matching_reference_calls=len(same),
			matching_references=",".join(sorted({s["source"] for s in same if s["source"]})),
			matched_interval_max_lod=max((s["score"] for s in same if s["score"] is not None), default=None),
			other_lineage_overlap=int(bool(other)),
			other_lineage_max_lod=max((s["score"] for s in other if s["score"] is not None), default=None),
			anchor_inside_this_candidate=int(contains_pos1(interval, c["anchor_pos_1based"]))
			if c["anchor_pos_1based"]
			else None,
			ibdmix_covers_anchor_in_same_individual=(
				int(any(contains_pos1((s["start"], s["end"]), c["anchor_pos_1based"]) for s in same))
				if c["anchor_pos_1based"]
				else None
			),
			completion_metadata_present=known,
			evidence_eligible=int(eligible),
			strict_coverage_threshold=threshold,
			strict_coverage_pass=(
				int(metric["best_single_fraction"] >= threshold)
				if metric["best_single_fraction"] is not None and eligible
				else None
			),
			callability="unknown_individual_callable_bases_not_in_normalized_summary",
			phase_corroboration="individual_only_IBDmix_phase_unknown",
			status=status,
		)
		out.append(r)
	return out


def collapse_exact(rows, thresholds=(0.25, 0.5, 0.8)):
	grouped = defaultdict(list)
	for r in rows:
		grouped[(*key(r), r["candidate_group"], r.get("locus_key", ""))].append(r)
	out = []
	for k, rr in grouped.items():
		samples = {r["sample_id"] for r in rr}
		eligible = [r for r in rr if r["evidence_eligible"] == 1 and r["best_single_fraction"] is not None]
		for t in thresholds:
			passing = {r["sample_id"] for r in eligible if r["best_single_fraction"] >= t}
			seen = {r["sample_id"] for r in eligible}
			out.append(
				dict(
					locus_key=k[4] or None,
					dataset_id=k[0],
					genome_build=k[1],
					locus_id=k[2],
					candidate_group=k[3],
					threshold=t,
					n_candidate_individuals=len(samples),
					n_candidate_copies=len(rr),
					n_individuals_with_eligible_comparison=len(seen),
					n_individuals_passing=len(passing),
					support_fraction=len(passing) / len(seen) if seen else None,
					comparison="single_matching_lineage_IBDmix_call_covers_candidate_fraction",
					phase="not_confirmed",
				)
			)
	return out


def prepare_review_main(argv=None):
	ap = argparse.ArgumentParser(description=__doc__)
	ap.add_argument("--normalize", required=True, type=Path)
	ap.add_argument("--output", required=True, type=Path)
	ap.add_argument("--database", type=Path, help="Optional existing GU SQLite, read-only")
	ap.add_argument(
		"--summary-only", action="store_true", help="Skip strict IBDmix comparison; dual-lead analysis remains enabled"
	)
	ap.add_argument(
		"--loci",
		type=Path,
		help="Optional CHR START END lead file; exact coordinate match, fourth column is display metadata only",
	)
	ap.add_argument("--loci-format", choices=["bed", "1based"], default="bed")
	ap.add_argument("--skip-tags", action="store_true", help="Explicitly skip dual-lead sequence analysis")
	ap.add_argument("--min-match-sites", type=int, default=10)
	ap.add_argument("--min-match-callable", type=float, default=0.8)
	ap.add_argument("--min-tag-call-rate", type=float, default=0.95)
	ap.add_argument("--min-tag-target-copies", type=int, default=2)
	ap.add_argument(
		"--strong-tag-r2",
		type=float,
		default=0.8,
		help="In-sample proxy quality label; does not change the regional introgression call",
	)
	ap.add_argument("--top-tags", type=int, default=20)
	ap.add_argument(
		"--tag-population-level",
		choices=["none", "super", "all"],
		default="super",
		help="Reselect tag within superpopulations or all groups; fixed-tag metrics are always exported for all known populations",
	)
	ap.add_argument("--coverage-threshold", type=float, default=0.8)
	ap.add_argument(
		"--ibdmix-coordinate-offset",
		type=int,
		choices=[-1, 0, 1],
		default=0,
		help="Explicit audit-only shift of BOTH IBDmix bounds; default preserves normalized coordinates",
	)
	args = ap.parse_args(argv)
	if not 0 < args.coverage_threshold <= 1:
		ap.error("--coverage-threshold must be in (0,1]")
	if args.min_match_sites < 1 or args.min_tag_target_copies < 1 or args.top_tags < 1:
		ap.error("Site, copy and top-tag counts must be positive")
	if not all(0 < v <= 1 for v in [args.min_match_callable, args.min_tag_call_rate, args.strong_tag_r2]):
		ap.error("Match/tag fractions and r2 thresholds must be in (0,1]")
	root = args.normalize.resolve()
	out = args.output.resolve()
	if not root.is_dir():
		ap.error("normalize directory does not exist")
	if out == root or root in out.parents:
		ap.error("output must be outside the normalized input tree")
	if out in root.parents:
		ap.error("output must not be an ancestor of normalized input")
	auth, files = read_authoritative(root)
	if args.loci:
		auth = apply_loci_file(auth, args.loci, args.loci_format)
		files.append(args.loci.resolve())

	def load(name, required=True):
		path = root / "summary" / (name + ".tsv.gz")
		if not path.exists():
			if required:
				raise ValueError(f"Missing required summary: {path}")
			return []
		files.append(path)
		return list(read_tsv(path))

	legacy = load("locus_evidence")
	support = load("locus_method_support")
	runs = load("method_runs", False)
	overview, population = prepare_summary(auth, legacy, support)
	pops = load("sample_populations", False)
	if args.skip_tags:
		tags = {
			k: []
			for k in ["summary", "comparisons", "ranked", "population", "by_population", "haplotypes", "gwas", "files"]
		}
		for r in auth:
			tags["summary"].append(
				dict(
					locus_key=coord_key(r),
					dataset_id=r["dataset_id"],
					genome_build=r["genome_build"],
					locus_id=r["locus_id"],
					input_lead_snp=r.get("_input_lead_override", r.get("named_anchor_id", r["locus_id"])),
					tag_status="skipped_by_user",
				)
			)
	else:
		tags = analyse_all(auth, pops, args)
		files += tags.pop("files")
	ti = {t["locus_key"]: t for t in tags["summary"]}
	for r, a in zip(
		overview,
		sorted(
			auth,
			key=lambda z: (
				z["dataset_id"],
				z["genome_build"],
				99 if prepare_review_chrom(z["chr"]) == "X" else prepare_review_integer(z["chr"], 98),
				prepare_review_integer(z.get("core_start"), 0),
			),
		),
	):
		r["locus_key"] = coord_key(a)
		r["input_lead_snp"] = a.get("_input_lead_override", a.get("named_anchor_id", a.get("name", a["locus_id"])))
		t = ti.get(r["locus_key"], {})
		a["_display_input_pos1"] = (
			t.get("input_pos_1based")
			if "_input_lead_override" in a and a["_input_lead_override"] != a.get("named_anchor_id", a["locus_id"])
			else prepare_review_integer(a.get("named_anchor_pos"))
		)
		for field in [
			"best_haplotype_id",
			"best_haplotype_tier",
			"best_haplotype_lineage",
			"best_tag_snp",
			"best_tag_pos_1based",
			"best_tag_allele",
			"best_tag_r2",
			"best_tag_quality",
			"input_tag_r2",
			"input_tag_allele",
			"input_vs_best_snp_r2",
			"same_input_and_best",
			"n_equivalent_best_tags",
			"tag_status",
			"tag_detail",
		]:
			r[field] = t.get(field)
		if "_input_lead_override" in a and a["_input_lead_override"] != a.get("named_anchor_id", a["locus_id"]):
			r["upstream_anchor_snp"] = a.get("named_anchor_id", a["locus_id"])
			r["anchor_assessment"] = "not_evaluable_input_changed"
			r["anchor_pos_1based"] = t.get("input_pos_1based")
			r["anchor_consensus_allele_index"] = None
			r["upstream_anchor_call"] = None
	coordinates = {key(r): r["locus_key"] for r in overview}
	for r in population:
		r["locus_key"] = coordinates.get(key(r))
	copies = []
	tracks = []
	strict = []
	status = {"status": "not_computed_summary_only"}
	if not args.summary_only:
		candidates, cfiles = candidate_copies(auth)
		files += cfiles
		if not candidates:
			status = {"status": "not_computed_per_haplotype_candidates_missing_or_empty"}
		else:
			tracks, status = selected_segments(root, args.database, candidates, args.ibdmix_coordinate_offset)
			if tracks is not None:
				copies = exact_support(candidates, tracks, runs, pops, args.coverage_threshold)
				strict = collapse_exact(copies, tuple(sorted({0.25, 0.5, 0.8, args.coverage_threshold})))
			else:
				tracks = []
	manifest = {
		"version": prepare_review_VERSION,
		"created_utc": datetime.now(timezone.utc).isoformat(),
		"normalize_root": str(root),
		"output_root": str(out),
		"input_files": [],
		"n_loci": len(overview),
		"strict_comparison": status,
		"dual_lead": {
			"status_counts": dict(__import__("collections").Counter(t.get("tag_status") for t in tags["summary"])),
			"analysis_unit": "dataset + genome build + CHR:START-END (0-based half-open)",
			"input_lead_role": "annotation/comparator only; never a target or rank constraint",
			"parameters": {
				k: getattr(args, k)
				for k in [
					"min_match_sites",
					"min_match_callable",
					"min_tag_call_rate",
					"min_tag_target_copies",
					"strong_tag_r2",
					"top_tags",
					"tag_population_level",
				]
			},
			"selection_rule": "evidence tier, diagnostic match proportion, contiguous diagnostic matches, raw identity, informative sites, copies; tag ranked by copy-level r2 then F1/call rate; no GWAS",
			"scope_limit": "Best among stored biallelic SNPs with phased A/C/G/T calls. Exact input indels use anchor_copies.tsv. Not a genome-wide or all-variant VCF rerun.",
		},
		"coordinate_note": "Candidate start/end are GU 0-based half-open; anchor is 1-based. IBDmix input coordinates are preserved by default. Verify native coordinate provenance before endpoint-specific SNP claims.",
		"inference_limits": [
			"No introgression probability is calculated.",
			"Legacy overlap is not tract coverage.",
			"IBDmix phase and individual callable-base denominator are unknown.",
			"Multiple archaic references are not independent methods.",
			"No VCF, phylogenetic tree, IBDmix, TRACE or AS3 caller is rerun.",
			"Best tag and all proxy quality metrics are in-sample, not external validation or GWAS association.",
			"Population tag reselection keeps the same global haplotype target.",
		],
	}
	for path in sorted(set(files)):
		data = path.read_bytes()
		manifest["input_files"].append(
			{
				"path": str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
				"size": len(data),
				"sha256": hashlib.sha256(data).hexdigest(),
				"git_blob_sha1": hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest(),
			}
		)
	out.mkdir(parents=True, exist_ok=True)
	write_tsv(out / "loci_review.tsv", overview)
	write_tsv(out / "population_concordance.tsv", population)
	write_tsv(out / "method_runs.tsv", runs)
	write_tsv(out / "candidate_copy_support.tsv", copies)
	write_tsv(out / "strict_support_summary.tsv", strict)
	write_tsv(out / "selected_ibdmix_calls.tsv", tracks)
	for label, name in [
		("summary", "locus_tag_summary"),
		("comparisons", "lead_comparison"),
		("ranked", "haplotype_tag_candidates"),
		("population", "tag_population_metrics"),
		("by_population", "best_tag_by_population"),
		("haplotypes", "best_haplotypes"),
		("gwas", "gwas_lookup"),
	]:
		write_tsv(out / (name + ".tsv"), tags[label])
	payload = {
		"dual_lead": tags,
		"overview": overview,
		"population": population,
		"methods": runs,
		"copies": copies,
		"strict": strict,
		"tracks": tracks,
		"manifest": manifest,
	}
	atomic_text(out / "review_data.json", json.dumps(payload, ensure_ascii=False, allow_nan=False, indent=2))
	atomic_text(out / "manifest.json", json.dumps(manifest, ensure_ascii=False, indent=2))
	template = Path(__file__).with_name("review_template.html").read_text(encoding="utf-8")
	safe_json = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c")
	atomic_text(out / "review.html", template.replace("/*__GU_DATA__*/", safe_json))
	print(
		json.dumps(
			{
				"loci": len(overview),
				"candidate_copies_compared": len(copies),
				"strict_status": status["status"],
				"output": str(out),
			},
			ensure_ascii=False,
			indent=2,
		)
	)
	return 0


def prepare_review_cli():
	try:
		raise SystemExit(prepare_review_main())
	except (OSError, ValueError, KeyError, csv.Error, sqlite3.Error) as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		raise SystemExit(2)


SUBCOMMANDS = {}


def main():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] in SUBCOMMANDS:
		command = sys.argv.pop(1)
		return SUBCOMMANDS[command]()
	return prepare_review_cli()


if __name__ == "__main__":
	main()
