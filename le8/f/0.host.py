#!/usr/bin/env python3
"""0.host utilities and command entry points."""
from __future__ import annotations
# C1_LOCKPAPER::2026-10-17.c1-lockpaper-v1
import sys

sys.dont_write_bytecode = True
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

if "le8_common" not in sys.modules:
	_io_spec = spec_from_file_location("le8_common", Path(__file__).with_name("0.common.py"))
	_io_module = module_from_spec(_io_spec)
	sys.modules[_io_spec.name] = _io_module
	try:
		_io_spec.loader.exec_module(_io_module)
	except BaseException:
		sys.modules.pop(_io_spec.name, None)
		raise
import argparse, ast, hashlib, json, os, re, runpy, shutil, subprocess, tempfile, time
import numpy as np
import pandas as pd

cohort_HERE = Path(__file__).resolve().parent


def cohort_parser():
	p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
	p.add_argument("--Y", dest="trait", default="cvd_cad")
	p.add_argument("--biom", choices=["prot", "met"], default="prot")
	p.add_argument("--ukb-phe", default=os.getenv("UKB_PHE", "/mnt/d/data/ukb"))
	for name in [
		"phe-file",
		"omics-file",
		"met-map",
		"diagnosis-col",
		"end-date",
		"covariates",
		"residualize",
		"categorical",
		"group-col",
		"group-file",
		"ukb-rel-file",
		"outer-roster",
	]:
		p.add_argument("--" + name)
	p.add_argument("--met-input", choices=["raw", "named"])
	for name, default in [
		("id-col", "eid"),
		("baseline-col", "date_attend"),
		("death-col", "date_death"),
		("lost-col", "date_lost"),
		("disease-evidence-col", ""),
		("healthy-date-cols", ""),
		("exclude-features", ""),
	]:
		p.add_argument("--" + name, default=default)
	p.add_argument("--transform", choices=["none", "log1p"])
	p.add_argument("--seed", type=int, default=2026)
	p.add_argument("--r-bin", default="Rscript")
	p.add_argument("--max-samples", type=int, default=0)
	p.add_argument("--demo-features", type=int, default=80)
	p.add_argument("--demo", action="store_true")
	p.add_argument("--include-prevalent", action="store_true",
		help="Prepare a separate build-family historical source; never change incident eligibility")
	p.add_argument("--history-age-col", default="age")
	p.add_argument("--history-sex-col", default="sex")
	p.add_argument("--history-treatment-cols", default="")
	return p


def cohort_configure(a):
	common = runpy.run_path(str(cohort_HERE / "0.common.py"), run_name="le8_config")
	settings, sources, _ = common["shared_analysis_settings"](a)
	for name, value in settings.items():
		setattr(a, name, value)
	a.shared_config_sources = sources
	a.group_source_column = a.group_col
	if a.group_file:
		a.group_col = ".le8_family"
	if a.demo and a.outer_roster:
		raise ValueError("Demo cannot use a real outer roster")
	a.covariates = (
		a.covariates if a.covariates is not None else (a.shared_covariates or "age,sex,tdi,PC1,PC2")
	)
	a.outcome_type, a.target_col = "survival", a.trait
	a.diagnosis_col = a.diagnosis_col or "fod_icd10_" + a.trait
	base = Path(a.ukb_phe)
	a.phe_file = a.phe_file or str(base / "Rdata/all.rds")
	# Use the same cleaned, named matrices as native C1. Raw assays are an
	# explicit alternative, never an implicit second cohort/preprocessing path.
	a.omics_file = a.omics_file or str(base / f"Rdata/{a.biom}.rds")
	a.met_map = a.met_map or str(base / "common/met.lst")
	a.met_input = a.met_input or "named"
	a.residualize = "" if a.residualize is None else a.residualize
	a.categorical = "sex" if a.categorical is None else a.categorical
	# Numeric labels imported from R factors are quantitative unless the user
	# explicitly declares the variable categorical. This version participates in
	# cache identity; old age-dummy models remain available for read-only audits.
	a.clinical_encoding_version = "numeric-labels-v2-20261021"
	a.excluded_technical_covariates = []
	for name in ("covariates", "residualize", "categorical"):
		kept = []
		for col in cohort_words(getattr(a, name)):
			if cohort_technical_category(col):
				if col not in a.excluded_technical_covariates:
					a.excluded_technical_covariates.append(col)
			else:
				kept.append(col)
		setattr(a, name, ",".join(kept))
	if a.excluded_technical_covariates:
		cohort_log("SKIP", "technical_covariates", ",".join(a.excluded_technical_covariates))
	a.transform = a.transform or (
		"log1p" if a.biom == "met" and a.met_input == "raw" and not a.demo else "none"
	)
	if a.max_samples < 0 or a.demo_features < 3:
		raise ValueError("Invalid sample limit or demo feature count")
	return a


def cohort_assign_roles(p, a):
	"""Deterministic, outcome-blind family assignment for the current protocol."""
	from abm import ROLES, tie_hash

	groups = cohort_family_groups(p, a)
	outer = cohort_shared_outer_roles(p, a)
	# Independent hash stream: outer test families occupy the lowest outer hashes.
	values = np.array([tie_hash(g, a.seed + 7919 if outer is not None else a.seed) / 2**64 for g in groups])
	edges = np.array([0, 0.50, 0.60, 0.68, 0.76, 0.84])
	if outer is not None:
		edges = edges / edges[-1]
	role = np.full(len(p), "test", dtype="U25")
	for j, name in enumerate(ROLES):
		selected = (values >= edges[j]) & (values < edges[j + 1])
		if outer is not None:
			selected &= outer == "training"
		role[selected] = name
	if set(role) != set(ROLES) | {"test"}:
		raise ValueError("Each development role and test must have eligible participants")
	check_group_roles(groups, role)
	return role


cohort_RDS_EXPORT_SCRIPT = '# Read-only bridge. No sourcing of phenotype pipelines or installing R packages.\nargs <- commandArgs(trailingOnly = TRUE)\nstopifnot(length(args) == 3L)\nx <- readRDS(args[1])\nif (!is.data.frame(x)) stop("Expected one data.frame in RDS")\ncols <- readLines(args[3], warn = FALSE)\nif (length(cols)) {\n  absent <- setdiff(cols, names(x))\n  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))\n  x <- x[, cols, drop = FALSE]\n}\nfor (nm in names(x)) {\n  if (inherits(x[[nm]], "Date") || inherits(x[[nm]], "POSIXt"))\n    x[[nm]] <- format(x[[nm]], "%Y-%m-%d")\n  if (inherits(x[[nm]], "integer64")) x[[nm]] <- as.character(x[[nm]])\n}\n# Binary transfer preserves subnormal numeric metadata used as category codes.\nsaveRDS(x, args[2], compress = FALSE)\n\n'


def cohort_words(value):
	return [s.strip() for s in (value or "").split(",") if s.strip()]


def cohort_technical_category(name):
	"""Technical location/batch categories excluded from this ABM specification."""
	return bool(set(re.split(r"[.\s_:-]+", name.casefold())) & {"plate", "well", "center", "centre"})


def cohort_json_safe(value):
	if isinstance(value, dict):
		return {str(k): cohort_json_safe(v) for k, v in value.items()}
	if isinstance(value, (list, tuple, np.ndarray)):
		return [cohort_json_safe(v) for v in value]
	if isinstance(value, (np.integer, np.bool_)):
		return value.item()
	if isinstance(value, (float, np.floating)):
		return float(value) if np.isfinite(value) else None
	if isinstance(value, Path):
		return str(value)
	return value


def cohort_dump(path, value):
	path = Path(path)
	temp = path.with_suffix(path.suffix + ".tmp")
	temp.write_text(
		json.dumps(cohort_json_safe(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n",
		encoding="utf-8",
	)
	temp.replace(path)


def cohort_log(action, step, detail=""):
	print(f"[ABM] {action} {step}" + (f" | {detail}" if detail else ""), flush=True)


def cohort_parsed_date(raw, name):
	if pd.api.types.is_numeric_dtype(raw) and raw.notna().any():
		raise ValueError(f"{name}: numeric dates are ambiguous; use ISO dates or R Date")
	converted = pd.to_datetime(raw, errors="coerce")
	bad = raw.notna() & raw.astype(str).str.strip().ne("") & converted.isna()
	if bad.any():
		raise ValueError(f"Invalid dates in {name}")
	return converted


def cohort_outcomes(p, a):
	p = p.copy()
	if a.outcome_type == "quantitative":
		raw = p[a.target_col]
		value = pd.to_numeric(raw, errors="coerce")
		if (raw.notna() & value.isna()).any():
			raise ValueError(f"Non-numeric quantitative outcome: {a.target_col}")
		p["target"] = value
		p["eligible"] = np.isfinite(value)
		return p, dict(
			joined=len(p),
			eligible=int(p.eligible.sum()),
			missing_target=int((~p.eligible).sum()),
			outcome_type="quantitative",
		)
	baseline = cohort_parsed_date(p[a.baseline_col], a.baseline_col)
	diagnosis = cohort_parsed_date(p[a.diagnosis_col], a.diagnosis_col)
	censor = pd.concat(
		[
			cohort_parsed_date(p[a.death_col], a.death_col),
			cohort_parsed_date(p[a.lost_col], a.lost_col),
			pd.Series(pd.Timestamp(a.end_date), index=p.index),
		],
		axis=1,
	).min(axis=1)
	prevalent = diagnosis.notna() & (diagnosis <= baseline)
	unknown = pd.Series(False, index=p.index)
	if a.disease_evidence_col:
		evidence = pd.to_numeric(p[a.disease_evidence_col], errors="raise")
		unknown = evidence.gt(0) & diagnosis.isna()
	valid = baseline.notna() & (censor >= baseline)
	healthy = pd.Series(True, index=p.index)
	for col in cohort_words(a.healthy_date_cols):
		d = cohort_parsed_date(p[col], col)
		healthy &= ~(d.notna() & (d <= baseline))
	event = valid & ~prevalent & diagnosis.notna() & (diagnosis > baseline) & (diagnosis <= censor)
	p["time"] = (censor.where(~event, diagnosis) - baseline).dt.days / 365.25
	p["event"] = event.astype(int)
	p["prevalent"] = prevalent
	p["censor_date"] = censor
	p["endpoint_valid"] = valid & ~prevalent & ~unknown
	p["eligible"] = p.endpoint_valid & healthy & (p.time > 0)
	return p, dict(
		joined=len(p),
		eligible=int(p.eligible.sum()),
		prevalent=int(prevalent.sum()),
		same_day=int((diagnosis == baseline).sum()),
		invalid_followup=int((~valid).sum()),
		unknown_diagnosis=int(unknown.sum()),
		other_baseline_disease=int((~healthy & ~prevalent).sum()),
		diagnosis_after_censor=int((diagnosis > censor).sum()),
		incident=int(p.loc[p.eligible, "event"].sum()),
		outcome_type="survival",
	)


def cohort_validate_ids(frame, id_col):
	if id_col not in frame:
		raise ValueError(f"Missing ID column: {id_col}")
	if frame.columns.duplicated().any():
		raise ValueError("Duplicate column names")
	if frame[id_col].isna().any():
		raise ValueError("Missing person IDs")
	raw = frame[id_col]
	if pd.api.types.is_numeric_dtype(raw):
		num = raw.to_numpy(float)
		if not np.isfinite(num).all() or (num != np.floor(num)).any():
			raise ValueError("Numeric IDs must be finite integers; supply strings otherwise")
		frame[id_col] = raw.astype("int64").astype(str)
	else:
		frame[id_col] = raw.astype(str).str.strip()
	if frame[id_col].eq("").any() or frame[id_col].duplicated().any():
		raise ValueError("Empty or duplicate person IDs; select one baseline visit first")
	return frame


def cohort_read_table(path, id_col="eid", columns=None, r_bin="Rscript", string_columns=()):
	path = Path(path)
	if path.suffix.lower() == ".rds":
		import pyreadr

		if shutil.which(r_bin):
			# R selects columns before serialization, avoiding a second full RDS in Python.
			with tempfile.TemporaryDirectory(prefix="abm_rds_") as tmp:
				out = Path(tmp) / "selected.rds"
				cols = Path(tmp) / "columns.txt"
				cols.write_text("\n".join(columns or []))
				bridge = Path(tmp) / "export_rds.R"
				bridge.write_text(cohort_RDS_EXPORT_SCRIPT, encoding="utf-8")
				subprocess.run([r_bin, str(bridge), str(path), str(out), str(cols)], check=True)
				frame = pyreadr.read_r(str(out))[None].reset_index(drop=True)
		else:
			objects = pyreadr.read_r(str(path))
			if len(objects) != 1 or not isinstance(next(iter(objects.values())), pd.DataFrame):
				raise ValueError(f"{path}: expected a single R data.frame")
			frame = next(iter(objects.values()))
		# R NA_real_ can arrive as a signaling NaN. Canonicalize missing payloads
		# before NumPy/sklearn reductions, preserving all finite values and dtypes.
		frame = frame.mask(frame.isna(), np.nan)
	elif path.suffix.lower() == ".parquet":
		frame = pd.read_parquet(path, columns=columns)
	else:
		sep = "," if ".csv" in path.name.lower() else "\t"
		# Reject duplicate raw headers before pandas can silently suffix them.
		header = pd.read_csv(path, sep=sep, header=None, nrows=1).iloc[0].astype(str)
		if header.duplicated().any():
			raise ValueError(f"Duplicate header in {path}")
		frame = pd.read_csv(
			path,
			sep=sep,
			dtype={c: str for c in (id_col, *string_columns)},
			usecols=columns,
			low_memory=False,
		)
	if columns:
		missing = set(columns) - set(frame)
		if missing:
			raise ValueError(f"{path}: missing columns {sorted(missing)}")
		frame = frame[columns].copy()
	return cohort_validate_ids(frame, id_col)


def cohort_numeric(frame, columns, context, dtype="float32"):
	values = []
	for name in columns:
		raw = frame[name]
		val = pd.to_numeric(raw, errors="coerce")
		bad = raw.notna() & val.isna()
		if bad.any():
			raise ValueError(f"{context}: nonnumeric values in {name}")
		values.append(val.to_numpy(dtype=dtype, na_value=np.nan))
	if not values:
		raise ValueError(f"{context}: no molecular columns")
	result = np.column_stack(values)
	result[~np.isfinite(result)] = np.nan
	return result


def cohort_molecular_input_qc(raw, biom, demo=False):
	"""Treat negative metabolite concentrations as missing; preserve signed NPX."""
	if biom == "met" and not demo:
		return np.where(raw < 0, np.nan, raw)
	return raw


def cohort_safe_expression(expression, frame):
	"""Only arithmetic over explicit columns; never eval spreadsheet expressions."""

	def calculate(node):
		if isinstance(node, ast.Name):
			if node.id not in frame:
				raise KeyError(node.id)
			return pd.to_numeric(frame[node.id], errors="raise").to_numpy(float)
		if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
			return float(node.value)
		if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
			value = calculate(node.operand)
			return value if isinstance(node.op, ast.UAdd) else -value
		if isinstance(node, ast.BinOp):
			left, right = calculate(node.left), calculate(node.right)
			with np.errstate(divide="ignore", invalid="ignore"):
				if isinstance(node.op, ast.Div):
					return np.divide(left, right)
				if isinstance(node.op, ast.Mult):
					return left * right
				if isinstance(node.op, ast.Add):
					return left + right
				if isinstance(node.op, ast.Sub):
					return left - right
		raise ValueError(f"Unsupported metabolite expression: {expression}")

	return calculate(ast.parse(expression, mode="eval").body)


def cohort_map_metabolites(frame, mapping, id_col):
	# Match ukb/f/phe.R: baseline _i0 fields and common/met.lst.
	repeated = [c for c in frame if re.search(r"_i[1-9]\d*$", str(c))]
	frame = frame.drop(columns=repeated).copy()
	frame.columns = [re.sub(r"_i0$", "", str(c)) for c in frame]
	if frame.columns.duplicated().any():
		raise ValueError("Metabolite columns collide after removing _i0")
	mapping = Path(mapping)
	if not mapping.is_file():
		raise FileNotFoundError(f"Raw metabolite mapping required: {mapping}")
	spec = pd.read_csv(
		mapping,
		sep="\t",
		header=None,
		usecols=[0, 1],
		names=["expression", "feature"],
		dtype=str,
	)
	# Support both the original headerless mapping and the current UKB mapping.
	if len(spec) and tuple(spec.iloc[0]) == ("data_field", "met_name"):
		spec = spec.iloc[1:].copy()
	if spec.isna().any().any() or spec.feature.duplicated().any():
		raise ValueError("met.lst must contain unique, nonempty feature names")
	mapped = {id_col: frame[id_col].to_numpy()}
	audit = []
	for row in spec.itertuples():
		expr, feature = row.expression.strip(), row.feature.strip()
		try:
			values = (
				pd.to_numeric(frame[expr], errors="raise").to_numpy(float)
				if expr in frame
				else cohort_safe_expression(expr, frame)
			)
			values = np.broadcast_to(values, (len(frame),)).copy()
			invalid = ~np.isfinite(values)
			values[invalid] = np.nan
			mapped[feature] = values.astype("float32")
			audit.append(
				dict(
					feature=feature,
					expression=expr,
					status="mapped",
					nonfinite_to_missing=int(invalid.sum()),
				)
			)
		except KeyError as exc:
			audit.append(
				dict(
					feature=feature,
					expression=expr,
					status="missing_source",
					reason=str(exc),
				)
			)
	result = pd.DataFrame(mapped, index=frame.index)
	if result.shape[1] < 4:
		raise ValueError("Fewer than three metabolites map; check raw baseline fields and met.lst")
	return result, pd.DataFrame(audit)


def cohort_phenotype_columns(a):
	cols = [a.id_col] + cohort_words(a.covariates) + cohort_words(a.residualize)
	if getattr(a, "include_prevalent", False):
		cols += [a.history_age_col, a.history_sex_col] + cohort_words(a.history_treatment_cols)
	if a.outcome_type == "survival":
		cols += [a.baseline_col, a.diagnosis_col, a.death_col, a.lost_col]
		cols += cohort_words(a.healthy_date_cols)
		if a.disease_evidence_col:
			cols += [a.disease_evidence_col]
	else:
		cols += [a.target_col]
	if a.group_source_column:
		cols += [a.group_source_column]
	return list(dict.fromkeys(cols))


def cohort_generate_demo(a):
	rng = np.random.default_rng(a.seed)
	n = a.max_samples or 1500
	nf = a.demo_features
	state = rng.integers(0, 3, n)
	latent = rng.normal(size=(n, 6))
	latent[:, 0] += (state - 1) * 4
	x = (latent @ rng.normal(size=(6, nf)) + rng.normal(size=(n, nf))).astype("float32")
	age = rng.uniform(40, 75, n)
	sex = rng.integers(0, 2, n)
	plate = rng.integers(1, 6, n)
	x += plate[:, None] * 0.06
	x[rng.random(x.shape) < 0.03] = np.nan
	risk = 0.03 * (age - 55) + 0.35 * latent[:, 1] + 0.4 * (state - 1) * latent[:, 2]
	times = rng.exponential(18 * np.exp(-risk))
	censor = rng.uniform(7, 12, n)
	base = pd.Timestamp("2010-01-01")
	p = pd.DataFrame(
		{
			a.id_col: [f"demo_{i}" for i in range(n)],
			"age": age,
			"sex": sex,
			"tdi": rng.normal(size=n),
			"PC1": rng.normal(size=n),
			"PC2": rng.normal(size=n),
			"center": rng.choice(["A", "B", "C"], n),
			f"{a.biom}.plate": plate,
			a.baseline_col: "2010-01-01",
			a.death_col: pd.Series([None] * n, dtype=object),
			a.lost_col: (base + pd.to_timedelta(censor * 365.25, unit="D")).strftime("%Y-%m-%d"),
			a.diagnosis_col: (base + pd.to_timedelta(np.minimum(times, 100) * 365.25, unit="D")).strftime(
				"%Y-%m-%d"
			),
		}
	)
	if a.outcome_type == "quantitative":
		p[a.target_col] = (
			165 + 7 * sex + 2.5 * latent[:, 1] + 2 * (state - 1) * latent[:, 2] + rng.normal(size=n)
		)
	if a.group_col:
		p[a.group_source_column or a.group_col] = (np.arange(n) // 2).astype(str)
	features = [f"F{j:04d}" for j in range(nf)]
	omics = pd.DataFrame(x, columns=features)
	omics.insert(0, a.id_col, p[a.id_col])
	return p, omics


def cohort_attach_shared_groups(p, a):
	p = cohort_validate_ids(p.copy(), a.id_col)
	if a.group_file:
		m = cohort_read_table(a.group_file, "eid", r_bin=a.r_bin, string_columns=("group",))
		if "group" not in m or m["group"].isna().any() or m["group"].astype(str).str.strip().eq("").any():
			raise ValueError("LE8_GROUP_FILE requires complete eid,group columns")
		g = p[a.id_col].map(m.set_index("eid")["group"].astype(str))
		if g.isna().any():
			raise ValueError("Family mapping does not cover every ABM input participant")
		if a.group_source_column:
			old = p[a.group_source_column]
			if old.isna().any() or not old.astype(str).eq(g).all():
				raise ValueError("Conflicting family definitions in group file and phenotype column")
		p[a.group_col] = g
	if a.group_col:
		g = p[a.group_col]
		if g.isna().any() or g.astype(str).str.strip().eq("").any():
			raise ValueError("Incomplete family column")
		p[a.group_col] = g.astype(str)
	return p


def cohort_shared_outer_roles(p, a):
	if not getattr(a, "outer_roster", ""):
		return None
	m = cohort_read_table(a.outer_roster, "eid", r_bin=a.r_bin)
	if "role" not in m or not set(m.role) <= {"training", "test"}:
		raise ValueError("Shared outer roster requires eid,role with training/test")
	role = p[a.id_col].map(m.set_index("eid").role)
	if role.isna().any() or set(role) != {"training", "test"}:
		raise ValueError("Outer roster must cover input participants and both roles")
	groups = cohort_family_groups(p, a)
	if pd.DataFrame(dict(group=groups, role=role.to_numpy())).groupby("group").role.nunique().max() > 1:
		raise ValueError("Families cross shared outer roles")
	return role.to_numpy()


def cohort_prepare(a, out):
	mapping_audit = None
	if a.demo:
		p, omics = cohort_generate_demo(a)
	else:
		p = cohort_read_table(
			a.phe_file,
			a.id_col,
			cohort_phenotype_columns(a),
			a.r_bin,
			string_columns=([a.group_source_column] if a.group_source_column else []),
		)
		omics = cohort_read_table(a.omics_file, a.id_col, r_bin=a.r_bin)
		if a.biom == "met" and a.met_input == "raw":
			omics, mapping_audit = cohort_map_metabolites(omics, a.met_map, a.id_col)
		if a.biom == "prot":
			omics.columns = [c if c == a.id_col else str(c).upper() for c in omics]
			if omics.columns.duplicated().any():
				raise ValueError("Protein names collide after uppercasing")
	missing = set(cohort_phenotype_columns(a)) - set(p)
	if missing:
		raise ValueError(f"Missing phenotype columns: {sorted(missing)}")
	omics = omics[omics[a.id_col].isin(p[a.id_col])].copy()
	if a.max_samples and len(omics) > a.max_samples:
		omics = omics.sample(a.max_samples, random_state=a.seed)
	if omics.empty:
		raise ValueError("No matched phenotype/omics IDs")
	p = p.set_index(a.id_col).loc[omics[a.id_col]].reset_index()
	p = cohort_attach_shared_groups(p, a)
	features = [
		c
		for c in omics
		if c != a.id_col and c not in cohort_words(a.exclude_features) and not cohort_technical_category(c)
	]
	forbidden = set(cohort_phenotype_columns(a)) | {"event", "time", "target", "split", "Y"}
	collision = {c.casefold() for c in features} & {c.casefold() for c in forbidden}
	if collision:
		raise ValueError("Metadata/outcome columns in omics matrix: " + ", ".join(sorted(collision)))
	x = cohort_numeric(omics, features, "omics", dtype="float64")
	negative = np.sum(x < 0, axis=0)
	# Named/preprocessed matrices keep their supplied scale, including signed
	# values. Concentration-domain QC belongs only to explicitly raw metabolites.
	raw_metabolites = a.biom == "met" and a.met_input == "raw" and not a.demo
	if raw_metabolites:
		x = cohort_molecular_input_qc(x, a.biom)
	negative_to_missing = int(negative.sum()) if raw_metabolites else 0
	domain_invalid = np.sum(x <= -1, axis=0)
	if a.transform == "log1p" and domain_invalid.any():
		bad = [f"{name} ({int(n)})" for name, n in zip(features, domain_invalid) if n]
		raise ValueError("log1p requires x > -1; invalid input features: " + ", ".join(bad))
	transform_audit = dict(
		transform=a.transform,
		negative_values=int(negative.sum()),
		negative_features={name: int(n) for name, n in zip(features, negative) if n},
		negative_to_missing_values=negative_to_missing,
		log1p_domain_invalid_values=int(domain_invalid.sum()),
	)
	cohort_log(
		"DONE",
		"input_domain",
		f"transform={a.transform}; negative_values={int(negative.sum())}; "
		f"negative_features={int(np.count_nonzero(negative))}; negative_to_missing={negative_to_missing}; "
		f"log1p_domain_invalid={int(domain_invalid.sum())}",
	)
	# Preserve missing observations and invalid concentrations for train-only imputation.
	p.to_csv(out / "phenotype.csv", index=False)
	np.save(out / "raw.npy", x)
	(out / "features.txt").write_text("\n".join(features) + "\n")
	if mapping_audit is not None:
		mapping_audit.to_csv(out / "metabolite_mapping.csv", index=False)
	upstream = "supplied matrix; upstream preprocessing not inferred"
	if str(a.omics_file).endswith(f"Rdata/{a.biom}.rds"):
		upstream = "LE8 cleaned RDS: upstream all-person QC/imputation may already have occurred"
	cohort_dump(
		out / "input_audit.json",
		dict(
			n=len(p),
			features=len(features),
			source="SYNTHETIC" if a.demo else upstream,
			met_input=a.met_input,
			transform_audit=transform_audit,
			missing_fraction=float(np.mean(~np.isfinite(x))),
		),
	)
	return p


def cohort_family_groups(p, a):
	values = p[a.group_col if a.group_col else a.id_col]
	if values.isna().any() or values.astype(str).str.strip().eq("").any():
		raise ValueError("Missing participant/family identifiers")
	return values.astype(str).to_numpy()


ROOT = Path(__file__).resolve().parents[1]


def host_config(args):
	import shlex

	values = ["--Y", args.Y, "--biom", args.biom, "--seed", str(args.seed), "--r-bin", args.r_bin]
	for key, flag in [
		("ukb_phe", "ukb-phe"),
		("ukb_rel_file", "ukb-rel-file"),
		("group_file", "group-file"),
		("group_col", "group-col"),
		("outer_roster", "outer-roster"),
		("end_date", "end-date"),
		("diagnosis_col", "diagnosis-col"),
		("shared_covariates", "covariates"),
	]:
		value = getattr(args, key, None)
		# An unset shared adjustment is serialized as "" by the public worker.
		# Keep the established clinical default instead of passing --covariates "".
		if value is not None and (key != "shared_covariates" or value):
			values += ["--" + flag, str(value)]
	# Scope and random seed are owned by the public invocation.
	extra = shlex.split(args.host_args)
	if getattr(args, "include_prevalent", False) and "--include-prevalent" not in extra:
		extra += ["--include-prevalent"]
	if any(t.split("=")[0] in {"--Y", "--biom", "--seed"} for t in extra):
		raise ValueError("Set outcome, layer and seed using the public options")
	a = cohort_parser().parse_args(values + extra)
	if not a.demo:
		from c1_lock_roster import configure_roster
		a.analysis_root = args.analysis_root
		configure_roster(a)
	a = cohort_configure(a)
	a.cores = args.cores
	return a


def check_group_roles(groups, roles):
	groups = np.asarray(groups, str)
	roles = np.asarray(roles, str)
	order = np.argsort(groups, kind="stable")
	if np.any((groups[order][1:] == groups[order][:-1]) & (roles[order][1:] != roles[order][:-1])):
		raise ValueError("Roles divide a family; check the supplied grouping and roster")


def clinical_matrix(p, columns, categorical, training):
	"""Build-only vocabulary, no label-derived preprocessing or rare-level pooling.

    Numeric imputation/scale and constant-column removal remain fold-local inside
    twin. Dummy vocabulary is schema metadata; held-out levels stay explicitly missing.
    Technical factors are fixed adjustment inputs, outside the molecular budget.
    """
	arrays = []
	names = []
	audit = []
	for col in columns:
		s = p[col]
		numeric = pd.to_numeric(s.astype("string") if isinstance(s.dtype, pd.CategoricalDtype) else s, errors="coerce")
		# Infer the encoding from build rows only; a held-out malformed label
		# cannot change a continuous predictor into a dummy-variable panel.
		numeric_labels = bool((s.isna() | numeric.notna())[training].all())
		explicit_category = col in categorical
		if explicit_category or not numeric_labels:
			values = s.astype("string")
			levels = sorted(values[training].dropna().unique())
			if not levels:
				raise ValueError("No training levels for " + col)
			valid = values.isin(levels).to_numpy()
			for level in levels:
				x = values.eq(level).fillna(False).to_numpy(dtype=float)
				x[~valid] = np.nan
				arrays.append(x)
				names.append("clinical::" + col + "=" + str(level))
			audit.append(
				dict(
					column=col,
					encoding="training-vocabulary one-hot",
					encoding_reason="explicit_categorical" if explicit_category else "nonnumeric_labels",
					source_dtype=str(s.dtype),
					levels=list(map(str, levels)),
					unknown_or_missing_N=int((~valid).sum()),
				)
			)
		else:
			arrays.append(numeric.to_numpy(dtype=float, na_value=np.nan))
			names.append("clinical::" + col)
			audit.append(dict(column=col, encoding="numeric", source_dtype=str(s.dtype),
				encoding_reason="build_numeric_or_numeric_labels", explicitly_categorical=False,
				nonnumeric_labels_treated_missing_N=int((s.notna() & numeric.isna()).sum())))
	if not arrays:
		raise ValueError("Clinical covariates required")
	return np.column_stack(arrays), names, audit


def prepare_prevalent_source(p, raw, features, a, eligible, roles, clinical_columns, clinical_names):
	"""Historical Yang stays outside the incident target, risk set and test labels.

	All source families must be in build under the same outcome-blind split. A
	prevalent person absent from an explicit outer roster is unavailable, not
	a new training case. No outcome, diagnosis date, duration or medication field
	is ever appended to the incident model as a direct predictor.
	"""
	from abm import ROLES, tie_hash
	groups = cohort_family_groups(p, a)
	ids = p[a.id_col].astype(str).to_numpy()
	incident_groups = groups[eligible]
	inherited = dict(zip(incident_groups, roles))
	values = np.asarray([tie_hash(g, a.seed + 7919 if a.outer_roster else a.seed) / 2**64 for g in groups])
	edges = np.asarray([0, .50, .60, .68, .76, .84])
	outer = np.full(len(p), "training", dtype="U12")
	if a.outer_roster:
		roster = cohort_read_table(a.outer_roster, "eid", r_bin=a.r_bin)
		if "role" not in roster or not set(roster.role) <= {"training", "test"}:
			raise ValueError("Historical source requires the same explicit training/test roster")
		outer = p[a.id_col].map(roster.set_index("eid").role).fillna("unassigned").to_numpy(str)
		edges = edges / edges[-1]
	all_roles = np.full(len(p), "test", dtype="U25")
	for j, name in enumerate(ROLES):
		mask = (values >= edges[j]) & (values < edges[j+1])
		if a.outer_roster:
			mask &= outer == "training"
		all_roles[mask] = name
	all_roles[outer == "unassigned"] = "unassigned"
	# A declared roster must agree with the existing incident split; otherwise
	# this is a lineage conflict rather than a reason to borrow the person.
	for j, group in enumerate(groups):
		if group in inherited and all_roles[j] != "unassigned" and all_roles[j] != inherited[group]:
			raise ValueError("Historical source family crosses the incident role assignment")
	if not np.array_equal(all_roles[eligible], roles):
		raise ValueError("Historical source changed the incident roster")
	prev = p.prevalent.to_numpy(bool)
	use = prev & (all_roles == "build")
	build_mask = np.zeros(len(p), bool)
	build_mask[np.flatnonzero(eligible)[roles == "build"]] = True
	C, names, _ = clinical_matrix(p, clinical_columns, cohort_words(a.categorical), build_mask)
	if names != clinical_names:
		raise ValueError("Incident and historical source clinical schemas disagree")
	M = np.log1p(raw) if a.transform == "log1p" else np.asarray(raw)
	age = pd.to_numeric(p[a.history_age_col], errors="coerce").to_numpy(float)
	sex = p[a.history_sex_col].astype("string").fillna("unknown").to_numpy(str)
	duration = (cohort_parsed_date(p[a.baseline_col], a.baseline_col) -
		cohort_parsed_date(p[a.diagnosis_col], a.diagnosis_col)).dt.days.to_numpy(float) / 365.25
	treatment = np.full(len(p), np.nan)
	columns = cohort_words(a.history_treatment_cols)
	if columns:
		v = p[columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
		if np.any(np.isfinite(v) & ~np.isin(v, [0, 1])):
			raise ValueError("Historical treatment columns must be explicit baseline binary 0/1 indicators")
		treatment[np.any(v == 1, axis=1)] = 1.0
		treatment[np.all(v == 0, axis=1)] = 0.0
	return dict(
		incident_metadata={"age": age[eligible], "sex": sex[eligible]},
		prevalent_X=M[use], prevalent_ids=ids[use], prevalent_groups=groups[use],
		prevalent_roles=all_roles[use],
		prevalent_metadata={"age": age[use], "sex": sex[use], "duration_years": duration[use],
			"treatment": treatment[use]},
		feature_names=list(features),
		audit={"all_prevalent_N": int(prev.sum()), "build_prevalent_N": int(use.sum()),
			"heldout_or_unassigned_prevalent_N": int(np.sum(prev & ~use)),
			"treatment_columns": columns, "treatment_adjustment_available": bool(columns),
			"source_roles": "build families only", "prevalent_labels_added_to_incident": False,
			"incident_roster_unchanged": True})


def prepare_host(args, work):
	a = host_config(args)
	print(
		f"[LE8] START ABM inputs: biom={a.biom}; omics={a.omics_file}; "
		f"met_input={a.met_input}; transform={a.transform}; "
		f'covariates={a.covariates}; extra_adjustment={a.residualize or "none"}',
		flush=True,
	)
	prepared = work / "host_inputs_LOCAL_ONLY"
	prepared.mkdir(parents=True, mode=0o700)
	p = cohort_prepare(a, prepared)
	raw = np.load(prepared / "raw.npy")
	features = (prepared / "features.txt").read_text().splitlines()
	p, outcome_audit = cohort_outcomes(p, a)
	eligible = p.eligible.to_numpy(bool)
	eligible_p = p.loc[eligible].reset_index(drop=True)
	roles = cohort_assign_roles(eligible_p, a)
	groups = cohort_family_groups(eligible_p, a)
	check_group_roles(groups, roles)
	ids = eligible_p[a.id_col].astype(str).to_numpy()
	columns = list(dict.fromkeys(cohort_words(a.covariates) + cohort_words(a.residualize)))
	C, cn, encoding = clinical_matrix(eligible_p, columns, cohort_words(a.categorical), roles == "build")
	x = raw[eligible]
	if a.transform == "log1p":
		x = np.log1p(x)
	# Deviation preserves prevalent participants and does not exclude future cases.
	from deviation import fixed_roles

	all_groups = cohort_family_groups(p, a)
	dr = fixed_roles(all_groups, a.seed)
	mapping = {
		"build": "train",
		"tune_model": "tune",
		"tune_gate": "tune",
		"calibration_fit": "calibrate",
		"calibration_audit": "calibrate",
		"test": "test",
	}
	inherited = dict(zip(groups, [mapping[r] for r in roles]))
	dr = np.asarray([inherited.get(g, r) for g, r in zip(all_groups, dr)])
	check_group_roles(all_groups, dr)
	dc, dcn, de = clinical_matrix(p, columns, cohort_words(a.categorical), dr == "train")
	baseline = cohort_parsed_date(p[a.baseline_col], a.baseline_col)
	known = p.prevalent.astype(float).to_numpy()
	known[baseline.isna().to_numpy()] = np.nan
	if a.disease_evidence_col:
		unknown = (pd.to_numeric(p[a.disease_evidence_col], errors="raise") > 0) & p[a.diagnosis_col].isna()
		known[unknown.to_numpy()] = np.nan
	endpoint = pd.DataFrame(
		{"eid": p[a.id_col].astype(str), "known_cad_at_index": known, "index_date": baseline}
	)
	if args.objective_file:
		from le8_common import id_table

		extra = id_table(args.objective_file).reset_index(drop=True)
		# A supplied baseline status takes precedence; future-event columns never do.
		fields = [
			c
			for c in extra
			if c
			in {
				"eid",
				"known_cad_at_index",
				"index_date",
				"objective_date",
				"objective_cad",
				"objective_modality",
				"intervening_cad",
				"intervening_treatment",
			}
		]
		for col in set(fields) & set(endpoint) - {"eid"}:
			endpoint = endpoint.drop(columns=col)
		endpoint = endpoint.merge(extra[fields], on="eid", how="left", validate="one_to_one")
	audit = dict(
		outcome=outcome_audit,
		split_source="explicit_outer_roster" if a.outer_roster else "fresh_outcome_blind_family_assignment",
		roster_sha256=hashlib.sha256(
			json.dumps(list(zip(ids.tolist(), roles.tolist())), separators=(",", ":")).encode()
		).hexdigest(),
		clinical_encoding=encoding,
		deviation_encoding=de,
		clinical_columns=columns,
		technical_adjustment="fixed covariate inputs; no full-cohort residualization",
		family_status="verified_supplied_groups" if a.group_col else "unavailable_individual_ID_groups_only",
		independent_family_support_verified=bool(a.group_col),
		validation_status="exploratory",
		input_config=vars(a),
	)
	from c1_lock_roster import RELEASE, eid_hash, roster_record, grouping_audit
	audit.update(lockpaper_release=RELEASE, ABM_test_eid_hash=eid_hash(ids[roles == "test"]),
	             outer_roster=a.outer_roster,
	             outer_roster_sha256=roster_record(a.outer_roster)["roster_sha256"] if a.outer_roster else "synthetic")
	if a.outer_roster:
		audit.update(grouping_audit(roster_record(a.outer_roster)))
	# This flag is emitted only after checking the prepared, eligible ABM roster.
	audit["family_disjoint_ABM_roles"] = bool(a.group_col)
	history = None
	if getattr(a, "include_prevalent", False):
		history = prepare_prevalent_source(p, raw, features, a, eligible, roles, columns, cn)
		audit["prevalent_source"] = history["audit"]
	return (
		dict(
			X=np.c_[x, C],
			time=eligible_p.time.to_numpy(),
			event=eligible_p.event.to_numpy(),
			eid=ids,
			family_group=groups,
			role=roles,
			feature_names=features + cn,
			clinical_feature_names=cn,
		),
		dict(
			M=raw,
			C=dc,
			covariate_names=dcn,
			features=features,
			endpoint=endpoint,
			eid=np.asarray(p[a.id_col], str),
			group=np.asarray(all_groups, str),
			role=np.asarray(dr, str),
			prevalent_aux=history,
		),
		audit,
	)


# 🚩 run host
# Run twin confidence and molecular deviation from current cohort inputs.
from dataclasses import replace
from types import SimpleNamespace
import fcntl
import uuid
from threadpoolctl import threadpool_limits
from abm import json_write
from abm import load_config
from le8_common import write_tables, load_pgs, merge_native_scores
import deviation as deviation
from c1_abm_compare import run_comparison
from c1_lock_roster import grouping_audit


def parser():
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("--Y", required=True)
	p.add_argument("--biom", choices=["prot", "met"], required=True)
	p.add_argument("--analysis-root", type=Path, required=True)
	p.add_argument("--config", required=True)
	p.add_argument("--seed", type=int, default=2026)
	p.add_argument("--cores", type=int, default=4)
	p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE")
	p.add_argument("--r-bin", default="Rscript")
	for flag in [
		"ukb-phe",
		"ukb-rel-file",
		"group-file",
		"group-col",
		"outer-roster",
		"end-date",
		"diagnosis-col",
		"shared-covariates",
		"deviation-config",
		"objective-file",
		"pgs-root",
		"pgs-mapping",
	]:
		p.add_argument("--" + flag)
	p.add_argument("--host-args", default="")
	p.add_argument("--remerge-pgs", action="store_true")
	p.add_argument("--preflight", action="store_true")
	return p


def publish(stage, target):
	"""Retain the previous complete run and publish only a successful replacement."""
	previous = None
	if target.exists():
		previous = target.parent / "_history" / (target.name + "_" + uuid.uuid4().hex[:12])
		previous.parent.mkdir(exist_ok=True)
		target.rename(previous)
	try:
		stage.rename(target)
	except BaseException:
		if previous:
			previous.rename(target)
		raise


def completed(target):
	marker = target / "COMPLETE.json"
	if not marker.is_file():
		return False
	saved = json.loads(marker.read_text())
	for name, digest in saved.get("files", {}).items():
		path = target / name
		if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
			return False
	return bool(saved.get("files"))


def check_input_reuse(target, current):
	"""Do not silently reuse a raw-assay fit under the cleaned-RDS contract."""
	if target.name == "abm_twin":
		from abm import VERSION
		marker = target / "fit/MODEL_FROZEN.json"
		if not marker.is_file() or json.loads(marker.read_text()).get("version") != VERSION:
			raise ValueError("C1 LOCKPAPER requires a fresh molecular-only fit; run c1_correlate,c1_abm --replace TRUE")
	contract = target / "host_contract.json"
	if not contract.is_file():
		return
	previous = json.loads(contract.read_text()).get("input_config", {})
	if not previous:
		return
	keys = ("phe_file", "omics_file", "met_input", "transform", "covariates", "residualize", "categorical", "clinical_encoding_version", "outer_roster", "group_file", "group_source_column", "seed")
	changed = [key for key in keys if previous.get(key) != getattr(current, key)]
	if changed:
		raise ValueError(
			"ABM input configuration differs from the saved fit: "
			+ ",".join(changed)
			+ ". Run c1_abm --replace TRUE to fit the new inputs and archive the previous results: "
			+ str(target)
		)


def finish(stage, audit):
	contract = stage / "host_contract.json"
	existing = json.loads(contract.read_text()) if contract.is_file() else {}
	json_write(contract, {**existing, **audit})
	# Only small integrity anchors; do not rehash every large local array on resume.
	files = {
		str(p.relative_to(stage)): hashlib.sha256(p.read_bytes()).hexdigest()
		for p in stage.rglob("*")
		if p.is_file() and (p.name in {"MODEL_FROZEN.json", "table_manifest.json", "host_contract.json"})
	}
	json_write(stage / "COMPLETE.json", dict(status="completed", files=files))


def pgs_path(args):
	if args.pgs_root:
		return Path(args.pgs_root)
	cached = args.analysis_root / "_pgs_fullprecision" / args.biom / (args.biom + ".pgs.fullprecision.gz")
	return cached if cached.is_file() else Path("/mnt/f/gwas") / args.biom / "pgs" / (args.biom + ".pgs.gz")


def run_deviation(args, data, stage, config, audit):
	source = pgs_path(args)
	if args.remerge_pgs:
		manifest = Path("/mnt/f/gwas") / args.biom / "pgs" / (args.biom + ".pgs.files.tsv")
		source = args.analysis_root / "_pgs_fullprecision" / args.biom / (args.biom + ".pgs.fullprecision.gz")
		if not source.exists():
			temp = source.with_name("_pending_" + uuid.uuid4().hex + ".gz")
			merge_native_scores(manifest, temp)
			temp.rename(source)
	if not source.exists():
		stage.mkdir()
		write_tables(
			stage / "aggregate",
			{
				"status": pd.DataFrame(
					[
						dict(
							status="PGS_input_unavailable",
							source=str(source),
							validation_status="not_validated_for_diagnosis",
						)
					]
				)
			},
			"deviation",
			args.biom,
			args.Y,
			"unavailable",
		)
		print("[LE8 twin] Deviation unavailable: missing PGS " + str(source), flush=True)
		finish(stage, audit)
		return
	names, G, pgs_audit = load_pgs(source, args.biom, data["eid"], data["features"], args.pgs_mapping)
	if len(names) < 2:
		raise ValueError("Fewer than two exactly matched molecular PGSs")
	indices = [data["features"].index(n) for n in names]
	M = data["M"][:, indices]
	role = data["role"]
	dev = role != "test"
	eligible, y, reason = deviation.objective_endpoint(data["endpoint"], config)
	known = pd.to_numeric(data["endpoint"].known_cad_at_index, errors="coerce").to_numpy()
	stage.mkdir()
	local = stage / "_LOCAL_ONLY"
	local.mkdir(mode=0o700)
	fields = dict(
		M=M,
		G=G,
		C=data["C"],
		eid=data["eid"],
		group=data["group"],
		role=role,
		reference=known == 0,
		y=y,
		endpoint_eligible=eligible,
	)
	schema = dict(feature_names=np.asarray(names), covariate_names=np.asarray(data["covariate_names"]))
	np.savez_compressed(local / "development.npz", **{k: v[dev] for k, v in fields.items()}, **schema)
	np.savez_compressed(
		local / "test_X.npz", **{k: fields[k][~dev] for k in ["M", "G", "C", "eid", "group"]}, **schema
	)
	pd.DataFrame(
		dict(
			eid=data["eid"][~dev], y=y[~dev], endpoint_eligible=eligible[~dev], exclusion_reason=reason[~dev]
		)
	).to_csv(local / "test_truth.tsv.gz", sep="\t", index=False)
	pd.DataFrame(dict(eid=data["eid"], group=data["group"], role=role)).to_csv(
		local / "roster.tsv.gz", sep="\t", index=False
	)
	from dataclasses import asdict

	json_write(stage / "config.json", asdict(config))
	flow = (
		pd.DataFrame(dict(role=role, reason=reason)).groupby(["role", "reason"]).size().reset_index(name="N")
	)
	write_tables(
		stage / "aggregate",
		{
			"pgs_quality_audit": pgs_audit,
			"grouping_audit": pd.DataFrame([grouping_audit(audit)]),
			"endpoint_flow": flow,
			"support_status": pd.DataFrame(
				[
					{
						"family_status": audit["family_status"],
						"validation_status": "exploratory",
						"objective_endpoint_N": int(eligible.sum()),
					}
				]
			),
		},
		"deviation",
		args.biom,
		args.Y,
		"host",
	)
	A = SimpleNamespace
	deviation.fit(A(data=local / "development.npz", config=stage / "config.json", output=stage / "fit"))
	model = stage / "fit/_LOCAL_ONLY/model.joblib"
	deviation.project(A(model=model, data=local / "test_X.npz", output=local / "prediction.joblib"))
	deviation.evaluate(
		A(
			model=model,
			predictions=local / "prediction.joblib",
			truth=local / "test_truth.tsv.gz",
			output=str(stage / "results"),
			anchors="L_VLDL_TG.pct,Total_TG,ApoB,PCSK9,LPA,GDF15",
		)
	)
	finish(stage, audit)


def main(args):
	from c1_lock_roster import RELEASE
	scope = args.analysis_root / args.Y / args.biom / "c1_abm"
	cfg = replace(load_config(args.config), trait=args.Y, layer=args.biom, seed=args.seed, cores=args.cores, final_multiroute_release=RELEASE)
	cfg.validate()
	dc = deviation.read_config(args.deviation_config or ROOT / f"config/deviation_{args.biom}.json")
	dc = replace(dc, trait=args.Y, layer=args.biom, seed=args.seed)
	dc.validate()
	# Explicit alternative endpoints have distinct directories and completion
	# identities; they cannot overwrite or resume an objective-plaque analysis.
	deviation_target = (
		scope / "molecular_deviation"
		if dc.task == "objective_plaque"
		else scope / dc.task / "molecular_deviation"
	)
	if args.Y != "cvd_cad":
		# Do not reinterpret another disease's prevalence as known CAD.
		deviation_supported = False
	else:
		deviation_supported = True
	if args.preflight:
		a = host_config(args)
		from c1_lock_roster import PRIMARY_SELECTOR, COVERAGES, PREDICTABILITY_TARGET_MODEL, sha
		print(f"[LE8 C1 LOCK] outer_roster={a.outer_roster}; SHA256={sha(a.outer_roster) if a.outer_roster else 'synthetic'}; "
		      f"release={RELEASE}; seed={cfg.seed}; primary_selector={PRIMARY_SELECTOR}; "
		      f"primary_coverage={cfg.primary_coverage:.2f}; coverage_grid={COVERAGES}; "
		      f"selector_bootstrap_requested={cfg.bootstrap}; predictability_target_model={PREDICTABILITY_TARGET_MODEL}; "
		      f"grouping_unit={'family' if a.group_col else 'individual_unverified'}; "
		      f"independent_family_support_verified={bool(a.group_col)}", flush=True)
		if not a.demo:
			for path in [
				a.phe_file,
				a.omics_file,
				*([a.met_map] if a.biom == "met" and a.met_input == "raw" else []),
			]:
				if not Path(path).is_file():
					raise FileNotFoundError(path)
		for path in [args.objective_file, args.pgs_mapping, args.group_file, args.outer_roster]:
			if path and not Path(path).is_file():
				raise FileNotFoundError(path)
		if args.replace != "TRUE":
			for target in (scope / "abm_twin", deviation_target):
				if completed(target):
					check_input_reuse(target, a)
		print(
			"[LE8 twin] Preflight passed; family support="
			+ ("supplied" if a.group_col else "unverified; individual IDs only")
			+ "; objective endpoint="
			+ ("supplied" if args.objective_file else "unavailable; reference calibration only"),
			flush=True,
		)
		print(
			f"[LE8] DONE ABM input preflight: biom={a.biom}; omics={a.omics_file}; "
			f"met_input={a.met_input}; transform={a.transform}; covariates={a.covariates}; "
			f'extra_adjustment={a.residualize or "none"}',
			flush=True,
		)
		return
	scope.mkdir(parents=True, exist_ok=True)
	lock_key = hashlib.sha256(str(scope.resolve()).encode()).hexdigest()[:20]
	with open("/tmp/le8-abm-" + lock_key + ".lock", "a") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		targets = {"abm_twin": scope / "abm_twin", "molecular_deviation": deviation_target}
		todo = {name for name, target in targets.items() if args.replace == "TRUE" or not completed(target)}
		if args.replace != "TRUE" and set(targets) - todo:
			current = host_config(args)
			for name in set(targets) - todo:
				check_input_reuse(targets[name], current)
		if not todo:
			print("[LE8 twin] SKIP: both twin fits completed; refreshing C1 comparison if needed", flush=True)
			with tempfile.TemporaryDirectory(prefix="_abm_staging_", dir=scope) as work:
				run_comparison(args, scope, Path(work), publish, completed, finish)
				refresh_performance(args, scope, Path(work))
				try:
					subprocess.run([sys.executable,str(Path(__file__).with_name("c1_article.py")),"--analysis-root",str(args.analysis_root),"--Y",args.Y],check=True)
				except Exception as exc: print("[LE8 C1 LOCK] manuscript reporting warning:",exc,file=sys.stderr)
			return
		work = Path(tempfile.mkdtemp(prefix="_abm_staging_", dir=scope))
		try:
			abm, data, audit = prepare_host(args, work)
			unit = "family" if audit["independent_family_support_verified"] else "individual_unverified"
			cfg = replace(cfg, grouping_unit=unit)
			dc = replace(dc, grouping_unit=unit)
			if "abm_twin" in todo:
				print("[LE8 twin] START twin fit/project/evaluate", flush=True)
				from abm import run_prepared_numeric_cohort

				stage = work / "abm_twin"
				result = run_prepared_numeric_cohort(**abm, output=stage, config=cfg)
				write_tables(
					stage / "aggregate",
					{
						"grouping_audit": pd.DataFrame([grouping_audit(audit)]),
						"support_status": pd.DataFrame(
							[
								{
									"family_status": audit["family_status"],
									"independent_family_support_verified": audit[
										"independent_family_support_verified"
									],
									"validation_status": "exploratory",
									"split_source": audit["split_source"],
								}
							]
						)
					},
					"ABM twin",
					args.biom,
					args.Y,
					"host",
				)
				finish(stage, audit)
				publish(stage, scope / "abm_twin")
				print("[LE8 twin] DONE twin", flush=True)
			del abm
			if "molecular_deviation" in todo:
				stage = work / "molecular_deviation"
				if deviation_supported:
					print("[LE8 twin] START molecular deviation", flush=True)
					run_deviation(args, data, stage, dc, audit)
				else:
					stage.mkdir()
					write_tables(
						stage / "aggregate",
						{"status": pd.DataFrame([dict(status="not_applicable_non_CAD_trait", trait=args.Y)])},
						"deviation",
						args.biom,
						args.Y,
						"not_applicable",
					)
					finish(stage, audit)
				deviation_target.parent.mkdir(parents=True, exist_ok=True)
				publish(stage, deviation_target)
				print("[LE8 twin] DONE deviation", flush=True)
			run_comparison(args, scope, work, publish, completed, finish)
			refresh_performance(args, scope, work)
			try:
				subprocess.run([sys.executable,str(Path(__file__).with_name("c1_article.py")),"--analysis-root",str(args.analysis_root),"--Y",args.Y],check=True)
			except Exception as exc: print("[LE8 C1 LOCK] manuscript reporting warning:",exc,file=sys.stderr)
		except BaseException as exc:
			# The public console recognizes this prefix before the outer
			# CalledProcessError hides the actual cause behind a long command.
			print(
				f"[LE8] FAIL c1_abm Y={args.Y} biom={args.biom}: {type(exc).__name__}: {exc}",
				file=sys.stderr,
				flush=True,
			)
			print(
				"[LE8 twin] Failed; earlier published results retained. Local diagnostics: " + str(work),
				file=sys.stderr,
				flush=True,
			)
			raise
		else:
			shutil.rmtree(work)


def refresh_performance(args, scope, work):
	from final import (
		performance_collect as collect,
		performance_write_snapshot as write_snapshot,
		performance_complete,
	)

	if args.replace != "TRUE" and performance_complete(scope / "performance"):
		print(
			f"[LE8] SKIP performance Y={args.Y} biom={args.biom} stage=c1_abm: saved report files validated",
			flush=True,
		)
		return
	stage = work / "performance"
	tables = collect(args.analysis_root, args.analysis_root, args.Y, args.biom, "c1_abm")
	write_snapshot(tables, stage, f"{args.Y} / {args.biom} | through c1_abm")
	publish(stage, scope / "performance")


if __name__ == "__main__":
	os.umask(0o077)
	args = parser().parse_args()
	with threadpool_limits(args.cores):
		main(args)
