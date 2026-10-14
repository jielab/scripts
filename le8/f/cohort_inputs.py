"""Current cohort loading, assay mapping and incident outcome definitions."""
from pathlib import Path
import argparse, ast, hashlib, json, os, re, runpy, shutil, subprocess, tempfile, time
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument('--Y', dest='trait', default='cvd_cad')
    p.add_argument('--biom', choices=['prot', 'met'], default='prot')
    p.add_argument('--ukb-phe', default=os.getenv('UKB_PHE', '/mnt/d/data/ukb/phe'))
    for name in ['phe-file', 'omics-file', 'met-map', 'diagnosis-col', 'end-date',
                 'covariates', 'residualize', 'categorical', 'group-col', 'group-file', 'outer-roster']:
        p.add_argument('--' + name)
    p.add_argument('--met-input', choices=['raw', 'named'])
    for name, default in [('id-col', 'eid'), ('baseline-col', 'date_attend'),
                          ('death-col', 'date_death'), ('lost-col', 'date_lost'),
                          ('disease-evidence-col', ''), ('healthy-date-cols', ''),
                          ('exclude-features', '')]:
        p.add_argument('--' + name, default=default)
    p.add_argument('--transform', choices=['none', 'log1p'])
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--r-bin', default='Rscript')
    p.add_argument('--max-samples', type=int, default=0)
    p.add_argument('--demo-features', type=int, default=80)
    p.add_argument('--demo', action='store_true')
    return p


def configure(a):
    common = runpy.run_path(str(HERE / '0.common.py'), run_name='le8_config')
    settings, sources, _ = common['shared_analysis_settings'](a)
    for name, value in settings.items():
        setattr(a, name, value)
    a.shared_config_sources = sources
    a.group_source_column = a.group_col
    if a.group_file:
        a.group_col = '.le8_family'
    if a.demo and a.outer_roster:
        raise ValueError('Demo cannot use a real outer roster')
    a.covariates = a.covariates if a.covariates is not None else (a.shared_covariates or 'age,sex,tdi,PC1,PC2')
    a.outcome_type, a.target_col = 'survival', a.trait
    a.diagnosis_col = a.diagnosis_col or 'fod_icd10_' + a.trait
    base = Path(a.ukb_phe)
    a.phe_file = a.phe_file or str(base / 'Rdata/all.rds')
    # Use the same cleaned, named matrices as native C1. Raw assays are an
    # explicit alternative, never an implicit second cohort/preprocessing path.
    a.omics_file = a.omics_file or str(base / f'Rdata/{a.biom}.rds')
    a.met_map = a.met_map or str(base / 'common/met.lst')
    a.met_input = a.met_input or 'named'
    a.residualize = '' if a.residualize is None else a.residualize
    a.categorical = 'sex' if a.categorical is None else a.categorical
    a.excluded_technical_covariates = []
    for name in ('covariates', 'residualize', 'categorical'):
        kept = []
        for col in words(getattr(a, name)):
            if technical_category(col):
                if col not in a.excluded_technical_covariates:
                    a.excluded_technical_covariates.append(col)
            else:
                kept.append(col)
        setattr(a, name, ','.join(kept))
    if a.excluded_technical_covariates:
        log('SKIP', 'technical_covariates', ','.join(a.excluded_technical_covariates))
    a.transform = a.transform or ('log1p' if a.biom == 'met' and a.met_input == 'raw' and not a.demo else 'none')
    if a.max_samples < 0 or a.demo_features < 3:
        raise ValueError('Invalid sample limit or demo feature count')
    return a


def assign_roles(p, a):
    """Deterministic, outcome-blind family assignment for the current protocol."""
    from abm_twin_core import ROLES, tie_hash
    groups = family_groups(p, a)
    outer = shared_outer_roles(p, a)
    values = np.array([tie_hash(g, a.seed) / 2**64 for g in groups])
    edges = np.array([0, .50, .60, .68, .76, .84])
    if outer is not None:
        edges = edges / edges[-1]
    role = np.full(len(p), 'test', dtype='U25')
    for j, name in enumerate(ROLES):
        selected = (values >= edges[j]) & (values < edges[j + 1])
        if outer is not None:
            selected &= outer == 'training'
        role[selected] = name
    if set(role) != set(ROLES) | {'test'}:
        raise ValueError('Each development role and test must have eligible participants')
    return role


RDS_EXPORT_SCRIPT = '# Read-only bridge. No sourcing of phenotype pipelines or installing R packages.\nargs <- commandArgs(trailingOnly = TRUE)\nstopifnot(length(args) == 3L)\nx <- readRDS(args[1])\nif (!is.data.frame(x)) stop("Expected one data.frame in RDS")\ncols <- readLines(args[3], warn = FALSE)\nif (length(cols)) {\n  absent <- setdiff(cols, names(x))\n  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))\n  x <- x[, cols, drop = FALSE]\n}\nfor (nm in names(x)) {\n  if (inherits(x[[nm]], "Date") || inherits(x[[nm]], "POSIXt"))\n    x[[nm]] <- format(x[[nm]], "%Y-%m-%d")\n  if (inherits(x[[nm]], "integer64")) x[[nm]] <- as.character(x[[nm]])\n}\n# Binary transfer preserves subnormal numeric metadata used as category codes.\nsaveRDS(x, args[2], compress = FALSE)\n\n'

def words(value):
	return [s.strip() for s in (value or "").split(",") if s.strip()]

def technical_category(name):
	"""Technical location/batch categories excluded from this ABM specification."""
	return bool(set(re.split(r'[.\s_:-]+', name.casefold())) & {'plate', 'well', 'center', 'centre'})

def json_safe(value):
	if isinstance(value, dict):
		return {str(k): json_safe(v) for k, v in value.items()}
	if isinstance(value, (list, tuple, np.ndarray)):
		return [json_safe(v) for v in value]
	if isinstance(value, (np.integer, np.bool_)):
		return value.item()
	if isinstance(value, (float, np.floating)):
		return float(value) if np.isfinite(value) else None
	if isinstance(value, Path):
		return str(value)
	return value

def dump(path, value):
	path = Path(path)
	temp = path.with_suffix(path.suffix + ".tmp")
	temp.write_text(
		json.dumps(json_safe(value), ensure_ascii=False, indent=2, allow_nan=False)
		+ "\n",
		encoding="utf-8",
	)
	temp.replace(path)

def log(action, step, detail=""):
	print(f"[ABM] {action} {step}" + (f" | {detail}" if detail else ""), flush=True)

def parsed_date(raw, name):
	if pd.api.types.is_numeric_dtype(raw) and raw.notna().any():
		raise ValueError(
			f"{name}: numeric dates are ambiguous; use ISO dates or R Date"
		)
	converted = pd.to_datetime(raw, errors="coerce")
	bad = raw.notna() & raw.astype(str).str.strip().ne("") & converted.isna()
	if bad.any():
		raise ValueError(f"Invalid dates in {name}")
	return converted

def outcomes(p, a):
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
	baseline = parsed_date(p[a.baseline_col], a.baseline_col)
	diagnosis = parsed_date(p[a.diagnosis_col], a.diagnosis_col)
	censor = pd.concat(
		[
			parsed_date(p[a.death_col], a.death_col),
			parsed_date(p[a.lost_col], a.lost_col),
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
	for col in words(a.healthy_date_cols):
		d = parsed_date(p[col], col)
		healthy &= ~(d.notna() & (d <= baseline))
	event = (
		valid
		& ~prevalent
		& diagnosis.notna()
		& (diagnosis > baseline)
		& (diagnosis <= censor)
	)
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

def validate_ids(frame, id_col):
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
			raise ValueError(
				"Numeric IDs must be finite integers; supply strings otherwise"
			)
		frame[id_col] = raw.astype("int64").astype(str)
	else:
		frame[id_col] = raw.astype(str).str.strip()
	if frame[id_col].eq("").any() or frame[id_col].duplicated().any():
		raise ValueError(
			"Empty or duplicate person IDs; select one baseline visit first"
		)
	return frame

def read_table(path, id_col="eid", columns=None, r_bin="Rscript", string_columns=()):
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
				bridge.write_text(RDS_EXPORT_SCRIPT, encoding="utf-8")
				subprocess.run(
					[r_bin, str(bridge), str(path), str(out), str(cols)], check=True
				)
				frame = pyreadr.read_r(str(out))[None].reset_index(drop=True)
		else:
			objects = pyreadr.read_r(str(path))
			if len(objects) != 1 or not isinstance(
				next(iter(objects.values())), pd.DataFrame
			):
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
			path, sep=sep, dtype={c: str for c in (id_col, *string_columns)}, usecols=columns, low_memory=False
		)
	if columns:
		missing = set(columns) - set(frame)
		if missing:
			raise ValueError(f"{path}: missing columns {sorted(missing)}")
		frame = frame[columns].copy()
	return validate_ids(frame, id_col)

def numeric(frame, columns, context, dtype="float32"):
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

def molecular_input_qc(raw, biom, demo=False):
	"""Treat negative metabolite concentrations as missing; preserve signed NPX."""
	if biom == "met" and not demo:
		return np.where(raw < 0, np.nan, raw)
	return raw

def safe_expression(expression, frame):
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

def map_metabolites(frame, mapping, id_col):
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
				else safe_expression(expr, frame)
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
		raise ValueError(
			"Fewer than three metabolites map; check raw baseline fields and met.lst"
		)
	return result, pd.DataFrame(audit)

def phenotype_columns(a):
	cols = [a.id_col] + words(a.covariates) + words(a.residualize)
	if a.outcome_type == "survival":
		cols += [a.baseline_col, a.diagnosis_col, a.death_col, a.lost_col]
		cols += words(a.healthy_date_cols)
		if a.disease_evidence_col:
			cols += [a.disease_evidence_col]
	else:
		cols += [a.target_col]
	if a.group_source_column:
		cols += [a.group_source_column]
	return list(dict.fromkeys(cols))

def generate_demo(a):
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
			a.lost_col: (base + pd.to_timedelta(censor * 365.25, unit="D")).strftime(
				"%Y-%m-%d"
			),
			a.diagnosis_col: (
				base + pd.to_timedelta(np.minimum(times, 100) * 365.25, unit="D")
			).strftime("%Y-%m-%d"),
		}
	)
	if a.outcome_type == "quantitative":
		p[a.target_col] = (
			165
			+ 7 * sex
			+ 2.5 * latent[:, 1]
			+ 2 * (state - 1) * latent[:, 2]
			+ rng.normal(size=n)
		)
	if a.group_col:
		p[a.group_source_column or a.group_col] = (np.arange(n) // 2).astype(str)
	features = [f"F{j:04d}" for j in range(nf)]
	omics = pd.DataFrame(x, columns=features)
	omics.insert(0, a.id_col, p[a.id_col])
	return p, omics

def attach_shared_groups(p, a):
	p = validate_ids(p.copy(), a.id_col)
	if a.group_file:
		m = read_table(a.group_file, "eid", r_bin=a.r_bin, string_columns=("group",))
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

def shared_outer_roles(p, a):
	if not getattr(a, "outer_roster", ""):
		return None
	m = read_table(a.outer_roster, "eid", r_bin=a.r_bin)
	if "role" not in m or not set(m.role) <= {"training", "test"}:
		raise ValueError("Shared outer roster requires eid,role with training/test")
	role = p[a.id_col].map(m.set_index("eid").role)
	if role.isna().any() or set(role) != {"training", "test"}:
		raise ValueError("Outer roster must cover input participants and both roles")
	groups = family_groups(p, a)
	if pd.DataFrame(dict(group=groups, role=role.to_numpy())).groupby("group").role.nunique().max() > 1:
		raise ValueError("Families cross shared outer roles")
	return role.to_numpy()

def prepare(a, out):
	mapping_audit = None
	if a.demo:
		p, omics = generate_demo(a)
	else:
		p = read_table(a.phe_file, a.id_col, phenotype_columns(a), a.r_bin, string_columns=([a.group_source_column] if a.group_source_column else []))
		omics = read_table(a.omics_file, a.id_col, r_bin=a.r_bin)
		if a.biom == "met" and a.met_input == "raw":
			omics, mapping_audit = map_metabolites(omics, a.met_map, a.id_col)
		if a.biom == "prot":
			omics.columns = [c if c == a.id_col else str(c).upper() for c in omics]
			if omics.columns.duplicated().any():
				raise ValueError("Protein names collide after uppercasing")
	missing = set(phenotype_columns(a)) - set(p)
	if missing:
		raise ValueError(f"Missing phenotype columns: {sorted(missing)}")
	omics = omics[omics[a.id_col].isin(p[a.id_col])].copy()
	if a.max_samples and len(omics) > a.max_samples:
		omics = omics.sample(a.max_samples, random_state=a.seed)
	if omics.empty:
		raise ValueError("No matched phenotype/omics IDs")
	p = p.set_index(a.id_col).loc[omics[a.id_col]].reset_index()
	p = attach_shared_groups(p, a)
	features = [
		c for c in omics if c != a.id_col and c not in words(a.exclude_features) and not technical_category(c)
	]
	forbidden = set(phenotype_columns(a)) | {"event", "time", "target", "split", "Y"}
	collision = {c.casefold() for c in features} & {c.casefold() for c in forbidden}
	if collision:
		raise ValueError(
			"Metadata/outcome columns in omics matrix: " + ", ".join(sorted(collision))
		)
	x = numeric(omics, features, "omics", dtype="float64")
	negative = np.sum(x < 0, axis=0)
	# Named/preprocessed matrices keep their supplied scale, including signed
	# values. Concentration-domain QC belongs only to explicitly raw metabolites.
	raw_metabolites = a.biom == "met" and a.met_input == "raw" and not a.demo
	if raw_metabolites:
		x = molecular_input_qc(x, a.biom)
	negative_to_missing = int(negative.sum()) if raw_metabolites else 0
	domain_invalid = np.sum(x <= -1, axis=0)
	if a.transform == "log1p" and domain_invalid.any():
		bad = [f"{name} ({int(n)})" for name, n in zip(features, domain_invalid) if n]
		raise ValueError("log1p requires x > -1; invalid input features: " + ", ".join(bad))
	transform_audit = dict(transform=a.transform, negative_values=int(negative.sum()),
		negative_features={name:int(n) for name,n in zip(features,negative) if n},
		negative_to_missing_values=negative_to_missing,
		log1p_domain_invalid_values=int(domain_invalid.sum()))
	log("DONE", "input_domain", f"transform={a.transform}; negative_values={int(negative.sum())}; "
		f"negative_features={int(np.count_nonzero(negative))}; negative_to_missing={negative_to_missing}; "
		f"log1p_domain_invalid={int(domain_invalid.sum())}")
	# Preserve missing observations and invalid concentrations for train-only imputation.
	p.to_csv(out / "phenotype.csv", index=False)
	np.save(out / "raw.npy", x)
	(out / "features.txt").write_text("\n".join(features) + "\n")
	if mapping_audit is not None:
		mapping_audit.to_csv(out / "metabolite_mapping.csv", index=False)
	upstream = "supplied matrix; upstream preprocessing not inferred"
	if str(a.omics_file).endswith(f"Rdata/{a.biom}.rds"):
		upstream = "LE8 cleaned RDS: upstream all-person QC/imputation may already have occurred"
	dump(
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

def family_groups(p, a):
    values = p[a.group_col if a.group_col else a.id_col]
    if values.isna().any() or values.astype(str).str.strip().eq("").any():
        raise ValueError("Missing participant/family identifiers")
    return values.astype(str).to_numpy()
