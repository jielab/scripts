"""Read LE8-layout inputs without executing the LE8 analysis environment."""
from pathlib import Path
import ast
import re
import shutil
import subprocess
import tempfile
import numpy as np
import pandas as pd
from common import words, dump

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
            raise ValueError("Numeric IDs must be finite integers; supply strings otherwise")
        frame[id_col] = raw.astype("int64").astype(str)
    else:
        frame[id_col] = raw.astype(str).str.strip()
    if frame[id_col].eq("").any() or frame[id_col].duplicated().any():
        raise ValueError("Empty or duplicate person IDs; select one baseline visit first")
    return frame

def read_table(path, id_col="eid", columns=None, r_bin="Rscript"):
    path = Path(path)
    if path.suffix.lower() == ".rds":
        if shutil.which(r_bin):
            # R selects columns before serialization, avoiding a second full RDS in Python.
            with tempfile.TemporaryDirectory(prefix="panome_rds_") as tmp:
                out = Path(tmp) / "selected.csv"
                cols = Path(tmp) / "columns.txt"
                cols.write_text("\n".join(columns or []))
                subprocess.run([r_bin, str(Path(__file__).with_name("export_rds.R")),
                                str(path), str(out), str(cols)], check=True)
                frame = pd.read_csv(out, dtype={id_col: str}, low_memory=False)
        else:
            import pyreadr
            objects = pyreadr.read_r(str(path))
            if len(objects) != 1 or not isinstance(next(iter(objects.values())), pd.DataFrame):
                raise ValueError(f"{path}: expected a single R data.frame")
            frame = next(iter(objects.values()))
    elif path.suffix.lower() == ".parquet":
        frame = pd.read_parquet(path, columns=columns)
    else:
        sep = "," if ".csv" in path.name.lower() else "\t"
        # Reject duplicate raw headers before pandas can silently suffix them.
        header = pd.read_csv(path, sep=sep, header=None, nrows=1).iloc[0].astype(str)
        if header.duplicated().any():
            raise ValueError(f"Duplicate header in {path}")
        frame = pd.read_csv(path, sep=sep, dtype={id_col: str},
                            usecols=columns, low_memory=False)
    if columns:
        missing = set(columns) - set(frame)
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        frame = frame[columns].copy()
    return validate_ids(frame, id_col)

def numeric(frame, columns, context):
    values = []
    for name in columns:
        raw = frame[name]
        val = pd.to_numeric(raw, errors="coerce")
        bad = raw.notna() & val.isna()
        if bad.any():
            raise ValueError(f"{context}: nonnumeric values in {name}")
        values.append(val.to_numpy(dtype="float32", na_value=np.nan))
    if not values:
        raise ValueError(f"{context}: no molecular columns")
    result = np.column_stack(values)
    result[~np.isfinite(result)] = np.nan
    return result

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
    spec = pd.read_csv(mapping, sep="\t", header=None, usecols=[0, 1],
                       names=["expression", "feature"], dtype=str)
    # Support both the original headerless mapping and the current UKB mapping.
    if len(spec) and tuple(spec.iloc[0]) == ("data_field", "met_name"):
        spec = spec.iloc[1:].copy()
    if spec.isna().any().any() or spec.feature.duplicated().any():
        raise ValueError("met.lst must contain unique, nonempty feature names")
    result = pd.DataFrame({id_col: frame[id_col]})
    audit = []
    for row in spec.itertuples():
        expr, feature = row.expression.strip(), row.feature.strip()
        try:
            values = (pd.to_numeric(frame[expr], errors="raise").to_numpy(float)
                      if expr in frame else safe_expression(expr, frame))
            values = np.broadcast_to(values, (len(frame),)).copy()
            invalid = ~np.isfinite(values)
            values[invalid] = np.nan
            result[feature] = values.astype("float32")
            audit.append(dict(feature=feature, expression=expr, status="mapped",
                              nonfinite_to_missing=int(invalid.sum())))
        except KeyError as exc:
            audit.append(dict(feature=feature, expression=expr, status="missing_source",
                              reason=str(exc)))
    if result.shape[1] < 4:
        raise ValueError("Fewer than three metabolites map; check raw baseline fields and met.lst")
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
    if a.group_col:
        cols += [a.group_col]
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
    x += plate[:, None] * .06
    x[rng.random(x.shape) < .03] = np.nan
    risk = .03 * (age - 55) + .35 * latent[:, 1] + .4 * (state - 1) * latent[:, 2]
    times = rng.exponential(18 * np.exp(-risk))
    censor = rng.uniform(7, 12, n)
    base = pd.Timestamp("2010-01-01")
    p = pd.DataFrame({a.id_col: [f"demo_{i}" for i in range(n)], "age": age, "sex": sex,
                      "tdi": rng.normal(size=n), "PC1": rng.normal(size=n),
                      "PC2": rng.normal(size=n), "center": rng.choice(["A", "B", "C"], n),
                      f"{a.biom}.plate": plate,
                      a.baseline_col: "2010-01-01",
                      a.death_col: pd.Series([None] * n, dtype=object),
                      a.lost_col: (base + pd.to_timedelta(censor * 365.25, unit="D")).strftime("%Y-%m-%d"),
                      a.diagnosis_col: (base + pd.to_timedelta(np.minimum(times, 100) * 365.25, unit="D")).strftime("%Y-%m-%d")})
    if a.outcome_type == "quantitative":
        p[a.target_col] = 165 + 7 * sex + 2.5 * latent[:, 1] + 2 * (state - 1) * latent[:, 2] + rng.normal(size=n)
    if a.group_col:
        p[a.group_col] = np.arange(n) // 2
    features = [f"F{j:04d}" for j in range(nf)]
    omics = pd.DataFrame(x, columns=features)
    omics.insert(0, a.id_col, p[a.id_col])
    return p, omics

def prepare(a, out):
    mapping_audit = None
    if a.demo:
        p, omics = generate_demo(a)
    else:
        p = read_table(a.phe_file, a.id_col, phenotype_columns(a), a.r_bin)
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
    features = [c for c in omics if c != a.id_col and c not in words(a.exclude_features)]
    forbidden = set(phenotype_columns(a)) | {"event", "time", "target", "split", "Y"}
    collision = {c.casefold() for c in features} & {c.casefold() for c in forbidden}
    if collision:
        raise ValueError("Metadata/outcome columns in omics matrix: " + ", ".join(sorted(collision)))
    x = numeric(omics, features, "omics")
    # Retain original missingness. Scale transforms are part of train-fitted preprocessing.
    p.to_csv(out / "phenotype.csv", index=False)
    np.save(out / "raw.npy", x)
    (out / "features.txt").write_text("\n".join(features) + "\n")
    if mapping_audit is not None:
        mapping_audit.to_csv(out / "metabolite_mapping.csv", index=False)
    upstream = "supplied matrix; upstream preprocessing not inferred"
    if str(a.omics_file).endswith(f"Rdata/{a.biom}.rds"):
        upstream = "LE8 cleaned RDS: upstream all-person QC/imputation may already have occurred"
    dump(out / "input_audit.json", dict(n=len(p), features=len(features),
         source="SYNTHETIC" if a.demo else upstream, met_input=a.met_input,
         missing_fraction=float(np.mean(~np.isfinite(x)))))
