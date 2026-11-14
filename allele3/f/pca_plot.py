#!/usr/bin/env python3
"""Consolidated pca_plot workflow; internal stages share one explicit module."""
from __future__ import annotations


# ========================================================================
# pca atlas
# ========================================================================
_description_pca_atlas = 'Re-render fixed PC1/PC2 coordinates with an explicit chronological third axis.\n\nNo PCA or genotype inference is performed. The source may be an existing\nallele3 HTML data payload or a shared coordinate table plus sample metadata.\nAll original points and dates remain in the downloadable audit, including\nancient people outside the primary display interval.'
#!/usr/bin/env python3
"""Re-render fixed PC1/PC2 coordinates with an explicit chronological third axis.

No PCA or genotype inference is performed. The source may be an existing
allele3 HTML data payload or a shared coordinate table plus sample metadata.
All original points and dates remain in the downloadable audit, including
ancient people outside the primary display interval.
"""

import argparse
from collections import Counter
from contextlib import contextmanager
import copy
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import html
import json
import math
from pathlib import Path
import re
import sys

VERSION = "pca-time-atlas-2.0"
EXPLORER_TEMPLATE = Path(__file__).with_name("pca_explorer.html")
EXPLORER_SCRIPT = Path(__file__).with_name("pca_explorer.js")
REGION_COLORS = {"C": "#2166AC", "E": "#D47728", "N": "#188A72", "SE": "#A65F98", "SW": "#B99B18"}
REGION_NAMES = {"C": "Central", "E": "East", "N": "North", "SE": "Southeast", "SW": "Southwest"}
LOG_SCALE_BP = 1000.0


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path):
    path = Path(path).expanduser().resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha256(path))


def read_source_html(path):
    text = Path(path).read_text(encoding="utf-8")
    found = re.findall(r'<script\b(?=[^>]*\bid=["\']allele3-data["\'])[^>]*>(.*?)</script\s*>', text, re.S | re.I)
    if len(found) != 1:
        raise ValueError("Source HTML must contain exactly one allele3-data JSON script; no JavaScript is executed")
    payload = json.loads(found[0])
    return payload, dict(source_type="existing_html_embedded_coordinates", inputs={"html": file_identity(path)},
                         plan_validation_status="embedded_provenance_unverified_plan",
                         action="lossless re-rendering of existing PC1/PC2; no new PCA fit")


def read_coordinate_source(coordinates, metadata, basis_label=None, *, simulated=False):
    from evolve import _coordinates, _pca_metadata
    coords, kind = _coordinates(coordinates)
    meta = _pca_metadata(metadata)
    included = [sid for sid, row in meta.items() if row["include"] == "1"]
    if set(included) - set(coords) or set(coords) - set(meta):
        raise ValueError("Coordinate/metadata identity mismatch: all included samples need coordinates; unknown IDs are rejected")
    points = []
    for sid in included:
        row = meta[sid]
        points.append(dict(sample_id=sid, pc1=coords[sid][0], pc2=coords[sid][1], age_bp=row["age"],
                           cohort=row["cohort"], region=row.get("region", row["stratum"]),
                           source=row.get("data_source", row["cohort"]), data_types=row.get("data_types", "."),
                           coverage=row.get("target_coverage", "."), iqs=row.get("iqs", ".")))
    payload = dict(points=points, coordinate_kind=kind, basis_label=basis_label or "supplied shared coordinates", simulated=bool(simulated))
    provenance = dict(source_type="coordinate_table_plus_metadata", inputs={"coordinates": file_identity(coordinates),
                      "metadata": file_identity(metadata)}, source_excluded_metadata_samples=len(meta) - len(included),
                      plan_validation_status="coordinate_table_provenance_unverified_plan",
                      action="display existing PC1/PC2; no new PCA fit")
    return payload, provenance


def validate_points(payload):
    if not isinstance(payload, dict) or not isinstance(payload.get("points"), list) or not payload["points"]:
        raise ValueError("A nonempty points list is required")
    seen = set()
    for p in payload["points"]:
        sid = p.get("sample_id")
        if not isinstance(sid, str) or not sid.strip() or sid in seen:
            raise ValueError("Sample IDs must be nonempty and unique")
        seen.add(sid)
        for key in ("pc1", "pc2", "age_bp"):
            value = p.get(key)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{sid}: {key} must be finite numeric data")
        if p["age_bp"] < 0 or p.get("cohort") not in {"ancient", "modern"}:
            raise ValueError(f"{sid}: require nonnegative BP and explicit ancient/modern cohort")
        if not isinstance(p.get("region", "unknown"), str):
            raise ValueError(f"{sid}: region must be a string label")


def validate_overlays(payload):
    """Retain supplied measurements; never infer them from frequency or PCA."""
    n = len(payload["points"])
    variants = payload.get("variants", [])
    features = payload.get("features", [])
    if not isinstance(variants, list) or not isinstance(features, list):
        raise ValueError("variants and features must be explicit lists")
    ids = set()
    for variant in variants:
        if not isinstance(variant, dict):
            raise ValueError("Every variant overlay must be an object")
        vid = variant.get("variant_id")
        if not isinstance(vid, str) or not vid.strip() or vid in ids:
            raise ValueError("Variant overlay IDs must be nonempty and unique")
        ids.add(vid)
        if variant.get("ref") not in tuple("ACGT") or variant.get("alt") not in tuple("ACGT") or variant["ref"] == variant["alt"]:
            raise ValueError("Variant overlays require distinct single canonical REF/ALT bases")
        if variant.get("count_kind") not in {"hardcall", "dosage"}:
            raise ValueError("Variant overlays must declare hardcall or dosage")
        if "values" not in variant:
            raise ValueError("Variant overlay values must be explicit: aligned values or null when unavailable")
        values = variant.get("values")
        if values is not None and (not isinstance(values, list) or len(values) != n):
            raise ValueError("Variant overlay values must align to every original point")
        for value in values or []:
            if value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 2:
                raise ValueError("Variant ALT dosage must be finite in [0,2], or null")
            if variant["count_kind"] == "hardcall" and value not in {0, 1, 2}:
                raise ValueError("Hardcall overlays accept only 0, 1, 2 or null")
    ids = set()
    for feature in features:
        if not isinstance(feature, dict):
            raise ValueError("Every sample feature must be an object")
        fid = feature.get("feature_id", "")
        if not isinstance(fid, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", fid) or fid in ids or fid == "sample_id":
            raise ValueError("Feature IDs must be unique simple names")
        ids.add(fid)
        if feature.get("kind") not in {"numeric", "categorical"}:
            raise ValueError("Sample feature kind must be numeric or categorical")
        for key in ("label", "units", "provenance"):
            if not isinstance(feature.get(key), str) or not feature[key].strip():
                raise ValueError("Sample feature " + key + " must be supplied")
        values = feature.get("values")
        if not isinstance(values, list) or len(values) != n:
            raise ValueError("Sample feature values must align to every original point")
        for value in values:
            if value is None:
                continue
            if feature["kind"] == "numeric":
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                    raise ValueError("Numeric sample features require finite numbers or null")
            elif not isinstance(value, str) or not value.strip():
                raise ValueError("Categorical sample features require strings or null")


def read_sample_features(path, schema_path, sample_ids):
    """Exact-ID TSV join; missing sample observations remain null.

    Schema: feature_id, label, kind (numeric/categorical), units, provenance.
    Data: sample_id plus exactly the declared feature columns. Numeric zero is
    observed zero; empty, '.', NA and null are unavailable, never imputed.
    """
    required = {"feature_id", "label", "kind", "units", "provenance"}
    def rows(source):
        opener = gzip.open if str(source).endswith(".gz") else open
        with opener(source, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError("Feature TSV requires nonempty unique headers")
            header = reader.fieldnames
            result = []
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError("Malformed feature TSV row")
                result.append(row)
        return header, result
    header, records = rows(schema_path)
    if not required <= set(header):
        raise ValueError("Feature schema requires: " + ",".join(sorted(required)))
    features = []
    for row in records:
        record = {k: row[k].strip() for k in required}
        record["values"] = [None] * len(sample_ids)
        features.append(record)
    if not features:
        raise ValueError("Feature schema is empty")
    validate_overlays({"points": [{}] * len(sample_ids), "features": features})
    if len(features) > 100:
        raise ValueError("At most 100 explicit sample feature tracks are supported")
    header, records = rows(path)
    if set(header) != {"sample_id", *(f["feature_id"] for f in features)}:
        raise ValueError("Sample feature columns must equal sample_id plus the declared schema")
    indexes = {sid: i for i, sid in enumerate(sample_ids)}
    if len(indexes) != len(sample_ids):
        raise ValueError("Source sample IDs must be unique")
    seen = set()
    missing = {"", ".", "NA", "NAN", "NULL", "NONE"}
    for row in records:
        sid = row["sample_id"]
        if sid not in indexes or sid in seen:
            raise ValueError("Feature sample IDs must be unique and present in the coordinate source: " + sid)
        seen.add(sid)
        i = indexes[sid]
        for feature in features:
            value = row[feature["feature_id"]].strip()
            if value.upper() in missing:
                continue
            if feature["kind"] == "numeric":
                try:
                    value = float(value)
                except ValueError as exc:
                    raise ValueError("Sample numeric feature is not numeric: " + feature["feature_id"]) from exc
                if not math.isfinite(value):
                    raise ValueError("Sample numeric feature must be finite")
            feature["values"][i] = value
    return features, dict(supplied_rows=len(seen), source_samples=len(sample_ids),
                         absent_sample_rows=len(sample_ids) - len(seen),
                         tracks=[dict(feature_id=f["feature_id"], kind=f["kind"],
                                      observed=sum(v is not None for v in f["values"]),
                                      unavailable=sum(v is None for v in f["values"])) for f in features],
                         interpretation="Supplied sample-level observations; provenance is caller-declared, not independently verified")


def time_transform(age_bp, scale="linear"):
    """Return display units only. Raw BP ages are never modified."""
    if age_bp < 0 or not math.isfinite(age_bp):
        raise ValueError("Finite nonnegative time required")
    if scale == "linear":
        return age_bp / 1000.0
    if scale == "log1p":
        return math.log1p(age_bp / LOG_SCALE_BP) / math.log(10.0)
    raise ValueError("time-scale must be linear or log1p")


def layer_index(point, layer_years):
    return -1 if point["cohort"] == "modern" else int(point["age_bp"] // layer_years)


def display_age(point, layer_years, view="layered"):
    if view == "continuous":
        return point["age_bp"]
    if point["cohort"] == "modern":
        return 0.0
    if view == "layered":
        return (layer_index(point, layer_years) + .5) * layer_years
    raise ValueError("time-view must be layered or continuous")


def visible_points(points, max_age_bp, scope="primary"):
    if scope not in {"primary", "all"}:
        raise ValueError("age scope must be primary or all")
    return [p for p in points if p["cohort"] == "modern" or scope == "all" or p["age_bp"] <= max_age_bp]


def shared_limits(points):
    """Identical numeric span and physical scale for both genetic axes."""
    x, y = [p["pc1"] for p in points], [p["pc2"] for p in points]
    span = max(max(x) - min(x), max(y) - min(y), 1e-8) * 1.12
    cx, cy = (min(x) + max(x)) / 2, (min(y) + max(y)) / 2
    return dict(x=[cx - span / 2, cx + span / 2], y=[cy - span / 2, cy + span / 2], shared_span=span)


def palette(points):
    import colorsys
    colors = dict(REGION_COLORS)
    for i, region in enumerate(sorted({p.get("region", "unknown") for p in points} - set(colors))):
        rgb = colorsys.hsv_to_rgb((.13 + i * .381966) % 1, .55, .65)
        colors[region] = "#" + "".join(f"{round(v * 255):02x}" for v in rgb)
    return colors


def layer_rows(points, layer_years, max_age_bp):
    records = []
    for scope in ("primary", "all"):
        counts = Counter(layer_index(p, layer_years) for p in visible_points(points, max_age_bp, scope))
        for index, n in sorted(counts.items()):
            start, end = (0, 0) if index == -1 else (index * layer_years, (index + 1) * layer_years)
            records.append(dict(scope=scope, layer_index=index, start_bp=start, end_bp=end,
                                display_bp=0 if index == -1 else (start + end) / 2, n=n,
                                cohort="modern" if index == -1 else "ancient",
                                label="Modern" if index == -1 else f"{start / 1000:g}–{end / 1000:g} ka"))
    return records


def write_tsv(path, fields, records):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as stream:
        output = csv.DictWriter(stream, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        output.writeheader()
        for record in records:
            output.writerow({k: record.get(k, ".") for k in fields})


def draw_static(data, out, scale, scope, basename):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    points = visible_points(data["points"], data["display"]["max_age_bp"], scope)
    width = data["display"]["layer_years"]
    layers = [r for r in data["layers"] if r["scope"] == scope]
    limits, colors = data["limits"], data["region_colors"]
    upper = max(width, math.ceil(max([p["age_bp"] for p in points] + [0]) / width) * width)
    if scope == "primary":
        upper = max(upper, data["display"]["max_age_bp"])
    ztop = time_transform(upper, scale)
    x, y, z = ([p["pc1"] for p in points], [p["pc2"] for p in points],
               [time_transform(display_age(p, width, "continuous"), scale) for p in points])
    ancient_n = sum(p["cohort"] == "ancient" for p in points)
    modern_n = len(points) - ancient_n
    fig = plt.figure(figsize=(12.0, 11.4), facecolor="white")
    ax = fig.add_axes([.02, .13, .73, .77], projection="3d")
    ax.set_proj_type("ortho")
    ax.set_box_aspect((1, 1, 1.60))
    ax.view_init(elev=19, azim=-54)
    ax.set(xlim=limits["x"], ylim=limits["y"], zlim=(-.025 * ztop, 1.035 * ztop),
           xlabel="PC1 · shared score", ylabel="PC2 · shared score", zlabel="Age (ka BP)")
    ax.xaxis.labelpad = 9; ax.yaxis.labelpad = 9; ax.zaxis.labelpad = 7
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor((1, 1, 1, 0))
        axis._axinfo["grid"]["color"] = (.77, .81, .84, .22)
        axis._axinfo["grid"]["linewidth"] = .4
    ax.tick_params(axis="both", labelsize=8, pad=1)
    lo_x, hi_x = limits["x"]; lo_y, hi_y = limits["y"]
    # A sparse frame shows each plane without an opaque sheet in front of data.
    for reference_bp in range(0, int(upper) + 1, width):
        zz = time_transform(reference_bp, scale)
        ax.plot([lo_x, hi_x, hi_x], [lo_y, lo_y, hi_y], [zz] * 3,
                color="#9CADB7", alpha=.42 if reference_bp == 0 else .20, lw=.65)
    for cohort, marker, size, alpha in (("ancient", "o", 2.8, .52), ("modern", "D", 5.5, .90)):
        idx = [i for i, p in enumerate(points) if p["cohort"] == cohort]
        if idx:
            ax.scatter(np.asarray(x)[idx], np.asarray(y)[idx], np.asarray(z)[idx],
                       c=[colors[points[i].get("region", "unknown")] for i in idx], s=size, alpha=alpha,
                       marker=marker, edgecolors="none", depthshade=False, rasterized=True)
    if scale == "linear":
        ticks = list(range(0, int(upper) + 1, 2000 if upper <= 12000 else 4000))
    else:
        ticks = [v for v in (0, 500, 1000, 2000, 4000, 6000, 10000, 14000, 18000, 20000) if v <= upper]
    ax.set_zticks([time_transform(v, scale) for v in ticks], [f"{v / 1000:g}" for v in ticks])
    fig.text(.07, .958, "Genetic structure through time", fontsize=22, weight="semibold", color="#173B51")
    interval = f"0–{data['display']['max_age_bp'] / 1000:g} ka BP" if scope == "primary" else "All recorded ages"
    fig.text(.07, .924, f"{interval}  ·  {ancient_n:,} ancient people + {modern_n:,} modern references", fontsize=12, color="#4A6271")
    side = fig.add_axes([.79, .23, .19, .53])
    side.set_axis_off()
    side.text(0, 1.04, "People per age interval", fontsize=10.5, weight="semibold", color="#173B51", transform=side.transAxes)
    ordered = list(reversed(layers))
    gap = min(.061, .93 / max(1, len(ordered)))
    for i, layer in enumerate(ordered):
        yy = .98 - i * gap
        side.text(0, yy, layer["label"], fontsize=8.2, color="#526775", transform=side.transAxes)
        side.text(.99, yy, f"{layer['n']:,}", ha="right", fontsize=8.5, color="#183F55", transform=side.transAxes)
        side.plot([0, .96 * layer["n"] / max(r["n"] for r in layers)],
                  [yy - .015] * 2, color="#AFC1CD", lw=1.2, transform=side.transAxes)
    legend = [Line2D([0], [0], marker="o", lw=0, color=color, markersize=5,
                     label=f"{region} · {REGION_NAMES[region]}" if region in REGION_NAMES else region)
              for region, color in colors.items() if any(p.get("region", "unknown") == region for p in points)]
    fig.legend(handles=legend, loc="upper left", bbox_to_anchor=(.79, .91), frameon=False, fontsize=9,
               title="Broad geographic region", title_fontsize=9.5)
    footer = (f"PC1 and PC2 use one fixed coordinate system and equal axis scales. Third axis = time, not PC3.\n"
              f"Points retain exact recorded dates; {width:,}-year frames are guides. Circles: ancient; diamonds: modern.\n")
    if scope == "primary":
        footer += f"{data['counts']['ancient_older_than_primary']:,} older ancient people are retained in the all-ages supplement and audit. "
    footer += "Linear time." if scale == "linear" else "Nonlinear display: log10(1 + age/1,000 BP); ticks show original ages."
    fig.text(.07, .055, footer, fontsize=9, color="#526775", linespacing=1.5)
    if data.get("simulated"):
        fig.text(.07, .014, "SIMULATED INPUT — not measured ancient coordinates", color="#A63C27", weight="bold")
    paths = []
    for ext in ("png", "pdf"):
        path = Path(out) / f"{basename}.{ext}"
        fig.savefig(path, dpi=200, facecolor="white", bbox_inches="tight", pad_inches=.15)
        paths.append(path.name)
    plt.close(fig)
    return paths


def explorer_template():
    """Load the local template that is embedded into each offline output."""
    return EXPLORER_TEMPLATE.read_text(encoding="utf-8")


def _write_manifest(destination, value):
    path = destination / "atlas_manifest.json"
    partial = destination / "atlas_manifest.json.partial"
    partial.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    partial.replace(path)


@contextmanager
def _atlas_write_transaction(destination, provenance):
    """Never leave a previous complete receipt over partially replaced figures."""
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / ".pca_atlas.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("Another atlas render is using this output directory") from None
        status = dict(version=VERSION, status="running", created_utc=datetime.now(timezone.utc).isoformat(),
                      provenance=provenance, code_sha256=sha256(__file__))
        _write_manifest(destination, status)
        try:
            yield
        except BaseException as exc:
            _write_manifest(destination, dict(status, status="failed", reason=str(exc)))
            raise


def render_atlas(payload, out, time_scale="linear", max_age_bp=10000, layer_years=1000, provenance=None, *, static=True, theme="dark", color_by="region", replace=False):
    validate_points(payload)
    validate_overlays(payload)
    if theme not in {"dark", "light"}:
        raise ValueError("theme must be dark or light")
    if color_by not in {"region", "layer", "age", "cohort", "source"}:
        raise ValueError("color-by must select a defined core track")
    if time_scale not in {"linear", "log1p"}:
        raise ValueError("time-scale must be linear or log1p")
    if isinstance(layer_years, bool) or not isinstance(layer_years, int) or layer_years < 1:
        raise ValueError("layer-years must be a positive integer")
    if isinstance(max_age_bp, bool) or not isinstance(max_age_bp, (int, float)) or not math.isfinite(max_age_bp) or max_age_bp <= 0:
        raise ValueError("max-age-bp must be finite and positive")
    destination = Path(out).expanduser().resolve()
    original = copy.deepcopy(payload["points"])
    counts = dict(total=len(original), ancient=sum(p["cohort"] == "ancient" for p in original),
                  modern=sum(p["cohort"] == "modern" for p in original),
                  ancient_primary=sum(p["cohort"] == "ancient" and p["age_bp"] <= max_age_bp for p in original),
                  ancient_older_than_primary=sum(p["cohort"] == "ancient" and p["age_bp"] > max_age_bp for p in original))
    provenance = copy.deepcopy(provenance or dict(source_type="supplied_payload", plan_validation_status="embedded_provenance_unverified_plan"))
    outputs = ["pca_time_atlas.html", "atlas_payload.json.gz", "original_points.tsv.gz", "sample_display_audit.tsv.gz", "layer_counts.tsv",
               "sample_feature_metadata.tsv", "sample_feature_values.tsv.gz", "variant_overlay_metadata.tsv"]
    if static:
        for name in ("pca_time_stack", "pca_time_stack_alternative", "pca_all_ages"):
            outputs += [name + ".png", name + ".pdf"]
    reserved = {(destination / name).resolve() for name in outputs + ["atlas_manifest.json", "atlas_manifest.json.partial"]}
    input_paths = {Path(v["path"]).resolve() for v in provenance.get("inputs", {}).values() if isinstance(v, dict) and v.get("path")}
    if reserved & input_paths:
        raise ValueError("Output would overwrite a source input; use a separate output directory")
    if any((destination / name).is_symlink() for name in outputs + ["atlas_manifest.json", "atlas_manifest.json.partial", ".pca_atlas.lock"]):
        raise ValueError("Atlas output files must not be symbolic links; use a separate output directory")
    data = dict(points=original, variants=copy.deepcopy(payload.get("variants", [])),
                features=copy.deepcopy(payload.get("features", [])), simulated=bool(payload.get("simulated", False)),
                coordinate_kind=payload.get("coordinate_kind", "unspecified_shared_coordinates"),
                basis_label=payload.get("basis_label", "unspecified source basis"), provenance=provenance,
                display=dict(time_scale=time_scale, max_age_bp=max_age_bp, layer_years=layer_years, default_view="continuous",
                             log1p_scale_bp=LOG_SCALE_BP, xy_offset=0, date_jitter=0,
                             theme=theme, point_size=2.7, opacity=.75, default_color=color_by,
                             display_downsampling=False, linked_projection=True),
                counts=counts, limits=shared_limits(original), region_colors=palette(original), region_names=REGION_NAMES,
                layers=layer_rows(original, layer_years, max_age_bp))
    receipt = destination / 'atlas_manifest.json'
    if not replace and receipt.is_file():
        saved = json.loads(receipt.read_text())
        if (saved.get('status') == 'complete' and saved.get('provenance') == provenance
                and saved.get('display') == data['display']
                and saved.get('code_sha256') == sha256(__file__)
                and saved.get('template_sha256') == sha256(EXPLORER_TEMPLATE)
                and saved.get('explorer_script_sha256') == sha256(EXPLORER_SCRIPT)
                and saved.get('source_points_sha256') == hashlib.sha256(json.dumps(original, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
                and saved.get('payload_sha256') == hashlib.sha256(json.dumps(data, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
                and set(saved.get('output_files', {})) == set(outputs)):
            if all((destination / name).is_file() and file_identity(destination / name) == spec
                   for name, spec in saved['output_files'].items()):
                print('[pca3d] SKIP verified atlas: ' + str(destination), file=sys.stderr, flush=True)
                return dict(saved, reused=True)
    with _atlas_write_transaction(destination, provenance):
        result = _render_atlas_outputs(data, destination, outputs, static)
        if original != payload["points"]:
            raise AssertionError("Rendering changed source points")
        _write_manifest(destination, result)
        return result


def _render_atlas_outputs(data, destination, outputs, static):
    from plotly.offline import get_plotlyjs
    original, counts, provenance = data["points"], data["counts"], data["provenance"]
    time_scale = data["display"]["time_scale"]
    max_age_bp, layer_years = data["display"]["max_age_bp"], data["display"]["layer_years"]
    with gzip.open(destination / "atlas_payload.json.gz", "wt", encoding="utf-8") as stream:
        json.dump(data, stream, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    fields = list(dict.fromkeys(k for p in original for k in p))
    write_tsv(destination / "original_points.tsv.gz", fields, original)
    audit = [dict(sample_id=p["sample_id"], cohort=p["cohort"], age_bp=p["age_bp"], pc1=p["pc1"], pc2=p["pc2"],
                  primary_display=int(p["cohort"] == "modern" or p["age_bp"] <= max_age_bp),
                  age_scope="modern_separate" if p["cohort"] == "modern" else "ancient_primary" if p["age_bp"] <= max_age_bp else "ancient_older_supplement",
                  layer_index=layer_index(p, layer_years), layer_display_bp=display_age(p, layer_years)) for p in original]
    write_tsv(destination / "sample_display_audit.tsv.gz", list(audit[0]), audit)
    write_tsv(destination / "layer_counts.tsv", list(data["layers"][0]), data["layers"])
    metadata_fields = ["feature_id", "label", "kind", "units", "provenance"]
    write_tsv(destination / "sample_feature_metadata.tsv", metadata_fields, data["features"])
    feature_rows = (dict(sample_id=p["sample_id"], **{f["feature_id"]: f["values"][i]
                    if f["values"][i] is not None else "NA" for f in data["features"]}) for i, p in enumerate(original))
    write_tsv(destination / "sample_feature_values.tsv.gz", ["sample_id", *(f["feature_id"] for f in data["features"])], feature_rows)
    variant_fields = ["variant_id", "chrom", "pos", "ref", "alt", "status", "count_kind", "dosage_evidence", "observed_sample_values"]
    write_tsv(destination / "variant_overlay_metadata.tsv", variant_fields,
              (dict(v, observed_sample_values=sum(x is not None for x in (v.get("values") or []))) for v in data["variants"]))
    encoded = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    # Replace placeholders only in the template, never in supplied sample IDs,
    # labels or bundled scripts that happen to contain placeholder-like text.
    replacements = {"__PLOTLY__": get_plotlyjs(), "__DATA__": encoded,
                    "__EXPLORER_JS__": EXPLORER_SCRIPT.read_text(encoding="utf-8")}
    document = re.sub(r"__PLOTLY__|__DATA__|__EXPLORER_JS__", lambda match: replacements[match.group()], explorer_template())
    (destination / "pca_time_atlas.html").write_text(document, encoding="utf-8")
    if static:
        draw_static(data, destination, time_scale, "primary", "pca_time_stack")
        draw_static(data, destination, "log1p" if time_scale == "linear" else "linear", "primary", "pca_time_stack_alternative")
        draw_static(data, destination, time_scale, "all", "pca_all_ages")
    manifest = dict(version=VERSION, status="complete", created_utc=datetime.now(timezone.utc).isoformat(),
                    payload_sha256=hashlib.sha256(json.dumps(data, sort_keys=True, ensure_ascii=False).encode()).hexdigest(),
                    counts=counts, display=data["display"], coordinate_kind=data["coordinate_kind"], basis_label=data["basis_label"],
                    simulated=data["simulated"], overlay_variants=len(data["variants"]), sample_features=len(data["features"]),
                    provenance=provenance, coordinates_unchanged=True, ages_unchanged=True, genetic_axis_equal_span=True,
                    recomputed_PCA=False, PC3_available=False, explained_variance_available=False,
                    source_points_sha256=hashlib.sha256(json.dumps(original, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest(),
                    code_sha256=sha256(__file__), template_sha256=sha256(EXPLORER_TEMPLATE),
                    explorer_script_sha256=sha256(EXPLORER_SCRIPT),
                    self_contained=True, network_requests_required=False, output_files={name: file_identity(destination / name) for name in outputs},
                    limitations=["Third axis is time, not PC3; genetic coordinates and dates are not jittered or re-estimated.",
                                 "Original PCA basis size, marker set, weights, eigenvalues and variance explained are not inferred from two coordinate columns.",
                                 "Modern comparison cohort is exactly what the source contains; no unavailable UKB participants are added.",
                                 "Layer counts reflect uneven sampling; a denser point cloud is not higher population abundance or genetic diversity.",
                                 "No trajectory lines, group-centroid ancestry paths or natural-selection conclusions are inferred.",
                                 "Variant overlays are the supplied ALT hardcall/dosage values; 0/1/2 never means three alleles. Missing stays null.",
                                 "Optional sample features are explicitly supplied observations with declared provenance; no window annotation or aggregate frequency is converted into individual phenotype data.",
                                 "Color and marker-size ranges are fixed over the complete supplied feature track. Filtering and time animation change display membership only.",
                                 "3D and linked 2D contain all currently filtered points; no stochastic or hidden display downsampling is performed."])
    return manifest


def parser():
    p = argparse.ArgumentParser(description=_description_pca_atlas)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--source-html", help="Existing allele3 HTML containing actual PC1/PC2/date payload")
    source.add_argument("--coordinates", help="Shared .sscore/.eigenvec; requires --metadata")
    p.add_argument("--metadata")
    p.add_argument("--sample-features", help="Optional exact sample_id TSV with measured/estimated per-sample features; requires --feature-schema")
    p.add_argument("--feature-schema", help="TSV: feature_id,label,kind(numeric/categorical),units,provenance")
    p.add_argument("--basis-label")
    p.add_argument("--simulated", action="store_true", help="Explicitly label supplied coordinates as simulated; source HTML simulation labels are always retained")
    p.add_argument("--out", required=True)
    p.add_argument('--replace', choices=('TRUE', 'FALSE'), default='FALSE')
    p.add_argument("--time-scale", choices=("linear", "log1p"), default="linear")
    p.add_argument("--max-age-bp", type=float, default=10000)
    p.add_argument("--layer-years", type=int, default=1000)
    p.add_argument("--theme", choices=("dark", "light"), default="dark")
    p.add_argument("--color-by", choices=("region", "layer", "age", "cohort", "source"), default="region")
    p.add_argument("--no-static", action="store_true", help="Generate the complete offline explorer and audits without additional Matplotlib figures")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.source_html:
        if args.metadata:
            raise ValueError("--metadata is for --coordinates; source HTML already contains its metadata")
        if args.basis_label:
            raise ValueError("--basis-label is for --coordinates; source HTML retains its original basis label")
        payload, provenance = read_source_html(args.source_html)
        if args.simulated:
            payload["simulated"] = True
    else:
        if not args.metadata:
            raise ValueError("--coordinates requires --metadata")
        payload, provenance = read_coordinate_source(args.coordinates, args.metadata, args.basis_label, simulated=args.simulated)
    if bool(args.sample_features) != bool(args.feature_schema):
        raise ValueError("--sample-features and --feature-schema must be provided together")
    if args.sample_features:
        if payload.get("features"):
            raise ValueError("Source HTML already contains sample features; use its existing tracks or provide a coordinate source")
        payload["features"], provenance["sample_feature_audit"] = read_sample_features(
            args.sample_features, args.feature_schema, [p["sample_id"] for p in payload["points"]])
        provenance["inputs"].update(sample_features=file_identity(args.sample_features), feature_schema=file_identity(args.feature_schema))
    result = render_atlas(payload, args.out, args.time_scale, args.max_age_bp, args.layer_years, provenance,
                          static=not args.no_static, theme=args.theme, color_by=args.color_by, replace=args.replace == 'TRUE')
    print(json.dumps(dict(status=result["status"], counts=result["counts"], out=str(Path(args.out).resolve()))))


def _entry_pca_atlas():
    main()


# ========================================================================
# pca trajectory
# ========================================================================
_description_pca_trajectory = 'Temporal structure in frozen shared projected PCA coordinates; never fits PCA.'
#!/usr/bin/env python3
"""Temporal structure in frozen shared projected PCA coordinates; never fits PCA."""

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
trajectory_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(trajectory_paths_root / "0f"))
from workflow_paths import configured_path
from collections import defaultdict
import numpy as np
from analysis_integrity import atomic_json, identity, transactional, code_identity
from pca import frame_read_scores as read_scores, frame_historical_identity as historical_identity, frame_coordinate_join as coordinate_join
from allele_measurement import (stats_number as number, stats_read_table as read_table, stats_write_table as write_table, stats_technology as technology, stats_age_bin as age_bin, stats_seed_for as seed_for,
                            stats_bootstrap_mean as bootstrap_mean, stats_fixed_composition as fixed_composition, stats_linear_model as linear_model, stats_adjust_fdr as adjust_fdr, STATS_NA as NA)

TRAJECTORY_LIMITATIONS = [
    'PC1/PC2/age BP are the original three axes; time is not PC3 and never enters the genotype PCA fit.',
    'PC changes can reflect migration, mixture and sampling. No selection, ancestral path, or evolutionary endpoint is inferred.',
    'Intervals resample unrelated individuals conditional on frozen coordinates, the supplied cohort and declared composition.',
    'Source is a measurement category, not an independent study/site cluster; no invented clusters or dating errors.',
    'PCA dispersion is not nucleotide diversity. Modern samples are a separate endpoint, never ancient age-zero observations.',
    'Shotgun labels do not establish independent, unimputed whole-genome validation. Genotype encoding needs source documentation.',
    'No missing-marker reprojection, new basis fitting, or independent read validation is performed by this coordinate module.',
]


def trajectory_read_unique(path):
    result = {}
    for row in read_table(path):
        sid = row.get('sample_id')
        if not sid or sid in result: raise ValueError('Missing/duplicate metadata sample_id: ' + str(sid))
        result[sid] = row
    if not result: raise ValueError('Empty metadata table')
    return result


def trajectory_training_members(args, contract):
    basis = Path(contract['basis_dir']) if contract.get('basis_dir') else None
    path = Path(args.training_ids) if args.training_ids else basis / 'basis.keep' if basis else None
    if not path or not path.is_file(): return set(), {'status': 'unavailable'}
    recorded = None
    if basis:
        plan = json.loads((basis / 'pca_plan.json').read_text())
        matched = [digest for name, digest in plan.get('plan_file_sha256', {}).items() if Path(name).name == path.name]
        if len(matched) == 1: recorded = {'sha256': matched[0]}
    receipt = historical_identity(path, args.provenance_manifest, recorded)
    ids = []
    with path.open() as handle:
        for line in handle:
            parts = line.split()
            if not parts or parts[0].startswith('#') or parts[-1] == 'IID': continue
            ids.append(parts[-1])
    if len(set(ids)) != len(ids): raise ValueError('Duplicate training ID')
    return set(ids), dict(receipt, status=receipt['historical_metadata_verification'])


def trajectory_load_inputs(args):
    requested = args.pcs or 2
    coords, columns, contract = read_scores(args.coordinates, requested, basis_dir=args.basis_dir,
        receipt=args.coordinate_receipt, require_verified=False)
    # Optional prespecified 5/10-PC supplements only when every column exists.
    if args.pcs is None:
        for k in (10, 5):
            try:
                extra, names, detail = read_scores(args.coordinates, k, basis_dir=args.basis_dir,
                    receipt=args.coordinate_receipt, require_verified=False)
            except ValueError as exc:
                if 'Missing consecutive requested PCs' in str(exc): continue
                raise
            coords, columns, contract = extra, names, detail; break
    contract.update(primary_pc_count=requested, supplemental_k=[k for k in (5, 10) if k <= len(columns)])
    meta = trajectory_read_unique(args.metadata)
    meta_receipt = historical_identity(args.metadata, args.provenance_manifest, contract.get('metadata_recorded'))
    training, train_receipt = trajectory_training_members(args, contract)
    supplement_info = None
    if args.supplement_metadata:
        supplemental = trajectory_read_unique(args.supplement_metadata)
        for sid, row in supplemental.items():
            if sid not in meta: raise ValueError('Supplement has unmatched sample_id: ' + sid)
            for key, value in row.items():
                if value in NA or key == 'sample_id': continue
                if key in meta[sid] and meta[sid][key] not in NA and meta[sid][key] != value:
                    raise ValueError('Supplement would rewrite existing metadata: ' + sid + ':' + key)
                meta[sid][key] = value
        supplement_info = identity(args.supplement_metadata, full=True)
    biological = Counter(r.get('individual_id') for r in meta.values() if r.get('individual_id') not in NA | {None})
    rows = []
    # Coordinate order is unchanged, with unmatched metadata explicitly appended.
    ids = list(coords) + [sid for sid in meta if sid not in coords]
    for sid in ids:
        raw = dict(meta.get(sid, {})); age = number(raw.get('age_bp')); cohort = raw.get('cohort', 'unknown')
        region = raw.get('region', raw.get('stratum', 'unknown'))
        if region in NA: region = 'unknown'
        source = raw.get('data_source', raw.get('source', 'unknown'))
        if source in NA: source = 'unknown'
        tech = technology(raw)
        reason = []
        if sid not in meta: reason.append('metadata_missing')
        if sid not in coords: reason.append('coordinates_missing')
        if raw.get('include') != '1': reason.append('include_not_explicit_1')
        if raw.get('unrelated') != '1': reason.append('related' if raw.get('unrelated') == '0' else 'relationship_unknown')
        if cohort not in {'ancient', 'modern'}: reason.append('cohort_unknown')
        if age is None or age < 0: reason.append('invalid_date')
        if region == 'unknown': reason.append('region_unknown')
        if biological.get(raw.get('individual_id'), 0) > 1: reason.append('duplicate_biological_individual')
        formal = not reason
        membership = 'training' if sid in training else 'projected' if train_receipt.get('status') == 'verified' else 'unknown'
        if train_receipt.get('status') != 'verified': membership = 'unknown'
        row = dict(raw, sample_id=sid, IID=sid, cohort=cohort, age_bp=age, original_age_bp=raw.get('age_bp'),
                   region=region, data_source=source, technology=tech, measurement_source=tech + '|' + source,
                   formal=formal, formal_exclusion=';'.join(reason), coordinate_present=int(sid in coords),
                   training_status=membership, individual_audit='known' if raw.get('individual_id') not in NA | {None} else 'not_provided',
                   basis_id=contract.get('basis_id'), coordinate_kind=contract['coordinate_kind'],
                   time_bin=age_bin(age, args.bin_years, args.max_age_bp) if cohort == 'ancient' else None,
                   date_uncertainty='not_provided', iqs_definition='as supplied; not inferred from coverage',
                   genotype_encoding='not_established_from_coordinates')
        for k in range(len(columns)): row[f'PC{k+1}'] = coords[sid][k] if sid in coords else None
        rows.append(row)
    # Retain FID, denominators and measurement fields from the exact original score.
    opener = gzip.open if str(args.coordinates).endswith('.gz') else open
    score_info = {}
    with opener(args.coordinates, 'rt') as handle:
        lines = (line for line in handle if line.strip() and not line.startswith('##'))
        names = next(lines).lstrip('#').split()
        for line in lines:
            record = dict(zip(names, line.split())); score_info[record['IID']] = record
    for row in rows:
        record = score_info.get(row['sample_id'], {})
        for key in ('FID', 'ALLELE_CT', 'DENOM', 'NAMED_ALLELE_DOSAGE_SUM', 'marker_count', 'call_fraction'):
            if key in record: row[key] = record[key]
    return rows, coords, columns, dict(contract, join=coordinate_join(list(meta.values()), coords),
        metadata_provenance=meta_receipt, training_provenance=train_receipt, supplement=supplement_info)


def trajectory_composition(samples, width, maximum):
    groups = defaultdict(list)
    for s in samples:
        epoch = age_bin(s['age_bp'], width, maximum) if s['cohort'] == 'ancient' else 'modern' if s['cohort'] == 'modern' else 'unknown'
        if epoch is None: epoch = 'older_or_invalid'
        groups[s['region'], epoch, s['cohort'], s['measurement_source']].append(s)
    result = []
    for (region, epoch, cohort, source), ss in sorted(groups.items(), key=lambda x: str(x[0])):
        row = dict(region=region, time_bin=epoch, cohort=cohort, measurement_source=source, n=len(ss),
                   n_formal=sum(s['formal'] for s in ss), n_related=sum(s.get('unrelated') == '0' for s in ss),
                   n_relationship_unknown=sum(s.get('unrelated') not in {'0','1'} for s in ss),
                   n_training=sum(s['training_status'] == 'training' for s in ss),
                   n_projected=sum(s['training_status'] == 'projected' for s in ss),
                   n_modern=sum(s['cohort'] == 'modern' for s in ss),
                   exclusion_counts=json.dumps(dict(Counter(r for s in ss for r in s['formal_exclusion'].split(';') if r)), sort_keys=True))
        for name in ('target_coverage', 'iqs', 'marker_count', 'call_fraction', 'ALLELE_CT', 'DENOM'):
            values = [number(s.get(name)) for s in ss]; values = [v for v in values if v is not None]
            row.update({name+'_n': len(values), name+'_median': float(np.median(values)) if values else None,
                        name+'_q25': float(np.quantile(values,.25)) if values else None, name+'_q75': float(np.quantile(values,.75)) if values else None})
        result.append(row)
    return result


def trajectory_scaling(samples, k):
    reference = [s for s in samples if s['coordinate_present'] and s['training_status'] == 'training']
    definition = 'verified_frozen_training_members'
    if len(reference) < 2:
        reference = [s for s in samples if s['formal'] and s['coordinate_present']]
        definition = 'fixed_entire_formal_reference_cohort; training_members_unavailable'
    if len(reference) < 2: raise ValueError('At least two fixed reference samples required for PC scaling')
    matrix = np.array([[s[f'PC{i+1}'] for i in range(k)] for s in reference])
    scale = matrix.std(axis=0, ddof=1)
    result = [dict(pc=f'PC{i+1}', mean=float(matrix[:,i].mean()), scale=float(scale[i]), n=len(reference),
                   reference_population=definition, sample_ids=';'.join(s['sample_id'] for s in reference),
                   status='usable' if scale[i] > 1e-15 else 'zero_reference_variance') for i in range(k)]
    scale[scale <= 1e-15] = np.nan
    return scale, result


def trajectory_trajectory_tables(samples, k, args, scales):
    centroids, dispersions, contrasts, support, cache = [], [], [], [], {}
    regions = sorted({s['region'] for s in samples if s['region'] != 'unknown'})
    epochs = range(math.ceil(args.max_age_bp / args.bin_years))
    for region in ['all', *regions]:
        eligible = [s for s in samples if s['formal'] and (region == 'all' or s['region'] == region)]
        modern = [s for s in eligible if s['cohort'] == 'modern']
        for epoch in [*epochs, 'modern']:
            ss = modern if epoch == 'modern' else [s for s in eligible if s['cohort'] == 'ancient' and s['time_bin'] == epoch]
            n = len(ss); ok = n >= args.min_samples
            base = dict(analysis_set='formal_include_unrelated', region=region, cohort='modern' if epoch == 'modern' else 'ancient',
                        time_bin=epoch, age_start_bp=None if epoch == 'modern' else epoch*args.bin_years,
                        age_end_bp=None if epoch == 'modern' else min((epoch+1)*args.bin_years,args.max_age_bp),
                        n=n, n_eff=n, status='estimated' if ok else 'not_estimable', reason='' if ok else 'below_min_samples',
                        weights_id='equal_individual', support_id='fixed_bins', seed=args.seed)
            support.append(dict(base, minimum=args.min_samples, modern_reference_n=len(modern),
                                modern_reference_status='available' if len(modern) >= args.min_samples else 'no_matched_modern_support'))
            values = np.array([[s[f'PC{i+1}'] for i in range(k)] for s in ss]).reshape(n,k)
            se, lo, hi, draws = bootstrap_mean(values,args.bootstrap,seed_for(args.seed,(region,epoch))) if ok else (None,None,None,None)
            if ok: cache[region,epoch] = (values.mean(axis=0),draws,n)
            for i in range(k):
                centroids.append(dict(base, measure='raw_mean', pc_or_k=f'PC{i+1}', estimate=float(values[:,i].mean()) if ok else None,
                    descriptive_mean=float(values[:,i].mean()) if n else None, robust_center=float(np.median(values[:,i])) if n else None,
                    individual_q25=float(np.quantile(values[:,i],.25)) if n else None, individual_q75=float(np.quantile(values[:,i],.75)) if n else None,
                    age_median=float(np.median([s['age_bp'] for s in ss])) if n else None,
                    actual_age_min=min((s['age_bp'] for s in ss),default=None),actual_age_max=max((s['age_bp'] for s in ss),default=None),
                    se=float(se[i]) if se is not None else None,ci_low=float(lo[i]) if lo is not None else None,ci_high=float(hi[i]) if hi is not None else None,
                    sample_ids=';'.join(s['sample_id'] for s in ss), uncertainty='mean_bootstrap_CI; individual_quartiles_are_distribution'))
            for dimension in sorted({2, k, *[v for v in (5,10) if v <= k]}):
                if dimension > k: continue
                valid = ok and np.isfinite(scales[:dimension]).all()
                z = values[:,:dimension]/scales[:dimension] if valid else None
                dispersion = float(np.mean(np.sum((z-z.mean(axis=0))**2,axis=1))) if valid else None
                dbase = dict(base,pc_or_k=dimension)
                if not valid: dbase.update(status='not_estimable',reason=base['reason'] or 'zero_reference_variance')
                dispersions.append(dict(dbase,measure='mean_squared_distance_to_centroid_in_fixed_reference_scale',estimate=dispersion))
                matched = valid and epoch != 'modern' and region != 'all' and len(modern) >= args.min_samples
                distance = None
                if matched:
                    endpoint = np.array([[s[f'PC{i+1}'] for i in range(dimension)] for s in modern]).mean(axis=0)
                    distance = float(np.linalg.norm((values[:,:dimension].mean(axis=0)-endpoint)/scales[:dimension]))
                dispersions.append(dict(dbase,measure='distance_to_region_matched_modern',estimate=distance,
                    status='estimated' if matched else 'not_estimable',reason='' if matched else 'modern_separate_or_no_region_matched_support', modern_reference_n=len(modern)))
        # Adjacent bins only; gaps are never bridged.
        for epoch in list(epochs)[1:]:
            left,right = cache.get((region,epoch-1)),cache.get((region,epoch))
            for i in range(k):
                ok = left is not None and right is not None
                draws = right[1][:,i]-left[1][:,i] if ok and right[1] is not None and left[1] is not None else None
                contrasts.append(dict(region=region,cohort='ancient',analysis_set='formal_include_unrelated',time_bin=epoch,
                    age_start_bp=(epoch-1)*args.bin_years,age_end_bp=min((epoch+1)*args.bin_years,args.max_age_bp),
                    measure='older_bin_minus_adjacent_younger_bin',pc_or_k=f'PC{i+1}',estimate=float(right[0][i]-left[0][i]) if ok else None,
                    ci_low=float(np.quantile(draws,.025)) if draws is not None else None,ci_high=float(np.quantile(draws,.975)) if draws is not None else None,
                    n=left[2]+right[2] if ok else 0,n_eff=left[2]+right[2] if ok else None,
                    status='estimated' if ok else 'not_estimable',reason='' if ok else 'adjacent_bin_support_gap',seed=args.seed))
    return centroids,dispersions,contrasts,support


def trajectory_standardized_tables(samples, k, args):
    primary = [s for s in samples if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None]
    result, weights_out, support = [], [], []
    definitions = [('fixed_region','all','region',primary),('fixed_tech_source','all','measurement_source',primary)]
    definitions += [('fixed_tech_source',r,'measurement_source',[s for s in primary if s['region']==r]) for r in sorted({s['region'] for s in primary})]
    for method,region,key,ss in definitions:
        if not ss: continue
        counts = Counter(s[key] for s in ss)
        weights = {g:n/len(ss) for g,n in sorted(counts.items())}
        support_id=f'{method}:{region}:0-{args.max_age_bp}:{args.bin_years}'
        for g,w in weights.items():
            weights_out.append(dict(weights_id=support_id,stratum=g,weight=w,n_reference=counts[g],reference_population='entire_formal_ancient_window',region=region, fixed_across_all_bins=True))
        estimated = fixed_composition([[s[f'PC{i+1}'] for i in range(k)] for s in ss],[s['time_bin'] for s in ss],[s[key] for s in ss],
            weights,range(math.ceil(args.max_age_bp/args.bin_years)),args.min_samples,0 if args.command == 'trajectory-plan' else args.bootstrap,args.seed)
        for e in estimated:
            base={key:value for key,value in e.items() if key not in {'estimate','ci_low','ci_high','se','stratum_n'}}
            base.update(region=region,cohort='ancient',analysis_set='formal_include_unrelated',measure=method,weights_id=support_id,support_id=support_id,
                age_start_bp=e['time_bin']*args.bin_years,age_end_bp=min((e['time_bin']+1)*args.bin_years,args.max_age_bp),seed=args.seed)
            support.append(dict(base,stratum_n=json.dumps(e['stratum_n'],sort_keys=True)))
            for i in range(k): result.append(dict(base,pc_or_k=f'PC{i+1}',**{key:float(e[key][i]) if e[key] is not None else None for key in ('estimate','se','ci_low','ci_high')}))
    return result,weights_out,support


def trajectory_models(samples, k, verified, args):
    primary = [s for s in samples if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None]
    results = []
    for region in ['all', *sorted({s['region'] for s in primary})]:
        subset = [s for s in primary if region=='all' or s['region']==region]
        for i in range(k):
            fit=linear_model(subset,[s[f'PC{i+1}'] for s in subset],minimum=args.min_samples) if verified else dict(status='not_estimable',reason='unverified_basis_or_metadata',n=len(subset),estimate=None,p=None,q=None)
            results.append(dict(fit,region=region,cohort='ancient',measure='PC_response_time_association',pc_or_k=f'PC{i+1}',
                analysis_set='formal_include_unrelated',age_start_bp=0,age_end_bp=args.max_age_bp,seed=args.seed,
                hypothesis_family='all_prespecified_PC_by_region_and_global_linear_models', genetic_PC_covariates='none_PC_is_response'))
    return adjust_fdr(results)


def trajectory_sensitivities(samples, k, args, verified):
    rows=[]
    for name,predicate in [('all_formal',lambda s:True),('coverage_ge_1',lambda s:(number(s.get('target_coverage')) or 0)>=1),
                           ('shotgun_label',lambda s:s['technology']=='shotgun_label'),('training_only',lambda s:s['training_status']=='training'),
                           ('projected_only',lambda s:s['training_status']=='projected')]:
        subset=[s for s in samples if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None and predicate(s)]
        for i in range(min(k,2)):
            fit=linear_model(subset,[s[f'PC{i+1}'] for s in subset],minimum=args.min_samples) if verified else dict(status='not_estimable',reason='unverified_basis_or_metadata',n=len(subset),estimate=None)
            rows.append(dict(fit,analysis_set=name,measure='technical_subset_time_slope',pc_or_k=f'PC{i+1}',region='all',cohort='ancient',seed=args.seed))
    for field in ('target_coverage','marker_count','call_fraction','ALLELE_CT','DENOM'):
        for i in range(min(k,2)):
            pairs=[(number(s.get(field)),s[f'PC{i+1}']) for s in samples if s['formal'] and number(s.get(field)) is not None]
            valid=len(pairs)>=3 and np.std(np.array(pairs)[:,0])>0 and np.std(np.array(pairs)[:,1])>0
            rows.append(dict(analysis_set='formal',measure='descriptive_PC_correlation_'+field,pc_or_k=f'PC{i+1}',n=len(pairs),
                            estimate=float(np.corrcoef(np.array(pairs).T)[0,1]) if valid else None,
                            status='descriptive' if valid else 'not_estimable',reason='' if valid else 'missing_or_constant_covariate'))
    for width,offset in [(500,0),(2000,0),(args.bin_years,args.bin_years/2)]:
        groups=defaultdict(list)
        for s in samples:
            if s['formal'] and s['cohort']=='ancient' and s['age_bp'] is not None and 0<=s['age_bp']<=args.max_age_bp:
                groups[s['region'], math.floor((s['age_bp']-offset)/width)].append(s)
        for (region,epoch),ss in sorted(groups.items()):
            for i in range(min(k,2)):
                rows.append(dict(analysis_set=f'bin_{width}_offset_{offset:g}',region=region,cohort='ancient',measure='raw_mean_boundary_sensitivity',
                    age_start_bp=max(0,offset+epoch*width),age_end_bp=min(args.max_age_bp,offset+(epoch+1)*width),pc_or_k=f'PC{i+1}',n=len(ss),
                    estimate=float(np.mean([s[f'PC{i+1}'] for s in ss])) if len(ss)>=args.min_samples else None,
                    status='descriptive' if len(ss)>=args.min_samples else 'not_estimable',reason='' if len(ss)>=args.min_samples else 'below_min_samples'))
    for name,need in [('missing_marker_reprojection','actual missing-marker masks and explicit separate projection'),('dosage_vs_hardcall','actual dosage-bearing genotype files'),('independent_read_validation','independent reads and study IDs'),('alternative_basis','explicit new basis fit and its receipts')]:
        rows.append(dict(analysis_set=name,status='not_run',reason=need,measure='conditional_extension'))
    return rows


def trajectory_date_sensitivity(samples, args):
    if not args.date_distributions and not args.date_interval_uniform:
        return [], {'status':'not_provided','point_ages_preserved':True}
    distributions=json.loads(Path(args.date_distributions).read_text()) if args.date_distributions else {}
    known={s['sample_id'] for s in samples}
    if set(distributions)-known: raise ValueError('Date distribution contains unmatched sample IDs')
    prepared={}; modern_ignored=[]
    for s in samples:
        sid=s['sample_id']
        if s['cohort']!='ancient':
            if sid in distributions: modern_ignored.append(sid)
            continue
        if sid in distributions:
            d=distributions[sid]; ages=np.asarray(d['ages_bp'],dtype=float); probability=np.asarray(d['probabilities'],dtype=float)
            if ages.ndim!=1 or not len(ages) or ages.shape!=probability.shape or not np.isfinite(ages).all() or np.any(ages<0) or not np.isfinite(probability).all() or np.any(probability<0) or not math.isclose(float(probability.sum()),1.,abs_tol=1e-8):
                raise ValueError('Invalid calibrated age distribution: '+sid)
            prepared[sid]=('supplied_distribution',ages,probability/probability.sum())
        elif args.date_interval_uniform:
            lo,hi=number(s.get('age_lower_bp')),number(s.get('age_upper_bp'))
            if lo is not None and hi is not None:
                if lo<0 or hi<lo: raise ValueError('Invalid supplied age interval: '+sid)
                prepared[sid]=('assumed_uniform_interval',lo,hi)
    for s in samples:
        if s['sample_id'] in prepared: s['date_uncertainty'] = prepared[s['sample_id']][0]
    rng=np.random.default_rng(args.seed); rows=[]
    for draw in range(args.date_draws):
        grouped=defaultdict(list)
        for s in samples:
            if not s['formal'] or s['cohort']!='ancient': continue
            d=prepared.get(s['sample_id']); age=s['age_bp']
            if d: age=float(rng.choice(d[1],p=d[2])) if d[0]=='supplied_distribution' else float(rng.uniform(d[1],d[2]))
            epoch=age_bin(age,args.bin_years,args.max_age_bp)
            if epoch is not None: grouped[s['region'],epoch].append(s)
        for (region,epoch),ss in sorted(grouped.items()):
            for i in (1,2): rows.append(dict(draw=draw,region=region,time_bin=epoch,pc_or_k=f'PC{i}',n=len(ss),
                estimate=float(np.mean([s[f'PC{i}'] for s in ss])) if len(ss)>=args.min_samples else None,
                status='date_sensitivity' if len(ss)>=args.min_samples else 'not_estimable',seed=args.seed))
    return rows,dict(status='sensitivity_only',n_supplied=len(prepared),methods=dict(Counter(d[0] for d in prepared.values())),
        modern_not_sampled=modern_ignored,draws=args.date_draws,seed=args.seed,point_ages_preserved=True,
        distribution_input=identity(args.date_distributions,full=True) if args.date_distributions else None)


def trajectory_display_sample(samples, limit, seed):
    """Stable stratified display subset; statistical functions never receive it."""
    if len(samples)<=limit: return samples
    groups=defaultdict(list)
    for s in samples: groups[s['cohort'],s['region'],s.get('time_bin')].append(s)
    rng=np.random.default_rng(seed); selected=[]
    for key,ss in sorted(groups.items(),key=lambda x:str(x[0])):
        n=max(1,int(limit*len(ss)/len(samples)))
        selected.extend(ss[i] for i in rng.choice(len(ss),min(n,len(ss)),replace=False))
    if len(selected)>limit: selected=[selected[i] for i in rng.choice(len(selected),limit,replace=False)]
    return selected


def trajectory_figures(samples, centroids, dispersion, standardized, sensitivity, composition_rows, out, args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figdir=out/'figures';figdir.mkdir(exist_ok=True)
    points=[s for s in samples if s['coordinate_present'] and s['age_bp'] is not None]
    shown=trajectory_display_sample(points,args.display_limit,args.seed)
    write_table(figdir/'F1_plot_points.tsv.gz',shown)
    write_table(figdir/'F1_full_points.tsv.gz',points)
    def save(fig,name,caption):
        fig.suptitle(caption,fontsize=14)
        fig.text(.01,.012,f'N coordinates={len(points):,}; shown={len(shown):,}. Statistical tables use full eligible data. See ../manifest.json and sample_audit.tsv.gz.',fontsize=8)
        fig.tight_layout(rect=(0,.04,1,.94))
        for ext in ('png','pdf'):fig.savefig(figdir/f'{name}.{ext}',dpi=180)
        plt.close(fig)
    fig=plt.figure(figsize=(16,9));ax=fig.add_subplot(121,projection='3d');b=fig.add_subplot(122)
    for cohort,color in [('ancient','#1676a3'),('modern','#d16b38')]:
        ss=[s for s in shown if s['cohort']==cohort]
        ax.scatter([s['PC1'] for s in ss],[s['PC2'] for s in ss],[s['age_bp'] for s in ss],s=4,alpha=.4,label=f'{cohort} n={len(ss):,}',color=color)
        b.scatter([s['PC1'] for s in ss],[s['PC2'] for s in ss],s=4,alpha=.4,label=cohort,color=color)
    ax.set(xlabel='PC1 projected AVG',ylabel='PC2 projected AVG',zlabel='Original age BP (time)');b.set(xlabel='PC1 projected AVG',ylabel='PC2 projected AVG');ax.legend();b.legend()
    save(fig,'F1_shared_coordinates','Fixed genetic coordinates with original age; modern endpoint shown separately')
    write_table(figdir/'F2_composition.tsv',composition_rows)
    fig,axes=plt.subplots(2,2,figsize=(16,9));regions=sorted({r['region'] for r in composition_rows});bins=range(math.ceil(args.max_age_bp/args.bin_years))
    matrix=np.array([[sum(r['n_formal'] for r in composition_rows if r['region']==reg and r['cohort']=='ancient' and r['time_bin']==ep) for ep in bins] for reg in regions])
    axes[0,0].imshow(matrix,aspect='auto');axes[0,0].set(yticks=range(len(regions)),yticklabels=regions,xlabel='Ancient time bin',title='Formal sample N by region and time')
    for source in sorted({s['technology'] for s in samples}):
        counts=[sum(s['formal'] and s['cohort']=='ancient' and s['time_bin']==ep and s['technology']==source for s in samples) for ep in bins]
        axes[0,1].plot(list(bins),counts,label=source)
    axes[0,1].set(xlabel='Ancient time bin',ylabel='N',title='Technology labels');axes[0,1].legend(fontsize=8)
    ss=[s for s in shown if number(s.get('target_coverage')) is not None and s['cohort']=='ancient']
    axes[1,0].scatter([s['age_bp'] for s in ss],[float(s['target_coverage']) for s in ss],s=5,alpha=.3);axes[1,0].set(xlabel='Age BP',ylabel='Reported target coverage')
    reasons=Counter(r for s in samples for r in s['formal_exclusion'].split(';') if r)
    axes[1,1].barh(list(reasons),list(reasons.values()));axes[1,1].set(title='Formal exclusions (reasons may overlap)',xlabel='N')
    save(fig,'F2_sample_composition','Sampling and measurement audit; unknown relationships excluded from formal estimates')
    write_table(figdir/'F3_centroids.tsv',centroids)
    fig,axes=plt.subplots(1,2,figsize=(16,9))
    for i,ax in enumerate(axes,1):
        for region in regions:
            rows=[r for r in centroids if r['region']==region and r['cohort']=='ancient' and r['pc_or_k']==f'PC{i}']
            x=[(r['age_start_bp']+r['age_end_bp'])/2 for r in rows];y=[r['estimate'] if r['estimate'] is not None else np.nan for r in rows]
            ax.plot(x,y,marker='o',ms=3,label=region)
            ax.fill_between(x,[r['ci_low'] if r['ci_low'] is not None else np.nan for r in rows],[r['ci_high'] if r['ci_high'] is not None else np.nan for r in rows],alpha=.12)
        ss=[s for s in shown if s['formal'] and s['cohort']=='ancient' and s['time_bin'] is not None]
        ax.scatter([s['age_bp'] for s in ss],[s[f'PC{i}'] for s in ss],color='gray',s=2,alpha=.12,label='Individual distribution')
        ax.set(xlabel='Age BP',ylabel=f'PC{i} projected AVG',title=f'PC{i}: mean and individual distribution');ax.legend(fontsize=8)
    save(fig,'F3_regional_trajectories','Within-region means; shaded intervals are bootstrap CIs for means (not individual ellipses)')
    write_table(figdir/'F4_standardized.tsv',standardized);write_table(figdir/'F4_dispersion.tsv',dispersion)
    fig,axes=plt.subplots(2,2,figsize=(16,9))
    for i in (1,2):
        ax=axes[0,i-1]
        for method in ('raw','fixed_region','fixed_tech_source'):
            rows=[r for r in (centroids if method=='raw' else standardized) if r['region']=='all' and r['cohort']=='ancient' and r['pc_or_k']==f'PC{i}' and (method=='raw' or r['measure']==method)]
            ax.plot([(r['age_start_bp']+r['age_end_bp'])/2 for r in rows],[r['estimate'] if r['estimate'] is not None else np.nan for r in rows],marker='o',label=method)
        ax.set(xlabel='Age BP',ylabel=f'PC{i} projected AVG');ax.legend()
    for ax,measure in zip(axes[1],('mean_squared_distance_to_centroid_in_fixed_reference_scale','distance_to_region_matched_modern')):
        for region in regions:
            rows=[r for r in dispersion if r['region']==region and r['cohort']=='ancient' and r['pc_or_k']==2 and r['measure']==measure]
            ax.plot([(r['age_start_bp']+r['age_end_bp'])/2 for r in rows],[r['estimate'] if r['estimate'] is not None else np.nan for r in rows],label=region)
        ax.set(xlabel='Age BP',ylabel='Fixed reference scale',title='Within-bin dispersion' if measure.startswith('mean') else 'Distance to region-matched modern (NA without support)');ax.legend(fontsize=8)
    save(fig,'F4_fixed_composition','Raw and fixed-composition estimands; missing support remains a gap')
    write_table(figdir/'F5_sensitivity.tsv',sensitivity)
    fig,axes=plt.subplots(1,2,figsize=(16,9))
    for i,ax in enumerate(axes,1):
        rows=[r for r in sensitivity if r.get('measure')=='technical_subset_time_slope' and r.get('pc_or_k')==f'PC{i}']
        for j,r in enumerate(rows):
            if r.get('estimate') is not None:
                ax.plot(r['estimate'],j,'o')
                if r.get('ci_low') is not None:ax.plot([r['ci_low'],r['ci_high']],[j,j],'-')
            else:ax.text(0,j,'NA: '+str(r.get('reason')),fontsize=8)
        ax.axvline(0,color='.7');ax.set(yticks=range(len(rows)),yticklabels=[r['analysis_set'] for r in rows],xlabel=f'PC{i} per 1000 years toward present',title='Prespecified technical subsets')
    save(fig,'F5_technical_sensitivity','Conditional sensitivity; shotgun labels are not independent WGS validation')
    return [str(p.relative_to(out)) for p in sorted(figdir.glob('*.png'))]


@transactional('manifest.json')
def trajectory_run(args):
    if args.min_samples<2 or args.bin_years<1 or args.max_age_bp<args.bin_years or args.bootstrap<2 or args.display_limit<1 or args.date_draws<1 or args.pcs is not None and args.pcs<2:
        raise ValueError('Require >=2 PCs/samples/bootstrap, positive bins, display limit and date draws')
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    samples,coords,columns,contract=trajectory_load_inputs(args);k=len(columns)
    verified=contract['status']=='verified' and contract['metadata_provenance']['historical_metadata_verification']=='verified'
    for s in samples:
        s['basis_id']=contract.get('basis_id');s['coordinate_kind']=contract['coordinate_kind']
    counts=dict(total_coordinates=len(coords),total_metadata=contract['join']['n_metadata'],joined=contract['join']['n_joined'],
                displayed_cohorts=dict(Counter(s['cohort'] for s in samples if s['coordinate_present'])),formal=sum(s['formal'] for s in samples))
    plan=dict(status='planned',command=args.command,coordinate_contract=contract,sample_counts=counts,parameters=vars(args),
              bins='left-closed right-open; last ancient bin includes max_age_bp; modern separate',
              weights='fixed over all formal ancient samples in prespecified window; no missing-stratum extrapolation',
              inference='PC response linear HC3, region/source/technology and available coverage; all PC x region/global BH family',
              bootstrap='unrelated individual, conditional on fixed coordinates/composition',date_uncertainty='not_provided_unless_explicit_inputs',
              planned_outputs=['sample_audit.tsv.gz','sample_composition.tsv','support_audit.tsv','reference_weights.tsv','pc_scaling.tsv',
                               'temporal_centroids.tsv','temporal_dispersion.tsv','temporal_contrasts.tsv','standardized_trajectories.tsv',
                               'temporal_models.tsv','sensitivity_summary.tsv','figures','report.html'],code=code_identity())
    atomic_json(out/'trajectory_plan.json',plan)
    write_table(out/'sample_audit.tsv.gz',samples)
    comp=trajectory_composition(samples,args.bin_years,args.max_age_bp);write_table(out/'sample_composition.tsv',comp)
    standardized,weights,stdsupport=trajectory_standardized_tables(samples,k,args)
    write_table(out/'reference_weights.tsv',weights)
    if args.command=='trajectory-plan':
        write_table(out/'support_audit.tsv',stdsupport)
        return dict(plan,analysis_execution='not_run_plan_only',limitations=TRAJECTORY_LIMITATIONS)
    scales,scalerows=trajectory_scaling(samples,k);write_table(out/'pc_scaling.tsv',scalerows)
    centroid,dispersion,contrasts,support=trajectory_trajectory_tables(samples,k,args,scales)
    fits=trajectory_models(samples,k,verified,args);sens=trajectory_sensitivities(samples,k,args,verified)
    date_rows,date_status=trajectory_date_sensitivity(samples,args);write_table(out/'date_sensitivity.tsv.gz',date_rows)
    write_table(out/'sample_audit.tsv.gz',samples)
    for name,records in [('temporal_centroids',centroid),('temporal_dispersion',dispersion),('temporal_contrasts',contrasts),
                          ('standardized_trajectories',standardized),('support_audit',support+stdsupport),('temporal_models',fits),('sensitivity_summary',sens)]:
        for r in records:r.update(basis_id=contract.get('basis_id'),coordinate_kind=contract['coordinate_kind'])
        write_table(out/(name+'.tsv'),records)
    pictures=[] if args.no_plots else trajectory_figures(samples,centroid,dispersion,standardized,sens,comp,out,args)
    text='<html lang="en"><meta charset="utf-8"><title>Fixed PCA temporal structure</title><body style="font:16px system-ui;max-width:1400px;margin:32px auto"><h1>Fixed PCA temporal structure</h1>'
    text+='<p>Raw, fixed-composition, and adjusted PC-response estimands are separate. HERV conditional associations require the independent time stage. Modern is a separate endpoint.</p>'
    text+='<pre>'+html.escape(json.dumps(counts,indent=2))+'</pre>'
    for picture in pictures:text+='<img style="width:100%" src="'+picture+'">'
    text+='<ul>'+''.join('<li>'+html.escape(s)+'</li>' for s in TRAJECTORY_LIMITATIONS)+'</ul>'
    text+='<p>Provenance: '+html.escape('verified' if verified else 'unverified; formal adjusted models unavailable')+'</p><a href="manifest.json">Manifest and output inventory</a></body></html>'
    (out/'report.html').write_text(text)
    return dict(status='complete' if verified else 'complete_descriptive_unverified',schema='pca-trajectory-1',sample_counts=counts,
                coordinate_contract=contract,parameters=vars(args),date_uncertainty=date_status,limitations=TRAJECTORY_LIMITATIONS,
                hypothesis_family='all_prespecified_PC_by_region_and_global_linear_models',n_models=len(fits),
                n_estimable_models=sum(r['status']=='estimated' for r in fits),genotype_data_read=False,
                original_points_and_files_modified=False,unexecuted=['new_reference_fits','missing_marker_reprojection','dosage_validation','independent_reads'],
                figures=pictures,interpretation='coordinate temporal structure; no selection or causal HERV evidence')


def trajectory_parser():
    p=argparse.ArgumentParser(description=_description_pca_trajectory)
    p.add_argument('command',choices=['trajectory-plan','trajectory'])
    p.add_argument('--coordinates','--eigenvec',dest='coordinates',required=True)
    p.add_argument('--metadata',required=True)
    p.add_argument('--out',default=configured_path('ANALYSIS_ROOT', 'allele3/pca3d/temporal'))
    for name in ('basis-dir','coordinate-receipt','provenance-manifest','training-ids','supplement-metadata','date-distributions'):p.add_argument('--'+name)
    p.add_argument('--pcs',type=int,default=None,help='Default PC1/2 primary plus available 5/10 supplements; explicit K requires every PC')
    for name,value in [('bin-years',1000),('max-age-bp',10000),('min-samples',20),('bootstrap',500),('seed',20261028),('display-limit',20000),('date-draws',100)]:p.add_argument('--'+name,type=int,default=value)
    p.add_argument('--date-interval-uniform',action='store_true',help='Explicit assumption using real age_lower_bp/age_upper_bp fields')
    p.add_argument('--no-plots',action='store_true')
    return p


def _entry_pca_trajectory():
    try:
        trajectory_result=trajectory_run(trajectory_parser().parse_args());print(json.dumps({k:trajectory_result.get(k) for k in ('status','sample_counts','n_models','n_estimable_models')}))
    except (ValueError,OSError,RuntimeError) as trajectory_exc:raise SystemExit('pca-trajectory: '+str(trajectory_exc))


# ========================================================================
# pca display reuse
# ========================================================================
_description_pca_display_reuse = 'Reuse the historical plain PCA display without regenerating its figures.'
"""Reuse the historical plain PCA display without regenerating its figures."""


def display_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def display_reuse_display(html, inputs, basis_label, *, layer_years, time_view, color_by, simulated):
    html = Path(html)
    manifest = Path(str(html) + '.json')
    atlas = html.parent / 'atlas/atlas_manifest.json'
    if not manifest.is_file() or not atlas.is_file():
        return None
    saved = json.loads(manifest.read_text())
    if (saved.get('inputs') != inputs or saved.get('variants') != 0
            or saved.get('basis_label') != str(basis_label)
            or saved.get('simulated') != bool(simulated)):
        return None
    display = saved.get('display', {})
    if any(display.get(k) != v for k, v in dict(layer_years=layer_years, time_view=time_view, color_by=color_by).items()):
        return None
    # The historical renderer did not hash HTML itself, but its atlas did.
    parent = json.loads(atlas.read_text())
    source = parent.get('provenance', {}).get('inputs', {}).get('html', {})
    if source.get('path') != str(html.resolve()) or not html.is_file():
        return None
    if source.get('sha256') != display_sha256(html):
        raise ValueError('Completed PCA HTML changed: ' + str(html))
    static = saved.get('static_figures', {})
    paths = [Path(p) for p in static.get('outputs', [])]
    if static.get('layer_table'):
        paths.append(Path(static['layer_table']))
    if not paths or any(not p.is_file() or p.stat().st_size == 0 for p in paths):
        return None
    print('[pca3d] SKIP existing PCA display and figures: ' + str(html), flush=True)
    return dict(saved, manifest_path=str(manifest), reused=True)


if __name__ == '__main__':
    routes = {'atlas': _entry_pca_atlas, 'trajectory': _entry_pca_trajectory}
    arguments = list(sys.argv[1:])
    action = arguments.pop(0) if arguments and arguments[0] in routes else 'atlas'
    if action is None:
        raise SystemExit('Choose a stage: ' + ', '.join(routes))
    sys.argv[1:] = arguments
    routes[action]()
