#!/usr/bin/env python3
"""Re-render fixed PC1/PC2 coordinates with an explicit chronological third axis.

No PCA or genotype inference is performed. The source may be an existing
allele3 HTML data payload or a shared coordinate table plus sample metadata.
All original points and dates remain in the downloadable audit, including
ancient people outside the primary display interval.
"""
from __future__ import annotations

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

VERSION = "pca-time-atlas-1.0"
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
    if point["cohort"] == "modern":
        return 0.0
    if view == "continuous":
        return point["age_bp"]
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
    ax.set_zticks([time_transform(v, scale) for v in ticks], ["0 · modern" if v == 0 else f"{v / 1000:g}" for v in ticks])
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
                     label=f"{region} · {REGION_NAMES.get(region, region)}")
              for region, color in colors.items() if any(p.get("region", "unknown") == region for p in points)]
    fig.legend(handles=legend, loc="upper left", bbox_to_anchor=(.79, .91), frameon=False, fontsize=9,
               title="Broad geographic region", title_fontsize=9.5)
    footer = (f"PC1 and PC2 use one fixed coordinate system and equal axis scales. Third axis = time, not PC3.\n"
              f"Points retain their exact recorded dates; {width:,}-year frames are visual guides. Modern references are at 0 BP.\n")
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


HTML_TEMPLATE = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Ancient DNA · fixed PCA through time</title><style>
*{box-sizing:border-box}body{margin:0;background:#f4f7f9;color:#183b50;font:15px system-ui,-apple-system,Segoe UI,sans-serif}header{padding:26px 32px 18px;background:white;border-bottom:1px solid #dce5ea}h1{margin:0 0 9px;font-size:28px;font-weight:650}header p{margin:5px 0;color:#556d7b}.layout{display:grid;grid-template-columns:minmax(480px,1fr) 292px;gap:18px;padding:20px}.card{background:white;border:1px solid #dce5ea;border-radius:10px;padding:17px}#plot{height:760px;min-height:560px}.controls h2{font-size:16px;margin:0 0 18px}.field{margin:15px 0}label{display:block;font-size:13px;font-weight:600;margin-bottom:6px}select,button{width:100%;padding:9px;border:1px solid #c6d4dd;background:#fff;border-radius:5px;color:#173b50;font:inherit}button{cursor:pointer}.note{font-size:12px;line-height:1.5;color:#587080}#counts{font-size:14px;line-height:1.6;background:#edf4f8;border-radius:6px;padding:12px;margin:18px 0}.legend{display:flex;gap:9px 16px;flex-wrap:wrap;padding:8px 12px;font-size:12px}.legend i{display:inline-block;width:9px;height:9px;margin-right:6px;border-radius:50%}.layers{max-height:350px;overflow:auto;margin-top:16px}.layers div{display:flex;justify-content:space-between;font-size:12px;padding:5px 0;border-bottom:1px solid #eef2f5}.details{padding:0 20px 22px}.details summary{cursor:pointer;font-weight:600;margin-bottom:10px}.details p{font-size:13px;line-height:1.6;max-width:1000px}.badge{color:#a23924;font-weight:700}#error{color:#a23924;padding:10px}@media(max-width:900px){.layout{grid-template-columns:1fr}.controls{display:grid;grid-template-columns:1fr 1fr;gap:0 16px}#plot{height:650px}.controls h2,#counts,.layers,.note{grid-column:1/-1}header{padding:20px}.layout{padding:12px}}
</style></head><body><header><h1>Genetic structure through time</h1><p>Shared PC1 and PC2 coordinates with a chronological third axis.</p><p id="subtitle"></p><p id="simulated" class="badge" hidden>SIMULATED INPUT — not observed ancient data</p></header>
<main class="layout"><section class="card"><div id="plot"></div><div id="error" hidden></div><div id="legend" class="legend"></div></section>
<aside class="card controls"><h2>One coordinate system</h2>
<div class="field"><label for="scope">Age coverage</label><select id="scope"><option value="primary">Primary window</option><option value="all">All original ages</option></select></div>
<div class="field"><label for="timeScale">Time spacing</label><select id="timeScale"><option value="linear">Linear · primary interpretation</option><option value="log1p">log1p · recent-time detail</option></select></div>
<div class="field"><label for="view">Age placement</label><select id="view"><option value="continuous">Exact recorded dates · primary</option><option value="layered">1,000-year bin midpoints · display only</option></select></div>
<div class="field"><label for="region">Broad region</label><select id="region"><option value="all">All regions</option></select></div>
<div id="counts"></div><button id="reset">Reset view and camera</button><p class="note">Modern people are separate diamonds at 0 BP. Region colors and PC coordinates stay fixed. Layer placement changes the display only; hovering always shows the original date.</p><div id="layers" class="layers"></div></aside></main>
<section class="details"><details class="card"><summary>Data, scale and interpretation</summary><p id="provenance"></p><p>No PCA is fitted by this viewer. No PC3, explained-variance percentage or basis sample count is inferred. Original dates and coordinates remain unchanged. A nonlinear time axis changes visual spacing, not the observations or statistical model. Samples share a display coordinate system; their movement across periods does not by itself establish migration, adaptation or selection.</p><p id="formula"></p><button id="download" style="max-width:260px">Download all original points (TSV)</button></details></section>
<script>__PLOTLY__</script><script id="allele3-data" type="application/json">__DATA__</script><script>
'use strict';const D=JSON.parse(document.getElementById('allele3-data').textContent),P=D.points,el=id=>document.getElementById(id),esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const width=D.display.layer_years,maxAge=D.display.max_age_bp,colors=D.region_colors;
el('timeScale').value=D.display.time_scale;el('view').value=D.display.default_view;el('view').options[1].text=width.toLocaleString()+'-year bin midpoints · display only';el('scope').options[0].text='0–'+(maxAge/1000)+' ka BP + modern';
el('subtitle').textContent=D.counts.ancient.toLocaleString()+' ancient + '+D.counts.modern.toLocaleString()+' modern people; '+D.counts.ancient_older_than_primary+' older ancient people retained outside the primary window.';el('simulated').hidden=!D.simulated;
for(const r of Object.keys(colors)){const o=document.createElement('option');o.value=r;o.text=r+' · '+(D.region_names[r]||r);el('region').appendChild(o)}
el('legend').innerHTML=Object.entries(colors).map(([r,c])=>'<span><i style="background:'+c+'"></i>'+esc(r+' · '+(D.region_names[r]||r))+'</span>').join('');
el('provenance').textContent='Source basis label: '+D.basis_label+'. Coordinate units: '+D.coordinate_kind+'. Provenance: '+D.provenance.plan_validation_status+'. This is a re-rendering of supplied fixed coordinates; original PCA plan, weights and eigenvalues are not revalidated here.';
el('formula').textContent='Linear: age/1,000. Optional log1p: log10(1 + age/1,000 BP), defined at modern age zero. Axis ticks and hover labels show original BP/ka. Genetic axes have identical fixed spans in both views.';
const transform=(a,scale)=>scale==='linear'?a/1000:Math.log1p(a/1000)/Math.LN10;
const index=p=>p.cohort==='modern'?-1:Math.floor(p.age_bp/width);
const label=k=>k<0?'Modern':(k*width/1000)+'–'+((k+1)*width/1000)+' ka';
const plane=p=>p.cohort==='modern'?0:(index(p)+.5)*width;
const camera={projection:{type:'orthographic'},eye:{x:1.45,y:-1.85,z:.84}};
let pending=null,drawing=false,resetCamera=true;
async function draw(traces,layout){pending={traces,layout};if(drawing)return;drawing=true;try{while(pending){const q=pending;pending=null;await Plotly.react('plot',q.traces,q.layout,{responsive:true,displaylogo:false,toImageButtonOptions:{format:'png',filename:'ancient_pca_time',scale:2}})}el('error').hidden=true}catch(e){el('error').hidden=false;el('error').textContent='Plot unavailable: '+e.message}finally{drawing=false}}
function update(){const scope=el('scope').value,scale=el('timeScale').value,view=el('view').value,region=el('region').value;
const pts=P.filter(p=>(p.cohort==='modern'||scope==='all'||p.age_bp<=maxAge)&&(region==='all'||p.region===region));
const scoped=P.filter(p=>p.cohort==='modern'||scope==='all'||p.age_bp<=maxAge);
const rawTop=scoped.reduce((m,p)=>Math.max(m,view==='layered'?(index(p)+1)*width:Math.ceil(p.age_bp/width)*width),Math.max(width,scope==='primary'?maxAge:0));const top=transform(rawTop,scale);
const layers=new Map();for(const p of pts)layers.set(index(p),(layers.get(index(p))||0)+1);
const traces=[];
const frames=view==='layered'?[...layers.keys()].map(k=>k<0?0:(k+.5)*width):Array.from({length:Math.floor(rawTop/width)+1},(_,i)=>i*width);
for(const a of frames){const z=transform(a,scale);traces.push({type:'scatter3d',mode:'lines',x:[D.limits.x[0],D.limits.x[1],D.limits.x[1]],y:[D.limits.y[0],D.limits.y[0],D.limits.y[1]],z:[z,z,z],line:{color:a===0?'#8EA4B1':'#D7E1E7',width:a===0?2:1},hoverinfo:'skip',showlegend:false})}
for(const cohort of ['ancient','modern']){const a=pts.filter(p=>p.cohort===cohort);if(!a.length)continue;traces.push({type:'scatter3d',mode:'markers',name:cohort==='ancient'?'Ancient':'Modern',x:a.map(p=>p.pc1),y:a.map(p=>p.pc2),z:a.map(p=>transform(view==='layered'?plane(p):p.cohort==='modern'?0:p.age_bp,scale)),marker:{size:cohort==='modern'?3.3:2.2,opacity:cohort==='modern'?.94:.64,symbol:cohort==='modern'?'diamond':'circle',color:a.map(p=>colors[p.region]||'#778899'),line:{width:0}},text:a.map(p=>esc(p.sample_id)+'<br>'+esc(p.cohort+' · '+p.region)+'<br>Original date: '+p.age_bp.toLocaleString()+' BP'+(view==='layered'?'<br>Display plane: '+label(index(p)):'')+'<br>PC1: '+p.pc1+'<br>PC2: '+p.pc2+'<br>'+esc(p.source||'')),hoverinfo:'text'})}
let ticks=scale==='linear'?Array.from({length:Math.floor(rawTop/(rawTop<=12000?2000:4000))+1},(_,i)=>i*(rawTop<=12000?2000:4000)):[0,500,1000,2000,4000,6000,10000,14000,18000,20000].filter(v=>v<=rawTop);
const scene={xaxis:{title:'PC1 · shared score',range:D.limits.x,autorange:false,showbackground:false,gridcolor:'#edf1f4'},yaxis:{title:'PC2 · shared score',range:D.limits.y,autorange:false,showbackground:false,gridcolor:'#edf1f4'},zaxis:{title:'Age (ka BP)'+(scale==='log1p'?' · nonlinear spacing':''),range:[-.025*top,1.035*top],autorange:false,showbackground:false,gridcolor:'#edf1f4',tickmode:'array',tickvals:ticks.map(v=>transform(v,scale)),ticktext:ticks.map(v=>v===0?'0 · modern':v/1000)},aspectmode:'manual',aspectratio:{x:1,y:1,z:1.6}};
if(resetCamera){scene.camera=camera;resetCamera=false}
draw(traces,{margin:{l:0,r:0,t:15,b:0},paper_bgcolor:'white',font:{family:'system-ui',color:'#2f5267'},showlegend:false,scene,uirevision:'fixed-pca-time-atlas'});
const anc=pts.filter(p=>p.cohort==='ancient').length;el('counts').textContent=pts.length.toLocaleString()+' displayed: '+anc.toLocaleString()+' ancient + '+(pts.length-anc).toLocaleString()+' modern. '+(D.counts.total-pts.length).toLocaleString()+' retained outside current filters.';
el('layers').innerHTML='<b>People per age bin</b>'+[...layers].sort((a,b)=>b[0]-a[0]).map(([k,n])=>'<div><span>'+label(k)+'</span><span>'+n.toLocaleString()+'</span></div>').join('');}
for(const id of ['scope','timeScale','view','region'])el(id).addEventListener('change',update);
el('reset').onclick=()=>{el('scope').value='primary';el('timeScale').value=D.display.time_scale;el('view').value=D.display.default_view;el('region').value='all';resetCamera=true;update()};
el('download').onclick=()=>{const names=new Set();for(const p of P)for(const k of Object.keys(p))names.add(k);const keys=[...names],clean=v=>String(v??'').replace(/[\t\r\n]/g,' ');const tsv=keys.join('\t')+'\n'+P.map(p=>keys.map(k=>clean(p[k])).join('\t')).join('\n')+'\n';const u=URL.createObjectURL(new Blob([tsv],{type:'text/tab-separated-values'})),a=document.createElement('a');a.href=u;a.download='original_pca_points.tsv';a.click();setTimeout(()=>URL.revokeObjectURL(u),1000)};
update();</script></body></html>'''


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


def render_atlas(payload, out, time_scale="linear", max_age_bp=10000, layer_years=1000, provenance=None, *, static=True):
    validate_points(payload)
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
    outputs = ["pca_time_atlas.html", "atlas_payload.json.gz", "original_points.tsv.gz", "sample_display_audit.tsv.gz", "layer_counts.tsv"]
    if static:
        for name in ("pca_time_stack", "pca_time_stack_alternative", "pca_all_ages"):
            outputs += [name + ".png", name + ".pdf"]
    reserved = {(destination / name).resolve() for name in outputs + ["atlas_manifest.json", "atlas_manifest.json.partial"]}
    input_paths = {Path(v["path"]).resolve() for v in provenance.get("inputs", {}).values() if isinstance(v, dict) and v.get("path")}
    if reserved & input_paths:
        raise ValueError("Output would overwrite a source input; use a separate output directory")
    if any((destination / name).is_symlink() for name in outputs + ["atlas_manifest.json", "atlas_manifest.json.partial", ".pca_atlas.lock"]):
        raise ValueError("Atlas output files must not be symbolic links; use a separate output directory")
    data = dict(points=original, variants=[], simulated=bool(payload.get("simulated", False)),
                coordinate_kind=payload.get("coordinate_kind", "unspecified_shared_coordinates"),
                basis_label=payload.get("basis_label", "unspecified source basis"), provenance=provenance,
                display=dict(time_scale=time_scale, max_age_bp=max_age_bp, layer_years=layer_years, default_view="continuous",
                             log1p_scale_bp=LOG_SCALE_BP, xy_offset=0, date_jitter=0),
                counts=counts, limits=shared_limits(original), region_colors=palette(original), region_names=REGION_NAMES,
                layers=layer_rows(original, layer_years, max_age_bp))
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
    encoded = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    document = HTML_TEMPLATE.replace("__PLOTLY__", get_plotlyjs()).replace("__DATA__", encoded)
    (destination / "pca_time_atlas.html").write_text(document, encoding="utf-8")
    if static:
        draw_static(data, destination, time_scale, "primary", "pca_time_stack")
        draw_static(data, destination, "log1p" if time_scale == "linear" else "linear", "primary", "pca_time_stack_alternative")
        draw_static(data, destination, time_scale, "all", "pca_all_ages")
    manifest = dict(version=VERSION, status="complete", created_utc=datetime.now(timezone.utc).isoformat(),
                    counts=counts, display=data["display"], coordinate_kind=data["coordinate_kind"], basis_label=data["basis_label"],
                    simulated=data["simulated"],
                    provenance=provenance, coordinates_unchanged=True, ages_unchanged=True, genetic_axis_equal_span=True,
                    recomputed_PCA=False, PC3_available=False, explained_variance_available=False,
                    source_points_sha256=hashlib.sha256(json.dumps(original, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest(),
                    code_sha256=sha256(__file__), output_files={name: file_identity(destination / name) for name in outputs},
                    limitations=["Third axis is time, not PC3; genetic coordinates and dates are not jittered or re-estimated.",
                                 "Original PCA basis size, marker set, weights, eigenvalues and variance explained are not inferred from two coordinate columns.",
                                 "Modern comparison cohort is exactly what the source contains; no unavailable UKB participants are added.",
                                 "Layer counts reflect uneven sampling; a denser point cloud is not higher population abundance or genetic diversity.",
                                 "No trajectory lines, group-centroid ancestry paths or natural-selection conclusions are inferred."])
    return manifest


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--source-html", help="Existing allele3 HTML containing actual PC1/PC2/date payload")
    source.add_argument("--coordinates", help="Shared .sscore/.eigenvec; requires --metadata")
    p.add_argument("--metadata")
    p.add_argument("--basis-label")
    p.add_argument("--simulated", action="store_true", help="Explicitly label supplied coordinates as simulated; source HTML simulation labels are always retained")
    p.add_argument("--out", required=True)
    p.add_argument("--time-scale", choices=("linear", "log1p"), default="linear")
    p.add_argument("--max-age-bp", type=float, default=10000)
    p.add_argument("--layer-years", type=int, default=1000)
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
    result = render_atlas(payload, args.out, args.time_scale, args.max_age_bp, args.layer_years, provenance)
    print(json.dumps(dict(status=result["status"], counts=result["counts"], out=str(Path(args.out).resolve()))))


if __name__ == "__main__":
    main()
