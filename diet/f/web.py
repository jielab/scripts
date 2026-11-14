#!/usr/bin/env python3
"""Build the offline, aggregate-only diet.html report."""
from __future__ import annotations
import argparse
import base64
from datetime import datetime,timezone
import html
import io
import json
import os
from pathlib import Path
import tempfile
from results import rows,number,collect_effects

BASE=Path(__file__).resolve().parents[1]

def thumbnail(path):
    if path.suffix.lower()!=".png" or not path.is_file():return ""
    try:
        from PIL import Image
        with Image.open(path) as img:
            img=img.convert("RGB");img.thumbnail((1100,850))
            b=io.BytesIO();img.save(b,format="JPEG",quality=85)
        return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()
    except ImportError:
        if path.stat().st_size>2_000_000:return ""
        return "data:image/png;base64,"+base64.b64encode(path.read_bytes()).decode()

def build_data(out):
    scores=[]
    for p in sorted((out/"disco").glob("*.scores.tsv")):scores.extend(rows(p))
    effects=collect_effects(out)
    bayes=rows(out/"bayes/effects.tsv")
    finals=rows(out/"final/index.tsv")
    stages=rows(out/"logs/stages.tsv")
    availability=[]
    for p in sorted((out/"disco").glob("*.availability.tsv")):availability.extend(rows(p))
    # Completed and explicitly partial manuscript images are embedded. Full native figures stay
    # linked and downloadable; no individual-level data are serialized.
    embedded={}
    for item in finals:
        if item.get("kind")=="manuscript_item" and item.get("status") in ("available","partial") and item.get("file"):
            path=(out/item["file"]).resolve()
            if not path.is_relative_to(out):continue
            uri=thumbnail(path)
            if uri:embedded[item["file"]]=uri
    return dict(schema_version=1,generated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        source_commit="971a61b44b84500ab46e2f5592478505a0eae50e",scores=scores,effects=effects,bayes=bayes,
        figures=finals,stages=stages,availability=availability,images=embedded)

def render(out,dest=None):
    out=Path(out).resolve();dest=Path(dest or out/"diet.html")
    data=build_data(out)
    # Keep embedding safe even when a source outcome/label contains HTML-like text.
    packed=json.dumps(data,ensure_ascii=False,allow_nan=False,separators=(",",":"))
    packed=packed.replace("&","\\u0026").replace("<","\\u003c").replace(">","\\u003e").replace("\u2028","\\u2028").replace("\u2029","\\u2029")
    template=(BASE/"f/web.html").read_text(encoding="utf-8")
    if template.count("__DIET_DATA__")!=1:raise ValueError("Invalid diet report template")
    document=template.replace("__DIET_DATA__",packed)
    dest.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile("w",encoding="utf-8",dir=dest.parent,delete=False) as h:
        h.write(document);tmp=Path(h.name)
    os.replace(tmp,dest)
    return dest

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();print(render(a.output))
if __name__=="__main__":main()
