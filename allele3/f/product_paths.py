"""Portable permanent-product roots and separate analysis runtime roots."""
from pathlib import Path
import os
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'0f'))
from workflow_paths import configured_path


RESEARCH_MODULES = {'avi_window_enrichment','herv_region_enrichment','archaic_tract_enrichment','alt_impact_similarity'}


def defaults():
    analysis = Path(os.environ.get('ALLELE_ANALYSIS_ROOT', configured_path('ANALYSIS_ROOT','allele3')))
    annot = Path(os.environ.get('ALLELE_ANNOTATION_ROOT', configured_path('DIR_ANNOT','allele3/38')))
    results = Path(os.environ.get('ALLELE_RESULT_ROOT', str(analysis/'main')))
    return dict(analysis=analysis, annotations=annot, results=results,
        summary=Path(os.environ.get('ALLELE_SUMMARY_STORE', str(analysis/'3.allele3_scan'))),
        detailed=Path(os.environ.get('ALLELE_SCAN_STORE', str(analysis/'3.allele3_scan/gnomad_detailed'))),
        avi=Path(os.environ.get('ALLELE_AVI_STORE', str(analysis/'4.avi_scan'))))


def binding(output):
    """Return runtime module + subpath for configured permanent destinations."""
    output = Path(output).expanduser().absolute()
    d = defaults()
    for key, module in (('summary','3.allele3_scan'), ('detailed','3.allele3_scan/gnomad_detailed'), ('avi','4.avi_scan')):
        root = d[key].expanduser().absolute()
        if output.is_relative_to(root):
            # Explicit legacy stores under analysis retain their existing layout.
            if root.is_relative_to(d['analysis'].expanduser().absolute()): return None
            return d['analysis']/module, output.relative_to(root)
    root = d['results'].expanduser().absolute()
    if output.is_relative_to(root):
        relative = output.relative_to(root)
        if len(relative.parts)>=2 and relative.parts[0] in RESEARCH_MODULES:
            return d['analysis'], relative
        if relative.parts and relative.parts[0] in {'hotspots','density','herv_enrichment','base_simulation','base_occupancy'} and not (len(relative.parts)>1 and relative.parts[1] in {'tmp','logs','jobs','qc'}):
            return root/relative.parts[0], Path(*relative.parts[1:])
    return None


def support(output, kind):
    pair = binding(output)
    return pair[0]/kind/pair[1] if pair else None


def permanent(output):
    return binding(output) is not None
