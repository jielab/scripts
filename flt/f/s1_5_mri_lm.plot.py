"""FLT anatomical MRI plots adapted from dose/scripts/s4_4_mri.plot.py.
Uses existing regression CSV; never refits models or changes FDR.
Atlas/package assets are shared read-only from the reference project.
Run directly with a Python environment containing the MRI plotting dependencies:
    ./flt.paper.sh fig4
Add --skip-anatomy for coefficient plots only; no PowerShell wrapper is needed.
"""
import argparse
import os
import sys
from pathlib import Path

ROOT = Path(os.environ["FLT_OUT"])

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tools-dir', type=Path, default=Path(os.environ['FLT_MRI_TOOLS']))
    p.add_argument('--results', type=Path, default=ROOT/'01lm/mri/mri_lm_results.csv')
    p.add_argument('--out', type=Path, default=ROOT/'01lm/mri/figures/anatomy')
    p.add_argument('--exposure', choices=['p20628','F-Salad_vegetables','p20638'])
    p.add_argument('--skip-anatomy', action='store_true', help='Only create tract/volume/WMH coefficient plots')
    args = p.parse_args()
    atlas = Path(os.getenv('FLT_MRI_ATLAS', args.tools_dir/'mri_atlases'))
    sys.path[:0] = [str(args.tools_dir/x) for x in
                   ['enigma_packages','ENIGMA','mri_plot_packages']]
    import pandas as pd
    from mri_plot_utils import reference_plots as plots
    r = pd.read_csv(args.results, dtype={'field_id':str})
    required = {'exposure','family','field_id','beta','LCI','HCI','FDR_family','status','model'}
    if not required.issubset(r):
        raise ValueError(f'Missing columns: {required-set(r)}')
    if r.duplicated(['exposure','field_id']).any():
        raise ValueError('Duplicate exposure/field pairs')
    foods = {'p20628':'Butter on bread','F-Salad_vegetables':'Total salad vegetable','p20638':'Chicken'}
    errors = []
    for exposure, food in foods.items():
        if args.exposure and exposure != args.exposure:
            continue
        d = r[r.exposure.eq(exposure)].copy()
        if len(d)!=260:
            raise ValueError(f'{exposure}: expected 260 planned MRI results')
        out = args.out/exposure
        out.mkdir(parents=True, exist_ok=True)
        plots.EXPOSURE = exposure
        print(f'Plotting {food}', flush=True)
        jobs = [('tract_forest', lambda: plots.tract_forest(d,out)),
                ('volume_WMH_forests', lambda: plots.supplementary_forests(d,out))]
        if not args.skip_anatomy:
            from mri_plot_utils import tracts_3d, volume_3d
            tracts_3d.EXPOSURE = volume_3d.EXPOSURE = exposure
            jobs += [('thickness_surface',lambda: plots.panel_a(d,out)),
                     ('tract_anatomy',lambda: tracts_3d.generate(d,out,atlas/'ukb_tracts27')),
                     ('volume_anatomy',lambda: volume_3d.generate(d,out,atlas))]
        status = []
        for name, fn in jobs:
            try:
                fn()
                status.append(f'{name}: OK')
            except Exception as exc:
                status.append(f'{name}: FAILED: {exc}')
                errors.append(f'{exposure}/{name}: {exc}')
            print(status[-1], flush=True)
        (out/'plot_status.txt').write_text('\n'.join(status),encoding='utf-8')
        (out/'README.txt').write_text(
            f'Exposure: {food} ({exposure}).\n'
            'beta: MRI SD per one FLT point. covs.all; ICV additionally for volumes/WMH.\n'
            'Stars/coloured anatomical regions use existing BH FDR within this FLT and full MRI family.\n'
            'Anatomy is a reference atlas, not individual segmentation or voxelwise significance.\n'
            'Colour ranges may differ between panels/exposures; consult each colour bar.\n'
            'Thickness surfaces show all estimates and FDR<0.05 separately.\n'
            'Tract and volume anatomy show FDR<0.05. Forests retain all estimates.\n'
            'WMH is shown as a coefficient, without inventing a spatial lesion map.\n',encoding='utf-8')
    if errors:
        raise RuntimeError('\n'.join(errors))

if __name__ == '__main__':
    main()
