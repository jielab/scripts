#!/usr/bin/env python3
"""Effective regional dimensionality of fitted PRS-CSx posterior signal.

Fixed beta_squared primary and abs_beta sensitivity; no model fitting. Pure
metrics use only posterior weights and GRID boundaries. The optional aggregate
performance join runs after those metrics and never defines regional mass.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from block_inputs import (read_config, resolve_block_files, build_manifest,
                          locate_weights, read_weights, assign_blocks, POPS)

MASS_DEFINITIONS = ('beta_squared', 'abs_beta')
SCHEMA = 'grid.regional_dimensionality.1'
TABLES = ('regional_dimensionality_by_population.tsv', 'regional_dimensionality_trait.tsv',
          'regional_top_blocks.tsv', 'regional_dimensionality_vs_gain.tsv')


def identity(path):
    p = Path(path).resolve()
    s = p.stat()
    return dict(resolved_path=str(p), size=s.st_size, mtime_ns=s.st_mtime_ns)


def concentration(masses):
    x = np.asarray(masses, dtype=float)
    if not len(x) or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('Expected finite nonnegative mass in a nonempty block map')
    total = x.sum()
    names = ('top1pct_share', 'top5pct_share', 'top10pct_share', 'hhi',
             'effective_blocks', 'effective_blocks_fraction', 'shannon_effective_blocks',
             'shannon_effective_blocks_fraction', 'gini')
    if total == 0:
        return dict.fromkeys(names, np.nan)
    p = x / total
    hhi = np.square(p).sum()
    shannon = np.exp(-(p[p > 0] * np.log(p[p > 0])).sum())
    result = dict(hhi=hhi, effective_blocks=1/hhi, effective_blocks_fraction=1/hhi/len(x),
                  shannon_effective_blocks=shannon, shannon_effective_blocks_fraction=shannon/len(x),
                  gini=np.sum((2*np.arange(1, len(x)+1)-len(x)-1)*np.sort(p))/len(x))
    for q, name in zip((.01, .05, .10), names[:3]):
        result[name] = np.sort(p)[::-1][:max(1, int(np.ceil(q*len(x))))].sum()
    return result


def population_metrics(weights, bed, definition):
    beta = weights.BETA.to_numpy(float)
    mass = np.square(beta) if definition == 'beta_squared' else np.abs(beta)
    if definition not in MASS_DEFINITIONS or not np.isfinite(mass).all():
        raise ValueError('Invalid mass definition or nonfinite mass')
    assigned = assign_blocks(weights, bed)
    good = assigned >= 0
    block_mass = np.bincount(assigned[good], weights=mass[good], minlength=len(bed))
    counts = np.bincount(assigned[good], minlength=len(bed))
    total = block_mass.sum()
    metrics = dict(n_blocks=len(bed), n_weight_snps_total=len(weights),
                   n_weight_snps_assigned=int(good.sum()),
                   assigned_snp_fraction=good.mean() if len(good) else np.nan,
                   assigned_mass_fraction=total/mass.sum() if mass.sum() else np.nan,
                   total_assigned_mass=total, status='completed' if total > 0 else 'zero_mass',
                   **concentration(block_mass))
    blocks = bed[['block_id', 'CHR', 'start', 'end']].copy()
    blocks['n_snps'] = counts
    blocks['block_mass'] = block_mass
    blocks = blocks.sort_values('block_mass', ascending=False, kind='stable')
    blocks['rank'] = np.arange(1, len(blocks)+1)
    blocks['block_mass_share'] = blocks.block_mass/total if total else np.nan
    blocks['cumulative_mass_share'] = blocks.block_mass_share.cumsum()
    return metrics, blocks.head(25)


def summarize(populations):
    rows = []
    for (trait, definition), group in populations.groupby(['trait', 'mass_definition'], sort=False):
        g = group[group.status.eq('completed')]
        row = dict(trait=trait, mass_definition=definition, populations_available=len(g))
        for metric in ('effective_blocks_fraction', 'top5pct_share'):
            for label, q in (('median', .5), ('q25', .25), ('q75', .75)):
                row[label+'_'+metric] = g[metric].quantile(q)
        for label in ('min', 'max'):
            row[label+'_effective_blocks_fraction'] = getattr(g.effective_blocks_fraction, label)()
        for metric in ('hhi', 'gini', 'shannon_effective_blocks_fraction'):
            row['median_'+metric] = g[metric].median()
        rows.append(row)
    return pd.DataFrame(rows)


def join_gain(summary, performance):
    rows = []
    for trait, g in summary.groupby('trait', sort=False):
        primary = g[g.mass_definition.eq('beta_squared')].iloc[0]
        sensitivity = g[g.mass_definition.eq('abs_beta')].iloc[0]
        loss = 'brier' if trait == 't2dm' else 'mse'
        row = dict(trait=trait, primary_loss=loss, baseline_loss=np.nan, candidate_loss=np.nan,
                   absolute_loss_reduction=np.nan, relative_primary_loss_reduction=np.nan,
                   interpretation_status='not_available')
        for metric in ('effective_blocks_fraction', 'top5pct_share', 'hhi', 'gini'):
            row['median_'+metric+'_beta2'] = primary['median_'+metric]
        row['median_effective_blocks_fraction_abs_beta'] = sensitivity.median_effective_blocks_fraction
        if {'trait', 'method', loss} <= set(performance):
            d = performance[performance.trait.eq(trait)]
            if 'stage' in d:
                d = d[d.stage.eq('test')]
            base = d[d.method.eq('csx-linear')]
            candidate = d[d.method.eq('csx-blocknet')]
            if len(base) > 1 or len(candidate) > 1:
                raise ValueError('Ambiguous primary aggregate rows: '+trait)
            if len(base) == len(candidate) == 1:
                b, c = float(base.iloc[0][loss]), float(candidate.iloc[0][loss])
                if np.isfinite([b,c]).all() and b > 0 and c >= 0:
                    row.update(baseline_loss=b, candidate_loss=c, absolute_loss_reduction=b-c,
                               relative_primary_loss_reduction=(b-c)/b, interpretation_status='descriptive_only')
        rows.append(row)
    out = pd.DataFrame(rows)
    valid = out.dropna(subset=['relative_primary_loss_reduction', 'median_effective_blocks_fraction_beta2'])
    # Require all three prespecified traits and strict inverse ordering. No fit or P-value.
    if set(valid.trait) == {'height', 'ldl', 't2dm'}:
        ordered = valid.sort_values('median_effective_blocks_fraction_beta2')
        consistent = (np.diff(ordered.median_effective_blocks_fraction_beta2) > 0).all() and (np.diff(ordered.relative_primary_loss_reduction) < 0).all()
        out['interpretation_status'] = 'consistent_with_hypothesis_descriptive' if consistent else 'regional_concentration_does_not_fully_explain_gain'
        a = out.sort_values('median_effective_blocks_fraction_beta2').trait.tolist()
        b = out.sort_values('median_effective_blocks_fraction_abs_beta').trait.tolist()
        if a != b:
            out['interpretation_status'] += ';sensitivity_order_disagrees'
    return out


def plot(populations, gain, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    traits = list(gain.trait)
    labels = {'height':'Height', 'ldl':'LDL', 't2dm':'T2DM'}
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    primary = populations[populations.mass_definition.eq('beta_squared')]
    for ax, metric, title in zip(axes[:2], ('effective_blocks_fraction', 'top5pct_share'),
                                 ('A  Effective regional dimensionality', 'B  Top 5% block mass share')):
        for j, pop in enumerate(POPS):
            y = [primary[(primary.trait.eq(t)) & (primary.population.eq(pop))][metric].iloc[0] for t in traits]
            ax.scatter(np.arange(len(traits))+(j-1.5)*.09, y, label=pop, s=25)
        for i,t in enumerate(traits):
            median = primary[primary.trait.eq(t)][metric].median()
            ax.plot([i-.2,i+.2], [median,median], color='black', linewidth=2)
        ax.set_xticks(range(len(traits)), [labels[t] for t in traits])
        ax.set_title(title, fontsize=10)
        ax.set_ylabel('Effective blocks / all blocks' if metric.startswith('effective') else 'Fraction of assigned posterior mass')
    axes[0].legend(fontsize=7, frameon=False)
    axes[2].set_title('C  Dimensionality versus GRID gain', fontsize=10)
    for row in gain.itertuples():
        x,y = row.median_effective_blocks_fraction_beta2,100*row.relative_primary_loss_reduction
        if np.isfinite([x,y]).all():
            axes[2].scatter(x,y,s=30)
            axes[2].annotate(labels[row.trait],(x,y),xytext=(4,5),textcoords='offset points',fontsize=8)
    if gain.relative_primary_loss_reduction.isna().all():
        axes[2].text(.5,.5,'Gain not available',ha='center',transform=axes[2].transAxes)
    axes[2].margins(.25)
    axes[2].set_xlabel('Median effective blocks / all blocks')
    axes[2].set_ylabel('Relative primary loss reduction (%)')
    fig.text(.5,.02,'Three-trait comparison is descriptive; regional dimensionality is calculated from PRS-CSx posterior weights\nindependently of phenotype outcomes. Primary mass: beta squared.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.12,1,1))
    for ext in ('png','pdf'):
        fig.savefig(out/('regional_dimensionality.'+ext),dpi=220)
    plt.close(fig)


def validate_tables(populations, summary, top, gain, traits):
    if len(populations) != len(traits)*len(POPS)*2 or populations.duplicated(['trait','population','mass_definition']).any():
        raise ValueError('Incomplete or duplicate population metrics')
    if len(summary) != len(traits)*2 or len(gain) != len(traits):
        raise ValueError('Incomplete trait summaries')
    valid = populations[populations.status.eq('completed')]
    for col in ('assigned_snp_fraction','assigned_mass_fraction','top1pct_share','top5pct_share','top10pct_share','effective_blocks_fraction','shannon_effective_blocks_fraction','gini'):
        if not valid[col].between(-1e-12,1+1e-12).all():
            raise ValueError('Invalid concentration range: '+col)
    if not np.isfinite(valid.effective_blocks).all() or (valid.effective_blocks < 1-1e-12).any():
        raise ValueError('Invalid effective number')
    if not (top.groupby(['trait','population','mass_definition']).size() <= 25).all():
        raise ValueError('Invalid top block table')


def completed_tables(root):
    """Read validated completed aggregates; never assign blocks in final."""
    out = Path(root)/'architecture'
    marker = out/'ARCHITECTURE_COMPLETE.json'
    if not marker.is_file():
        return None
    receipt = json.loads(marker.read_text())
    if receipt.get('schema') != SCHEMA or receipt.get('status') != 'completed' or receipt.get('mass_definitions') != list(MASS_DEFINITIONS):
        raise ValueError('Invalid architecture completion')
    for record in receipt['inputs'] + list(receipt['output_identities'].values()):
        if identity(record['resolved_path']) != record:
            raise ValueError('Architecture input/output changed: '+Path(record['resolved_path']).name)
    tables = [pd.read_csv(out/name, sep='\t') for name in TABLES]
    validate_tables(*tables, receipt['traits'])
    return tables


def paper_table(root, traits):
    columns = ['trait','primary_mass_definition','median_effective_blocks_fraction',
               'q25_effective_blocks_fraction','q75_effective_blocks_fraction','median_top5pct_share',
               'median_hhi','median_gini','primary_loss','relative_primary_loss_reduction','interpretation']
    tables = completed_tables(root)
    if tables is None:
        return pd.DataFrame(columns=columns)
    _, summary, _, gain = tables
    d = summary[summary.mass_definition.eq('beta_squared') & summary.trait.isin(traits)].copy()
    d = d.rename(columns={'mass_definition':'primary_mass_definition'})
    d = d.merge(gain[['trait','primary_loss','relative_primary_loss_reduction','interpretation_status']],on='trait',validate='one_to_one')
    return d.rename(columns={'interpretation_status':'interpretation'})[columns]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-root',default='/mnt/d/analysis/grid')
    p.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config/abm.blocks.json'))
    p.add_argument('--traits',default='height,ldl,t2dm')
    p.add_argument('--dry-run',action='store_true')
    a = p.parse_args(argv)
    traits = [t.strip() for t in a.traits.split(',')]
    if not traits or len(set(traits)) != len(traits) or not set(traits) <= {'height','ldl','t2dm'}:
        p.error('Use unique traits from height,ldl,t2dm')
    if a.dry_run:
        print('Architecture: posterior weights + frozen LD maps; aggregate outputs only')
        return
    config = read_config(a.config)
    inputs = [identity(a.config)]
    layout = resolve_block_files(config['ld_block_dir'], config['ld_block_files'], lightweight=True)
    inputs += [identity(f) for f in layout['files']]
    manifest, _, _ = build_manifest(layout, config)
    out = Path(a.analysis_root)/'architecture'
    out.mkdir(parents=True,exist_ok=True)
    marker = out/'ARCHITECTURE_COMPLETE.json'
    marker.unlink(missing_ok=True)
    rows, top = [], []
    for trait in traits:
        inputs.append(identity(Path(config['csx_analysis_root'])/trait/'provenance.json'))
        files, _ = locate_weights(config, trait, lightweight=True)
        for i,pop in enumerate(POPS):
            inputs.append(identity(files[pop]))
            weights = read_weights(files[pop])
            bed = manifest[manifest.source_population.eq(pop if layout['mode']=='population_multiview' else 'ALL')].copy()
            for definition in MASS_DEFINITIONS:
                metrics, blocks = population_metrics(weights, bed, definition)
                metadata = dict(trait=trait,population=pop,mass_definition=definition)
                rows.append(dict(**metadata,genome_build=37,ld_layout=layout['layout_id'],
                                 block_file=layout['files'][i if len(layout['files'])==4 else 0].name,
                                 weight_file_basename=files[pop].name,**metrics))
                top.append(blocks.assign(**metadata))
            print(f'[architecture] {trait}/{pop}: {len(weights)} valid SNPs, {len(bed)} blocks',flush=True)
    populations = pd.DataFrame(rows)
    summary = summarize(populations)
    top = pd.concat(top,ignore_index=True)
    # Only now read the optional, already aggregated primary performance table.
    perf = Path(a.analysis_root)/'final/paper_main_performance.tsv'
    performance = pd.read_csv(perf,sep='\t') if perf.is_file() else pd.DataFrame()
    gain = join_gain(summary,performance)
    tables = (populations, summary, top, gain)
    validate_tables(*tables,traits)
    for record in inputs:
        if identity(record['resolved_path']) != record:
            raise ValueError('Input changed during analysis')
    for name,table in zip(TABLES,tables):
        temp=out/(name+'.tmp')
        table.to_csv(temp,sep='\t',index=False)
        temp.replace(out/name)
    plot(populations,gain,out)
    # Validate serialized TSVs as well, before publishing completion.
    validate_tables(*[pd.read_csv(out/name,sep='\t') for name in TABLES],traits)
    outputs = [*TABLES,'regional_dimensionality.png','regional_dimensionality.pdf']
    receipt = dict(schema=SCHEMA,traits=traits,populations=list(POPS),mass_definitions=list(MASS_DEFINITIONS),
                   primary_mass_definition='beta_squared',genome_build=37,block_layout_id=layout['layout_id'],
                   inputs=inputs,weight_basenames=populations.weight_file_basename.unique().tolist(),
                   block_filenames=[f.name for f in layout['files']],output_filenames=outputs,
                   output_identities={name:identity(out/name) for name in outputs},
                   gain_source_identity=identity(perf) if perf.is_file() else None,status='completed')
    temp=out/'ARCHITECTURE_COMPLETE.json.tmp'
    temp.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    temp.replace(marker)
    print(gain.to_string(index=False),flush=True)


if __name__ == '__main__':
    main()
