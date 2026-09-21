#!/usr/bin/env python3
"""Synthetic end-to-end checks, including independent OOF metric verification.
Run after activating the GRID environment: python3 f/tests/test_yeval.py
Temporary inputs/outputs are removed automatically.
"""
import os
from html.parser import HTMLParser
from pathlib import Path
import subprocess
import tempfile
from urllib.parse import urlsplit, unquote
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


class ReportLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.files = set()

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('src', 'href') and value:
                url = urlsplit(value)
                if not url.scheme and not url.netloc and url.path:
                    self.files.add(unquote(url.path))


def check_report(folder):
    links = ReportLinks()
    links.feed((folder/'report.html').read_text())
    assert len(links.files) == 9
    for name in links.files:
        assert (folder/name).is_file(), f'Broken report link: {name}'
        assert (folder/name).stat().st_size > 0
    expected = links.files | {'report.html', 'eval.log', 'run.lock'}
    assert {f.name for f in folder.iterdir()} == expected


def main():
    rng = np.random.default_rng(91)
    with tempfile.TemporaryDirectory(prefix='yeval-test-') as tmp:
        p = Path(tmp); n = 1200
        groups = ['EUR', 'AFR', 'EAS', 'SAS']; pop = np.repeat(groups, 300)
        iid = np.arange(1, n + 1).astype(str); scores = rng.normal(size=(n, 4))
        centers = rng.normal(0, 4, (4, 10)); pcs = rng.normal(size=(n, 10)) + np.repeat(centers, 300, axis=0)
        x = pd.DataFrame(dict(eid=iid, genetic_ancestry=pop, age=rng.normal(55, 6, n),
                              sex=rng.integers(0, 2, n), PC1=pcs[:, 0], PC2=pcs[:, 1]))
        x['y'] = .3 * scores[:, 0] + .5 * scores[:, 1] + rng.normal(size=n) + .03 * x.age
        x['d'] = rng.binomial(1, 1 / (1 + np.exp(-(.8 * scores[:, 0] + .3 * scores[:, 1] - .4))))
        x['e'] = rng.binomial(1, .6, n); x['time'] = rng.exponential(5, n) + .1
        x.to_csv(p / 'phe.tsv', sep='\t', index=False)
        csx = pd.DataFrame(dict(eid=iid, **{f'csx.{z}':scores[:, j] for j, z in enumerate(['AFR','EAS','EUR','SAS'])},
                               **{'csx.auto':scores.mean(axis=1), 'csx.meta':scores.mean(axis=1) + rng.normal(0,.1,n)}))
        csx.to_csv(p / 'csx.tsv', sep='\t', index=False)
        pd.DataFrame(dict(eid=iid, disco=scores[:,0]*.3+scores[:,1]*.5)).sample(frac=1,random_state=1).to_csv(p/'disco.tsv',sep='\t',index=False)
        pd.DataFrame(dict(eid=iid, **{f'pt.{z}':scores[:,j]+rng.normal(size=n) for j,z in enumerate(['AFR','EAS','EUR','SAS'])})).to_csv(p/'pt.tsv',sep='\t',index=False)
        pd.DataFrame(dict(eid=iid, **{f'PC{j+1}':pcs[:,j] for j in range(10)})).to_csv(p/'pca.tsv',sep='\t',index=False)
        pd.DataFrame(dict(POP=groups, **{f'PC{j+1}':centers[:,j] for j in range(10)})).to_csv(p/'centers.tsv',sep='\t',index=False)
        common = ['--pheno-file',str(p/'phe.tsv'),'--pgs-file',str(p/'csx.tsv'),
                  '--disco-file',str(p/'disco.tsv'),'--pt-file',str(p/'pt.tsv'),
                  '--pca-file',str(p/'pca.tsv'),'--med-file',str(p/'centers.tsv'),
                  '--remove','','--out-root',str(p/'out'),'--bootstrap','20','--min-n','30']
        for kind, extra in [('ct',['--phenotype-col','y']),
                            ('dt',['--phenotype-col','d','--prevalence','EUR=0.05,AFR=0.1,EAS=0.08,SAS=0.16']),
                            ('t2e',['--event-col','e','--time-col','time'])]:
            cmd = [str(ROOT/'Yeval.sh'),'--trait',kind,'--type',kind,*common,*extra]
            result = subprocess.run(cmd,capture_output=True,text=True)
            assert result.returncode == 0, result.stdout + result.stderr
            folder = p/'out'/kind
            perf = pd.read_csv(folder/'performance.tsv',sep='\t')
            assert len(perf) == 36 and np.isfinite(perf.estimate).all()
            check_report(folder)
            assert (folder/'plots.pdf').stat().st_size > 1000
            assert not (folder/'csx').exists()
            if kind == 'ct':
                # Reconstruct the first ancestry's seeded folds independently;
                # no production prediction/fold exports are needed for this check.
                r_env = {k:v for k,v in os.environ.items() if k not in ('R_ENVIRON_USER','R_LIBS_USER')}
                assigned = subprocess.check_output(['/usr/bin/Rscript','-e',
                    'set.seed(20260904);cat(sample(rep(1:5,length.out=300)),sep="\\n")'],env=r_env,text=True)
                input_data = x.merge(csx,on='eid')
                for name in ('pt', 'disco'):
                    input_data = input_data.merge(pd.read_csv(p/f'{name}.tsv',sep='\t',dtype={'eid':str}),on='eid')
                input_data = input_data[input_data.genetic_ancestry=='EUR'].sort_values('eid').reset_index(drop=True)
                input_data['fold'] = np.fromstring(assigned,sep='\n',dtype=int)
                base_columns = ['age','sex','PC1','PC2']
                models = {'PT':['pt.EUR'], 'PRS-CS-multi':['csx.auto'],
                          'PRS-CSX':[f'csx.{z}' for z in ['AFR','EAS','EUR','SAS']],
                          'DiscoDivas':['disco'], 'CSx-meta':['csx.meta']}
                models.update({f'csx.{z}':[f'csx.{z}'] for z in groups})
                def predict(columns):
                    values = np.empty(len(input_data))
                    for fold in range(1,6):
                        train = input_data[input_data.fold!=fold]; test = input_data[input_data.fold==fold]
                        design = lambda z: np.column_stack([np.ones(len(z)),z[columns].to_numpy()])
                        b = np.linalg.lstsq(design(train),train.y.to_numpy(),rcond=None)[0]
                        values[test.index] = design(test) @ b
                    return values
                baseline = predict(base_columns)
                for method, columns in models.items():
                    prediction = predict(base_columns + columns)
                    expected = 1-((input_data.y-prediction)**2).sum()/((input_data.y-baseline)**2).sum()
                    actual = perf.loc[(perf.target=='EUR')&(perf.method==method),'estimate'].item()
                    assert abs(actual-expected)<1e-10, (method,actual,expected)
                # A check-only run must leave an existing report and its downloads intact.
                saved = {f:f.read_bytes() for f in folder.iterdir() if f.name!='eval.log'}
                checked = subprocess.run(cmd+['--check'],capture_output=True,text=True)
                assert checked.returncode == 0, checked.stdout + checked.stderr
                for f, content in saved.items():
                    assert f.read_bytes() == content
                check_report(folder)
                # Full reruns remove legacy files, while preserving user-added files.
                for name in ['command.sh','comparison.pdf','combined_scores.pdf','paired_improvement.pdf',
                             'genetic_landscape.pdf','distance_performance.pdf','distance_performance.tsv',
                             'fold_coefficients.tsv','folds.tsv.gz','genetic_distance.tsv.gz','manifest.tsv',
                             'methods.tsv','paired_comparison.tsv','predictions.tsv.gz','prevalence.tsv',
                             'skipped.tsv','pt.commands.jsonl','pt.log','pt.matches.tsv','pt.plink.log','pt.variants.tsv']:
                    (folder/name).write_text('legacy')
                custom = folder/'user-notes.tsv'; custom.write_text('keep me')
                rerun = subprocess.run(cmd,capture_output=True,text=True)
                assert rerun.returncode == 0, rerun.stdout + rerun.stderr
                assert custom.read_text() == 'keep me'; custom.unlink()
                check_report(folder)
                pd.testing.assert_frame_equal(perf,pd.read_csv(folder/'performance.tsv',sep='\t'))
            if kind=='dt':
                assert set(perf.metric)=={'OOF_liability_R2'}
                assert set(perf.K)=={.05,.1,.08,.16}
                assert perf.AUC.between(0,1).all()
                assert 'Prevalence assumptions</h2><table>' in (folder/'report.html').read_text()
            if kind=='t2e':
                assert set(perf.metric)=={'OOF_Harrell_C'} and perf.estimate.between(0,1).all()
        csx.drop(columns=['csx.auto']).to_csv(p/'csx.tsv',sep='\t',index=False)
        bad = subprocess.run([str(ROOT/'Yeval.sh'),'--trait','missing','--type','ct',*common,'--phenotype-col','y','--check'],capture_output=True,text=True)
        assert bad.returncode != 0 and 'Missing CSx columns: csx.auto' in bad.stderr+bad.stdout
        partial = subprocess.run([str(ROOT/'Yeval.sh'),'--trait','partial','--type','ct',*common,
                                  '--phenotype-col','y','--allow-missing-scores'],capture_output=True,text=True)
        assert partial.returncode == 0, partial.stdout + partial.stderr
        check_report(p/'out'/'partial')
        report = (p/'out'/'partial'/'report.html').read_text()
        assert 'PARTIAL REPORT' in report and 'Skipped evaluations' in report and 'PRS-CS-multi' in report
        print('PASS: continuous/binary/survival, independent OOF reconstruction, report links, minimal outputs, legacy cleanup, check-only preservation and partial reports.')


if __name__=='__main__':
    main()
