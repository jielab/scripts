"""Integrated ABM must use the requested cleaned cohort and adjustment schema."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import cohort_inputs as host
from host_inputs import prepare_host
from run_host import check_input_reuse, parser


@pytest.fixture(autouse=True)
def clean_shared_environment(monkeypatch):
    for key in ('LE8_VARS_ADJ', 'LE8_GROUP_FILE', 'LE8_GROUP_COLUMN', 'PGS_GROUP_COLUMN',
                'LE8_OUTER_ROSTER', 'DATE_FOLLOW_END', 'LE8_Y_DATE'):
        monkeypatch.delenv(key, raising=False)


@pytest.mark.parametrize('layer', ['prot', 'met'])
def test_defaults_use_cleaned_named_rds_without_technical_categories(tmp_path, layer):
    a = host.configure(host.parser().parse_args(['--biom', layer, '--ukb-phe', str(tmp_path)]))
    assert a.omics_file == str(tmp_path / f'Rdata/{layer}.rds')
    assert a.met_input == 'named' and a.transform == 'none'
    assert a.covariates == 'age,sex,tdi,PC1,PC2'
    assert a.categorical == 'sex' and a.residualize == ''
    assert not any(host.technical_category(c) for c in host.phenotype_columns(a))


def test_shared_and_explicit_technical_categories_are_excluded(monkeypatch, capsys):
    monkeypatch.setenv('LE8_VARS_ADJ', 'age,sex,center,met.plate,prot.well,PC1')
    a = host.configure(host.parser().parse_args([
        '--biom', 'met', '--residualize', 'assay_plate,assay_well',
        '--categorical', 'sex,assessment_center,met.plate,well']))
    assert a.covariates == 'age,sex,PC1'
    assert a.residualize == '' and a.categorical == 'sex'
    assert not any(host.technical_category(c) for c in host.phenotype_columns(a))
    assert 'technical_covariates' in capsys.readouterr().out
    b = host.configure(host.parser().parse_args(['--covariates', 'age,center,met.plate,well']))
    assert b.covariates == 'age'
    assert not host.technical_category('wellbeing_score')


def test_raw_metabolite_mapping_requires_explicit_selection(tmp_path):
    a = host.configure(host.parser().parse_args([
        '--biom', 'met', '--omics-file', str(tmp_path/'raw.tsv'), '--met-input', 'raw']))
    assert a.met_input == 'raw' and a.transform == 'log1p'
    assert a.residualize == '' and a.categorical == 'sex'


@pytest.mark.parametrize('layer', ['prot', 'met'])
def test_integrated_preparation_preserves_named_values_and_excludes_metadata(tmp_path, monkeypatch, layer):
    cfg = host.configure(host.parser().parse_args([
        '--biom', layer, '--ukb-phe', str(tmp_path), '--max-samples', '1200', '--demo-features', '6']))
    phenotype, molecules = host.generate_demo(cfg)
    phenotype = phenotype.drop(columns=['center', f'{layer}.plate'])
    molecules.iloc[0, 1] = -2.5  # Named matrices may be supplied on signed scales.
    original = molecules.drop(columns='eid').to_numpy(float).copy()
    for col in ('center', f'{layer}.plate', f'{layer}.well'):
        molecules[col] = np.arange(len(molecules))
    reads = []

    def read(path, id_col='eid', columns=None, r_bin=None, **kwargs):
        reads.append(str(path))
        if Path(path) == tmp_path/'Rdata/all.rds':
            assert set(columns) <= set(phenotype)
            return phenotype[columns].copy()
        assert Path(path) == tmp_path/f'Rdata/{layer}.rds'
        return molecules.copy()

    monkeypatch.setattr(host, 'read_table', read)
    def no_mapping(*args):
        raise AssertionError('Cleaned RDS must not be mapped as raw fields')
    monkeypatch.setattr(host, 'map_metabolites', no_mapping)
    args = parser().parse_args([
        '--Y', 'cvd_cad', '--biom', layer, '--analysis-root', str(tmp_path/'results'),
        '--config', '{}', '--ukb-phe', str(tmp_path)])
    abm, deviation, audit = prepare_host(args, tmp_path/'work')
    np.testing.assert_array_equal(deviation['M'], original)
    positions = {eid: i for i, eid in enumerate(molecules.eid)}
    np.testing.assert_array_equal(abm['X'][:, :6], original[[positions[eid] for eid in abm['eid']]])
    assert not any(host.technical_category(c) for c in abm['feature_names'])
    assert not any(host.technical_category(c) for c in deviation['covariate_names'])
    assert len(abm['clinical_feature_names']) == 6  # Four numeric columns plus two sex levels.
    assert reads == [str(tmp_path/'Rdata/all.rds'), str(tmp_path/f'Rdata/{layer}.rds')]
    assert audit['input_config']['transform'] == 'none'


def test_old_input_contract_requires_explicit_refit_without_mutation(tmp_path):
    a = host.configure(host.parser().parse_args(['--biom', 'met']))
    contract = tmp_path/'host_contract.json'
    previous = dict(vars(a), omics_file='/raw/met.tab.gz', met_input='raw',
                    transform='log1p', residualize='met.plate', categorical='sex,center,met.plate')
    saved = json.dumps({'input_config': previous})
    contract.write_text(saved)
    with pytest.raises(ValueError, match='c1_abm --replace TRUE'):
        check_input_reuse(tmp_path, a)
    assert contract.read_text() == saved
    contract.write_text(json.dumps({'input_config': vars(a)}))
    check_input_reuse(tmp_path, a)


def test_preflight_rejects_old_contract_before_analysis(tmp_path):
    from run_host import main, finish
    inputs = tmp_path/'Rdata'
    inputs.mkdir()
    for name in ('all.rds', 'met.rds'):
        (inputs/name).touch()  # Preflight checks paths without loading cohort data.
    config = tmp_path/'abm.json'
    config.write_text('{}')
    args = parser().parse_args([
        '--Y', 'cvd_cad', '--biom', 'met', '--analysis-root', str(tmp_path/'results'),
        '--ukb-phe', str(tmp_path), '--config', str(config), '--preflight'])
    a = host.configure(host.parser().parse_args(['--biom', 'met', '--ukb-phe', str(tmp_path)]))
    target = args.analysis_root/'cvd_cad/met/c1_abm/abm_twin'
    target.mkdir(parents=True)
    finish(target, {'input_config': dict(vars(a), covariates='age,sex,center')})
    before = (target/'COMPLETE.json').read_bytes()
    with pytest.raises(ValueError, match='covariates'):
        main(args)
    args.replace = 'TRUE'
    main(args)
    assert (target/'COMPLETE.json').read_bytes() == before
