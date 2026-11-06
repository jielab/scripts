"""ABM restart orchestration, separate from the statistical implementation."""
from dataclasses import asdict
from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from resume import Checkpoints, choose_run, run_lock, sha, validate_manifests


def canonical(value):
    return json.loads(json.dumps(value, default=str))


def contract_for(api, args, cfg, routes, raw, query, truth, extra, endpoint):
    return dict(format=1, config=canonical(asdict(cfg)), routes=canonical(asdict(routes)),
                data=joblib.hash((raw, query, truth, extra)), endpoint=endpoint,
                bootstrap=args.bootstrap,
                comparators=sha(args.comparators) if args.comparators else None,
                require_comparators=bool(args.require_comparators or routes.publication_mode),
                model_sources={name: sha(api.ROOT / 'f' / name) for name in ('twin_core.py', 'abm_routes.py')})


def validate_legacy(api, dest, contract, raw, query, truth, meta):
    protocol = json.loads((dest / 'protocol.json').read_text())
    source = json.loads((dest / 'SOURCE_FROZEN.json').read_text())
    if (protocol['config'] != contract['config'] or protocol['abm_routes_config'] != contract['routes']
            or protocol['endpoint_definition'] != contract['endpoint']
            or protocol.get('input', {}).get('input_sha256') != meta.get('input_sha256')):
        raise ValueError('Legacy run input/configuration mismatch; use --new-run')
    for name, digest in contract['model_sources'].items():
        if source['files_sha256'].get('f/' + name) != digest:
            raise ValueError('Legacy fitted-model source mismatch: ' + name)
    for name, current in [('development.npz', raw), ('test_X.npz', query)]:
        with np.load(dest / '_LOCAL_ONLY' / name, allow_pickle=False) as saved:
            if set(saved.files) != set(current):
                raise ValueError('Legacy array schema mismatch: ' + name)
            for key, value in current.items():
                a, b = saved[key], np.asarray(value)
                same = np.array_equal(a, b, equal_nan=True) if a.dtype.kind in 'fci' else np.array_equal(a, b)
                if not same:
                    raise ValueError('Legacy input changed: ' + name + ':' + key)
    saved = pd.read_csv(dest / '_LOCAL_ONLY/test_truth.tsv.gz', sep='\t', dtype={'eid': str})
    pd.testing.assert_frame_equal(saved.reset_index(drop=True), truth.reset_index(drop=True),
                                  check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    model_path = dest / 'fit/_LOCAL_ONLY/model.joblib'
    if not model_path.is_file():
        raise ValueError('Legacy run has no saved fit to resume; use --new-run')
    model = joblib.load(model_path)
    marker = json.loads((dest / 'fit/MODEL_FROZEN.json').read_text())
    if api.state_hash(model) != marker['model_hash']:
        raise ValueError('Legacy frozen model hash mismatch')
    prediction_path = dest / '_LOCAL_ONLY/predictions.joblib'
    if prediction_path.is_file():
        p = joblib.load(prediction_path)
        validate_prediction(api, model, p, query['ids'])
    return model_path, prediction_path


def validate_prediction(api, model, prediction, ids):
    if prediction['frozen_hash'] != api.state_hash(model) or not np.array_equal(prediction['eid'], ids):
        raise ValueError('Saved prediction/model/roster mismatch')


def run_trait(api, args, routes, context, comparator_check):
    trait,cfg,d,groups,meta,raw,query,truth,flow,extra,endpoint,prsformer_reuse = context
    root = Path(args.out_root) / 'twin' / trait
    contract = contract_for(api, args, cfg, routes, raw, query, truth, extra, endpoint)
    contract['comparator_check'] = canonical(comparator_check)
    if args.comparators:
        manifest = Path(args.comparators)
        contract['comparator_files'] = {item['path']: sha(manifest.parent/item['path']) for item in json.loads(manifest.read_text())['methods'] if (manifest.parent/item['path']).is_file()}
    with run_lock(root):
        dest = choose_run(root, args.run_id, getattr(args, 'new_run', False))
        receipt = dest / 'RESUME.json'
        legacy = dest.exists() and (not receipt.is_file() or (
            'fit' not in json.loads(receipt.read_text()).get('stages', {}) and
            (dest/'fit/_LOCAL_ONLY/model.joblib').is_file()))
        adopted = validate_legacy(api, dest, contract, raw, query, truth, meta) if legacy else None
        cp = Checkpoints(dest, contract)
        if cp.valid('completed_outputs') and (dest / 'COMPLETE.json').is_file():
            api.json_write(root/'CURRENT.json', dict(run=dest.name, model_hash=json.loads((dest/'COMPLETE.json').read_text())['model_hash']))
            print(f'[GRID twin] SKIP {trait}: completed run and all outputs verified ({dest.name})', flush=True)
            return
        if adopted:
            # Historical runs have table checksums and a frozen model digest.
            # Adopt a fully intact completion once, without repeating evaluation.
            required = [dest / x for x in ('fit/aggregate/table_manifest.json', 'audit/table_manifest.json',
                        'test/table_manifest.json', 'test/status.json', 'test/figures/same_subset.png',
                        'test/figures/effect_decomposition.png', 'test/figures/coverage.png')]
            complete = dest / 'COMPLETE.json'
            if (complete.is_file() and adopted[1].is_file() and args.bootstrap is None and not args.comparators
                    and all(p.is_file() for p in required) and validate_manifests(dest)):
                marker = json.loads(complete.read_text())
                frozen = json.loads((dest/'fit/MODEL_FROZEN.json').read_text())
                if marker.get('status') != 'completed' or marker.get('model_hash') != frozen['model_hash'] or marker.get('test_roster_hash') != api.id_hash(query['ids']):
                    raise ValueError('Legacy completion marker mismatch')
                cp.adopt('fit', adopted[0])
                cp.adopt('prediction', adopted[1])
                cp.mark('completed_outputs', [p for p in dest.rglob('*') if p.is_file() and p.name not in ('RESUME.json', 'FAILED.json')])
                api.json_write(root/'CURRENT.json', dict(run=dest.name,model_hash=marker['model_hash']))
                print(f'[GRID twin] SKIP {trait}: verified legacy completed run ({dest.name})', flush=True)
                return
            cp.adopt('fit', adopted[0])
            if adopted[1].is_file(): cp.adopt('prediction', adopted[1])
            # Reuse legacy aggregate evaluation tables when their checksums pass.
            # Missing PNG/XLSX exports do not justify repeating bootstrap inference.
            if args.bootstrap is None and not args.comparators:
                for section in ('audit', 'test'):
                    folder = dest / section
                    if validate_manifests(folder):
                        manifest = json.loads((folder/'table_manifest.json').read_text())
                        tables = {row['name']: pd.read_csv(folder/row['path'], sep='\t') for row in manifest}
                        if section == 'audit':
                            cp.mark('audit_report', [p for p in folder.rglob('*') if p.is_file()])
                        elif (folder/'status.json').is_file():
                            cp.object('test_evaluation', lambda: (tables,json.loads((folder/'status.json').read_text())))
            # Legacy completion markers did not hash figures. Rebuild reports once
            # from frozen predictions; subsequent runs validate every report file.
        local = dest / '_LOCAL_ONLY'; local.mkdir(exist_ok=True)
        if not (dest / 'protocol.json').exists():
            api.json_write(dest / 'protocol.json', dict(config=asdict(cfg),abm_routes_config=asdict(routes),
                abm_version=api.GRID_ABM_VERSION,prsformer_reuse=prsformer_reuse,comparator_check=comparator_check,
                input=meta,flow=flow,endpoint_definition=endpoint,test_roster_changed=False,validation_scope=cfg.validation_scope))
            api.write_source_provenance(dest, api.ROOT, routes, dict(test_roster_hash=api.id_hash(query['ids']),
                development_roster_hash=api.id_hash(raw['ids']),endpoint_definition=endpoint,
                covariate_definition=flow['clinical_columns'],grouping_status=cfg.grouping_status,
                upstream_gwas_overlap=cfg.upstream_gwas_overlap,validation_scope=cfg.validation_scope))
        try:
            def fit_core():
                print(f'[GRID twin] FIT {trait}: core model', flush=True)
                return api.TwinModel(cfg).fit(**raw)
            def fit_routes():
                core = cp.object('core_fit', fit_core)
                print(f'[GRID twin] FIT {trait}: ABM routes', flush=True)
                return api.fit_abm_routes(core, raw, routes)
            model = cp.object('fit', fit_routes)
            cp.report('fit_report', dest / 'fit', lambda: api.save_fit(model, dest / 'fit'))
            def audit_prediction():
                audit = raw['roles'] == 'audit'
                p = model.project(raw['X'][audit],raw['ids'][audit],raw['groups'][audit],raw['features'],audit=True)
                p = api.apply_abm_routes(model,p,raw['X'][audit],raw['ids'][audit],raw['groups'][audit])
                api.freeze_matched_controls(model,p,extra['audit_truth'].ancestry.to_numpy())
                return p
            if not cp.valid('audit_report'):
                pa = cp.object('audit_prediction', audit_prediction)
                validate_prediction(api, model, pa, raw['ids'][raw['roles'] == 'audit'])
                def audit_evaluation():
                    tables,status = api.evaluate(model,pa,extra['audit_truth'],bootstrap=args.bootstrap)
                    tables.update(api.evaluate_abm_routes(model,pa,extra['audit_truth'],routes,bootstrap=args.bootstrap))
                    return tables,status
                at,_ = cp.object('audit_evaluation', audit_evaluation)
                cp.report('audit_report', dest / 'audit', lambda: api.write_aggregate(dest/'audit',at,
                    dict(trait=trait,module='Twin audit',run_id=pa['frozen_hash'])))
            def predict():
                p = model.project(**query)
                p = api.apply_abm_routes(model,p,query['X'],query['ids'],query['groups'])
                api.freeze_matched_controls(model,p,truth.ancestry.to_numpy())
                return p
            p = cp.object('prediction', predict)
            validate_prediction(api, model, p, query['ids'])
            if not all(f'{routes.predictability_selector}:{q:g}' in p['masks'] for q in routes.coverages):
                raise ValueError('Incomplete frozen primary coverage masks')
            def comparators():
                ext,availability = api.benchmark_manifest(args.comparators,p['eid'],truth,trait,endpoint)
                disco,disco_model,da = api.calibrate_disco(extra,model); availability.append(da)
                if disco is not None: ext['Disco_reference'] = disco
                if args.analysis_root and contract['require_comparators']:
                    cov = meta.get('prepared_settings',{}).get('covariates','')
                    cov = cov.split(',') if isinstance(cov,str) and cov else list(groups.get('covariates',[]))
                    pp,stat = api.verified_prsformer(args.analysis_root,trait,d,truth,cfg,cov,endpoint)
                    availability.append(stat)
                    if pp is not None: ext['PRSformer'] = pp
                elif args.analysis_root:
                    availability.append(dict(method='PRSformer',status='comparison_deferred_independent_ABM_run'))
                return ext,availability,disco_model
            ext,availability,disco_model = cp.object('comparators', comparators)
            def evaluate_test():
                tables,status = api.evaluate(model,p,truth,external=ext,bootstrap=args.bootstrap)
                tables.update(api.evaluate_abm_routes(model,p,truth,routes,bootstrap=args.bootstrap))
                tables['comparator_provenance'] = pd.DataFrame(availability)
                return tables,status
            tables,status = cp.object('test_evaluation', evaluate_test)
            def report():
                api.write_aggregate(dest/'test',tables,dict(trait=trait,module='Twin test',run_id=p['frozen_hash']))
                api.json_write(dest/'test/status.json',status)
                api.draw_report(tables,dest/'test/figures')
            cp.report('test_report', dest/'test', report)
            # Preserve the existing downstream artifact paths, independently of reports.
            def exports():
                np.savez_compressed(local/'development.npz',**raw)
                np.savez_compressed(local/'test_X.npz',**query)
                truth.to_csv(local/'test_truth.tsv.gz',sep='\t',index=False,float_format='%.17g')
                pd.DataFrame(dict(eid=raw['ids'],group=raw['groups'],role=raw['roles'])).to_csv(local/'development_roles.tsv.gz',sep='\t',index=False)
                joblib.dump(p,local/'predictions.joblib',compress=3)
                api.export_projection(p,local/'confidence.tsv.gz')
                if ext: joblib.dump(dict(eid=p['eid'],predictions=ext),local/'comparator_predictions.joblib',compress=3)
                if disco_model is not None: joblib.dump(disco_model,local/'disco_calibration.joblib',compress=3)
            if not cp.valid('exports'):
                exports()
                cp.mark('exports', [local/n for n in ('development.npz','test_X.npz','test_truth.tsv.gz',
                    'development_roles.tsv.gz','predictions.joblib','confidence.tsv.gz')] +
                    ([local/'comparator_predictions.joblib'] if ext else []) +
                    ([local/'disco_calibration.joblib'] if disco_model is not None else []))
            api.json_write(dest/'COMPLETE.json',dict(status='completed',model_hash=p['frozen_hash'],test_roster_hash=api.id_hash(query['ids']),
                clinical_release=False,abm_version=api.GRID_ABM_VERSION,source_provenance='SOURCE_FROZEN.json',
                frozen_roster_reused=bool(prsformer_reuse.get('reused',False)),prsformer_reused='PRSformer' in ext,
                resume_receipt='RESUME.json'))
            cp.mark('completed_outputs', [q for q in dest.rglob('*') if q.is_file() and q.name not in ('RESUME.json', 'FAILED.json')])
            api.json_write(root/'CURRENT.json',dict(run=dest.name,model_hash=p['frozen_hash']))
            (dest/'FAILED.json').unlink(missing_ok=True)
            print(f'[GRID twin] COMPLETE {trait}: {dest}',flush=True)
        except BaseException:
            api.json_write(dest/'FAILED.json',dict(status='failed_resumable_run',resume_receipt='RESUME.json'))
            raise
