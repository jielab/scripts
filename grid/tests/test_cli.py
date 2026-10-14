from pathlib import Path
import sys,subprocess,json
import numpy as np,pytest
ROOT=Path(__file__).resolve().parents[1]


def run_cli(*args):return subprocess.run([sys.executable,str(ROOT/'f/twin_entry.py'),*args],capture_output=True,text=True)


def test_help_runs_without_genotypes_or_R():
    r=run_cli('--help');assert r.returncode==0 and 'genetic-profile twin' in r.stdout


def test_public_abm_dry_run():
    r=run_cli('abm','--traits','height','--dry-run');assert r.returncode==0
    assert 'twin_cli.py run' in r.stdout and 'height' in r.stdout


def test_sequence_rejects_shiny_not_last():
    r=run_cli('shiny,abm','--dry-run');assert r.returncode!=0


def test_sequence_dry_run_no_model_or_R():
    r=run_cli('abm,final,shiny','--traits','height','--dry-run');assert r.returncode==0
    assert 'browser_index.py' in r.stdout and 'app.R' in r.stdout


def test_no_silent_old_arguments():
    r=run_cli('abm','--abm-backend','selective_attention');assert r.returncode!=0


def test_explicit_projection_rejects_outcomes_before_predict(tmp_path):
    # The direct parser keeps no pathway for hidden test label fields.
    sys.path.insert(0,str(ROOT/'f'))
    from twin_cli import npz
    f=tmp_path/'bad.npz';np.savez(f,X=np.zeros((2,2)),ids=np.array(['a','b']),groups=np.array(['a','b']),features=np.array(['x','z']),y=np.array([0,1]))
    assert set(npz(f))!={'X','ids','groups','features'}


def test_legacy_logger_uses_preserved_entry(tmp_path,monkeypatch):
    sys.path.insert(0,str(ROOT/'f'))
    import twin_entry
    root=tmp_path/'grid';root.mkdir();(root/'grid.upstream.sh').touch()
    logger=tmp_path/'0f/console_run.py';logger.parent.mkdir();logger.touch()
    monkeypatch.setattr(twin_entry,'ROOT',root)
    monkeypatch.delenv('SCRIPT_CONSOLE_ACTIVE',raising=False)
    monkeypatch.delenv('SCRIPT_VERBOSE',raising=False)
    calls=[]
    monkeypatch.setattr(twin_entry.subprocess,'run',lambda cmd,**kwargs:calls.append((cmd,kwargs)))
    twin_entry.main(['legacy-abm','--check'])
    cmd,options=calls[0]
    assert cmd==[sys.executable,str(logger),'--script',str(root/'grid.upstream.sh'),'--','abm','--check']
    assert 'SCRIPT_LOG_DIR' in options['env']


def test_legacy_without_logger_cannot_reenter_new_dispatch(tmp_path,monkeypatch):
    sys.path.insert(0,str(ROOT/'f'))
    import twin_entry
    (tmp_path/'grid.upstream.sh').touch();monkeypatch.setattr(twin_entry,'ROOT',tmp_path)
    calls=[]
    monkeypatch.setattr(twin_entry.subprocess,'run',lambda cmd,**kwargs:calls.append((cmd,kwargs)))
    twin_entry.main(['legacy-final'])
    assert calls[0][0][2]=='final'
    assert calls[0][1]['env']['SCRIPT_CONSOLE_ACTIVE']=='1'
