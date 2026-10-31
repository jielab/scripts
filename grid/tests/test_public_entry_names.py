from pathlib import Path
import shutil
import subprocess
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1]


def test_public_files_and_references():
    names=['1.pca.sh','2.csx.sh','3.disco.sh','4.prsformer.sh','f/1.pca.R','f/1.pca.py','f/2.csx.py','f/3.disco.R']
    assert all((ROOT/n).is_file() for n in names)
    old=[f'f/{n}.{method}.{ext}' for n,method,ext in [(0,'pca','R'),(0,'pca','py'),(1,'csx','py'),(2,'disco','R')]]
    assert not any((ROOT/n).exists() for n in old)
    stale=[f'{n}.{method}.sh' for n,method in [(0,'pca'),(1,'csx'),(2,'disco'),(3,'prsformer')]]+old
    for p in ROOT.rglob('*'):
        if p.is_file() and p.suffix in {'.py','.sh','.R','.md'}:
            text=p.read_text(encoding='utf-8-sig')
            assert not any(s in text for s in stale),str(p)


@pytest.mark.parametrize('script',['1.pca.sh','2.csx.sh','3.disco.sh','4.prsformer.sh','grid.sh'])
def test_help(script):
    if not shutil.which('bash'):pytest.skip('bash unavailable')
    subprocess.run(['bash',str(ROOT/script),'--help'],check=True,capture_output=True)


def test_compile():
    for name in ['1.pca.py','2.csx.py','regional_dimensionality.py']:
        compile((ROOT/'f'/name).read_text(encoding='utf-8-sig'),name,'exec')


@pytest.mark.parametrize('name',['1.pca.R','3.disco.R'])
def test_r_parse(name):
    if not shutil.which('Rscript'):pytest.skip('R unavailable')
    subprocess.run(['Rscript','-e',f"parse(file='{(ROOT/'f'/name).as_posix()}')"],check=True,capture_output=True)
