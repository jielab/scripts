"""Install/verify software on a compute node, without preparing raw data."""
import argparse
from pathlib import Path
from hpc_common import load_config, require_compute, CODE_ROOT
from hpc_lsf import write_cmd, submit_cmd

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--submit', action='store_true')
p.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
a = p.parse_args()
cfg = load_config()
if a.worker:
    require_compute()
    from hpc_install import install
    from hpc_kunpeng_index import install as install_kunpeng
    install(cfg)
    install_kunpeng(cfg)
    print('SOFTWARE_INSTALLATION_VERIFIED')
else:
    path = write_cmd(cfg, Path(cfg['work_dir']) / 'install.cmd', 'setup',
                     ['bash', str(CODE_ROOT / 'install.sh'), '--worker'], 'metag.install')
    print(path)
    if a.submit:
        print(submit_cmd(cfg, path))
