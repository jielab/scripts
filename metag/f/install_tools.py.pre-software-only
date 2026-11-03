"""Install/verify software on a compute node, without preparing raw data."""
import argparse
from pathlib import Path
from hpc_common import load_config, require_compute, CODE_ROOT
from hpc_lsf import write_cmd, submit_cmd
from hpc_cli import submit_bool

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--submit', type=submit_bool, default=False, metavar='TRUE|FALSE')
p.add_argument('--config', help='HPC YAML configuration (default hpc.yml)')
p.add_argument('--root', help='Override root for relative configured paths')
p.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
a = p.parse_args()
cfg = load_config(a.config, a.root)
if a.worker:
    require_compute()
    from hpc_install import install
    from hpc_kunpeng_index import install as install_kunpeng
    install(cfg)
    install_kunpeng(cfg)
    print('SOFTWARE_INSTALLATION_VERIFIED')
else:
    command = ['bash', str(CODE_ROOT / 'install.sh'), '--worker']
    if a.config:
        command += ['--config', str(Path(a.config).resolve())]
    if a.root:
        command += ['--root', str(Path(a.root).resolve())]
    path = write_cmd(cfg, Path(cfg['work_dir']) / 'install.cmd', 'setup',
                     command, 'metag.install')
    print(path)
    if a.submit:
        print(submit_cmd(cfg, path))
