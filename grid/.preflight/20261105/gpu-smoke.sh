#!/usr/bin/env bash
set -uo pipefail
cd /work/sph-huangj/scripts/grid
source ../0f/workflow_env.sh
workflow_init grid
export PATH="$GRID_CONDA_ENV/bin:$PATH" PYTHONPATH="$PWD/f:$PWD/../0f" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
AUDIT=/work/sph-huangj/scripts/grid/.preflight/20261105
hostname
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
python3 - <<'PY'
import torch,natten
print('torch',torch.__version__,'CUDA',torch.version.cuda,'NATTEN',natten.__version__)
print('gpu',torch.cuda.get_device_name(), 'capability',torch.cuda.get_device_capability())
x=torch.randn(8,8,device='cuda',requires_grad=True)
(x @ x.T).sum().backward()
print('CUDA matrix forward/backward PASS')
PY
matrix_rc=$?
python3 f/prsformer.py runtime --genotypes /dev/null --data /dev/null --variants /dev/null --upstream-dir /work/sph-huangj/software/PRSformer --out-dir "$AUDIT/gpu-model" --device cuda --attention neighborhood --amp fp16 --kernel-size 385 --threads 2
prsformer_rc=$?
python3 - <<'PY'
import numpy as np
from count_abm.core import count_predict
r = count_predict(np.zeros((8, 2)), np.zeros((2, 2)), donor_ids=[str(i) for i in range(8)], query_ids=["q1", "q2"], donor_groups=[str(i) for i in range(8)], query_groups=["a", "b"], donor_y=np.arange(8.), task="continuous", radius=1., min_count=1, min_group_ess=1, device="cuda")
assert np.allclose(r.table.p_population.to_numpy(), 3.5)
print("GRID CUDA tiny prediction PASS")
PY
printf 'MATRIX_EXIT=%s\nPRSFORMER_EXIT=%s\n' "$matrix_rc" "$prsformer_rc" > "$AUDIT/gpu-exits.txt"
exit "$prsformer_rc"
