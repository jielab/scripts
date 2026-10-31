#!/usr/bin/env python3
"""Apply the small, reproducible adapter patch to a dedicated PhyML checkout."""
import pathlib, subprocess, sys
root=pathlib.Path(sys.argv[1])
commit='42bb4b52f479938a3427dfe4357b861f0e306686'
if subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()!=commit:
    raise SystemExit('PhyML source revision differs from the validated revision')
# Always regenerate from the pinned original; never accumulate patches.
for name in ('lk.c','main.c'):
    original=subprocess.check_output(['git','-C',str(root),'show',f'{commit}:src/{name}'],text=True)
    if name=='lk.c':
        original=original.replace('#include "lk.h"','#include "lk.h"\n#ifdef GU_CUDA\nint GU_Update_Partial(t_tree*,t_edge*,t_node*);\n#endif',1)
        needle='  if(d->tax) return;\n\n  \n#ifdef BEAGLE'
        assert original.count(needle)==1
        original=original.replace(needle,'  if(d->tax) return;\n#ifdef GU_CUDA\n  GU_Update_Partial(tree,b,d);\n  return;\n#endif\n#ifdef BEAGLE',1)
    else:
        needle='int main(int argc, char **argv)\n{'
        original=original.replace(needle,needle+'\n#ifdef GU_CUDA\n  extern void gu_cuda_check(void);\n  gu_cuda_check();\n  if(argc==2 && !strcmp(argv[1],"--gpu-check")) return 0;\n#endif',1)
    (root/'src'/name).write_text(original)
