"""Exercise the actual BED -> producer -> consumer contract with tiny real files."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import rdata
import block_inputs as b
from abm_core.grid_inputs import locate_blocks


def fixture(root, names, n=40):
    config=json.loads((Path(b.__file__).parents[1]/'config/abm.blocks.json').read_text())
    for key,folder in [('ld_block_dir','bed'),('prsformer_cache','genotype'),('csx_analysis_root','weights'),
                       ('score_root','scores'),('cache_root','cache'),('output_root','output')]:
        config[key]=str(root/folder);(root/folder).mkdir()
    config.update(ld_block_files=names,score_batch_size=13,reconstruction_sample_n=32)
    for name in names.split(','):(root/'bed'/name).write_text('1\t0\t100\n1\t100\t200\n')
    cohort=pd.DataFrame(dict(eid=[f'{i:06d}' for i in range(n)],family_id=[f'g{i//2}' for i in range(n)]))
    cohort.to_csv(root/'genotype/data.tsv.gz',sep='\t',index=False)
    variants=pd.DataFrame(dict(SNP=['rs1','rs2'],CHR=[1,1],BP=[20,150],REF=['A','C'],ALT=['G','T']))
    variants.to_csv(root/'genotype/variants.tsv.gz',sep='\t',index=False)
    geno=np.random.default_rng(37).uniform(0,2,(n,2)).astype('float32')
    np.save(root/'genotype/genotypes.npy',geno);(root/'genotype/prepare.json').write_text('{}')
    (root/'weights/height').mkdir();models={};scores=cohort[['eid']].copy()
    for i,pop in enumerate(b.POPS):
        weight=root/'weights'/f'{pop}.tsv'
        pd.DataFrame(dict(SNP=variants.SNP,CHR=variants.CHR,BP=variants.BP,A1=variants.ALT,
                          A2=variants.REF,BETA=[.125*(i+1),.25*(i+1)])).to_csv(weight,sep='\t',index=False)
        models['csx.'+pop]=dict(weights_file=str(weight),weights_sha256=b.sha256(weight))
        scores['csx.'+pop]=geno@np.array([.125*(i+1),.25*(i+1)])
    (root/'weights/height/provenance.json').write_text(json.dumps(dict(models=models)))
    (root/'scores/height').mkdir();rdata.write_rds(root/'scores/height/1csx.scores.rds',scores)
    path=root/'config.json';path.write_text(json.dumps(config))
    return path,config


class ProducerTests(unittest.TestCase):
    def test_failed_score_generation_removes_partial_memmap(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);path,cfg=fixture(root,'ALL.37.bed')
            layout=b.resolve_block_files(cfg['ld_block_dir'],'ALL.37.bed')
            base,score,_,_,receipt=b.cache_paths(cfg,layout,'height')
            with patch('torch.zeros',side_effect=RuntimeError('injected allocation failure')):
                with self.assertRaisesRegex(RuntimeError,'injected allocation'):
                    b.prepare_trait(path,'height',device='cpu')
            self.assertFalse(score.with_name(score.name+'.tmp').exists())
            self.assertFalse(score.exists())
            self.assertFalse(receipt.exists())

    def test_failed_memmap_allocation_removes_partial_file(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);path,cfg=fixture(root,'ALL.37.bed')
            layout=b.resolve_block_files(cfg['ld_block_dir'],'ALL.37.bed')
            _,score,_,_,receipt=b.cache_paths(cfg,layout,'height')
            def fail_allocation(target, **kwargs):
                Path(target).write_bytes(b'partial NPY allocation')
                raise OSError('injected disk full')
            with patch.object(np.lib.format,'open_memmap',side_effect=fail_allocation):
                with self.assertRaisesRegex(OSError,'injected disk full'):
                    b.prepare_trait(path,'height',device='cpu')
            self.assertFalse(score.with_name(score.name+'.tmp').exists())
            self.assertFalse(receipt.exists())

    def test_both_layouts_actual_producer_binding_and_consumer(self):
        for names in ['ALL.37.bed',','.join(p+'.37.bed' for p in b.POPS)]:
            with self.subTest(layout=names),tempfile.TemporaryDirectory() as td:
                root=Path(td);path,cfg=fixture(root,names)
                layout=b.resolve_block_files(cfg['ld_block_dir'],names)
                tokens,channels,audit=b.build_manifest(layout,cfg)
                self.assertNotIn('aliases',audit)
                self.assertEqual(channels,list(b.POPS) if names.startswith('ALL') else ['scalar'])
                self.assertEqual(tokens.token_index.tolist(),list(range(len(tokens))))
                _,base,receipt=b.prepare_trait(path,'height',device='cpu')
                outer=dict(block_config=str(path),execution={'device':'cpu'})
                arr,cohort,tokens,saved,arrpath=locate_blocks(outer,'height',allow_build=False)
                self.assertEqual(saved['array_binding']['axes'],['individual','token','channel'])
                self.assertEqual(saved['array_binding']['dtype'],'float32')
                original=np.array(arr);del arr
                # Select two rows omitted from the consumer's independent 32-row probe.
                probe=set(np.random.default_rng(20261020).choice(40,32,replace=False))
                a,c=[i for i in range(40) if i not in probe][:2]
                changed=original.copy();changed[[a,c]]=changed[[c,a]];np.save(arrpath,changed)
                with self.assertRaisesRegex(ValueError,'bytes changed'):locate_blocks(outer,'height',allow_build=False)
                np.save(arrpath,original);changed=original.copy();changed[-1,0,0]+=.125;np.save(arrpath,changed)
                with self.assertRaisesRegex(ValueError,'bytes changed'):locate_blocks(outer,'height',allow_build=False)
                np.save(arrpath,original)
                tokenpath=base/'tokens.tsv.gz';tokens.iloc[::-1].to_csv(tokenpath,sep='\t',index=False)
                with self.assertRaisesRegex(ValueError,'token'):locate_blocks(outer,'height',allow_build=False)
                tokens.to_csv(tokenpath,sep='\t',index=False)
                cohort.iloc[::-1].to_csv(base/'cohort.tsv.gz',sep='\t',index=False)
                with self.assertRaisesRegex(ValueError,'row identities'):locate_blocks(outer,'height',allow_build=False)
                cohort.to_csv(base/'cohort.tsv.gz',sep='\t',index=False)
                receipt.pop('array_binding');(base/'BLOCKS_COMPLETE.json').write_text(json.dumps(receipt,default=str))
                with self.assertRaisesRegex(ValueError,'trustworthy'):locate_blocks(outer,'height',allow_build=False)
                self.assertNotIn('array_binding',json.loads((base/'BLOCKS_COMPLETE.json').read_text()))

if __name__=='__main__':unittest.main()
