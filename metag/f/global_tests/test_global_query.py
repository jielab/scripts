"""Global streaming/query validation and original-read split integration.

The small native producer here tests supervision and the real public Kraken 2
five-column protocol. It is not a biological classifier; separate native tests
must validate this protocol against the pinned upstream executable.
"""
import gzip
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import global_query as query
import global_split
from library import file_state


PRODUCER = r'''#!/usr/bin/env python3
import json, os, pathlib, sys, time
if '-h' in sys.argv:
    for flag in ('H','t','o','p','T','g','Q','S'):
        print('  -'+flag+' option')
    raise SystemExit(0)
mode = os.environ.get('GLOBAL_TEST_MODE','normal')
trace = os.environ.get('GLOBAL_TEST_TRACE')
paired = '-S' in sys.argv
rows=[]
while True:
    head=sys.stdin.buffer.readline()
    if not head: break
    seq=sys.stdin.buffer.readline().rstrip(b'\r\n')
    plus=sys.stdin.buffer.readline()
    qual=sys.stdin.buffer.readline().rstrip(b'\r\n')
    if not head.startswith(b'@') or plus!=b'+\n' or len(seq)!=len(qual):
        raise SystemExit(8)
    ordinal=int(head[1:])
    size=str(len(seq))
    if paired:
        h2=sys.stdin.buffer.readline()
        s2=sys.stdin.buffer.readline().rstrip(b'\r\n')
        p2=sys.stdin.buffer.readline()
        q2=sys.stdin.buffer.readline().rstrip(b'\r\n')
        if h2!=head or p2!=b'+\n' or len(s2)!=len(q2): raise SystemExit(9)
        size+='|'+str(len(s2))
    tid=13 if ordinal%3==1 else (10 if ordinal%3==2 else 0)
    rows.append(('C' if tid else 'U',str(ordinal),str(tid),size,'13:1 |:| 13:1' if paired else '13:1'))
if trace:
    with open(trace,'a') as h: h.write(json.dumps({'argv':sys.argv,'n':len(rows),'pid':os.getpid()})+'\n')
if mode=='hang':
    while True: time.sleep(0.1)
if mode=='duplicate' and rows: rows.append(rows[0])
if mode=='missing' and rows: rows.pop()
if mode=='out_of_range' and rows:
    x=list(rows[0]);x[1]=str(len(rows)+999);rows[0]=tuple(x)
if mode=='bad_length' and rows:
    x=list(rows[0]);x[3]='99|99' if paired else '99';rows[0]=tuple(x)
if mode=='bad_taxid' and rows:
    x=list(rows[0]);x[0]='C';x[2]='9999';rows[0]=tuple(x)
for row in reversed(rows): print('\t'.join(row))
raise SystemExit(7 if mode=='nonzero' else 0)
'''


def record(name, bases=b'ACGT', quality=b'!#$%', newline=b'\n'):
    return b'@'+name+newline+bases+newline+b'+'+newline+quality+newline


class GlobalQueryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.db=self.root/'db'
        self.tax=self.db/'taxonomy'
        self.tax.mkdir(parents=True)
        entries=[(1,1,'no rank','root'),(2,1,'superkingdom','Bacteria'),
                 (10,2,'genus','Example'),(11,10,'species','Example alpha'),
                 (13,11,'subspecies','Example strain')]
        (self.tax/'nodes.dmp').write_text(''.join(f'{i}\t|\t{p}\t|\t{r}\t|\n' for i,p,r,n in entries))
        (self.tax/'names.dmp').write_text(''.join(f'{i}\t|\t{n}\t|\t\t|\tscientific name\t|\n' for i,p,r,n in entries))
        for name in ('hash.k2d','taxo.k2d','opts.k2d'):
            (self.db/name).write_text('test fixture\n')
        self.snapshot={'identity':'fixture-global-db','path':str(self.db),'taxonomy_path':str(self.tax),
                       'outputs':{str(self.db/name):file_state(self.db/name)
                                  for name in ('hash.k2d','taxo.k2d','opts.k2d')}}
        self.native=self.root/'classify'
        self.native.write_text(PRODUCER)
        self.native.chmod(0o755)
        self.work=self.root/'work'
        self.args=SimpleNamespace(threads=2,confidence=0.1,minimum_hit_groups=2,
                                  minimum_base_quality=0,min_free_gib=0,memory_reserve_gib=0,
                                  work_root=self.root/'runs_global',keep_work=False,benchmark_fragments=2)
        self.trace=self.root/'trace.jsonl'
        self.environment=patch.dict(os.environ,{'GLOBAL_TEST_TRACE':str(self.trace),'GLOBAL_TEST_MODE':'normal'})
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def sample(self, paired=True, count=3):
        original=[]
        paths=[]
        for mate in (1,2) if paired else (0,):
            rows=[record((b'duplicate/'+str(mate).encode()) if paired else b'duplicate',
                         newline=b'\r\n' if i==0 else b'\n') for i in range(count)]
            original.append(rows)
            path=self.root/f'R{mate}.fastq.gz'
            path.write_bytes(gzip.compress(b''.join(rows)))
            paths.append(path)
        return {'name':'Ready__RUN__kraken2','dataset':'Ready','run':'RUN','r1':paths[0],
                'r2':paths[1] if paired else None,'target':self.root/'output__kraken2'},original

    def test_out_of_order_paired_ids_normalize_strain_and_preserve_count(self):
        sample,_=self.sample()
        result=query.query_once(sample,self.snapshot,self.native,self.args,self.work)
        self.assertEqual((result['fragments'],result['reads']),(3,6))
        self.assertEqual(struct.unpack('<QQQ',(self.work/'taxids.bin').read_bytes()),(11,10,0))
        self.assertEqual(result['validated_scope'],'entire_input')
        self.assertIn('-S',result['command'])
        self.assertNotIn('-P',result['command'])
        self.assertFalse((self.work/'lengths.partial.bin').exists())

    def test_single_end_input_uses_no_paired_flag(self):
        sample,_=self.sample(paired=False)
        result=query.query_once(sample,self.snapshot,self.native,self.args,self.work)
        self.assertFalse(result['paired'])
        self.assertEqual(result['reads'],3)
        self.assertNotIn('-S',result['command'])

    def test_duplicate_missing_out_of_range_and_bad_lengths_never_publish(self):
        sample,_=self.sample()
        for mode in ('duplicate','missing','out_of_range','bad_length','bad_taxid','nonzero'):
            with self.subTest(mode=mode),patch.dict(os.environ,{'GLOBAL_TEST_MODE':mode}):
                work=self.root/mode
                with self.assertRaises(RuntimeError):
                    query.query_once(sample,self.snapshot,self.native,self.args,work)
                self.assertFalse((work/'taxids.bin').exists())

    def test_truncated_gzip_cannot_publish_even_when_native_exits_zero(self):
        sample,_=self.sample()
        sample['r1'].write_bytes(sample['r1'].read_bytes()[:-8])
        with self.assertRaisesRegex(RuntimeError,'FASTQ producer failed'):
            query.query_once(sample,self.snapshot,self.native,self.args,self.work)
        self.assertFalse((self.work/'taxids.bin').exists())

    def test_pair_mismatch_cannot_publish(self):
        sample,_=self.sample()
        sample['r2'].write_bytes(gzip.compress(record(b'unrelated/2')*3))
        with self.assertRaisesRegex(RuntimeError,'FASTQ producer failed'):
            query.query_once(sample,self.snapshot,self.native,self.args,self.work)
        self.assertFalse((self.work/'taxids.bin').exists())

    def test_sigterm_handler_stops_native_group_without_publishing(self):
        sample,_=self.sample()
        payload={"sample":{k:str(v) if isinstance(v,Path) else v for k,v in sample.items()},
                 "snapshot":self.snapshot,"executable":str(self.native),"work":str(self.work),
                 "args":vars(self.args).copy()}
        payload['args']['work_root']=str(payload['args']['work_root'])
        request=self.root/'request.json'
        request.write_text(json.dumps(payload))
        code='''import json, signal, sys
from types import SimpleNamespace
from pathlib import Path
import global_query
def stopped(sig, frame): raise KeyboardInterrupt
signal.signal(signal.SIGTERM, stopped)
p=json.loads(Path(sys.argv[1]).read_text())
try:
 global_query.query_once(p['sample'],p['snapshot'],p['executable'],SimpleNamespace(**p['args']),p['work'])
except KeyboardInterrupt:
 raise SystemExit(130)
'''
        env={**os.environ,'GLOBAL_TEST_MODE':'hang',
             'PYTHONPATH':str(Path(query.__file__).parent)}
        proc=subprocess.Popen([sys.executable,'-c',code,str(request)],env=env,
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            deadline=time.monotonic()+10
            while not self.trace.exists() and time.monotonic()<deadline:
                time.sleep(0.02)
            self.assertTrue(self.trace.exists(),'native producer did not start')
            pid=json.loads(self.trace.read_text().splitlines()[0])['pid']
            proc.terminate()
            stdout,stderr=proc.communicate(timeout=15)
            self.assertEqual(proc.returncode,130,(stdout,stderr))
            with self.assertRaises(ProcessLookupError):os.kill(pid,0)
            self.assertFalse((self.work/'taxids.bin').exists())
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.communicate()

    def test_benchmark_queries_prefix_without_publishing_species_bins(self):
        sample,_=self.sample()
        with patch.object(query,'_memory_check',return_value={}):
            folder=query.benchmark(sample,self.snapshot,self.native,self.args)
        receipt=json.loads((folder/'benchmark.json').read_text())
        self.assertEqual(receipt['query']['fragments'],2)
        self.assertFalse(receipt['species_output_created'])
        self.assertFalse(sample['target'].exists())
        self.assertFalse((folder/'taxids.bin').exists())

    def test_completed_query_resumes_failed_split_then_skips_completed_sample(self):
        sample,raw=self.sample()
        with patch.object(query,'_memory_check',return_value={}), \
                patch.object(global_split,'split',side_effect=RuntimeError('injected split interruption')):
            with self.assertRaisesRegex(RuntimeError,'injected'):
                query.run_sample(sample,self.snapshot,self.native,self.args)
        self.assertEqual(len(self.trace.read_text().splitlines()),1)
        with patch.object(query,'_memory_check',side_effect=AssertionError('query should be reused')):
            output=query.run_sample(sample,self.snapshot,self.native,self.args)
        self.assertEqual(len(self.trace.read_text().splitlines()),1)
        for mate in (1,2):
            self.assertEqual(gzip.decompress((output/'bacteria/11_Example_alpha'/f'PE.R{mate}.fastq.gz').read_bytes()),raw[mate-1][0])
        receipt=json.loads((output/'split.done').read_text())
        self.assertEqual((receipt['fragments'],receipt['reads']),(3,6))
        self.assertIn('Kraken 2',receipt['classifier'])
        self.assertTrue(all(s['fragments']>=1 for s in receipt['rows']))
        again=query.run_sample(sample,self.snapshot,self.native,self.args)
        self.assertEqual(again,output)
        self.assertEqual(len(self.trace.read_text().splitlines()),1)


class OrdinalStoreTests(unittest.TestCase):
    def test_grows_bounded_arrays_and_accepts_reversed_output(self):
        with tempfile.TemporaryDirectory() as directory:
            store=query.OrdinalStore(directory,chunk=2)
            for i in range(5):self.assertEqual(store.add_input(i+1,i+2),i+1)
            self.assertEqual(store.capacity,6)
            for i in reversed(range(5)):store.add_output(i+1,i+11,(i+1,i+2))
            target=Path(directory)/'taxids.bin'
            store.publish(target)
            self.assertEqual(struct.unpack('<QQQQQ',target.read_bytes()),(11,12,13,14,15))


if __name__=='__main__':unittest.main()
