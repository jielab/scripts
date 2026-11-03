#!/usr/bin/env python3
"""Exercise actual R publication using a synthetic table; no cohort data."""
from pathlib import Path
import importlib.util
import json
import tempfile

spec=importlib.util.spec_from_file_location('gu_common',Path(__file__).resolve().parents[1]/'0.common.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp);(root/'tiny.tsv').write_text('sample_id\tLOD\n0001\t5.0\n')
 job=root/'spec.json'
 job.write_text(json.dumps([{'method':'ibdmix','source':str(root),'files':[{'name':'tiny.tsv'}],
                            'figures':[],'destination':str(root/'tiny.xlsx')}]))
 m.gu_result_r('pack',job)
 assert (root/'tiny.xlsx').stat().st_size>100
 import openpyxl
 wb=openpyxl.load_workbook(root/'tiny.xlsx')
 assert wb.worksheets[0].cell(2,1).value=='0001'
 print('GU_PUBLICATION_PASS: actual R/Python workbook; leading-zero sample ID retained')
