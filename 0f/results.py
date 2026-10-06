#!/usr/bin/env python3
"""Read exact workbook exports and exchange named RDS tables through /tmp."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import sys
import zipfile


sys.dont_write_bytecode = True


# 🚩 R table storage
def run_r(*arguments):
	environment = os.environ.copy()
	rscript = environment.get("RESULTS_RSCRIPT", "/usr/bin/Rscript")
	if Path(rscript).resolve() == Path("/usr/bin/Rscript").resolve():
		for name in ("R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_ENVIRON_USER"):
			environment.pop(name, None)
	with tempfile.TemporaryDirectory(prefix = "result-code-", dir = "/tmp") as directory:
		script = Path(directory) / "results.R"
		script.write_bytes(Path(__file__).with_suffix(".R").read_bytes())
		return subprocess.run(
			[rscript, "--vanilla", str(script), *map(str, arguments)],
			check = True, capture_output = True, text = True, env = environment,
		).stdout.strip()


def write_rds(table, path, metadata = None, float_format = "%.17g"):
	path = Path(path)
	if path.suffix.lower() != ".rds":
		raise ValueError("Reusable raw data or fitted objects require a descriptive .rds filename")
	with tempfile.TemporaryDirectory(prefix = "result-write-", dir = "/tmp") as directory:
		if isinstance(table, dict):
			source = Path(directory) / "model.json"
			source.write_text(json.dumps(table), encoding = "utf-8")
			run_r("import-object", source, path)
			return path
		source = Path(directory) / (path.stem + ".tsv.gz")
		table.to_csv(source, sep = "\t", index = False, compression = "gzip", float_format = float_format)
		arguments = ["import-table", source, path]
		if metadata is not None:
			meta = Path(directory) / "metadata.json"
			meta.write_text(json.dumps(metadata), encoding = "utf-8")
			arguments.append(meta)
		run_r(*arguments)
	return path


def rds_metadata(path):
	return json.loads(run_r("metadata", path))


def write_workbook(tables, path):
	with tempfile.TemporaryDirectory(prefix = "result-book-", dir = "/tmp") as directory:
		sources = {}
		for name, table in tables.items():
			if table.empty:
				continue
			source = Path(directory) / (name + ".tsv")
			if source.parent != Path(directory):
				raise ValueError("Worksheet names cannot contain directory separators")
			table.to_csv(source, sep = "\t", index = False, float_format = "%.17g")
			sources[name] = str(source)
		spec = Path(directory) / "workbook.json"
		spec.write_text(json.dumps({"tables": sources}), encoding = "utf-8")
		run_r("workbook", spec, path)
	return Path(path)


def materialize_rds(path):
	path = Path(path).resolve()
	stat = path.stat()
	key = hashlib.sha256(f"{path}:{stat.st_size}:{stat.st_mtime_ns}".encode()).hexdigest()
	directory = Path("/tmp/analysis-results") / key
	directory.mkdir(parents = True, exist_ok = True, mode = 0o700)
	with (directory / "lock").open("a") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		marker = directory / "table-name"
		if marker.is_file():
			table = directory / marker.read_text()
			if table.is_file():
				return table
		table = Path(run_r("export-table", path, directory))
		if table.parent != directory or not table.is_file():
			raise ValueError("RDS exchange escaped its temporary directory")
		marker.write_text(table.name)
	return table


# 🚩 Exact aggregate workbook sources
def workbook_sources(path):
	with zipfile.ZipFile(path) as archive:
		if "results/manifest.json" not in archive.namelist():
			return {}
		manifest = json.loads(archive.read("results/manifest.json"))
		if manifest.get("format") != "analysis-tables-v1":
			raise ValueError(f"Unsupported result workbook: {path}")
		return manifest.get("files", {})


def materialize_workbook(path, table=None):
 path = Path(path)
 with zipfile.ZipFile(path) as archive:
  manifest = json.loads(archive.read('results/manifest.json'))
  candidates = manifest.get('tables', [])
  if table is not None:candidates = [x for x in candidates if x['name'] == table or x['source'] == table]
  if len(candidates) != 1:
   raise ValueError('Select one named table from this workbook: ' + str(path))
  name = candidates[0]['source']
  entry = manifest['files'][name]
  directory = Path('/tmp/analysis-results') / entry['sha256']
  directory.mkdir(parents=True, exist_ok=True)
  target = directory / Path(name).name
  if not target.is_file():
   fd,temporary = tempfile.mkstemp(dir=directory);os.close(fd)
   digest = hashlib.sha256()
   try:
    with archive.open(entry['part']) as source, open(temporary,'wb') as output:
     for block in iter(lambda:source.read(4*1024*1024),b''):output.write(block);digest.update(block)
    if digest.hexdigest() != entry['sha256']:raise ValueError('Workbook source checksum mismatch')
    os.replace(temporary,target)
   finally:Path(temporary).unlink(missing_ok=True)
 return target


def resolve_table(path):
	path = Path(path)
	if path.is_file():
		return materialize_rds(path) if path.suffix.lower() == ".rds" else materialize_workbook(path) if path.suffix.lower() == ".xlsx" else path
	if path.suffix.lower() == ".rds" and path.with_suffix(".xlsx").is_file():
		return materialize_workbook(path.with_suffix(".xlsx"))
	stem = path.name.removesuffix(".gz").removesuffix(".tsv").removesuffix(".csv")
	private = path.with_name(stem + ".rds")
	if private.is_file():
		return materialize_rds(private)
	for workbook in sorted(path.parent.glob("*.xlsx")):
		entry = workbook_sources(workbook).get(path.name)
		if entry is None:
			continue
		part = f"results/exports/{entry['sha256']}.bin"
		if entry["part"] != part or len(entry["sha256"]) != 64:
			raise ValueError("Invalid workbook export path")
		directory = Path("/tmp/analysis-results") / entry["sha256"]
		directory.mkdir(parents = True, exist_ok = True, mode = 0o700)
		target = directory / path.name
		if not target.is_file():
			with zipfile.ZipFile(workbook) as archive:
				data = archive.read(part)
			if hashlib.sha256(data).hexdigest() != entry["sha256"]:
				raise ValueError(f"Workbook export checksum mismatch: {workbook}")
			with tempfile.NamedTemporaryFile(dir = directory, delete = False) as stream:
				stream.write(data)
				name = stream.name
			os.replace(name, target)
		return target
	return path


def read_table(path, **options):
	import pandas as pd

	options.setdefault("sep", "\t")
	options.setdefault("float_precision", "round_trip")
	return pd.read_csv(resolve_table(path), **options)


# Large review tables are streamed into real Excel worksheets using only the
# standard library. Exact source tables remain inside the workbook for reuse.
def stream_workbook(specification, destination):
	import csv
	import gzip
	import math
	import re
	import shutil
	from xml.sax.saxutils import escape, quoteattr
	from xml.etree import ElementTree as ET
	from contextlib import ExitStack
	csv.field_size_limit(1024 * 1024 * 1024)

	spec = json.loads(Path(specification).read_text())
	destination = Path(destination)
	destination.parent.mkdir(parents=True, exist_ok=True)
	invalid = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f]')
	sheets = []
	used = set()
	manifest = {'format': 'analysis-tables-v1', 'files': {}, 'tables': [], 'metadata': spec.get('metadata')}
	long_rows = []
	strings = []
	string_index = {}
	def column_name(number):
		name = ''
		while number:
			number,remainder = divmod(number-1,26);name = chr(65+remainder)+name
		return name
	location = {'table':'','row':0,'column':''}
	long_sheet = '_long_text'
	while any(t['name'] == long_sheet for t in spec['tables']):long_sheet += '_'
	def string_cell(value):
		if len(value) > 32767:
			for part,start in enumerate(range(0,len(value),30000),1):
				long_rows.append([location['table'],location['row'],location['column'],part,value[start:start+30000]])
			value = f"Full text in {long_sheet}; source row {location['row']}, column {location['column']}"
		value = re.sub(r'_x[0-9A-Fa-f]{4}_', lambda m:'_x005F_'+m.group()[1:], value)
		value = invalid.sub(lambda m:f'_x{ord(m.group()):04X}_', value)
		if value not in string_index:
			string_index[value] = len(strings);strings.append(value)
		return '<c t="s"><v>'+str(string_index[value])+'</v></c>'
	def safe_name(name, part):
		base = re.sub(r'[\[\]:*?/\\]', '_', name).strip("'") or 'results'
		suffix = '' if part == 1 else '_' + str(part)
		candidate = base[:31-len(suffix)] + suffix
		number = 1
		while candidate.lower() in used:
			number += 1
			tail = '_' + str(number)
			candidate = (base[:31-len(suffix)-len(tail)] + suffix + tail)
		used.add(candidate.lower())
		return candidate
	fd, temporary = tempfile.mkstemp(prefix='result-workbook-', suffix='.xlsx', dir='/tmp')
	os.close(fd)
	long_path = Path(temporary + '.long_text.tsv')
	original_count = len(spec['tables'])
	try:
		with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive, ExitStack() as handles:
			for table_index,table in enumerate(spec['tables']):
				location['table'] = table['name']
				path = Path(table['path'])
				opener = gzip.open if path.suffix == '.gz' else open
				with opener(path, 'rt', newline='', encoding='utf-8-sig') as stream:
					reader = csv.reader(stream, delimiter=table.get('separator', '\t'))
					header = next(reader)
					columns = [column_name(i+1) for i in range(len(header))]
					if len(header) > 16384:
						raise ValueError('Too many Excel columns')
					types = table.get('types', ['character'] * len(header))
					if len(types) != len(header):raise ValueError('Column type/header mismatch')
					part = 0; count = 0; target = None; sheet_rows = 0
					def start_sheet():
						nonlocal part, target, sheet_rows
						part += 1
						sheets.append({'name': safe_name(table['name'],part), 'rows': 0, 'columns': len(header)})
						target = handles.enter_context(archive.open(f'xl/worksheets/sheet{len(sheets)}.xml', 'w', force_zip64=True))
						target.write(('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews><sheetData><row r="1">'+''.join(string_cell(x).replace('<c ',f'<c r="{column}1" ',1) for column,x in zip(columns,header))+'</row>').encode())
						sheet_rows = 0
					def close_sheet():
						target.write(b'</sheetData></worksheet>');target.close()
						sheets[-1]['rows'] = sheet_rows
					start_sheet()
					for row in reader:
						if len(row) != len(header):raise ValueError(f'Ragged result table: {path}')
						if sheet_rows == 1048575:close_sheet();start_sheet()
						cells = []
						location['row'] = count + 1
						for column,value,kind in zip(header,row,types):
							location['column'] = column
							if value == '':cells.append('<c/>')
							elif kind == 'numeric' and math.isfinite(float(value)):cells.append('<c><v>'+escape(value)+'</v></c>')
							elif kind == 'logical' and value in ('TRUE','FALSE'):cells.append('<c t="b"><v>'+str(int(value == 'TRUE'))+'</v></c>')
							else:cells.append(string_cell(value))
						cells = [cell.replace('<c',f'<c r="{column}{sheet_rows+2}"',1) for column,cell in zip(columns,cells)]
						sheet_rows += 1;count += 1
						target.write(('<row r="'+str(sheet_rows+1)+'">'+''.join(cells)+'</row>').encode())
					close_sheet()
				manifest['tables'].append({'name':table['name'],'source':path.name,'rows':count,'columns':len(header),'types':types})
				if table_index == original_count - 1 and long_rows:
					with long_path.open('w',newline='') as stream:
						writer = csv.writer(stream,delimiter='\t');writer.writerow(['source_table','source_row','column','part','text']);writer.writerows(long_rows)
					spec['tables'].append({'name':long_sheet,'path':str(long_path),'types':['character','numeric','character','numeric','character']})
					long_rows.clear()
			for source in dict.fromkeys([t['path'] for t in spec['tables']] + spec.get('sources', [])):
				path = Path(source)
				hash_value = hashlib.sha256()
				with path.open('rb') as stream:
					for block in iter(lambda:stream.read(4*1024*1024),b''):hash_value.update(block)
				sha = hash_value.hexdigest()
				part = f'results/exports/{sha}.bin'
				if path.name in manifest['files']:
					if manifest['files'][path.name]['sha256'] != sha:raise ValueError('Conflicting workbook sources: '+path.name)
					continue
				if part not in archive.namelist():archive.write(path, part)
				manifest['files'][path.name] = {'part':part,'sha256':sha}
			manifest['writer'] = 'shared-strings-v2'
			archive.writestr('results/manifest.json', json.dumps(manifest))
			with archive.open('xl/sharedStrings.xml','w',force_zip64=True) as target:
				target.write(b'<?xml version="1.0" encoding="UTF-8"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">')
				for value in strings:target.write(('<si><t xml:space="preserve">'+escape(value)+'</t></si>').encode())
				target.write(b'</sst>')
			ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
			rel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
			archive.writestr('xl/workbook.xml', '<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="'+ns+'" xmlns:r="'+rel+'"><sheets>'+''.join('<sheet name='+quoteattr(s['name'])+' sheetId="'+str(i)+'" r:id="rId'+str(i)+'"/>' for i,s in enumerate(sheets,1))+'</sheets></workbook>')
			archive.writestr('xl/_rels/workbook.xml.rels', '<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'+''.join(f'<Relationship Id="rId{i}" Type="{rel}/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(sheets)+1))+f'<Relationship Id="rId{len(sheets)+1}" Type="{rel}/sharedStrings" Target="sharedStrings.xml"/>'+'</Relationships>')
			archive.writestr('_rels/.rels', f'<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="{rel}/officeDocument" Target="xl/workbook.xml"/></Relationships>')
			archive.writestr('[Content_Types].xml', '<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Default Extension="json" ContentType="application/json"/><Default Extension="bin" ContentType="application/octet-stream"/><Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'+''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(sheets)+1))+'</Types>')
		with zipfile.ZipFile(temporary) as archive:
			if archive.testzip() is not None:raise ValueError('Workbook CRC validation failed')
			for i,sheet in enumerate(sheets,1):
				rows = 0
				with archive.open(f'xl/worksheets/sheet{i}.xml') as stream:
					for event,element in ET.iterparse(stream,events=('end',)):
						if element.tag.endswith('}row'):rows += 1
						element.clear()
				if rows != sheet['rows']+1:raise ValueError('Workbook row verification failed')
		fd,stage = tempfile.mkstemp(prefix='.'+destination.name+'.part.',dir=destination.parent);os.close(fd)
		try:
			shutil.copyfile(temporary,stage)
			with open(stage,'rb') as stream:os.fsync(stream.fileno())
			os.replace(stage,destination)
		finally:Path(stage).unlink(missing_ok=True)
	finally:
		Path(temporary).unlink(missing_ok=True)
		long_path.unlink(missing_ok=True)
	return manifest


def upgrade_workbook(path):
 """Repackage an early streaming workbook for Excel and R interoperability."""
 import math
 from xml.etree import ElementTree as ET
 path = Path(path)
 with tempfile.TemporaryDirectory(prefix='result-upgrade-',dir='/tmp') as directory:
  root = Path(directory)
  with zipfile.ZipFile(path) as archive:
   manifest = json.loads(archive.read('results/manifest.json'))
   if manifest.get('writer') == 'shared-strings-v2':return False
   spec = {'tables': [], 'sources': [], 'metadata': manifest.get('metadata')}
   sheet = 1
   for table in manifest['tables']:
    parts = max(1, math.ceil(table['rows'] / 1048575))
    if table['source'].endswith('.long_text.tsv'):
     sheet += parts;continue
    types = table.get('types', [None] * table['columns'])
    for index in range(sheet,sheet+parts):
     if all(types):break
     with archive.open(f'xl/worksheets/sheet{index}.xml') as source:
      for event,row in ET.iterparse(source, events=('end',)):
       if not row.tag.endswith('}row'):continue
       if row.get('r') != '1':
        for i,cell in enumerate(row):
         if i >= len(types) or types[i] or not len(cell):continue
         types[i] = 'character' if cell.get('t') in ('s','inlineStr','str') else 'logical' if cell.get('t') == 'b' else 'numeric'
       row.clear()
       if all(types):break
    types = [t or 'character' for t in types]
    name = Path(table['source']).name
    entry = manifest['files'][table['source']]
    target = root / name
    digest = hashlib.sha256()
    with archive.open(entry['part']) as source,target.open('wb') as output:
     for block in iter(lambda:source.read(4*1024*1024),b''):digest.update(block);output.write(block)
    if digest.hexdigest()!=entry['sha256']:raise ValueError('Stored workbook source checksum mismatch')
    spec['tables'].append({'name':table['name'],'path':str(target),'types':types})
    sheet += parts
  job = root / 'spec.json';job.write_text(json.dumps(spec))
  stream_workbook(job,path)
 return True


if __name__ == '__main__':
	if len(sys.argv) == 4 and sys.argv[1] == 'stream-workbook':
		stream_workbook(sys.argv[2],sys.argv[3])
	else:
		raise SystemExit('Usage: results.py stream-workbook SPEC.json OUTPUT.xlsx')
