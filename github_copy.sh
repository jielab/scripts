#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi


# 🚩 github_copy
# Build a bounded sharing copy from final figures and summary results.
set -Eeuo pipefail
command -v python3 >/dev/null 2>&1 || {
	echo 'ERROR: python3 is required' >&2
	exit 1
}
exec python3 - "$@" <<'PY'
import argparse
import codecs
import collections
from contextlib import contextmanager
import csv
import decimal
import fcntl
import gzip
import hashlib
from html.parser import HTMLParser
import io
import itertools
import json
import os
from pathlib import Path
import posixpath
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
import zipfile
from xml.etree import ElementTree as ET

SOURCE = Path(os.environ.get('GITHUB_COPY_SOURCE_ROOT', '/mnt/d/analysis'))
DEST = Path(os.environ.get('GITHUB_COPY_DEST_ROOT', '/mnt/d/github/analysis'))
MANIFEST_SUFFIX = '.files.lst'
CHECK_NAMES = ('size', 'type', 'id')
CHECKS = frozenset(CHECK_NAMES)
MAX_FILE_BYTES = 10_000_000
MAX_PROJECT_BYTES = 1_000_000_000
MAX_EXPANDED_BYTES = 512 * 1024 * 1024
MAX_XML_BYTES = 256 * 1024 * 1024
INTERMEDIATE_DIRECTORIES = frozenset({
	'raw', 'mask', 'masks', 'cache', 'caches', 'extract', 'infer', 'input', 'inputs',
	'tmp', 'temp', 'normalize', 'migration', 'work', 'scratch', 'shards', 'chunks',
	'checkpoint', 'checkpoints', 'models', 'weights', 'logs', 'log', 'debug',
	'test', 'tests', 'testing', 'validation', 'pilot', 'backup', 'backups', 'benchmark', 'benchmarks',
	'deployment', 'intermediate', 'intermediates', 'reference', 'references',
	'annotations', 'pgen', 'genotypes', 'metadata', 'join', 'panel', 'candidates',
	'bootstraps', 'permutation', 'permutations', 'folds', 'replicates',
})

# Conservatively reject seven-digit values in ID columns, including mixed
# or partially missing columns. Coordinates without an ID header are not IDs.
SEVEN_DIGIT_ID = re.compile(r'[0-9]{7}(?:\.0+)?')
NUMERIC_ID = re.compile(r'[0-9]+(?:\.[0-9]+)?[eE][+-]?[0-9]+')
EXCEL_ESCAPE = re.compile(r'_x([0-9a-fA-F]{4})_')
BAD_CONTROL = re.compile(r'[\x00-\x08\x0b\x0e-\x1f]')
RASTER = {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.tif', '.tiff', '.bmp'}
RESULT_SUFFIXES = {'.png', '.jpg', '.jpeg', '.webp', '.xlsx', '.pdf', '.svg'}
TEXT_RESULT_SUFFIXES = {'.csv', '.tsv', '.txt', '.md', '.html', '.htm', '.bed'}
SHEET_NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


class Rejected(Exception):
	def __init__(self, reason, locations=()):
		super().__init__(reason)
		self.locations = locations


class RecoveryFailed(Exception):
	pass


class SourceChanged(Exception):
	pass


@contextmanager
def staging_directory():
	stage = Path(tempfile.mkdtemp(prefix='.github_copy.stage.', dir='/tmp'))
	keep = False
	try:
		yield stage
	except RecoveryFailed:
		# Preserve the only recovery copy if restoring a destination fails.
		keep = True
		raise
	finally:
		if not keep:
			try:
				shutil.rmtree(stage)
			except OSError:
				print(f'WARNING: could not remove temporary directory: {stage}', file=sys.stderr)


def is_seven_digit_id(value):
	text = str(value).strip()
	if SEVEN_DIGIT_ID.fullmatch(text):
		return True
	if len(text) <= 64 and NUMERIC_ID.fullmatch(text):
		try:
			number = decimal.Decimal(text)
			return 1_000_000 <= number <= 9_999_999 and number == number.to_integral_value()
		except decimal.InvalidOperation:
			pass
	return False


class IdColumns:
	def __init__(self, header):
		self.candidates = {i for i, name in enumerate(header) if 'id' in str(name).lower()}

	def add(self, values):
		# Discover repeated headers and headers below a title row too.
		self.candidates.update(
			i for i, value in enumerate(values)
			if isinstance(value, str) and re.fullmatch(r'[\w #.-]*id[\w #.-]*', value, re.IGNORECASE)
		)
		if any(i < len(values) and is_seven_digit_id(values[i]) for i in self.candidates):
			raise Rejected('ukb_id_column')

	def finish(self):
		pass


def check_value(value):
	# JSON tables may be arrays of records, arrays with a header row, or
	# objects of column arrays. Isolated scalar values do not define a column.
	if isinstance(value, dict):
		for name, items in value.items():
			if 'id' in str(name).lower() and isinstance(items, (str, int, float, decimal.Decimal)) and is_seven_digit_id(items):
				raise Rejected('ukb_id_column')
			if (
				'id' in str(name).lower() and isinstance(items, (list, tuple))
				and any(is_seven_digit_id(item) for item in items)
			):
				raise Rejected('ukb_id_column')
			check_value(items)
	elif isinstance(value, (list, tuple)):
		if value and isinstance(value[0], dict):
			header = list(value[0])
			columns = IdColumns(header)
			for row in value:
				columns.add([row.get(name) for name in header] if isinstance(row, dict) else [])
			columns.finish()
		elif value and isinstance(value[0], (list, tuple)):
			columns = IdColumns(value[0])
			for row in value[1:]:
				columns.add(row if isinstance(row, (list, tuple)) else [])
			columns.finish()
		for item in value:
			check_value(item)


def relative_parts(relative):
	if (
		not relative
		or relative.startswith('/')
		or any(c in relative for c in '\\\t\r\n')
		or any(ord(c) < 32 for c in relative)
	):
		raise Rejected('unsafe_path')
	parts = relative.split('/')
	if len(parts) < 2 or any(p in {'', '.', '..'} for p in parts):
		raise Rejected('unsafe_path')
	return parts


def file_stat(root, relative):
	parts = relative_parts(relative)
	current = root
	for part in parts:
		current = current / part
		if current.is_symlink():
			raise Rejected('symlink')
	details = current.stat()
	if not stat.S_ISREG(details.st_mode):
		raise Rejected('not_a_regular_file')
	return details


def version(details):
	return details.st_size, details.st_mtime_ns, details.st_ctime_ns, details.st_ino


# 🚩 Output exclusions
def name_words(name):
	return set(re.split(r'[^a-z0-9]+', name.lower()))


def excluded_result_directory(relative):
	parts = relative.split('/')
	for part in parts[1:]:
		if part.startswith('.') or name_words(part) & INTERMEDIATE_DIRECTORIES:
			return 'intermediate_directory'
	return ''


def excluded_result(relative):
	path = Path(relative)
	reason = excluded_result_directory(path.parent.as_posix())
	if reason:
		return reason
	if path.name.startswith('.'):
		return 'hidden_file'
	# File names do not establish whether a table/figure is shareable.
	# In particular, performance.tsv, scores.xlsx and test_*.png must reach
	# the content checker instead of being rejected by keyword heuristics.
	if 'type' not in CHECKS or path.suffix.lower() in RESULT_SUFFIXES | TEXT_RESULT_SUFFIXES:
		return ''
	return 'not_a_result_file'


def excluded_artifact(path):
	if 'type' not in CHECKS:
		return ''
	name = path.name.lower()
	if name.startswith('~$'):
		return 'office_lock_file'
	if name.endswith('.gz'):
		name = name[:-3]
	# RDS may contain individual IDs and is opaque to this scanner. Keep
	# it out of scan lists even in confirmed non-UKB projects.
	if '.rds' in [s.lower() for s in path.suffixes]:
		return 'rds_file'
	if name.endswith('.log'):
		return 'log_file'
	if name.endswith('.done'):
		return 'done_file'
	if name in {'plots.pdf', 'rplots.pdf'}:
		return 'unnamed_plot'
	return ''


def file_policy(path, details):
	if 'size' in CHECKS:
		if details.st_size > MAX_FILE_BYTES:
			raise Rejected('over_file_size_limit')
		if details.st_size == 0:
			raise Rejected('empty_file')
	reason = excluded_artifact(path)
	if reason:
		raise Rejected(reason)


def text_handle(raw):
	prefix = raw.read(4)
	raw.seek(0)
	if prefix.startswith((codecs.BOM_UTF32_LE, codecs.BOM_UTF32_BE)):
		encoding = 'utf-32'
	elif prefix.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
		encoding = 'utf-16'
	else:
		encoding = 'utf-8-sig'
	return io.TextIOWrapper(raw, encoding=encoding, errors='strict', newline='')


class MarkupChecker(HTMLParser):
	"""Apply the same column rule to rendered HTML tables."""

	def __init__(self):
		super().__init__(convert_charrefs=True)
		self.tables = []
		self.visible_text = []
		self.ignored_text = 0

	def handle_starttag(self, tag, attrs):
		if tag in {'script', 'iframe', 'object', 'embed', 'image', 'foreignobject'}:
			raise Rejected('uninspectable_active_or_embedded_content')
		for name, value in attrs:
			if name.startswith('on') or (value and value.strip().lower().startswith(('data:', 'javascript:'))):
				raise Rejected('uninspectable_active_or_embedded_content')
		if tag in {'style', 'title'}:
			self.ignored_text += 1
		if tag == 'table':
			self.tables.append({'columns': None, 'row': None, 'cell': None})
		elif self.tables:
			table = self.tables[-1]
			if tag == 'tr':
				table['row'] = []
			elif tag in {'td', 'th'} and table['row'] is not None:
				table['cell'] = []

	def handle_endtag(self, tag):
		if tag in {'style', 'title'}:
			self.ignored_text = max(0, self.ignored_text - 1)
		if not self.tables:
			return
		table = self.tables[-1]
		if tag in {'td', 'th'} and table['cell'] is not None:
			table['row'].append(''.join(table['cell']))
			table['cell'] = None
		elif tag == 'tr' and table['row'] is not None:
			if table['columns'] is None:
				table['columns'] = IdColumns(table['row'])
			else:
				table['columns'].add(table['row'])
			table['row'] = None
		elif tag == 'table':
			if table['columns'] is not None:
				table['columns'].finish()
			self.tables.pop()

	def handle_data(self, data):
		if self.tables and self.tables[-1]['cell'] is not None:
			self.tables[-1]['cell'].append(data)
		elif not self.tables and not self.ignored_text:
			self.visible_text.extend(data.splitlines())

	def close(self):
		super().close()
		check_table(line.split() for line in self.visible_text if line.strip())


def checked_lines(handle):
	for line in handle:
		if BAD_CONTROL.search(line):
			raise Rejected('uninspectable_binary_content')
		yield line


def check_table(rows):
	rows = (row for row in rows if row)
	columns = IdColumns(next(rows, []))
	for row in rows:
		columns.add(row)
	columns.finish()


def check_text(raw, suffix):
	with text_handle(raw) as handle:
		if suffix == '.json':
			check_value(json.load(handle, object_pairs_hook=unique_json_object))
			return
		if suffix == '.jsonl':
			check_value([
				json.loads(line, object_pairs_hook=unique_json_object)
				for line in handle if line.strip()
			])
			return
		lines = checked_lines(handle)
		if suffix in {'.html', '.htm', '.svg', '.xml'}:
			markup = MarkupChecker()
			for line in lines:
				markup.feed(line)
			markup.close()
			return
		csv.field_size_limit(MAX_XML_BYTES)
		if suffix in {'.csv', '.tsv'}:
			check_table(csv.reader(lines, delimiter='\t' if suffix == '.tsv' else ',', strict=True))
			return
		if suffix == '.md':
			check_table(line.strip().strip('|').split('|') if '|' in line else line.split() for line in lines)
			return
		# Read plain text result tables and text exports embedded in XLSX.
		first = next((line for line in lines if line.strip()), '')
		lines = itertools.chain([first], lines)
		delimiter = '\t' if '\t' in first else ',' if ',' in first else None
		if delimiter:
			check_table(csv.reader(lines, delimiter=delimiter))
		else:
			check_table(line.split() for line in lines)


def check_printer_settings(data):
	# Known Excel/openxlsx printer metadata is not a table of participants.
	if len(data) > 1024 * 1024:
		raise Rejected('xlsx_printer_settings_too_large')


def unique_json_object(pairs):
	result = {}
	for key, value in pairs:
		if key in result:
			raise Rejected('duplicate_json_key')
		result[key] = value
	return result


def check_workbook_exports(book):
	checked = set()
	formats = {'le8': 'le8-workbook-v1', 'results': 'analysis-tables-v1'}
	for prefix, expected in formats.items():
		manifest_name = prefix + '/manifest.json'
		if manifest_name not in book.namelist():
			continue
		manifest = json.loads(book.read(manifest_name), object_pairs_hook=unique_json_object)
		if not isinstance(manifest, dict) or manifest.get('format') != expected:
			raise Rejected('xlsx_invalid_export_manifest')
		files = manifest.get('files')
		if not isinstance(files, dict):
			raise Rejected('xlsx_invalid_export_manifest')
		check_value({key: value for key, value in manifest.items() if key != 'files'})
		for source, entry in files.items():
			check_value(source)
			if (
				Path(source).name != source or '\\' in source
				or not isinstance(entry, dict)
			):
				raise Rejected('xlsx_invalid_export_manifest')
			sha = entry.get('sha256')
			if not isinstance(sha, str) or not re.fullmatch(r'[a-f0-9]{64}', sha):
				raise Rejected('xlsx_invalid_export_manifest')
			part = entry.get('part')
			if part != f'{prefix}/exports/{sha}.bin' or part not in book.namelist():
				raise Rejected('xlsx_invalid_export_manifest')
			# A verified byte count describes storage, not an individual.
			if 'size' in entry and (
				type(entry['size']) is not int or entry['size'] != book.getinfo(part).file_size
			):
				raise Rejected('xlsx_export_size_mismatch')
			check_value({key: value for key, value in entry.items() if key != 'size'})
			suffix = Path(source).suffix.lower()
			inner_suffix = Path(Path(source).stem).suffix.lower() if suffix == '.gz' else suffix
			if inner_suffix not in {'.csv', '.tsv', '.json', '.jsonl', '.txt'}:
				raise Rejected('xlsx_uninspectable_embedded_content')
			data = book.read(part)
			if hashlib.sha256(data).hexdigest() != sha:
				raise Rejected('xlsx_export_checksum_mismatch')
			try:
				with io.BytesIO(data) as raw:
					if suffix == '.gz':
						check_gzip(raw, inner_suffix)
					else:
						check_text(raw, suffix)
			except Rejected as exc:
				if str(exc) == 'ukb_id_column':
					raise Rejected('ukb_id_column', [('embedded table', source)]) from None
				raise
			checked.add(part)
		checked.add(manifest_name)
	return checked


def excel_text(value):
	return EXCEL_ESCAPE.sub(lambda match: chr(int(match.group(1), 16)), value)


def workbook_strings(book):
	if 'xl/sharedStrings.xml' not in book.namelist():
		return []
	strings = []
	with book.open('xl/sharedStrings.xml') as handle:
		for _, element in ET.iterparse(handle):
			if element.tag == SHEET_NS + 'si':
				strings.append(excel_text(''.join(t.text or '' for t in element.iter(SHEET_NS + 't'))))
				element.clear()
	return strings


def cell_value(cell, strings):
	if cell.get('t') == 's':
		index = int(cell.findtext(SHEET_NS + 'v'))
		if not 0 <= index < len(strings):
			raise Rejected('xlsx_invalid_shared_string')
		return strings[index]
	if cell.get('t') == 'inlineStr':
		return excel_text(''.join(t.text or '' for t in cell.iter(SHEET_NS + 't')))
	if cell.get('t') not in {None, 'n', 'str'}:
		# Boolean/error/date cells are populated non-ID values, not missing rows.
		return False
	return excel_text(cell.findtext(SHEET_NS + 'v', ''))


def check_worksheet(book, name, strings):
	columns, addresses = None, []
	with book.open(name) as handle:
		in_row = False
		for event, element in ET.iterparse(handle, events=('start', 'end')):
			if event == 'start':
				if element.tag == SHEET_NS + 'row':
					in_row = True
				continue
			if in_row and element.tag != SHEET_NS + 'row':
				continue
			if in_row:
				values = {}
				for index, cell in enumerate(element.findall(SHEET_NS + 'c')):
					address = re.sub(r'[0-9]+$', '', cell.get('r', '')) or index
					if address in values:
						raise Rejected('xlsx_duplicate_cell')
					values[address] = cell_value(cell, strings)
				# Formatting-only blank rows do not add data records. A missing
				# ID in a populated row disqualifies that column, as requested.
				if any(value not in {None, ''} for value in values.values()):
					if columns is None:
						addresses = list(values)
						columns = IdColumns(list(values.values()))
					else:
						addresses.extend(address for address in values if address not in addresses)
						columns.add([values.get(address) for address in addresses])
			element.clear()
			in_row = False
	if columns is not None:
		columns.finish()


def worksheet_names(book):
	names = book.namelist()
	if 'xl/workbook.xml' not in names or 'xl/_rels/workbook.xml.rels' not in names:
		return {}
	relationships = {}
	for item in ET.fromstring(book.read('xl/_rels/workbook.xml.rels')):
		target = item.get('Target', '')
		if item.get('TargetMode') != 'External':
			relationships[item.get('Id')] = posixpath.normpath(
				target.lstrip('/') if target.startswith('/') else posixpath.join('xl', target)
			)
	labels = {}
	for sheet in ET.fromstring(book.read('xl/workbook.xml')).iter(SHEET_NS + 'sheet'):
		relation = sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
		if relation in relationships:
			labels[relationships[relation]] = sheet.get('name', relationships[relation])
	return labels


def check_xlsx(path):
	with zipfile.ZipFile(path) as book:
		infos = book.infolist()
		if sum(item.file_size for item in infos) > MAX_EXPANDED_BYTES or any(
			item.file_size > MAX_XML_BYTES for item in infos
		):
			raise Rejected('xlsx_expansion_limit')
		names = book.namelist()
		if len(names) != len(set(names)):
			raise Rejected('xlsx_duplicate_archive_member')
		if not any(re.fullmatch(r'xl/worksheets/[^/]+\.xml', n) for n in names):
			raise Rejected('xlsx_content_check_failed')
		strings = workbook_strings(book)
		labels = worksheet_names(book)
		worksheets = [name for name in names if re.fullmatch(r'xl/worksheets/[^/]+\.xml', name)]
		blocked = []
		for name in worksheets:
			try:
				check_worksheet(book, name, strings)
			except Rejected as exc:
				if str(exc) != 'ukb_id_column':
					raise
				blocked.append(('worksheet', labels.get(name, name)))
		if blocked:
			raise Rejected('ukb_id_column', blocked)
		exports = check_workbook_exports(book)
		# Check every worksheet, including hidden ones, and exact exports.
		# XML metadata and unused strings alone do not define an ID column.
		for name in names:
			if name.endswith('/') or name in exports or name in worksheets or name == 'xl/sharedStrings.xml':
				continue
			if name.lower().endswith(('.xml', '.rels', '.vml')):
				with book.open(name) as handle:
					for _, element in ET.iterparse(handle):
						element.clear()
			elif Path(name).suffix.lower() in RASTER:
				with book.open(name) as handle:
					check_image(io.BytesIO(handle.read()))
			elif re.fullmatch(r'xl/printerSettings/printerSettings[0-9]+\.bin', name):
				check_printer_settings(book.read(name))
			else:
				raise Rejected('xlsx_uninspectable_embedded_content')
	return ''


def check_gzip(raw, suffix):
	with tempfile.TemporaryFile(dir='/tmp') as expanded:
		with gzip.GzipFile(fileobj=raw, mode='rb') as handle:
			while True:
				chunk = handle.read(1024 * 1024)
				if not chunk:
					break
				if expanded.tell() + len(chunk) > MAX_EXPANDED_BYTES:
					raise Rejected('gzip_expansion_limit')
				expanded.write(chunk)
		expanded.seek(0)
		check_text(expanded, suffix)


def check_image(path):
	# Check metadata and validate the actual format, without inspecting pixels
	# or requesting visual review.
	from PIL import Image

	with Image.open(path) as image:
		image.verify()
	if hasattr(path, 'seek'):
		path.seek(0)
	with Image.open(path) as image:
		for value in image.info.values():
			if isinstance(value, bytes):
				value = value.decode('utf-8', errors='ignore')
			check_value(value)
		check_value(dict(image.getexif()))
	return ''


def check_pdf(path):
	with path.open('rb') as handle:
		if handle.read(5) != b'%PDF-':
			raise Rejected('pdf_content_check_failed')
	if shutil.which('pdftotext'):
		command = ['pdftotext', str(path), '-']
	elif shutil.which('gs'):
		command = [
			'gs',
			'-q',
			'-dSAFER',
			'-dBATCH',
			'-dNOPAUSE',
			'-sDEVICE=txtwrite',
			'-sOutputFile=-',
			'-f',
			str(path),
		]
	else:
		raise Rejected('pdf_text_extractor_unavailable')
	with tempfile.TemporaryFile(dir='/tmp') as output:
		process = subprocess.run(
			command, stdout=output, stderr=subprocess.DEVNULL, timeout=60
		)
		if process.returncode or output.tell() > MAX_EXPANDED_BYTES:
			raise Rejected('pdf_content_check_failed')
		output.seek(0)
		check_text(output, '.txt')
	return 'pdf_visual_content_needs_manual_review'


def content_check(path):
	if 'id' not in CHECKS:
		return ''
	suffix = path.suffix.lower()
	if suffix == '.xlsx':
		return check_xlsx(path)
	if suffix in RASTER:
		return check_image(path)
	if suffix == '.pdf':
		return check_pdf(path)
	if 'type' in CHECKS and suffix not in TEXT_RESULT_SUFFIXES | {'.svg'}:
		raise Rejected('not_a_result_file')
	with path.open('rb') as handle:
		# Without the type filter, inspect additional text formats as well.
		# Opaque/binary content still fails closed when ID checking is enabled.
		if suffix == '.gz':
			check_gzip(handle, Path(path.stem).suffix.lower())
		else:
			check_text(handle, suffix)
	return ''


def inspect(root, relative):
	try:
		before = file_stat(root, relative)
		path = root / relative
		file_policy(path, before)
		reason = excluded_result(relative)
		if reason:
			raise Rejected(reason)
		review = content_check(path)
		if version(before) != version(file_stat(root, relative)):
			raise SourceChanged('source changed during content check; rerun: ' + relative)
		return before.st_size, review, version(before)
	except (Rejected, SourceChanged, FileNotFoundError):
		raise
	except ModuleNotFoundError as exc:
		raise RuntimeError('required content checker is unavailable: ' + str(exc.name)) from None
	except Exception as exc:
		# Never include file contents or parser exception messages in diagnostics.
		raise Rejected('content_check_failed_' + type(exc).__name__) from None


def project_list(arguments):
	projects = []
	for argument in arguments:
		projects.extend(p.strip() for p in argument.split(',') if p.strip())
	projects = list(dict.fromkeys(projects))
	if not projects:
		raise ValueError('provide an explicit project list')
	for project in projects:
		if (
			'/' in project
			or '\\' in project
			or project.startswith('.')
			or project.lower().endswith(MANIFEST_SUFFIX)
			or any(ord(c) < 32 for c in project)
		):
			raise ValueError('invalid project name')
		path = SOURCE / project
		if path.is_symlink() or not path.is_dir():
			raise ValueError(
				'selected project is not a real source directory: ' + project
			)
	return projects


def folder_list(arguments, projects):
	if not arguments:
		return []
	folders = list(dict.fromkeys(
		folder.strip() for argument in arguments for folder in argument.split(',') if folder.strip()
	))
	if not folders:
		raise ValueError('--folder requires a non-empty folder list')
	for folder in folders:
		parts = folder.split('/')
		if (
			any(not part or part.startswith('.') for part in parts)
			or '\\' in folder or any(ord(c) < 32 for c in folder)
		):
			raise ValueError('--folder must contain relative subdirectory names without . or ..')
		for project in projects:
			path = SOURCE / project
			for part in parts:
				path /= part
				if path.is_symlink() or not path.is_dir():
					raise ValueError('selected folder is not a real source directory: ' + project + '/' + folder)
	# Selecting an ancestor already covers its descendants. Avoid scanning or
	# replacing overlapping scopes twice, regardless of argument order.
	return [folder for folder in folders if not any(
		folder.startswith(other + '/') for other in folders if other != folder
	)]


def selected_scopes(projects, folders):
	return [project + '/' + folder for project in projects for folder in folders] if folders else projects


def prepare_roots():
	global SOURCE, DEST
	for path in (SOURCE, DEST):
		if not path.is_absolute() or path.is_symlink() or path == Path('/'):
			raise ValueError('roots must be absolute, non-symlink directories, not /')
	if not SOURCE.is_dir():
		raise ValueError('source root does not exist')
	SOURCE, DEST = (p.resolve() for p in (SOURCE, DEST))
	if SOURCE == DEST or SOURCE in DEST.parents or DEST in SOURCE.parents:
		raise ValueError('source and destination must be separate directories')
	if DEST.exists() and not DEST.is_dir():
		raise ValueError('destination root is not a directory')


def digest_file(path):
	digest = hashlib.sha256()
	with path.open('rb') as handle:
		for chunk in iter(lambda: handle.read(1024 * 1024), b''):
			digest.update(chunk)
	return digest.hexdigest()


def check_target(project):
	path = DEST / project
	current = DEST
	for part in Path(project).parts:
		current /= part
		if current.is_symlink() or (current.exists() and not current.is_dir()):
			raise ValueError('unsafe destination directory: ' + project)
		if (current / '.git').exists() or (current / '.git').is_symlink():
			raise ValueError('destination directory contains .git: ' + project)
	if not path.exists():
		return
	# Never replace a Git checkout, including one nested in an old copy.
	def walk_error(error):
		raise error
	for folder, directories, filenames in os.walk(path, followlinks=False, onerror=walk_error):
		if '.git' in directories or '.git' in filenames:
			raise ValueError('destination project contains .git: ' + project)


# 🚩 Select final results, shallow directories first

def result_paths(project, excluded, folders=()):
	queue = collections.deque(SOURCE / scope for scope in selected_scopes([project], folders))
	preferred = SOURCE / project / 'final'
	prefer_final = project == 'allele3' and preferred.is_dir() and not preferred.is_symlink()
	if prefer_final:
		print('  allele3: prefer final/ over matching old run directories; retain independent results.', flush=True)
	while queue:
		folder = queue.popleft()
		with os.scandir(folder) as entries:
			entries = sorted(entries, key=lambda e: (e.name.casefold(), e.name))
		children = []
		for entry in entries:
			path = Path(entry.path)
			relative = path.relative_to(SOURCE).as_posix()
			if entry.is_symlink():
				excluded['symlink'] += 1
			elif entry.is_dir(follow_symlinks=False):
				if (
					prefer_final and folder == SOURCE / project and path != preferred
					and (preferred / entry.name).is_dir()
				):
					excluded['earlier_run_directory'] += 1
				elif excluded_result_directory(relative):
					excluded['intermediate_directory'] += 1
				else:
					children.append(path)
			elif entry.is_file(follow_symlinks=False):
				yield relative
			else:
				excluded['not_a_regular_file'] += 1
		queue.extend(children)


def scan(projects, folders=()):
	print('Checks: ' + ','.join(name for name in CHECK_NAMES if name in CHECKS), flush=True)
	selected = {}
	for project in projects:
		print('Scanning ' + project + ' for shareable results...', flush=True)
		if folders:
			print('  Folder scope: ' + ','.join(folders), flush=True)
		excluded = collections.Counter()
		candidates, details = [], {}
		last_progress = time.monotonic()
		for relative in result_paths(project, excluded, folders):
			path = SOURCE / relative
			try:
				reason = excluded_artifact(path) or excluded_result(relative)
				if reason:
					raise Rejected(reason)
				info = file_stat(SOURCE, relative)
				file_policy(path, info)
			except Rejected as error:
				excluded[str(error)] += 1
				continue
			candidates.append(relative)
			details[relative] = info.st_size
		used, kept, seen, approved_png = 0, 0, set(), set()
		# At equal depth/name, try PNG before the corresponding PDF.
		candidates.sort(key=lambda r: (
			len(Path(r).parts), str(Path(r).with_suffix('')).casefold(),
			Path(r).suffix.lower() == '.pdf', r.casefold(), r,
		))
		for relative in candidates:
			path = SOURCE / relative
			stem = str(Path(relative).with_suffix(''))
			try:
				if path.suffix.lower() == '.pdf' and stem in approved_png:
					raise Rejected('pdf_already_shared_as_png')
				if 'size' in CHECKS and used + details[relative] > MAX_PROJECT_BYTES:
					print(f'  Budget skip: {relative} ({details[relative] / 1e6:.2f} MB)', flush=True)
					raise Rejected('project_size_limit')
				size, review, fingerprint = inspect(SOURCE, relative)
				sha = digest_file(path)
				if version(file_stat(SOURCE, relative)) != fingerprint:
					raise ValueError('source changed during scan: ' + relative)
				if sha in seen:
					raise Rejected('duplicate_content')
				if 'size' in CHECKS and used + size > MAX_PROJECT_BYTES:
					raise Rejected('project_size_limit')
			except Rejected as error:
				excluded[str(error)] += 1
				if str(error) == 'ukb_id_column':
					print('  Excluded ID-bearing table: ' + relative, flush=True)
				continue
			selected[relative] = (size, fingerprint, sha)
			seen.add(sha)
			used += size
			kept += 1
			if path.suffix.lower() == '.png':
				approved_png.add(stem)
			if review:
				print(f'  Visual review required: {relative} ({review})', flush=True)
			if time.monotonic() - last_progress >= 15:
				print(f'  Selected {kept} files, {used / 1e6:.2f} MB.', flush=True)
				last_progress = time.monotonic()
		budget = f' / {MAX_PROJECT_BYTES / 1e6:g}' if 'size' in CHECKS else ''
		print(f'  Selected: {kept} files, {used / 1e6:.2f}{budget} MB')
		for reason, count in sorted(excluded.items()):
			print(f'  Excluded {reason}: {count}')
		for relative in sorted((r for r in selected if r.split('/')[0] == project), key=lambda r: -selected[r][0])[:5]:
			print(f'  Largest: {selected[relative][0] / 1e6:.2f} MB  {relative}')
	return selected


def manifest_text(project, selected, folders=()):
	checks = ','.join(name for name in CHECK_NAMES if name in CHECKS)
	limits = f'{MAX_FILE_BYTES} bytes/file; {MAX_PROJECT_BYTES} bytes/project' if 'size' in CHECKS else 'disabled'
	types = 'Supported figures, workbooks and text tables' if 'type' in CHECKS else 'File type filter disabled'
	ids = 'ID-column checks include missing/mixed rows; no OCR of image pixels' if 'id' in CHECKS else 'ID-column checks disabled'
	return (
		'# Generated by github_copy.sh; selected project only.\n'
		f'# Source: {SOURCE}\n'
		'# Project scope: ' + project + '\n'
		'# Folder scope: ' + (','.join(folders) if folders else '(whole projects)') + '\n'
		f'# Checks: {checks}\n'
		f'# Limits: {limits}.\n'
		'# Priority: project root first, then increasing directory depth/name.\n'
		f'# {types}; no filename keyword filters.\n'
		f'# {ids}.\n'
		+ ''.join(relative + '\n' for relative in selected if relative.split('/', 1)[0] == project)
	)


# 🚩 Transactional replacement of selected sharing directories

@contextmanager
def destination_lock():
	key = hashlib.sha256(str(DEST).encode()).hexdigest()
	with open('/tmp/github_copy.' + key + '.lock', 'a') as handle:
		try:
			fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
		except BlockingIOError:
			raise ValueError('another github_copy process is using this destination') from None
		yield


def copy_files(projects, selected, folders=()):
	scopes = selected_scopes(projects, folders)
	for project in scopes:
		check_target(project)
	manifests = {project: DEST / (project + MANIFEST_SUFFIX) for project in projects}
	for manifest in manifests.values():
		if manifest.is_symlink() or (manifest.exists() and not manifest.is_file()):
			raise ValueError('unsafe manifest destination: ' + manifest.name)
	with staging_directory() as stage:
		payload = stage / 'payload'
		for project in scopes:
			(payload / project).mkdir(parents=True, exist_ok=True)
		for relative, (size, fingerprint, sha) in selected.items():
			if version(file_stat(SOURCE, relative)) != fingerprint:
				raise ValueError('source changed since scan; rerun: ' + relative)
			target = payload / relative
			target.parent.mkdir(parents=True, exist_ok=True)
			shutil.copy2(SOURCE / relative, target, follow_symlinks=False)
			if target.is_symlink() or target.stat().st_size != size or digest_file(target) != sha:
				raise ValueError('staged file differs from inspected source: ' + relative)
			if version(file_stat(SOURCE, relative)) != fingerprint:
				raise ValueError('source changed while staging: ' + relative)
		# The old copy is renamed, not duplicated. A short-lived sibling
		# transaction directory shares the destination filesystem, so project
		# renames and rollback are atomic. It is outside the sharing tree.
		DEST.mkdir(parents=True, exist_ok=True)
		transaction = Path(tempfile.mkdtemp(prefix='.github_copy.', dir=DEST.parent))
		keep = False
		try:
			incoming, backup = transaction / 'incoming', transaction / 'backup'
			shutil.copytree(payload, incoming, copy_function=shutil.copy2)
			for relative, (size, fingerprint, sha) in selected.items():
				if file_stat(incoming, relative).st_size != size or digest_file(incoming / relative) != sha:
					raise ValueError('destination staging verification failed: ' + relative)
			backup.mkdir()
			new_manifests = transaction / 'manifests'
			new_manifests.mkdir()
			for project, manifest in manifests.items():
				(new_manifests / manifest.name).write_text(
					manifest_text(project, selected, folders), encoding='utf-8', newline='\r\n'
				)
				if manifest.exists():
					shutil.copy2(manifest, backup / manifest.name)
			changed, changed_manifests = [], []
			try:
				for project in scopes:
					check_target(project)
					target = DEST / project
					target.parent.mkdir(parents=True, exist_ok=True)
					(backup / project).parent.mkdir(parents=True, exist_ok=True)
					changed.append(project)
					if target.exists():
						target.rename(backup / project)
					(incoming / project).rename(target)
				for manifest in manifests.values():
					changed_manifests.append(manifest)
					os.replace(new_manifests / manifest.name, manifest)
			except BaseException:
				failures = []
				for manifest in reversed(changed_manifests):
					try:
						if (backup / manifest.name).exists():
							os.replace(backup / manifest.name, manifest)
						else:
							manifest.unlink(missing_ok=True)
					except BaseException:
						failures.append(manifest.name)
				for project in reversed(changed):
					try:
						target = DEST / project
						if not (incoming / project).exists() and target.exists():
							target.rename(incoming / project)
						if (backup / project).exists():
							(backup / project).rename(target)
					except BaseException:
						failures.append(project)
				if failures:
					keep = True
					raise RecoveryFailed(f'rollback failed for {failures}; recovery copy retained at {transaction}') from None
				raise
		finally:
			if not keep:
				try:
					shutil.rmtree(transaction)
				except OSError:
					print(f'WARNING: remove leftover transaction directory: {transaction}', file=sys.stderr)
	print(f'COPY OK: {len(selected)} files, {sum(item[0] for item in selected.values()) / 1e6:.2f} MB.')
	print('Destination: ' + str(DEST))
	for manifest in manifests.values():
		print('Manifest: ' + str(manifest))
	print('Selected destination ' + ('folders' if folders else 'projects') + ' now contain only this selection; stale copies removed.')


def check_list(argument):
	names = [name.strip().lower() for name in argument.split(',')]
	if any(name not in CHECK_NAMES for name in names):
		raise argparse.ArgumentTypeError('--check requires a non-empty comma-separated list of size,type,id')
	return frozenset(names)


def main():
	global SOURCE, CHECKS
	arguments = sys.argv[1:]
	if not arguments or arguments[0] == 'help' or any(a in {'-h', '--help'} for a in arguments):
		print("""Modules:
    PROJECT[,PROJECT...]  Select sharing results, then rebuild selected destination projects

Main options:
    --dir DIR  Source root [GITHUB_COPY_SOURCE_ROOT or /mnt/d/analysis]
    --dat LIST  Project folders, comma-separated; positional names also accepted
    --folder LIST  Only these subfolders within each project, comma-separated
    --check LIST  Checks to enable, comma-separated [size,type,id]
                  size: file/project limits and empty files; type: file type filter
                  id: inspect content for seven-digit values in ID columns
    --dry-run  Read-only selection preview; no destination changes
    GITHUB_COPY_DEST_ROOT  Destination [/mnt/d/github/analysis]
    With size: 10 MB/file, 1 GB/project (decimal units); root files take priority
    With type: figures, XLSX and text tables/reports; exclude logs, RDS and unsupported types
    With id: uninspectable content is excluded; parser resource limits remain active
    Hidden files and intermediate directories are excluded regardless of --check
    With --folder, replace only selected subfolders; otherwise replace whole projects
    Write DEST/PROJECT.files.lst for each selected project; paths include PROJECT/
    -h, --help  Help; details: 0.README.md

Example:
    ./github_copy.sh --dir /mnt/d/analysis --dat grid --folder abm,csx,disco,final
    ./github_copy.sh --dat allele3 --dry-run
    ./github_copy.sh --dat allele3 --check size,type --dry-run
    ./github_copy.sh --dir /mnt/f/analysis/ --dat gu,ld
    ./github_copy.sh --dat allele3""")
		return
	parser = argparse.ArgumentParser(prog='github_copy.sh', add_help=False, allow_abbrev=False)
	parser.add_argument('--dir', metavar='DIR')
	parser.add_argument('--dat', metavar='LIST', action='append', default=[])
	parser.add_argument('--folder', metavar='LIST', action='append', default=[])
	parser.add_argument('--check', metavar='LIST', type=check_list, default=frozenset(CHECK_NAMES))
	parser.add_argument('--dry-run', action='store_true')
	parser.add_argument('projects', nargs='*', metavar='PROJECT')
	options = parser.parse_intermixed_args(arguments)
	CHECKS = options.check
	if options.dir is not None:
		if not options.dir.strip():
			raise ValueError('--dir requires a non-empty source directory')
		SOURCE = Path(options.dir)
	prepare_roots()
	projects = project_list(options.dat + options.projects)
	folders = folder_list(options.folder, projects)
	if options.dry_run:
		selected = scan(projects, folders)
		print(f'DRY RUN: {len(selected)} files, {sum(item[0] for item in selected.values()) / 1e6:.2f} MB; destination unchanged.')
		return
	with destination_lock():
		selected = scan(projects, folders)
		copy_files(projects, selected, folders)


def interrupted(signum, frame):
	raise InterruptedError('interrupted; destination transaction will be rolled back')


signal.signal(signal.SIGTERM, interrupted)
signal.signal(signal.SIGHUP, interrupted)
try:
	main()
except KeyboardInterrupt:
	print('ERROR: interrupted', file=sys.stderr)
	sys.exit(130)
except Exception as error:
	# Never print participant values or raw parser diagnostics.
	print('ERROR: ' + str(error), file=sys.stderr)
	sys.exit(1)
PY
