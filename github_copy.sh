#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi


# 🚩 github_copy
# Scan analysis results for UKB participant IDs before publishing to GitHub.
# Scan first, then copy the approved list automatically.
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
import stat
import subprocess
import sys
import tempfile
import time
import zipfile
from xml.etree import ElementTree as ET

SOURCE = Path(os.environ.get('GITHUB_COPY_SOURCE_ROOT', '/mnt/d/analysis'))
DEST = Path(os.environ.get('GITHUB_COPY_DEST_ROOT', '/mnt/d/github/analysis'))
MANIFEST = 'github_files.lst'
# 50 decimal MB; files exactly at the limit are allowed.
MAX_FILE_BYTES = 50_000_000
MAX_EXPANDED_BYTES = 512 * 1024 * 1024
MAX_XML_BYTES = 256 * 1024 * 1024
# User-confirmed projects without UKB participants. Public 1KG individuals
# and their formal results may be published without content/ID inspection.
NON_UKB_PROJECTS = frozenset({'gu'})
MAX_NON_UKB_FILES = 10_000
INTERMEDIATE_DIRECTORIES = frozenset({
	'raw', 'mask', 'masks', 'cache', 'extract', 'infer', 'inputs',
	'tmp', 'temp', 'normalize', 'storage-migration',
})

# A UKB ID column must satisfy BOTH user-defined conditions:
# its name contains "id" (case-insensitive), and every data row is a
# seven-digit integer. Numeric Excel/CSV exports may use .0 or scientific form.
SEVEN_DIGIT_ID = re.compile(r'[0-9]{7}(?:\.0+)?')
NUMERIC_ID = re.compile(r'[0-9]+(?:\.[0-9]+)?[eE][+-]?[0-9]+')
EXCEL_ESCAPE = re.compile(r'_x([0-9a-fA-F]{4})_')
BAD_CONTROL = re.compile(r'[\x00-\x08\x0b\x0e-\x1f]')
RASTER = {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.tif', '.tiff', '.bmp'}
RESULT_SUFFIXES = RASTER | {'.xlsx', '.xlsm', '.pdf', '.svg'}
SHEET_NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


class Rejected(Exception):
	def __init__(self, reason, locations=()):
		super().__init__(reason)
		self.locations = locations


class RecoveryFailed(Exception):
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
		self.rows = 0

	def add(self, values):
		self.rows += 1
		self.candidates = {
			i for i in self.candidates if i < len(values) and is_seven_digit_id(values[i])
		}

	def finish(self):
		if self.rows and self.candidates:
			raise Rejected('ukb_id_column')


def check_value(value):
	# JSON tables may be arrays of records, arrays with a header row, or
	# objects of column arrays. Isolated scalar values do not define a column.
	if isinstance(value, dict):
		for name, items in value.items():
			if (
				'id' in str(name).lower() and isinstance(items, (list, tuple))
				and items and all(is_seven_digit_id(item) for item in items)
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
def excluded_result_directory(relative):
	parts = relative.split('/')
	if parts[0] in NON_UKB_PROJECTS and any(
		part.lower() in INTERMEDIATE_DIRECTORIES for part in parts[1:]
	):
		return 'intermediate_directory'
	return ''


def excluded_result(relative):
	parts = relative.split('/')
	if parts[0] not in NON_UKB_PROJECTS:
		return ''
	reason = excluded_result_directory('/'.join(parts[:-1]))
	if reason:
		return reason
	path = Path(relative)
	if path.suffix.lower() in RESULT_SUFFIXES:
		return ''
	return 'not_a_result_file'


def check_file_count(files):
	counts = collections.Counter(p.split('/')[0] for p in files)
	for project in sorted(NON_UKB_PROJECTS):
		if counts[project] >= MAX_NON_UKB_FILES:
			raise ValueError(
				f'{project}: {counts[project]} selected files reaches the '
				f'{MAX_NON_UKB_FILES}-file limit; narrow the result selection first'
			)


def excluded_artifact(path):
	name = path.name.lower()
	if name.startswith('~$'):
		return 'office_lock_file'
	if name.endswith('.gz'):
		name = name[:-3]
	# RDS may contain individual IDs and is opaque to this scanner. Keep
	# it out of scan lists even in confirmed non-UKB projects.
	if name.endswith('.rds'):
		return 'rds_file'
	if name.endswith('.log'):
		return 'log_file'
	if name.endswith('.done'):
		return 'done_file'
	if name in {'plots.pdf', 'rplots.pdf'}:
		return 'unnamed_plot'
	return ''


def file_policy(path, details):
	if details.st_size > MAX_FILE_BYTES:
		raise Rejected('over_50_MB')
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

	def handle_starttag(self, tag, attrs):
		if tag == 'table':
			self.tables.append({'columns': None, 'row': None, 'cell': None})
		elif self.tables:
			table = self.tables[-1]
			if tag == 'tr':
				table['row'] = []
			elif tag in {'td', 'th'} and table['row'] is not None:
				table['cell'] = []

	def handle_endtag(self, tag):
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
		# Headered text tables also occur in .psam/.sscore/.out and unknown
		# suffixes. Headerless data and free text cannot satisfy both criteria.
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


def content_check(path, check_ids=True):
	suffix = path.suffix.lower()
	if not check_ids:
		# Eligibility relies on the explicit non-UKB project confirmation,
		# not on treating public genomic coordinates as participant IDs.
		# This also admits workbooks with embedded exact result tables.
		if suffix == '.pdf':
			return 'pdf_visual_content_needs_manual_review'
		return ''
	if suffix in {'.xlsx', '.xlsm'}:
		return check_xlsx(path)
	if suffix in RASTER:
		return check_image(path)
	if suffix == '.pdf':
		return check_pdf(path)
	if suffix == '.gz':
		# Unpack outside the repository, impose an expansion bound, and run
		# exactly the same detector on compressed and uncompressed text.
		inner_suffix = Path(path.stem).suffix.lower()
		if inner_suffix == '.rds':
			raise Rejected('rds_file')
		with path.open('rb') as raw:
			check_gzip(raw, inner_suffix)
		return ''
	# Text sniffing instead of a narrow extension whitelist: .psam, .fam,
	# .sscore, .out, .R, custom suffixes and extensionless files all work.
	with path.open('rb') as handle:
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
		check_ids = relative.split('/')[0] not in NON_UKB_PROJECTS
		review = content_check(path, check_ids)
		if version(before) != version(file_stat(root, relative)):
			raise Rejected('file_changed_during_check')
		return before.st_size, review
	except Rejected:
		raise
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
			or project in {'.', '..'}
			or any(ord(c) < 32 for c in project)
		):
			raise ValueError('invalid project name')
		path = SOURCE / project
		if path.is_symlink() or not path.is_dir():
			raise ValueError(
				'selected project is not a real source directory: ' + project
			)
	return projects


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
	DEST.mkdir(parents=True, exist_ok=True)


def write_manifest(path, text):
	if path.is_symlink() or (path.exists() and not path.is_file()):
		raise ValueError('unsafe manifest destination')
	with staging_directory() as stage:
		temporary, backup = stage / 'manifest', stage / 'backup'
		# Match the Windows Git checkout convention for this generated list.
		temporary.write_text(text, encoding='utf-8', newline='\r\n')
		os.chmod(temporary, 0o644)
		existed = path.exists()
		if existed:
			shutil.copy2(path, backup)
		try:
			shutil.copy2(temporary, path)
		except BaseException:
			try:
				if existed:
					shutil.copy2(backup, path)
				else:
					path.unlink(missing_ok=True)
			except BaseException:
				raise RecoveryFailed(
					f'manifest restoration failed; recovery files retained at {stage}'
				) from None
			raise


def print_summary(projects, approved, excluded, reviews, totals):
	summary = {
		'projects': projects,
		'approved_files': len(approved),
		'approved_bytes': sum(t['bytes'] for t in totals.values()),
		'manual_review_files': len(reviews),
		'excluded_reasons': dict(collections.Counter(r for _, r in excluded)),
		'by_project': totals,
		'max_file_bytes': MAX_FILE_BYTES,
	}
	for project in projects:
		t = totals[project]
		print(
			f'  {project}: {t["files"]} files, {t["bytes"] / 1048576:.2f} MiB, {t["review"]} visual reviews'
		)
	print(
		f'  Total: {summary["approved_files"]} files, {summary["approved_bytes"] / 1048576:.2f} MiB'
	)
	for reason, count in sorted(summary['excluded_reasons'].items()):
		print(f'  Excluded {reason}: {count}')
	for path, reason in reviews:
		print(f'  Manual review {path}: {reason}')


def totals_for(projects):
	return {p: {'files': 0, 'bytes': 0, 'review': 0} for p in projects}


# 🚩 Scan analysis results
def scan(projects):
	approved, excluded, reviews = [], [], []
	totals = totals_for(projects)
	checked = 0
	last_progress = time.monotonic()
	for project in projects:
		print('Scanning ' + project + ' recursively...', flush=True)
		if project in NON_UKB_PROJECTS:
			print(
				'  User-confirmed non-UKB project: content/ID inspection skipped; '
				'formal workbooks and figures only; RDS excluded. Size/path/artifact rules still apply.',
				flush=True,
			)

		def walk_error(error):
			raise error

		for folder, directories, filenames in os.walk(
			SOURCE / project, onerror=walk_error, followlinks=False
		):
			keep = []
			for name in sorted(directories):
				path = Path(folder) / name
				relative = path.relative_to(SOURCE).as_posix()
				reason = excluded_result_directory(relative)
				if path.is_symlink():
					excluded.append(
						(relative + '/', 'symlink')
					)
				elif reason:
					excluded.append((relative + '/', reason))
				else:
					keep.append(name)
			directories[:] = keep
			for name in sorted(filenames):
				relative = (Path(folder) / name).relative_to(SOURCE).as_posix()
				checked += 1
				try:
					size, review = inspect(SOURCE, relative)
				except Rejected as exc:
					excluded.append((relative, str(exc)))
					for kind, label in exc.locations:
						print(
							f'WARNING: {relative!r}: {kind} {label!r} contains a UKB ID column '
							'(name includes id; all data rows are seven-digit integers). '
							'File excluded; review after scan.',
							flush=True,
						)
				else:
					approved.append(relative)
					totals[project]['files'] += 1
					totals[project]['bytes'] += size
					if project in NON_UKB_PROJECTS and totals[project]['files'] >= MAX_NON_UKB_FILES:
						check_file_count(approved)
					if review:
						reviews.append((relative, review))
						totals[project]['review'] += 1
				if time.monotonic() - last_progress > 15:
					print(
						f'  Inspected {checked} files; admitted {len(approved)}.',
						flush=True,
					)
					last_progress = time.monotonic()
	approved.sort()
	check_file_count(approved)
	header = (
		'# Generated by github_copy.sh. Paths are relative to '
		+ str(SOURCE)
		+ '.\n'
		'# Project scope: ' + ','.join(projects) + '\n'
		'# Exclude .rds, .log, .done, unnamed plots.pdf/Rplots.pdf and files over 50 MB in every project.\n'
		'# UKB ID columns: name contains id (case-insensitive) AND every data row is a 7-digit integer.\n'
		'# Isolated numbers and filenames do not identify UKB IDs. Unreadable content remains excluded.\n'
		'# Workbook source exports are checksum-verified and content-scanned; Office lock files are excluded.\n'
		'# Non-UKB projects (content/ID inspection skipped; formal results only): '
		+ ','.join(p for p in projects if p in NON_UKB_PROJECTS) + '\n'
		'# Non-UKB results: XLSX/XLSM and figures only; intermediate directories excluded.\n'
		f'# Non-UKB selection must contain fewer than {MAX_NON_UKB_FILES} files per project.\n'
		'# RDS/log/done exclusions include gzip files.\n'
		'# Image pixels are not inspected and images do not require manual review.\n'
		'# After scanning, listed files are copied automatically without rechecks or confirmation.\n'
	)
	print(f'  Inspected {checked} files recursively.')
	print_summary(projects, approved, excluded, reviews, totals)
	write_manifest(DEST / MANIFEST, header + ''.join(p + '\n' for p in approved))
	print(
		f'SCAN OK: wrote {len(approved)} approved paths to {DEST / MANIFEST}',
		flush=True,
	)
	return approved


# 🚩 Copy approved results
def read_manifest(arguments):
	path = DEST / MANIFEST
	if path.is_symlink() or not path.is_file():
		raise ValueError('manifest missing or unsafe; rerun with the project list')
	lines = path.read_text(encoding='utf-8').splitlines()
	scopes = [
		line[len('# Project scope: ') :]
		for line in lines
		if line.startswith('# Project scope: ')
	]
	if len(scopes) != 1:
		raise ValueError(
			'manifest must contain exactly one project scope; rerun with the project list'
		)
	scope = project_list(scopes)
	projects = project_list(arguments) if arguments else scope
	if not set(projects) <= set(scope):
		raise ValueError('copy project is outside the manifest scope')
	files, seen = [], set()
	for line in lines:
		if not line or line.startswith('#'):
			continue
		parts = relative_parts(line)
		if parts[0] not in scope or line in seen:
			raise ValueError('invalid scope or duplicate manifest entry')
		seen.add(line)
		if parts[0] in projects:
			files.append(line)
	return projects, files


def check_target(project):
	path = DEST / project
	if path.is_symlink() or (path.exists() and not path.is_dir()):
		raise ValueError('unsafe destination project: ' + project)


def copy_files(arguments):
	projects, files = read_manifest(arguments)
	totals = totals_for(projects)
	staged_files, missing = [], []
	for project in projects:
		check_target(project)
	# The generated scan list is the sole selection. Do not reapply size,
	# artifact or content/ID checks here. Stage copies and backups in /tmp
	# so a failed copy can be rolled back without rebuilding whole projects.
	with staging_directory() as stage:
		payload, backup = stage / 'payload', stage / 'backup'
		payload.mkdir()
		backup.mkdir()
		last_progress = time.monotonic()
		for index, relative in enumerate(files, 1):
			try:
				file_stat(SOURCE, relative)
			except FileNotFoundError:
				missing.append(relative)
			else:
				target = payload / relative
				target.parent.mkdir(parents=True, exist_ok=True)
				shutil.copy2(SOURCE / relative, target, follow_symlinks=False)
				staged_files.append(relative)
				project = relative.split('/')[0]
				totals[project]['files'] += 1
				totals[project]['bytes'] += target.stat().st_size
			if time.monotonic() - last_progress > 15:
				print(f'  Prepared {len(staged_files)} copies from {index}/{len(files)} list entries.', flush=True)
				last_progress = time.monotonic()
		for relative in staged_files:
			try:
				file_stat(DEST, relative)
			except FileNotFoundError:
				continue
			saved = backup / relative
			saved.parent.mkdir(parents=True, exist_ok=True)
			shutil.copy2(DEST / relative, saved, follow_symlinks=False)
		changed = []
		try:
			for relative in staged_files:
				target = DEST / relative
				target.parent.mkdir(parents=True, exist_ok=True)
				changed.append(relative)
				shutil.copy2(payload / relative, target, follow_symlinks=False)
		except BaseException:
			failed = []
			for relative in reversed(changed):
				try:
					target, saved = DEST / relative, backup / relative
					if saved.exists():
						shutil.copy2(saved, target, follow_symlinks=False)
					else:
						target.unlink(missing_ok=True)
				except BaseException:
					failed.append(relative)
			if failed:
				raise RecoveryFailed(
					'restoration failed for ' + ', '.join(failed)
					+ f'; recovery files retained at {stage}'
				) from None
			raise
	for project, total in totals.items():
		print(f'  {project}: {total["files"]} files copied, {total["bytes"] / 1048576:.2f} MiB')
	if missing:
		print(f'WARNING: skipped {len(missing)} missing source files; manifest unchanged.')
		for relative in missing[:5]:
			print('  Missing source: ' + relative)
		if len(missing) > 5:
			print(f'  ... and {len(missing) - 5} more missing source files.')
	print(f'COPY OK: {len(staged_files)} files copied from the list; no content recheck.')


def main():
	global SOURCE
	arguments = sys.argv[1:]
	if not arguments or arguments[0] == 'help' or any(a in {'-h', '--help'} for a in arguments):
		print(f"""Modules:
    PROJECT[,PROJECT...]  Scan results, write github_files.lst and copy approved files

Main options:
    --dir DIR  Source root [GITHUB_COPY_SOURCE_ROOT or /mnt/d/analysis]
    --dat LIST  Project folders, comma-separated; positional names also accepted
    GITHUB_COPY_DEST_ROOT    Destination [/mnt/d/github/analysis]
    -h, --help  Show help
    Copy overwrites listed destinations; no Git commit/push. Rules: 0.README.md

Example:
    ./github_copy.sh --dir /mnt/f/analysis/ --dat gu,ld
    ./github_copy.sh le8 grid gu
    ./github_copy.sh gu""")
		return
	parser = argparse.ArgumentParser(prog='github_copy.sh', add_help=False, allow_abbrev=False)
	parser.add_argument('--dir', metavar='DIR')
	parser.add_argument('--dat', metavar='LIST', action='append', default=[])
	parser.add_argument('projects', nargs='*', metavar='PROJECT')
	options = parser.parse_intermixed_args(arguments)
	if options.dir is not None:
		if not options.dir.strip():
			raise ValueError('--dir requires a non-empty source directory')
		SOURCE = Path(options.dir)
	arguments = options.dat + options.projects
	if not any(p.strip() for argument in arguments for p in argument.split(',')):
		raise ValueError('provide an explicit project list with --dat or positional names')
	prepare_roots()
	projects = project_list(arguments)
	if not scan(projects):
		print('No approved files; nothing to copy.')
	else:
		copy_files(projects)


try:
	main()
except Exception as error:
	# Parser exceptions are redacted by inspect; no participant values printed.
	print('ERROR: ' + str(error), file=sys.stderr)
	sys.exit(1)
PY
