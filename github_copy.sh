#!/usr/bin/env bash


# 🚩 github_copy
# Scan analysis results for UKB participant IDs before publishing to GitHub.
# The copy list contains ONLY approved files. scan never copies analysis files.
set -Eeuo pipefail
command -v python3 >/dev/null 2>&1 || {
	echo 'ERROR: python3 is required' >&2
	exit 1
}
exec python3 - "$@" <<'PY'
import codecs
import collections
import decimal
import gzip
import html
from html.parser import HTMLParser
import io
import json
import os
from pathlib import Path
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
# 10 decimal MB; files exactly at the limit are allowed.
MAX_FILE_BYTES = 10_000_000
MAX_EXPANDED_BYTES = 512 * 1024 * 1024
MAX_XML_BYTES = 256 * 1024 * 1024

# Check VALUES regardless of header, extension, project, or row position.
# Boundaries avoid matching rs1234567, longer integers and the fractional part of
# 0.1234567. Also catch integers exported as 2545852.0 or 2.545852e+06.
# Seven-digit counts/coordinates are conservatively excluded too.
ID = re.compile(r'(?<![A-Za-z0-9.])[0-9]{7}(?:\.0+)?(?![A-Za-z0-9]|\.[0-9])')
SCIENTIFIC = re.compile(
	r'(?<![A-Za-z0-9.])[0-9]+(?:\.[0-9]+)?[eE][+-]?[0-9]+(?![A-Za-z0-9.])'
)
ESCAPED_CHAR = re.compile(r'\\(?:u[0-9a-fA-F]{4}|U[0-9a-fA-F]{8}|x[0-9a-fA-F]{2})')
BAD_CONTROL = re.compile(r'[\x00-\x08\x0b\x0e-\x1f]')
RASTER = {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.tif', '.tiff', '.bmp'}


class Rejected(Exception):
	pass


def has_id(value):
	text = str(value)
	if '\\' in text:
		text = ESCAPED_CHAR.sub(
			lambda m: (
				chr(int(m.group()[2:], 16))
				if int(m.group()[2:], 16) <= 0x10FFFF
				else m.group()
			),
			text,
		)
	if ID.search(text):
		return True
	for match in SCIENTIFIC.finditer(text):
		token = match.group()
		if len(token) > 64:
			continue
		try:
			number = decimal.Decimal(token)
			if (
				1_000_000 <= number <= 9_999_999
				and number == number.to_integral_value()
			):
				return True
		except decimal.InvalidOperation:
			continue
	return False


def check_value(value):
	if isinstance(value, dict):
		for key, item in value.items():
			check_value(key)
			check_value(item)
	elif isinstance(value, (list, tuple)):
		for item in value:
			check_value(item)
	elif value is not None and has_id(value):
		raise Rejected('ukb_id')


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
	if has_id(relative):
		raise Rejected('ukb_id_in_path')
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
def excluded_artifact(path):
	name = path.name.lower()
	if name.endswith('.gz'):
		name = name[:-3]
	if name.endswith('.log'):
		return 'log_file'
	if name.endswith('.done'):
		return 'done_file'
	if name in {'plots.pdf', 'rplots.pdf'}:
		return 'unnamed_plot'
	return ''


def file_policy(path, details):
	reason = excluded_artifact(path)
	if reason:
		raise Rejected(reason)
	# Rejected RDS/oversized files never enter the copy list, so reading
	# their participant records would add no safety and waste scan time.
	if path.suffix.lower() == '.rds':
		raise Rejected('rds_file')
	if details.st_size > MAX_FILE_BYTES:
		raise Rejected('over_10_MB')


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
	"""Check rendered text too: inline spans may split a visible identifier."""

	breaks = {'p', 'div', 'tr', 'td', 'th', 'br', 'hr', 'li', 'text', 'script', 'style'}

	def __init__(self):
		super().__init__(convert_charrefs=True)
		self.tail = ''

	def handle_starttag(self, tag, attrs):
		check_value(attrs)
		if tag in self.breaks:
			self.tail = ''

	def handle_endtag(self, tag):
		if tag in self.breaks:
			self.tail = ''

	def handle_data(self, data):
		combined = self.tail + data
		check_value(combined)
		self.tail = combined[-256:]


def check_text(raw, suffix):
	with text_handle(raw) as handle:
		if suffix == '.json':
			# Decode escaped strings and dictionary keys as well as numeric IDs.
			# Preserve duplicate keys; a later value must not hide an earlier ID.
			check_value(json.load(handle, object_pairs_hook=list))
			return
		if suffix == '.jsonl':
			for line in handle:
				if line.strip():
					check_value(json.loads(line, object_pairs_hook=list))
			return
		# Full-file scan, including headerless/unknown-extension files, all
		# columns, multiline fields, logs, scripts, and embedded HTML/JS data.
		markup = (
			MarkupChecker() if suffix in {'.html', '.htm', '.svg', '.xml'} else None
		)
		for line in handle:
			if BAD_CONTROL.search(line):
				raise Rejected('uninspectable_binary_content')
			check_value(html.unescape(line) if markup else line)
			if markup:
				markup.feed(line)
		if markup:
			markup.close()


def check_printer_settings(data):
	# Excel/openxlsx printer metadata may be binary or a hexadecimal dump.
	# Inspect its readable strings too; other embedded binary objects stay blocked.
	if len(data) > 1024 * 1024:
		raise Rejected('xlsx_printer_settings_too_large')
	check_value(data.decode('latin-1'))
	if re.fullmatch(rb'(?:[0-9a-fA-F]{2}\s*)+', data.strip()):
		data = bytes.fromhex(data.decode('ascii'))
	check_value(data.decode('latin-1').replace('\x00', ''))
	for encoding in ('utf-16-le', 'utf-16-be'):
		for offset in (0, 1):
			check_value(data[offset:].decode(encoding, errors='ignore'))


def check_xml_attributes(element):
	attributes = element.attrib
	if element.tag.startswith(
		'{http://schemas.openxmlformats.org/drawingml/2006/main}'
	):
		# DrawingML expresses angles in 1/60000 degree and lengths in EMUs.
		# These typed layout attributes (e.g. a shadow direction) are not cell
		# data. Still inspect labels/names and every text node, chart and cell.
		geometry = {
			'ang',
			'dir',
			'rot',
			'stAng',
			'endAng',
			'fov',
			'lat',
			'lon',
			'rev',
			'x',
			'y',
			'cx',
			'cy',
			'w',
			'h',
			'dist',
			'blurRad',
			'l',
			'r',
			't',
			'b',
			'sx',
			'sy',
			'kx',
			'ky',
		}
		attributes = {
			key: value
			for key, value in attributes.items()
			if key not in geometry or not re.fullmatch(r'[+-]?[0-9]+', value)
		}
	check_value(attributes)


def check_xlsx(path):
	visual = False
	with zipfile.ZipFile(path) as book:
		infos = book.infolist()
		if sum(item.file_size for item in infos) > MAX_EXPANDED_BYTES or any(
			item.file_size > MAX_XML_BYTES for item in infos
		):
			raise Rejected('xlsx_expansion_limit')
		names = book.namelist()
		if not any(re.fullmatch(r'xl/worksheets/[^/]+\.xml', n) for n in names):
			raise Rejected('xlsx_content_check_failed')
		# Inspect all XML, including unused shared strings, hidden sheets,
		# comments, formula caches and relationships, not just visible cells.
		for name in names:
			check_value(name)
			if name.endswith('/'):
				continue
			if name.lower().endswith(('.xml', '.rels', '.vml')):
				with book.open(name) as handle:
					rich_text = None
					for event, element in ET.iterparse(handle, events=('start', 'end')):
						tag = element.tag.rsplit('}', 1)[-1]
						if event == 'start':
							if tag in {'si', 'is'}:
								rich_text = []
							continue
						check_value(element.text)
						check_value(element.tail)
						check_xml_attributes(element)
						if rich_text is not None and tag == 't':
							rich_text.append(element.text or '')
						if tag in {'si', 'is'} and rich_text is not None:
							check_value(''.join(rich_text))
							rich_text = None
						element.clear()
			elif Path(name).suffix.lower() in RASTER:
				with book.open(name) as handle:
					check_image(io.BytesIO(handle.read()))
				visual = True
			elif re.fullmatch(r'xl/printerSettings/printerSettings[0-9]+\.bin', name):
				check_printer_settings(book.read(name))
			else:
				# Embedded objects, macros and opaque payloads cannot bypass
				# the detector by living inside an otherwise readable workbook.
				raise Rejected('xlsx_uninspectable_embedded_content')
	return 'xlsx_embedded_images' if visual else ''


def check_image(path):
	# Pixels need manual review, as requested by the user. Check metadata and
	# validate the actual format so renaming a data file .png cannot bypass us.
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
	return 'image_pixels_need_manual_review'


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
	with tempfile.TemporaryFile() as output:
		process = subprocess.run(
			command, stdout=output, stderr=subprocess.DEVNULL, timeout=60
		)
		if process.returncode or output.tell() > MAX_EXPANDED_BYTES:
			raise Rejected('pdf_content_check_failed')
		output.seek(0)
		check_text(output, '.txt')
	return 'pdf_visual_content_needs_manual_review'


def content_check(path):
	suffix = path.suffix.lower()
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
		with tempfile.TemporaryFile() as expanded:
			with gzip.open(path, 'rb') as handle:
				while True:
					chunk = handle.read(1024 * 1024)
					if not chunk:
						break
					if expanded.tell() + len(chunk) > MAX_EXPANDED_BYTES:
						raise Rejected('gzip_expansion_limit')
					expanded.write(chunk)
			expanded.seek(0)
			check_text(expanded, inner_suffix)
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
		review = content_check(path)
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


def atomic_text(path, text):
	descriptor, temporary = tempfile.mkstemp(
		prefix='.' + path.name + '.', dir=path.parent
	)
	try:
		with os.fdopen(descriptor, 'w', encoding='utf-8', newline='\n') as handle:
			handle.write(text)
		os.chmod(temporary, 0o644)
		os.replace(temporary, path)
	finally:
		if os.path.exists(temporary):
			os.unlink(temporary)


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

		def walk_error(error):
			raise error

		for folder, directories, filenames in os.walk(
			SOURCE / project, onerror=walk_error, followlinks=False
		):
			keep = []
			for name in sorted(directories):
				path = Path(folder) / name
				if path.is_symlink():
					excluded.append(
						(path.relative_to(SOURCE).as_posix() + '/', 'symlink')
					)
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
				else:
					approved.append(relative)
					totals[project]['files'] += 1
					totals[project]['bytes'] += size
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
	header = (
		'# Generated by github_copy.sh scan. Paths are relative to '
		+ str(SOURCE)
		+ '.\n'
		'# Project scope: ' + ','.join(projects) + '\n'
		'# Exclude .rds, .log, .done, unnamed plots.pdf/Rplots.pdf, files over 10 MB, and UKB IDs (7 digits).\n'
		'# Log/done exclusions include gzip files; unreadable content is excluded.\n'
		'# Images retained by request; visual review list is printed to the terminal.\n'
		'# sync rechecks an exact staged snapshot before changing destination projects.\n'
	)
	print_summary(projects, approved, excluded, reviews, totals)
	atomic_text(DEST / MANIFEST, header + ''.join(p + '\n' for p in approved))
	print(
		f'SCAN OK: wrote {len(approved)} approved paths to {DEST / MANIFEST}',
		flush=True,
	)


# 🚩 Sync approved results
def read_manifest(arguments):
	path = DEST / MANIFEST
	if path.is_symlink() or not path.is_file():
		raise ValueError('manifest missing or unsafe; run scan first')
	lines = path.read_text(encoding='utf-8').splitlines()
	scopes = [
		line[len('# Project scope: ') :]
		for line in lines
		if line.startswith('# Project scope: ')
	]
	if len(scopes) != 1:
		raise ValueError(
			'manifest must contain exactly one project scope; run scan again'
		)
	scope = project_list(scopes)
	projects = project_list(arguments) if arguments else scope
	if not set(projects) <= set(scope):
		raise ValueError('sync project is outside the manifest scope')
	files, excluded, seen = [], [], set()
	for line in lines:
		if not line or line.startswith('#'):
			continue
		parts = relative_parts(line)
		if parts[0] not in scope or line in seen:
			raise ValueError('invalid scope or duplicate manifest entry')
		seen.add(line)
		if parts[0] in projects:
			reason = excluded_artifact(Path(line))
			if reason:
				excluded.append((line, reason))
			else:
				files.append(line)
	return projects, files, excluded


def check_target(project):
	path = DEST / project
	if path.is_symlink() or (path.exists() and not path.is_dir()):
		raise ValueError('unsafe destination project: ' + project)


def sync(arguments):
	projects, files, excluded = read_manifest(arguments)
	totals = totals_for(projects)
	reviews = []
	represented = sorted({p.split('/')[0] for p in files})
	for project in projects:
		check_target(project)
	# Stage beside, never inside, the publication directory. Same filesystem
	# allows checked bytes to be renamed into place instead of being recopied.
	with tempfile.TemporaryDirectory(
		prefix='.github_copy.stage.', dir=DEST.parent
	) as directory:
		stage = Path(directory)
		payload, backup = stage / 'payload', stage / 'backup'
		payload.mkdir()
		backup.mkdir()
		last_progress = time.monotonic()
		for index, relative in enumerate(files, 1):
			before = file_stat(SOURCE, relative)
			file_policy(SOURCE / relative, before)
			target = payload / relative
			target.parent.mkdir(parents=True, exist_ok=True)
			shutil.copy2(SOURCE / relative, target, follow_symlinks=False)
			if version(before) != version(file_stat(SOURCE, relative)):
				raise Rejected('source_changed_during_copy')
			size, review = inspect(payload, relative)
			project = relative.split('/')[0]
			totals[project]['files'] += 1
			totals[project]['bytes'] += size
			if review:
				reviews.append((relative, review))
				totals[project]['review'] += 1
			if time.monotonic() - last_progress > 15:
				print(f'  Staged and rechecked {index}/{len(files)} files.', flush=True)
				last_progress = time.monotonic()
		# Complete every check before touching existing destination projects.
		# Retain originals until every selected replacement has succeeded.
		backed_up, installed, removed = [], [], []
		try:
			for project in represented:
				check_target(project)
				target = DEST / project
				if target.exists():
					os.replace(target, backup / project)
					backed_up.append(project)
				os.replace(payload / project, target)
				installed.append(project)
			# A project with no approved entries keeps its results, but old
			# logs, completion markers and unnamed PDFs must still be removed.
			for project in sorted(set(projects) - set(represented)):
				check_target(project)
				for folder, directories, filenames in os.walk(DEST / project):
					directories[:] = [
						name
						for name in directories
						if not (Path(folder) / name).is_symlink()
					]
					for name in filenames:
						target = Path(folder) / name
						if excluded_artifact(target):
							saved = backup / target.relative_to(DEST)
							saved.parent.mkdir(parents=True, exist_ok=True)
							os.replace(target, saved)
							removed.append((target, saved))
		except BaseException:
			for target, saved in reversed(removed):
				os.replace(saved, target)
			for project in reversed(installed):
				os.replace(DEST / project, payload / project)
			for project in reversed(backed_up):
				os.replace(backup / project, DEST / project)
			raise
	print_summary(projects, files, excluded, reviews, totals)
	print(f'SYNC OK: {len(represented)} projects rebuilt, {len(files)} files copied.')
	untouched = set(projects) - set(represented)
	if untouched:
		print(
			'No approved entries; retained existing results and removed '
			f'{len(removed)} excluded artifacts: ' + ','.join(sorted(untouched))
		)


def main():
	arguments = sys.argv[1:]
	if not arguments or arguments[0] in {'-h', '--help', 'help'}:
		print("""Usage:
  ./github_copy.sh scan le8,grid,abm,abm_TF
  ./github_copy.sh sync [le8,grid,abm,abm_TF]

scan recursively builds github_files.lst without copying analysis files.
Exclude .rds, .log, .done, unnamed plots.pdf/Rplots.pdf, files over 10 MB
(10,000,000 bytes), and files containing UKB IDs. Log/done exclusions include gzip.
UKB IDs of 7 digits are checked in paths and content, regardless of header or row.
Ambiguous 7-digit integers are conservatively excluded as possible IDs.
Other file names and directories remain eligible for content inspection.
Readable text, Excel, PDF and gzip text are inspected; unreadable content is excluded.
Images are retained with a manual-review list printed to the terminal (Pillow required).
Every admitted file is <=10 MB; there is no file-count or directory-depth cap.
sync rechecks staged contents, then rebuilds only projects with selected entries.
Excluded logs, completion markers and unnamed PDFs are skipped even in old lists.
Projects without approved entries retain existing results except those artifacts.
No Git commit or push is performed. PDF inspection uses pdftotext or Ghostscript.
Roots: GITHUB_COPY_SOURCE_ROOT, GITHUB_COPY_DEST_ROOT.
Summaries and manual-review lists are printed to the terminal; no report files are written.""")
		return
	if arguments[0] not in {'scan', 'sync'}:
		raise ValueError('unknown command; use scan, sync or --help')
	prepare_roots()
	if arguments[0] == 'scan':
		scan(project_list(arguments[1:]))
	else:
		sync(arguments[1:])


try:
	main()
except Exception as error:
	# Parser exceptions are redacted by inspect; no participant values printed.
	print('ERROR: ' + str(error), file=sys.stderr)
	sys.exit(1)
PY
