#!/usr/bin/env python3
"""按已有的权威 MD5 检查下载，并局部修复混入二进制文件的连接错误文本。

默认下载缺失文件、检查 MD5 并尝试局部修补；--check-only 只检查。
*.download_check.jsonl 合并结果与修补字节备份；*.download_check.fail.tsv
保存尚未解决的失败及原因。已有文件修补失败时，不自动重下整个文件。
需要 Python 3.9+，仅使用标准库。相对路径均相对于 --dir。
"""
import argparse
import base64
import csv
import fcntl
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
import time
from contextlib import ExitStack
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote, urlsplit
from urllib.error import HTTPError
from urllib.request import ProxyHandler, Request, build_opener


# 🚩 输入与记录
ERROR_TEXT = re.compile(rb'<Error><Code>ConnectionClosedException</Code>.{0,4096}?</Error>', re.DOTALL)
CHUNK_SIZE = 8 * 1024 * 1024


def resolve(base, value):
	return (base / value).resolve()


def read_jobs(base, url_file, md5_file, folder):
	checksums = {}
	for number, line in enumerate(md5_file.read_text(encoding='utf-8-sig').splitlines(), 1):
		if not line.strip() or line.startswith('#'):
			continue
		match = re.fullmatch(r'([0-9a-fA-F]{32}) [ *](.+)', line)
		if not match:
			raise ValueError(f'{md5_file.name}:{number}: 不是标准 md5sum 清单')
		digest, name = match.groups()
		# 清单可用“文件夹/文件名”或单独文件名；不跨文件夹猜测匹配。
		path = resolve(folder if len(Path(name).parts) == 1 else base, name)
		if path in checksums and checksums[path] != digest.lower():
			raise ValueError(f'MD5 清单有冲突: {name}')
		checksums[path] = digest.lower()
	entries = []
	for number, line in enumerate(url_file.read_text(encoding='utf-8-sig').splitlines(), 1):
		if not line.strip() or line.lstrip().startswith('#'):
			continue
		if line[0].isspace():
			if not entries or '=' not in line:
				raise ValueError(f'{url_file.name}:{number}: 无法识别 aria2 选项')
			key, value = line.strip().split('=', 1)
			if key not in {'out', 'checksum'}:
				raise ValueError(f'{url_file.name}:{number}: 不支持选项 {key}；请使用普通 URL 清单')
			entries[-1][key] = value
			continue
		urls = line.split()  # aria2 同一行的镜像地址仅使用第一个。
		if any(urlsplit(url).scheme not in {'http', 'https', 'ftp'} for url in urls):
			raise ValueError(f'{url_file.name}:{number}: 不是 HTTP(S)/FTP URL')
		entries.append({'url': urls[0]})
	metadata = {}
	metadata_file = base / 'ena-md5.metadata.tsv'
	if metadata_file.exists():
		with metadata_file.open(encoding='utf-8', newline='') as handle:
			for row in csv.DictReader(handle, delimiter='\t'):
				metadata[row['url']] = row
	jobs, seen = [], set()
	for entry in entries:
		name = entry.get('out', unquote(urlsplit(entry['url']).path.rsplit('/', 1)[-1]))
		if not name or Path(name).is_absolute() or '..' in Path(name).parts:
			raise ValueError(f'不安全或缺失的文件名: {name!r}')
		path = resolve(folder, name)
		if folder not in path.parents or path in seen:
			raise ValueError(f'重复文件或路径超出 --folder: {name}')
		seen.add(path)
		if path not in checksums:
			raise ValueError(f'MD5 清单没有该文件: {path}')
		digest = checksums[path]
		if 'checksum' in entry and entry['checksum'].lower() != f'md5={digest}':
			raise ValueError(f'URL 清单和 MD5 清单的校验值不一致: {name}')
		row = metadata.get(entry['url'])
		if row and row['md5'].lower() != digest:
			raise ValueError(f'ENA 元数据和 MD5 清单不一致: {name}')
		jobs.append({'path': path, 'url': entry['url'], 'md5': digest,
			'bytes': int(row['bytes']) if row else None})
	if not jobs:
		raise ValueError('URL 清单为空')
	return jobs


def record(journal, event, **fields):
	journal.write(json.dumps({'time': datetime.now().astimezone().isoformat(),
		'event': event, **fields}, ensure_ascii=False) + '\n')
	journal.flush()
	os.fsync(journal.fileno())


def read_failures(path):
	if not path.exists():
		return {}
	with path.open(encoding='utf-8', newline='') as handle:
		reader = csv.DictReader(handle, delimiter='\t')
		if reader.fieldnames != ['file', 'url', 'reason']:
			raise ValueError(f'失败清单格式不符，拒绝覆盖: {path}')
		return {row['file']: row for row in reader}


def save_failures(path, failures):
	# 在同一目录原子替换；中断不会先清空上一份失败清单。
	temporary = None
	try:
		with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='',
				dir=path.parent, prefix='.' + path.name + '.', suffix='.tmp', delete=False) as handle:
			temporary = Path(handle.name)
			writer = csv.DictWriter(handle, fieldnames=['file', 'url', 'reason'], delimiter='\t', lineterminator='\n')
			writer.writeheader()
			writer.writerows(failures.values())
			handle.flush()
			os.fsync(handle.fileno())
		os.replace(temporary, path)
	finally:
		if temporary is not None:
			temporary.unlink(missing_ok=True)


# 🚩 缺失文件下载：只在校验成功后发布正式文件
def download_missing(job, opener, timeout):
	path = job['path']
	part = path.with_name(path.name + '.download_check.part')
	if part.is_symlink():
		raise RuntimeError('下载临时文件是软链接，拒绝写入')
	part.parent.mkdir(parents=True, exist_ok=True)
	ensure_idle(part)
	print('  文件缺失，下载至临时文件；MD5 通过后再保存为正式文件…', flush=True)
	for attempt in range(3):
		offset = part.stat().st_size if part.exists() else 0
		if part.exists() and job['bytes'] is not None and offset >= job['bytes']:
			return part  # 完整或超长临时文件交给 MD5 流程判断，不能凭长度当成成功。
		is_http = urlsplit(job['url']).scheme in {'http', 'https'}
		headers = {'Accept-Encoding': 'identity', 'Cache-Control': 'no-cache'}
		if offset and is_http:
			headers['Range'] = f'bytes={offset}-'
		try:
			with opener.open(Request(job['url'], headers=headers), timeout=timeout) as response:
				length = response.headers.get('Content-Length')
				length = int(length) if length is not None else None
				if is_http and response.status == 206:
					match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
					if not offset or not match or int(match[1]) != offset or int(match[2]) != int(match[3]) - 1:
						raise RuntimeError('续传响应的字节范围不正确')
					total = int(match[3])
					if length is not None and length != total - offset:
						raise RuntimeError('续传响应的长度不正确')
					mode = 'ab'
				else:
					if is_http and response.status != 200:
						raise RuntimeError(f'下载 HTTP 状态异常: {response.status}')
					# 服务器不支持续传或使用 FTP 时，从头写本程序的临时文件。
					offset, total, mode = 0, length, 'wb'
				if job['bytes'] is not None and total is not None and total != job['bytes']:
					raise RuntimeError('服务器文件长度与已有 ENA 元数据不一致')
				written, last_progress = offset, time.monotonic()
				with part.open(mode) as handle:
					while data := response.read(CHUNK_SIZE):
						handle.write(data)
						written += len(data)
						if time.monotonic() - last_progress >= 30:
							print(f'  已下载 {written / 2**20:.1f} MiB', flush=True)
							last_progress = time.monotonic()
					handle.flush()
					os.fsync(handle.fileno())
				if (total is not None and written != total) or (job['bytes'] is not None and written != job['bytes']):
					raise RuntimeError('下载长度不足；临时文件保留供续传')
			return part
		except Exception as error:
			# 无元数据时，完整临时文件的续传可能返回 416，仍须做完整 MD5。
			if isinstance(error, HTTPError) and error.code == 416 and offset:
				match = re.fullmatch(r'bytes \*/(\d+)', error.headers.get('Content-Range', ''))
				if match and int(match[1]) == offset and job['bytes'] in {None, offset}:
					return part
			if attempt == 2:
				raise RuntimeError(f'下载未完成: {error}；临时文件供下次继续') from error
			print(f'  下载重试 {attempt + 1}/2: {error}', flush=True)
			time.sleep(2 * (attempt + 1))


# 🚩 完整校验与分段修补
def file_state(path):
	info = path.stat()
	if not stat.S_ISREG(info.st_mode):
		raise ValueError('不是普通文件')
	return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def ensure_idle(path):
	if path.with_name(path.name + '.aria2').exists():
		raise RuntimeError('存在 .aria2 控制文件；跳过，避免干扰未完成的下载')
	# 检查本机可见的写入者，兼顾没有 .aria2 标记的 wget。
	for process in Path('/proc').glob('[0-9]*'):
		if process.name == str(os.getpid()):
			continue
		try:
			for descriptor in (process / 'fd').iterdir():
				try:
					if os.readlink(descriptor) != str(path):
						continue
					info = (process / 'fdinfo' / descriptor.name).read_text()
					flags = re.search(r'^flags:\s+([0-7]+)', info, re.M)
					if flags and int(flags[1], 8) & os.O_ACCMODE:
						raise RuntimeError(f'PID {process.name} 正在写入此文件，跳过')
				except (OSError, ValueError):
					continue
		except OSError:
			continue


def scan(path):
	before = file_state(path)
	digest, errors, seen = hashlib.md5(), [], set()
	offset, tail = 0, b''
	with path.open('rb') as handle:
		while data := handle.read(CHUNK_SIZE):
			digest.update(data)
			buffer = tail + data
			for match in ERROR_TEXT.finditer(buffer):
				position = offset - len(tail) + match.start()
				if position not in seen:
					seen.add(position)
					errors.append((position, match[0]))
			tail = buffer[-8192:]
			offset += len(data)
	if file_state(path) != before:
		raise RuntimeError('扫描期间文件发生变化，停止处理')
	return digest.hexdigest(), sorted(errors), before


def fetch_range(opener, url, start, end, total, timeout):
	if urlsplit(url).scheme not in {'http', 'https'}:
		raise RuntimeError('局部修复需要 HTTP(S) URL；FTP 文件只做校验')
	for attempt in range(3):
		try:
			request = Request(url, headers={'Range': f'bytes={start}-{end}',
				'Accept-Encoding': 'identity', 'Cache-Control': 'no-cache'})
			with opener.open(request, timeout=timeout) as response:
				if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{total}':
					raise RuntimeError('服务器返回的字节范围或文件总长度不符合预期')
				data = response.read(end - start + 2)
			if len(data) != end - start + 1 or b'<Error>' in data:
				raise RuntimeError('远程分段不完整或含错误文本')
			return data
		except Exception:
			if attempt == 2:
				raise
			time.sleep(2 * (attempt + 1))


def repair(job, errors, before, opener, timeout, journal):
	path, total = job['path'], before[2]
	# 合并相邻错误的请求范围，使邻近的错误文本不会被误认为正常的保护字节。
	groups = []
	for offset, old in errors:
		start, end = max(0, offset - 64), min(total - 1, offset + len(old) + 63)
		if groups and start <= groups[-1][1] + 1:
			groups[-1][1] = end
			groups[-1][2].append((offset, old))
		else:
			groups.append([start, end, [(offset, old)]])
	patches = []
	for start, end, group in groups:
		remote = fetch_range(opener, job['url'], start, end, total, timeout)
		with path.open('rb') as handle:
			handle.seek(start)
			local = handle.read(end - start + 1)
		candidate = bytearray(local)
		for offset, old in group:
			left, right = offset - start, offset - start + len(old)
			if local[left:right] != old:
				raise RuntimeError('原始坏字节发生变化，停止修补')
			new = remote[left:right]
			candidate[left:right] = new
			patches.append((offset, old, new))
		if bytes(candidate) != remote:
			raise RuntimeError('错误文本之外的数据也不一致，不能安全局部修复')
	ensure_idle(path)
	if file_state(path) != before:
		raise RuntimeError('获取分段期间文件发生变化，停止修补')
	record(journal, 'before_patch', file=str(path), target_file=str(job.get('target', path)), expected_md5=job['md5'],
		patches=[{'offset': offset, 'old_base64': base64.b64encode(old).decode(),
			'new_base64': base64.b64encode(new).decode()} for offset, old, new in patches])
	with path.open('r+b') as handle:
		for offset, old, new in patches:
			handle.seek(offset)
			if handle.read(len(old)) != old:
				raise RuntimeError('写入前原始字节发生变化；已备份，停止处理')
			handle.seek(offset)
			handle.write(new)
		handle.flush()
		os.fsync(handle.fileno())
	return sum(len(old) for _, old, _ in patches)


# 🚩 命令入口
def main():
	parser = argparse.ArgumentParser(
		description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog=(
			'示例（在 /mnt/d/scripts 中运行）：\n'
			'  # 下载缺失文件，检查并尝试修复\n'
			'  ./download_check.py --dir /mnt/f/metag --url 全球污水.url --folder 全球污水\n'
			'  # 只检查，不修补\n'
			'  ./download_check.py --dir /mnt/f/metag --url 长江.url --folder 长江 --check-only\n\n'
			'  # 只检查并尝试修复指定文件\n'
			'  ./download_check.py --dir /mnt/f/metag --url 长江.url --folder 长江 --only SRR24791598_2.fastq.gz'
		))
	parser.add_argument('--dir', required=True, help='数据所在父目录，例如 /mnt/f/metag')
	parser.add_argument('--url', required=True, help='URL 清单，例如 长江.url')
	parser.add_argument('--folder', required=True, help='下载文件夹，例如 长江')
	parser.add_argument('--md5', help='权威 MD5 清单，默认与 URL 文件同名、后缀 .md5')
	parser.add_argument('--check-only', action='store_true', help='仅校验，不下载、不联网、不修补')
	parser.add_argument('--only', action='append', help='只处理指定文件名，可重复')
	parser.add_argument('--use-proxy', action='store_true', help='使用环境变量中的代理；默认直连')
	parser.add_argument('--timeout', type=float, default=30, help='网络请求的超时秒数（默认 30）')
	parser.add_argument('--plan', action='store_true', help='仅验证清单和路径，显示任务规模，不读取数据、不生成记录')
	args = parser.parse_args()
	base = Path(args.dir).resolve()
	url_file, folder = resolve(base, args.url), resolve(base, args.folder)
	md5_file = resolve(base, args.md5) if args.md5 else url_file.with_suffix('.md5')
	if not base.is_dir() or (folder.exists() and not folder.is_dir()) or args.timeout <= 0:
		parser.error('--dir 必须是已有目录，--folder 不能是文件，--timeout 必须大于 0')
	jobs = read_jobs(base, url_file, md5_file, folder)
	if args.only:
		selected = set(args.only)
		jobs = [job for job in jobs if job['path'].name in selected]
		if {job['path'].name for job in jobs} != selected:
			parser.error('--only 包含不在 URL 清单中的文件')
	journal_path = url_file.with_name(url_file.stem + '.download_check.jsonl')
	failure_path = url_file.with_name(url_file.stem + '.download_check.fail.tsv')
	print(f'文件数: {len(jobs)}；目录: {folder}\nMD5: {md5_file}\n记录与字节备份: {journal_path}\n失败清单: {failure_path}', flush=True)
	if args.plan:
		return 0
	if not folder.exists():
		if args.check_only:
			parser.error('--check-only 指定的文件夹不存在')
		folder.mkdir(parents=True)
	with ExitStack() as stack:
		# 锁目录本身，不生成 .lock 文件；同时尊重仍在运行的旧修复程序。
		directory_fd = os.open(folder, os.O_RDONLY | os.O_DIRECTORY)
		stack.callback(os.close, directory_fd)
		try:
			fcntl.flock(directory_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
			for legacy in {url_file.with_suffix('.repair.lock'), base / (folder.name + '.repair.lock')}:
				if legacy.exists():
					lock = stack.enter_context(legacy.open('r'))
					fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
		except BlockingIOError:
			raise RuntimeError('此文件夹已有检查或旧 repair 任务运行；请等它完成后再运行')
		journal = stack.enter_context(journal_path.open('a', encoding='utf-8'))
		fcntl.flock(journal.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
		failures = read_failures(failure_path)
		if not failure_path.exists():
			save_failures(failure_path, failures)
		opener = build_opener() if args.use_proxy else build_opener(ProxyHandler({}))
		passed, fixed, downloaded, failed, patched_bytes = 0, 0, 0, [], 0
		record(journal, 'start', folder=str(folder), url_file=str(url_file), md5_file=str(md5_file),
			files=len(jobs), check_only=args.check_only)
		for number, job in enumerate(jobs, 1):
			path = job['path']
			print(f'[{number}/{len(jobs)}] {path.name}', flush=True)
			actual, changed, staged, reason = None, 0, False, None
			try:
				ensure_idle(path)
				work = job
				if not path.exists():
					if args.check_only:
						raise RuntimeError('文件不存在（仅检查，不下载）')
					part = download_missing(job, opener, args.timeout)
					work = {**job, 'path': part, 'target': path}
					staged = True
				actual, errors, before = scan(work['path'])
				if actual != job['md5']:
					if args.check_only:
						raise RuntimeError(f'MD5 不符，找到 {len(errors)} 处连接错误文本（仅检查）')
					if job['bytes'] is not None and before[2] != job['bytes']:
						raise RuntimeError(f'大小不符: {before[2]} / {job["bytes"]} 字节；需要续传或重下')
					if not errors:
						raise RuntimeError('MD5 不符，未发现可定位的错误文本；需要进一步检查或重下')
					print(f'  发现 {len(errors)} 处错误文本，正在取回对应字节…', flush=True)
					changed = repair(work, errors, before, opener, args.timeout, journal)
					patched_bytes += changed
					actual, _, before = scan(work['path'])
					if actual != job['md5']:
						raise RuntimeError('局部修补后 MD5 仍不符；另有损坏，不能认定已修好；原始字节已备份')
				if staged:
					ensure_idle(path)
					ensure_idle(work['path'])
					if file_state(work['path']) != before:
						raise RuntimeError('下载临时文件在校验后发生变化，停止保存')
					# 同目录硬链接提供原子、不覆盖的发布；如果目标出现，保留临时文件并报错。
					os.link(work['path'], path)
					work['path'].unlink()
					downloaded += 1
				passed += 1
				fixed += bool(changed)
				print('  MD5 OK' + (f'，修补 {changed} 字节' if changed else ''), flush=True)
			except Exception as error:
				reason = str(error)
				failed.append(str(path))
				print(f'  FAILED: {error}', flush=True)
			record(journal, 'check', file=str(path), url=job['url'], status='FAILED' if reason is not None else 'OK',
				expected_md5=job['md5'], actual_md5=actual, reason=reason, replaced_bytes=changed, downloaded=staged and reason is None)
			if reason is None:
				failures.pop(str(path), None)
			else:
				failures[str(path)] = {'file': str(path), 'url': job['url'], 'reason': reason}
			save_failures(failure_path, failures)
		record(journal, 'summary', total=len(jobs), passed=passed, fixed=fixed,
			downloaded=downloaded, failed=failed, replaced_bytes=patched_bytes)
		print(f'完成：MD5 通过 {passed}/{len(jobs)}，本次下载 {downloaded}，修复 {fixed}，仍需处理 {len(failed)}，修补 {patched_bytes} 字节。', flush=True)
		print(f'失败清单现有 {len(failures)} 条: {failure_path}', flush=True)
		return 1 if failed else 0


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print('\n已中断；已有记录和字节备份保留，原命令可重新运行。', file=sys.stderr)
		sys.exit(130)
	except (OSError, ValueError, RuntimeError) as error:
		print(f'错误: {error}', file=sys.stderr)
		sys.exit(2)
