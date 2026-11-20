#!/usr/bin/env python3
"""Encrypt/decrypt selected ID columns in plain or gzip TSV files.

Default: deterministic FF1/AES encryption of seven-digit IDs to five A-Z letters.
The same key bytes and ID always give the same result, independent of file,
row order, or other IDs. No .table is created or needed in the default mode.
Encryption writes a sensitive plaintext original_id/encrypted_id .link beside
file-in: foo.tsv.gz -> foo.tsv.link; foo.psam -> foo.psam.link.
--first-part TRUE (default) preserves everything from the first hyphen onward.
Empty cells, negative integer markers such as -1, and one header line are preserved; --skip-lines 0 means no header.
Use --first-part FALSE for whole-cell processing. Output directories are created.
Requires pycryptodome or cryptography >= 37. No CSV quoting is interpreted.
Legacy mapping/text/alnum7 formats remain available only when explicitly selected.
FF1 ciphertext does not authenticate data or detect use of the wrong key.
"""
import argparse
import csv
import binascii
import hashlib
import json
import base64
import secrets
import gzip
import os
from pathlib import Path
import tempfile
from functools import lru_cache

ALPHABET = b'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'


class FF1:
    """NIST SP 800-38G FF1, fixed length, radix 10 or 36, empty tweak.

    AES is supplied by a crypto library; FF1 is implemented here and checked
    against NIST known-answer vectors. This implementation is not audited.
    """

    def __init__(self, key, radix=10, length=7):
        if radix not in (10, 36) or not 7 <= length <= 10:
            raise ValueError('unsupported FF1 domain')
        try:
            from Crypto.Cipher import AES
            self.aes = AES.new(key, AES.MODE_ECB).encrypt
        except ImportError:
            try:
                from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
            except ImportError:
                raise ImportError('install dependency: python3 -m pip install pycryptodome') from None
            self.aes = Cipher(algorithms.AES(key), modes.ECB()).encryptor().update
        self.radix, self.length = radix, length
        self.u = length // 2
        self.v = length - self.u
        self.b = ((radix ** self.v - 1).bit_length() + 7) // 8
        self.d = 4 * ((self.b + 3) // 4) + 4
        assert self.d <= 16
        param = bytes([1, 2, 1]) + radix.to_bytes(3, 'big') + bytes([10, self.u])
        param += length.to_bytes(4, 'big') + bytes(4)
        self.prf_start = self.aes(param)

    def transform(self, value, reverse=False):
        left, right = divmod(value, self.radix ** self.v)
        for i in (range(9, -1, -1) if reverse else range(10)):
            m = self.u if i % 2 == 0 else self.v
            q = bytes(15 - self.b) + bytes([i]) + (left if reverse else right).to_bytes(self.b, 'big')
            block = self.aes(bytes(a ^ b for a, b in zip(self.prf_start, q)))
            y = int.from_bytes(block[:self.d], 'big')
            if reverse:
                left, right = (right - y) % (self.radix ** m), left
            else:
                left, right = right, (left + y) % (self.radix ** m)
        return left * self.radix ** self.v + right


class LetterCipher:
    def __init__(self, secret, id_format='digits7'):
        if not secret:
            raise ValueError('key file is empty')
        self.radix = 10 if id_format == 'digits7' else 36
        self.width = 5 if self.radix == 10 else 8
        self.domain = self.radix ** 7
        key = hashlib.scrypt(secret, salt=b'column-crypt-ff1-v2/' + id_format.encode(),
                             n=16384, r=8, p=1, dklen=32)
        self.ff1 = FF1(key, self.radix)

    @lru_cache(maxsize=65536)
    def encrypt(self, value):
        if not value:
            return b''
        alphabet = ALPHABET[:self.radix]
        if len(value) != 7 or any(c not in alphabet for c in value):
            raise ValueError('expected seven digits' if self.radix == 10 else
                             'expected seven characters from 0-9A-Z')
        number = 0
        for c in value:
            number = number * self.radix + alphabet.index(c)
        number = self.ff1.transform(number)
        result = bytearray(b'A' * self.width)
        for i in range(self.width - 1, -1, -1):
            number, digit = divmod(number, 26)
            result[i] += digit
        return bytes(result)

    @lru_cache(maxsize=65536)
    def decrypt(self, value):
        if not value:
            return b''
        if len(value) != self.width or any(not 65 <= c <= 90 for c in value):
            raise ValueError(f'expected {self.width} uppercase letters')
        number = 0
        for c in value:
            number = number * 26 + c - 65
        if number >= self.domain:
            raise ValueError('letter code is outside this format')
        number = self.ff1.transform(number, reverse=True)
        result = bytearray(7)
        for i in range(6, -1, -1):
            number, digit = divmod(number, self.radix)
            result[i] = ALPHABET[digit]
        return bytes(result)


class TextCipher:
    AAD = b'column-crypt-siv-v3'
    TO_LETTERS = bytes.maketrans(b'0123456789abcdef', b'ABCDEFGHIJKLMNOP')
    TO_HEX = bytes.maketrans(b'ABCDEFGHIJKLMNOP', b'0123456789abcdef')

    def __init__(self, secret):
        if not secret:
            raise ValueError('key file is empty')
        self.key = hashlib.scrypt(secret, salt=self.AAD, n=16384, r=8, p=1, dklen=64)
        try:
            from Crypto.Cipher import AES
            self.AES = AES
            self.siv = None
        except ImportError:
            try:
                from cryptography.hazmat.primitives.ciphers.aead import AESSIV
                from cryptography.exceptions import InvalidTag
            except ImportError:
                raise ImportError('install dependency: python3 -m pip install pycryptodome') from None
            self.siv = AESSIV(self.key)
            self.invalid_tag = InvalidTag

    @lru_cache(maxsize=65536)
    def encrypt(self, value):
        if not value:
            return b''
        if self.siv is not None:
            data = self.siv.encrypt(value, [self.AAD])
        else:
            cipher = self.AES.new(self.key, self.AES.MODE_SIV)
            cipher.update(self.AAD)
            ciphertext, tag = cipher.encrypt_and_digest(value)
            data = tag + ciphertext
        return binascii.hexlify(data).translate(self.TO_LETTERS)

    @lru_cache(maxsize=65536)
    def decrypt(self, value):
        if not value:
            return b''
        if len(value) < 32 or len(value) % 2 or any(not 65 <= c <= 80 for c in value):
            raise ValueError('invalid AES-SIV letter encoding; expected an even number of A-P letters (at least 32)')
        data = binascii.unhexlify(value.translate(self.TO_HEX))
        if self.siv is not None:
            try:
                return self.siv.decrypt(data, [self.AAD])
            except self.invalid_tag:
                raise ValueError('decryption authentication failed: wrong key or modified ID') from None
        cipher = self.AES.new(self.key, self.AES.MODE_SIV)
        cipher.update(self.AAD)
        try:
            return cipher.decrypt_and_verify(data[16:], data[:16])
        except ValueError:
            raise ValueError('decryption authentication failed: wrong key or modified ID') from None


TABLE_HEADER = b'# crypt.py encrypted mapping v1; wrapped at 80 characters'


def unpack_table(data):
    if data.startswith(TABLE_HEADER + b'\n') or data.startswith(TABLE_HEADER + b'\r\n'):
        return b''.join(data.splitlines()[1:])
    return data  # Original single-line format remains supported.


def pack_table(data):
    return TABLE_HEADER + b'\n' + b'\n'.join(
        data[i:i + 80] for i in range(0, len(data), 80)) + b'\n'


class TableCipher(TextCipher):
    AAD = b'column-crypt-mapping-v1' 


def table_path(path):
    if path.suffix.lower() == '.gz':
        path = path.with_suffix('')
    return Path(str(path) + '.table')


class MappingCipher:
    def __init__(self, secret, source=None):
        self.codec = TableCipher(secret)
        self.forward, self.reverse = {}, {}
        if source is not None:
            try:
                obj = json.loads(self.codec.decrypt(unpack_table(source.read_bytes())))
                if obj['version'] != 1:
                    raise ValueError('unsupported table version')
                for code, encoded in obj['entries']:
                    value = base64.b64decode(encoded, validate=True)
                    token = code.encode('ascii')
                    if (len(token) != 5 or any(not 65 <= c <= 90 for c in token)
                            or not value or token in self.reverse or value in self.forward):
                        raise ValueError('invalid or duplicate table entry')
                    self.forward[value], self.reverse[token] = token, value
            except (ValueError, KeyError, TypeError, UnicodeError) as exc:
                raise ValueError(f'cannot read mapping table {source}: {exc}') from None

    def encrypt(self, value):
        if not value:
            return value
        if value in self.forward:
            return self.forward[value]
        if len(self.reverse) >= 26 ** 5:
            raise ValueError('mapping table exhausted (26**5 distinct IDs)')
        while True:
            n = secrets.randbelow(26 ** 5)
            code = bytearray(b'AAAAA')
            for i in range(4, -1, -1):
                n, digit = divmod(n, 26)
                code[i] += digit
            token = bytes(code)
            if token not in self.reverse:
                self.forward[value], self.reverse[token] = token, value
                return token

    def decrypt(self, value):
        if not value:
            return value
        try:
            return self.reverse[value]
        except KeyError:
            raise ValueError('ID not found in mapping table; check input and table') from None

    def save(self, target):
        obj = {'version': 1, 'entries': [
            [code.decode('ascii'), base64.b64encode(value).decode('ascii')]
            for code, value in self.reverse.items()]}
        data = self.codec.encrypt(json.dumps(obj, separators=(',', ':')).encode())
        fd, name = tempfile.mkstemp(prefix='.' + target.name + '.', dir=target.parent)
        try:
            with os.fdopen(fd, 'wb') as out:
                out.write(pack_table(data))
                out.flush()
                os.fsync(out.fileno())
            os.replace(name, target)
        finally:
            if os.path.exists(name):
                os.unlink(name)


def save_link(target, entries):
    fd, name = tempfile.mkstemp(prefix='.' + target.name + '.', dir=target.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as out:
            writer = csv.writer(out, delimiter='\t', lineterminator='\n')
            writer.writerow(['original_id', 'encrypted_id'])
            for original, encrypted in entries.items():
                writer.writerow([original.decode('utf-8'), encrypted.decode('ascii')])
            out.flush()
            os.fsync(out.fileno())
        os.replace(name, target)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  ./crypt.py encrypt --key secret.key --file-in raw0/chr1.psam --column 1,2 --file-out raw/chr1.psam
  ./crypt.py decrypt --key secret.key --file-in raw/chr1.psam --column 1,2 --file-out restored/chr1.psam

Default digits7 mode is deterministic, reversible and collision-free over seven
ASCII digits (including leading zeroes). The five-letter encoding is lossless.
Keep the key and plaintext .link files private. Neither .table nor .link is
required for decryption. All files can run independently with the same key.
--id-format mapping explicitly selects the old random-table format; its old
outputs need their original .table to decrypt and must not be mixed with digits7.
Use --skip-lines 0 for no header; --force to replace an existing output.
""")
    parser.add_argument('mode', choices=['encrypt', 'decrypt'])
    parser.add_argument('--key', required=True, type=Path)
    parser.add_argument('--file-in', required=True, type=Path)
    parser.add_argument('--file-out', required=True, type=Path)
    parser.add_argument('--column', required=True, help='comma-separated 1-based column numbers, e.g. 2,3,5')
    parser.add_argument('--first-part', type=str.upper, choices=['TRUE', 'FALSE'], default='TRUE',
                        help='transform only the prefix before the first hyphen; preserve suffixes (default: TRUE)')
    parser.add_argument('--delimiter', default='\t', help=r'default: tab; accepts \t')
    parser.add_argument('--skip-lines', type=int, default=1, help='initial lines to preserve (default: 1; use 0 for no header)')
    parser.add_argument('--id-format', choices=['mapping', 'text', 'digits7', 'alnum7'], default='digits7',
                        help='digits7: deterministic five letters without tables (default); other modes: legacy')
    parser.add_argument('--table', type=Path, help='existing encrypted mapping table to reuse; automatic by default')
    parser.add_argument('--table-out', type=Path, help='encrypted mapping output path (default: beside file-in, remove .gz and append .table)')
    parser.add_argument('--linked-id-columns', type=int, nargs='+', default=[],
                        help='columns containing the selected ID followed by a hyphen; transform that prefix too')
    parser.add_argument('--force', action='store_true', help='replace an existing output')
    args = parser.parse_args()
    temp = None
    lock = None
    try:
        first_part = args.first_part == 'TRUE'
        if first_part and args.linked_id_columns:
            raise ValueError('--first-part and --linked-id-columns cannot be combined')
        columns = [int(c) for c in args.column.split(',')]
        if len(set(columns)) != len(columns):
            raise ValueError('duplicate column numbers')
        args.column = columns[0]
        if min(columns) < 1 or args.skip_lines < 0:
            raise ValueError('column must be positive; skip-lines must be nonnegative')
        sep = ('\t' if args.delimiter == r'\t' else args.delimiter).encode()
        if not sep or b'\n' in sep or b'\r' in sep:
            raise ValueError('delimiter must be nonempty and cannot contain newlines')
        if any(c in sep for c in b'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-=:'):
            raise ValueError('delimiter conflicts with encrypted token alphabet')
        src, dst = args.file_in, args.file_out
        if dst.resolve() in (src.resolve(), args.key.resolve()):
            raise ValueError('output must differ from input and key')
        if dst.exists() and not args.force:
            raise ValueError('output exists; use --force to replace it')
        if args.table_out and (args.id_format != 'mapping' or args.mode != 'encrypt'):
            raise ValueError('--table-out requires mapping encryption')
        if args.linked_id_columns and len(columns) != 1:
            raise ValueError('linked-id-columns requires one primary column')
        if any(c < 1 or c == args.column for c in args.linked_id_columns):
            raise ValueError('linked ID columns must be positive and different from column')
        link_output = table_path(src).with_suffix('.link') if args.mode == 'encrypt' else None
        if link_output is not None and link_output.resolve() in (src.resolve(), dst.resolve(), args.key.resolve()):
            raise ValueError('link path conflicts with input, output or key')
        link_entries = {}
        secret = args.key.read_bytes()
        dst.parent.mkdir(parents=True, exist_ok=True)
        mapping_output = None
        if args.id_format == 'mapping':
            mapping_output = (args.table_out or table_path(src)) if args.mode == 'encrypt' else None
            source = args.table or (mapping_output if args.mode == 'encrypt' else table_path(src))
            if mapping_output is not None:
                if mapping_output.resolve() in (src.resolve(), dst.resolve(), args.key.resolve()):
                    raise ValueError('mapping path conflicts with input, output or key')
                if link_output is not None and mapping_output.resolve() == link_output.resolve():
                    raise ValueError('table and link paths must differ')
                mapping_output.parent.mkdir(parents=True, exist_ok=True)
                candidate = Path(str(mapping_output) + '.lock')
                fd = os.open(candidate, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
                os.close(fd)
                lock = candidate
                if mapping_output.exists() and source.resolve() != mapping_output.resolve():
                    raise ValueError('output table already exists; reuse it by omitting --table')
            if not source.exists() and (args.table or args.mode == 'decrypt'):
                raise ValueError(f'mapping table does not exist: {source}')
            if source.exists() and source.resolve() in (src.resolve(), dst.resolve(), args.key.resolve()):
                raise ValueError('mapping path conflicts with input, output or key')
            cipher = MappingCipher(secret, source if source.exists() else None)
        else:
            if args.table:
                raise ValueError('--table requires --id-format mapping')
            cipher = TextCipher(secret) if args.id_format == 'text' else LetterCipher(secret, args.id_format)
        fd, temp = tempfile.mkstemp(prefix='.' + dst.name + '.', dir=dst.parent)
        os.close(fd)
        opener_in = gzip.open if src.suffix.lower() == '.gz' else open
        opener_out = gzip.open if dst.suffix.lower() == '.gz' else open
        with opener_in(src, 'rb') as fin, opener_out(temp, 'wb') as fout:
            for number, line in enumerate(fin, 1):
                if number <= args.skip_lines or line in (b'\n', b'\r\n'):
                    fout.write(line)
                    continue
                end = b'\r\n' if line.endswith(b'\r\n') else (b'\n' if line.endswith(b'\n') else b'')
                fields = (line[:-len(end)] if end else line).split(sep)
                for column_index in columns:
                    if len(fields) < column_index:
                        raise ValueError(f'line {number}: missing column {column_index}')
                    original = value = fields[column_index - 1]
                    # Negative integer sentinels (e.g. -1) are not UKB IDs.
                    # Preserve the complete field in both directions; do not
                    # split its leading minus or include it in the ID link.
                    if (args.id_format == 'digits7' and original.startswith(b'-')
                            and original[1:].isdigit()):
                        continue
                    suffix = b''
                    if first_part:
                        value, separator, tail = original.partition(b'-')
                        suffix = separator + tail
                        if original and not value:
                            raise ValueError(f'line {number}: column {column_index} has empty ID prefix')
                    mapping_original = value
                    try:
                        if link_output is not None and value in link_entries:
                            value = link_entries[value]
                        else:
                            value = getattr(cipher, args.mode)(value)
                    except ValueError as exc:
                        raise ValueError(f'line {number}: {exc}') from None
                    for column in args.linked_id_columns:
                        if column > len(fields):
                            raise ValueError(f'line {number}: missing linked column {column}')
                        linked = fields[column - 1]
                        if not linked:
                            continue
                        if not original or not linked.startswith(original + b'-'):
                            raise ValueError(f'line {number}: linked column {column} does not start with selected ID and hyphen')
                        fields[column - 1] = value + linked[len(original):]
                    if link_output is not None and mapping_original:
                        link_entries[mapping_original] = value
                    fields[column_index - 1] = value + suffix
                fout.write(sep.join(fields) + end)
        # Publish the table first: a failed output publish can leave extra entries,
        # but never deletes old entries or leaves published output without its table.
        if mapping_output is not None:
            cipher.save(mapping_output)
        if link_output is not None:
            save_link(link_output, link_entries)
        if args.force:
            os.replace(temp, dst)
        else:
            os.link(temp, dst)
            os.unlink(temp)
        temp = None
    except (OSError, EOFError, ValueError, ImportError) as exc:
        parser.exit(1, f'Error: {exc}\n')
    finally:
        if temp is not None:
            os.unlink(temp)
        if lock is not None:
            lock.unlink()


if __name__ == '__main__':
    main()
