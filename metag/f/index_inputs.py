#!/usr/bin/env python3
"""Stream completed genome batches or legacy merged libraries to one index."""
import argparse
import hashlib
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile

from reference_layout import DATA_DIRS
from reference_batches import INPUTS, library_inputs


UPSTREAM_SHA256 = '8d10163a6baae9dd7e20266d88e051432a430afa9fd0f9b4add8de0df1c85a62'
UPSTREAM_LIST = """function list_sequence_files() {
  find library/ '(' -name '*.fna' -o -name '*.faa' ')' -print0
}
"""


def one_library(db, library):
    if library not in DATA_DIRS:
        raise RuntimeError(f'Unknown library: {library}')
    folder = Path(db) / 'library' / library
    if (folder / INPUTS).exists():
        return library_inputs(folder)
    for name in ('library.fna', 'prelim_map.txt'):
        path = folder / name
        info = path.stat()
        if not stat.S_ISREG(info.st_mode) or not info.st_size:
            raise RuntimeError(f'Missing/empty library product: {path}')
    return [folder / 'library.fna']


def merged_inputs(db):
    """Validate the whole configured input set before emitting any paths."""
    db = Path(db)
    state = db / 'metag_build_state'
    entries = [line.partition('=')[2].split()
               for line in (state / 'config.txt').read_text().splitlines()
               if line.startswith('libraries=')]
    if len(entries) != 1 or not entries[0]:
        raise RuntimeError('Missing/ambiguous libraries in metag build configuration')
    libraries = entries[0]
    if len(set(libraries)) != len(libraries) or any(x not in DATA_DIRS for x in libraries):
        raise RuntimeError('Duplicate or unknown library in metag build configuration')
    paths = []
    total = 0
    for library in libraries:
        if not (state / f'{library}.complete').is_file():
            raise RuntimeError(f'Library is not complete: {library}')
        for path in one_library(db, library):
            total += path.stat().st_size
            paths.append(path)
    return paths, total


def index_hook():
    helper = shlex.quote(str(Path(__file__).resolve()))
    python = shlex.quote(sys.executable)
    return f"""function list_sequence_files() {{
  # metag: both passes use the same validated, ordered genome/merged inputs.
  if [[ -f metag_build_state/config.txt ]]; then
    if [[ -n "${{KRAKEN2_PROTEIN_DB:-}}" ]]; then
      echo 'metag libraries require a nucleotide database' >&2
      return 1
    fi
    {python} -B {helper} list --dir-db .
  else
    find library/ '(' -name '*.fna' -o -name '*.faa' ')' -print0
  fi
}}
"""


def build_source(source):
    """Recognize upstream and the previous managed hook, without changing either."""
    pattern = r'^function list_sequence_files\(\) \{\n.*?^\}\n'
    matches = list(re.finditer(pattern, source, re.MULTILINE | re.DOTALL))
    if len(matches) != 1:
        raise RuntimeError('Unrecognized Kraken2 input-list function')
    original = re.sub(pattern, lambda _: UPSTREAM_LIST, source,
                      count=1, flags=re.MULTILINE | re.DOTALL)
    if hashlib.sha256(original.encode()).hexdigest() != UPSTREAM_SHA256:
        raise RuntimeError('Unrecognized Kraken2 build script; requires the checked 2.17.1 source')
    patched = original.replace(UPSTREAM_LIST, index_hook(), 1)
    capacity_line = ('estimate=$(list_sequence_files | xargs -0 cat | estimate_capacity '
                     '-k $KRAKEN2_KMER_LEN -l $KRAKEN2_MINIMIZER_LEN '
                     '-S $KRAKEN2_SEED_TEMPLATE -p $KRAKEN2_THREAD_CT $KRAKEN2XFLAG )')
    if patched.count(capacity_line) != 1:
        raise RuntimeError('Unrecognized Kraken2 capacity command')
    helper = shlex.quote(str(Path(__file__).with_name('index_capacity.py').resolve()))
    replacement = (f'estimate=$({shlex.quote(sys.executable)} -B {helper} --db . -- estimate_capacity '
                   '-k $KRAKEN2_KMER_LEN -l $KRAKEN2_MINIMIZER_LEN '
                   '-S $KRAKEN2_SEED_TEMPLATE -p $KRAKEN2_THREAD_CT $KRAKEN2XFLAG )')
    patched = patched.replace(capacity_line,
                              'if [[ -f metag_build_state/config.txt ]]; then\n'
                              + replacement + '\nelse\n' + capacity_line + '\nfi', 1)
    subprocess.run(['bash', '-n'], input=patched, text=True, check=True)
    return patched


def private_build(arguments):
    command = shutil.which('kraken2-build')
    if not command:
        raise RuntimeError('kraken2-build is not on PATH')
    installed = Path(command).resolve().parent
    wrapper = (installed / 'kraken2-build').read_text()
    declaration = r'^my \$KRAKEN2_DIR = [^\n]+;\n'
    if len(re.findall(declaration, wrapper, re.MULTILINE)) != 1:
        raise RuntimeError('Unrecognized Kraken2 build wrapper')
    wrapper = re.sub(declaration, "my $KRAKEN2_DIR = $ENV{\"METAG_PRIVATE_KRAKEN_DIR\"};\n",
                     wrapper, count=1, flags=re.MULTILINE)
    source = build_source((installed / 'build_kraken2_db.sh').read_text())
    # All writes are private. Shared executables are read through symlinks only.
    with tempfile.TemporaryDirectory(prefix='metag-build-tools-') as directory:
        private = Path(directory)
        for path in installed.iterdir():
            if path.name not in ('kraken2-build', 'build_kraken2_db.sh'):
                (private / path.name).symlink_to(path)
        for name, contents in [('kraken2-build', wrapper), ('build_kraken2_db.sh', source)]:
            path = private / name
            path.write_text(contents)
            path.chmod(0o700)
        print('[metag] START reference/build: one global index; direct genome inputs; '
              'shared Kraken2 installation unchanged', flush=True)
        result = subprocess.run([str(private / 'kraken2-build'), *arguments],
                                env=dict(os.environ, METAG_PRIVATE_KRAKEN_DIR=str(private)))
        return result.returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    listing = sub.add_parser('list', help='Emit validated genome/merged FASTA paths, NUL separated')
    listing.add_argument('--dir-db', type=Path, required=True)
    checking = sub.add_parser('check-library')
    checking.add_argument('--dir-db', type=Path, required=True)
    checking.add_argument('--library', required=True)
    build = sub.add_parser('build', help='Build with private tools; forward Kraken2 build options')
    build.add_argument('options', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    try:
        if args.command == 'build':
            options = args.options[1:] if args.options[:1] == ['--'] else args.options
            return private_build(options)
        elif args.command == 'check-library':
            one_library(args.dir_db, args.library)
        else:
            paths, total = merged_inputs(args.dir_db)
            print(f'[metag] Index input: {len(paths)} ordered FASTA files, '
                  f'{total / 2**30:.2f} GiB; no duplicate genome/merged copies.',
                  file=sys.stderr)
            sys.stdout.buffer.write(b''.join(os.fsencode(path) + b'\0' for path in paths))
        return 0
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f'[metag] Index input error: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
