"""Build metag_classify without changing the pinned upstream executables.

Public API: ``install(cfg) -> receipt``. The executable is
``Path(cfg['tools_dir']) / 'metag_union/bin/metag_classify'``.
No network access, package-manager install, or production input is required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hpc_common import CODE_ROOT, atomic_json, locked, sha256, software_sources, utcnow

MODULE_ROOT = Path(__file__).resolve().parent
SCHEMA = 1

LOAD_INDEX = r'''IndexData *load_index(Options &opts) {
  cerr << "Loading compatible category indexes (CPU mmap)...";
  if (opts.daemon_mode)
    errx(EX_USAGE, "metag_classify does not support daemon mode");
  opts.use_memory_mapping = true;
  const auto paths = metag_union::ValidateInputs(opts.taxonomy_filename, opts.options_filename);
  IndexOptions idx_opts = {};
  ifstream options(opts.options_filename, std::ios::binary);
  options.read(reinterpret_cast<char *>(&idx_opts), sizeof(idx_opts));
  if (!options) errx(EX_DATAERR, "metag union: unable to read index options");
  if (!idx_opts.dna_db)
    errx(EX_USAGE, "metag_classify accepts DNA category indexes only");
  opts.use_translated_search = false;
  IndexData *index_data = new IndexData();
  index_data->options = idx_opts;
  index_data->taxonomy = new Taxonomy(opts.taxonomy_filename, true);
  index_data->cht = new metag_union::UnionHashTable(paths, *index_data->taxonomy);
  cerr << " " << paths.size() << " categories; one FASTQ scan; global minimizer LCA." << endl;
  return index_data;
}

'''

TARGET = r'''

# Metag extension: separate binary; upstream classify is left intact.
metag_classify.o: metag_classify.cc metag_union.h kraken2_data.h kv_store.h taxonomy.h seqreader.h mmscanner.h compact_hash.h aa_translate.h reports.h utilities.h readcounts.h

metag_classify: metag_classify.o reports.o hyperloglogplus.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o aa_translate.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ metag_classify.o reports.o hyperloglogplus.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o aa_translate.o utilities.o $(LDFLAGS)
'''


def _replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"Pinned Kraken source does not match extension anchor: {old[:70]!r}")
    return text.replace(old, new, 1)


def patched_classifier(original: str) -> str:
    """Change index loading only; retain original ClassifySequence verbatim."""
    result = _replace_once(original, '#include "compact_hash.h"',
                           '#include "compact_hash.h"\n#include "metag_union.h"')
    start = result.index('IndexData *load_index(Options &opts) {')
    end = result.index('void classify(Options &opts, IndexData *index_data) {', start)
    result = result[:start] + LOAD_INDEX + result[end:]
    result = _replace_once(result, 'int main(int argc, char **argv) {',
                           'int main(int argc, char **argv) {\n'
                           '  metag_union::ExtractArguments(argc, argv);')
    result = _replace_once(result, '  if (opts.daemon_mode) {\n    ClassifyDaemon(opts);',
                           '  if (opts.daemon_mode) {\n'
                           '    errx(EX_USAGE, "metag_classify does not support daemon mode");\n'
                           '    ClassifyDaemon(opts);')
    result = _replace_once(result, '  ParseCommandLine(argc, argv, opts);\n  if (opts.daemon_mode)',
                           '  ParseCommandLine(argc, argv, opts);\n'
                           '  if (opts.quick_mode) errx(EX_USAGE, "metag_classify requires full evidence; -q is unsupported");\n'
                           '  if (opts.daemon_mode)')
    result = _replace_once(result, 'Usage: classify [options] <fasta/fastq file(s)>',
                           'Usage: metag_classify --metag-hashes LISTFILE [options] <fasta/fastq file(s)>')
    result = _replace_once(result, '  -q               Quick mode',
                           '  -q               Unsupported: metag union requires complete evidence')
    # This anchor is also a fail-closed upstream version check.
    marker = 'taxid_t ClassifySequence(Sequence &dna, Sequence &dna2, ostringstream &koss,'
    definition = original.rindex(marker)
    original_classification = original[definition:original.index('void AddHitlistString(', definition)]
    if original_classification not in result:
        raise RuntimeError('Extension unexpectedly changed upstream sequence classification')
    return result


def _extract(archive: Path, destination: Path) -> None:
    """Extract the bundled source, rejecting links and path traversal."""
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            p = Path(member.name)
            if p.is_absolute() or ".." in p.parts or not p.parts or p.parts[0] != "kraken":
                raise RuntimeError(f"Unsafe bundled archive member: {member.name}")
            if not (member.isfile() or member.isdir()):
                raise RuntimeError("Bundled source must contain only files and directories")
            target = destination.joinpath(*p.parts[1:])
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with bundle.extractfile(member) as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output)
                target.chmod(member.mode & 0o777)


def install(cfg: dict) -> dict:
    for program in ("g++", "make"):
        if not shutil.which(program):
            raise RuntimeError(f"Missing {program}; load the site's compiler module")
    source = software_sources(CODE_ROOT)["kraken"]
    archive = CODE_ROOT / "f" / "vendor" / source["archive"]
    if sha256(archive) != source["sha256"]:
        raise RuntimeError("Bundled Kraken source checksum mismatch")
    extension = hashlib.sha256()
    for path in (Path(__file__).resolve(), MODULE_ROOT / "metag_union.h"):
        extension.update(path.name.encode() + b"\0" + path.read_bytes())
    owner_data = {"schema": SCHEMA, "source_archive_sha256": source["sha256"],
                  "extension_sha256": extension.hexdigest()}
    prefix = Path(cfg["tools_dir"]).expanduser().absolute() / "metag_union"
    if prefix.is_symlink():
        raise RuntimeError(f"Refusing symlink at native extension installation: {prefix}")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with locked(prefix.parent / "metag_union.install.lock"):
        owner = prefix / ".metag-union-owner.json"
        if prefix.exists() and not owner.is_file():
            raise RuntimeError(f"Unowned native extension directory: {prefix}")
        prefix.mkdir(exist_ok=True)
        if owner.exists():
            if json.loads(owner.read_text()) != owner_data:
                raise RuntimeError("Different native extension owns this tools directory; choose a new tools_dir")
        else:
            atomic_json(owner, owner_data)
        receipt = prefix / "install.json"
        executable = prefix / "bin/metag_classify"
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if (saved.get("identity") != owner_data or
                    saved.get("files", {}).get("metag_classify") != sha256(executable) or
                    not os.access(executable, os.X_OK)):
                raise RuntimeError("Published metag_classify installation changed")
            print(f"REUSE verified CPU union classifier: {executable}", flush=True)
            return saved
        tree = prefix / "source"
        if not tree.exists():
            stage = prefix / "source.partial"
            if stage.exists():
                shutil.rmtree(stage)
            stage.mkdir()
            _extract(archive, stage)
            src = stage / "src"
            original = (src / "classify.cc").read_text()
            (src / "metag_classify.cc").write_text(patched_classifier(original))
            shutil.copyfile(MODULE_ROOT / "metag_union.h", src / "metag_union.h")
            with (src / "Makefile").open("a") as makefile:
                makefile.write(TARGET)
            stage.replace(tree)
        src = tree / "src"
        # An interrupted build may resume, but source changes never pass silently.
        if ((src / "metag_classify.cc").read_text() != patched_classifier((src / "classify.cc").read_text()) or
                sha256(src / "metag_union.h") != sha256(MODULE_ROOT / "metag_union.h")):
            raise RuntimeError("Unpublished native extension source changed")
        environment = dict(os.environ)
        environment.pop("KRAKEN2_SKIP_FOPENMP", None)
        environment.pop("MAKEFLAGS", None)
        jobs = max(1, int(cfg.get("setup_cores", 8)))
        log = prefix / "install.log"
        print(f"Compile CPU minimizer union classifier; log={log}", flush=True)
        with log.open("ab") as handle:
            completed = subprocess.run(["make", f"-j{jobs}", "metag_classify"], cwd=src,
                                       env=environment, stdout=handle, stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f"metag_classify compiler exited {completed.returncode}; inspect {log}")
        executable.parent.mkdir(exist_ok=True)
        temporary = executable.with_suffix(".partial")
        shutil.copyfile(src / "metag_classify", temporary)
        temporary.chmod(0o755)
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, executable)
        subprocess.run([str(executable), "-h"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        result = {"schema": SCHEMA, "identity": owner_data,
                  "source_commit": source["source_commit"],
                  "installed_at": utcnow(), "executable": str(executable),
                  "files": {"metag_classify": sha256(executable)},
                  "mode": "CPU mmap; per-minimizer union LCA; one FASTQ scan"}
        atomic_json(receipt, result)
        print(f"PASS CPU union classifier: {executable}", flush=True)
        return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tools-dir", required=True, type=Path)
    parser.add_argument("--jobs", default=8, type=int)
    args = parser.parse_args()
    print(json.dumps(install({"tools_dir": str(args.tools_dir), "setup_cores": args.jobs}), indent=2))


if __name__ == "__main__":
    main()
