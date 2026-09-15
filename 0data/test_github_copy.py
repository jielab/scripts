"""Isolated publishing safety regression test; uses synthetic data only."""
import os
from pathlib import Path
import subprocess
import tempfile
import zipfile


SCRIPT = Path(__file__).resolve().with_name("github_copy.sh")


def main():
    with tempfile.TemporaryDirectory(prefix="github-copy-test-") as directory:
        root = Path(directory)
        source, dest = root / "source", root / "dest"
        relative = Path("le8/cvd_cad/prot/c2_cause")
        folder = source / relative
        folder.mkdir(parents=True)
        fixtures = {
            "summary.csv": "protein,beta,p\nABC,0.2,0.01\n",
            "result.RDS": "blocked",
            "eid.csv": "eid,value\n42,0.1\n",
            "iid.tsv": "IID\tvalue\nsynthetic-subject\t0.1\n",
            "header_only.csv": "EiD,value\n",
            "late.csv": "note,value\n" + "summary,ok\n" * 300 + "IID,value\nx,1\n",
            "wide.csv": ",".join(["field"] * 300 + ["eid"]) + "\n",
            "large.csv": "note,value\n" + ("x" * 1024 + ",ok\n") * 4100 + "IID,value\n",
            "generic.csv": "id,value\n1000001,a\n1000002,b\n1000003,c\n",
            "quoted.csv": 'protein,note\nABC,"line one\nline two"\n',
            "broken.xlsx": "invalid workbook",
        }
        for name, value in fixtures.items():
            (folder / name).write_text(value, encoding="utf-8")
        (folder / "utf16.csv").write_text("IID,value\n42,1\n", encoding="utf-16")
        (folder / "binary.csv").write_bytes(b"\x00\x01\xff")
        for subdir in ("mrlink2", "cache", "too_deep"):
            (folder / subdir).mkdir()
            (folder / subdir / "excluded.csv").write_text("a,b\n1,2\n")
        with zipfile.ZipFile(folder / "late_sheet.xlsx", "w") as workbook:
            for index in range(1, 52):
                value = "IID" if index == 51 else "summary"
                workbook.writestr(
                    f"xl/worksheets/sheet{index:03}.xml",
                    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                    '<sheetData><row r="301"><c r="KO301" t="inlineStr">'
                    f"<is><t>{value}</t></is></c></row></sheetData></worksheet>",
                )
        env = dict(os.environ, GITHUB_COPY_SOURCE_ROOT=str(source),
                   GITHUB_COPY_DEST_ROOT=str(dest))

        def run(action, success=True):
            result = subprocess.run(["bash", str(SCRIPT), action, "le8"], env=env,
                                    text=True, capture_output=True)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return result

        result = run("scan")
        manifest = dest / "github_files.lst"
        def approved_paths():
            return {line for line in manifest.read_text().splitlines()
                    if line and not line.startswith("#")}
        assert approved_paths() == {str(relative / name) for name in ("summary.csv", "quoted.csv")}
        assert "1000001" not in result.stdout + result.stderr
        assert not (dest / "le8").exists(), "scan must not copy files"
        before = manifest.read_bytes()
        result = run("sync")
        assert "UKB CHECK" not in result.stdout + result.stderr
        assert manifest.read_bytes() == before, "sync must not rewrite the manifest"
        assert {str(path.relative_to(dest)) for path in (dest / "le8").rglob("*") if path.is_file()} == approved_paths()
        # sync performs copying only; source changes require a fresh scan.
        (folder / "summary.csv").write_text("eid,value\n42,1\n")
        run("sync")
        assert (dest / relative / "summary.csv").read_bytes() == (folder / "summary.csv").read_bytes()
        (folder / "quoted.csv").write_text("IID\n")
        run("scan")
        assert not approved_paths()
        before = manifest.read_bytes()
        run("sync")
        assert manifest.read_bytes() == before
        assert (dest / relative / "summary.csv").exists(), "sync must not delete existing files"
        assert (dest / relative / "quoted.csv").exists(), "empty allowlist must not clear destination"
        assert not list(dest.glob(".github*")), "temporary files were not cleaned up"
        assert (folder / "result.RDS").is_file(), "source files must never be removed"
    print("PASS: depth, blocked types, full text/Excel ID detection, scan content exclusion, copy-only sync, cleanup")


if __name__ == "__main__":
    main()
