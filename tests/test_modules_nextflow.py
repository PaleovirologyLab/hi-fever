import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_BIN = shutil.which("nextflow")


@unittest.skipIf(NEXTFLOW_BIN is None, "nextflow is not installed in PATH")
class TestNextflowModules(unittest.TestCase):
    def run_nf(self, script: Path, params: dict, run_dir: Path):
        cmd = [
            NEXTFLOW_BIN,
            "run",
            str(script),
            "-ansi-log",
            "false",
            "-work-dir",
            str(run_dir / "work"),
        ]
        for key, value in params.items():
            cmd.extend([f"--{key}", str(value)])

        result = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            self.fail(
                "Nextflow test workflow failed\n"
                f"Command: {' '.join(cmd)}\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )

    def test_parse_ftp_module(self):
        script = REPO_ROOT / "tests" / "nf" / "parse_ftp_test.nf"
        ftp_input = REPO_ROOT / "tests" / "fixtures" / "module" / "ftp_list.txt"

        with tempfile.TemporaryDirectory(prefix="hi-fever-parse-ftp-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            self.run_nf(
                script,
                {
                    "ftp_input": ftp_input,
                    "outdir": outdir,
                },
                run_dir,
            )

            manifest = outdir / "parsed_manifest.txt"
            self.assertTrue(manifest.exists(), "parsed_manifest.txt was not created")

            names = {line.strip() for line in manifest.read_text(encoding="utf-8").splitlines() if line.strip()}
            self.assertEqual(names, {"A.fa.gz.ftp.txt", "B.fa.gz.ftp.txt"})

    def test_concatenate_publish_tables_module(self):
        script = REPO_ROOT / "tests" / "nf" / "concat_tables_test.nf"
        fixture_dir = REPO_ROOT / "tests" / "fixtures" / "module"

        with tempfile.TemporaryDirectory(prefix="hi-fever-concat-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            self.run_nf(
                script,
                {
                    "input_glob": fixture_dir / "table_*.tsv",
                    "table_name": "merged.tsv",
                    "outdir": outdir,
                },
                run_dir,
            )

            merged = outdir / "sql" / "merged.tsv"
            self.assertTrue(merged.exists(), "merged.tsv was not published")
            self.assertEqual(merged.read_text(encoding="utf-8"), "header1\na1\nb2\n")


if __name__ == "__main__":
    unittest.main()
