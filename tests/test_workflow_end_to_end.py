import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_BIN = shutil.which("nextflow")


@unittest.skipIf(NEXTFLOW_BIN is None, "nextflow is not installed in PATH")
class TestEndToEndWorkflow(unittest.TestCase):
    def test_end_to_end_stub_run_local_custom_reciprocal(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-e2e-") as tmpdir:
            outdir = Path(tmpdir) / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            cmd = [
                NEXTFLOW_BIN,
                "run",
                str(REPO_ROOT / "workflows" / "hi-fever.nf"),
                "-ansi-log",
                "false",
                "-stub-run",
                "-entry",
                "HIFEVER",
                "--data_path",
                str(fixtures),
                "--assembly_mode",
                "local",
                "--assembly_file",
                "assembly.fna",
                "--query_file_aa",
                "query.fa",
                "--custom_reciprocal",
                "--custom_reciprocal_db",
                "reciprocal.fa",
                "--allow_missing_taxonomy",
                "true",
                "--outdir",
                str(outdir),
            ]

            result = subprocess.run(
                cmd,
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
                timeout=120,
            )
            if result.returncode != 0:
                self.fail(
                    "End-to-end stub run failed\n"
                    f"Command: {' '.join(cmd)}\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )


if __name__ == "__main__":
    unittest.main()
