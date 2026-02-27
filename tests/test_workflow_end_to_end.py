import shutil
import subprocess
import tempfile
import unittest
from csv import DictReader
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_BIN = shutil.which("nextflow")


@unittest.skipIf(NEXTFLOW_BIN is None, "nextflow is not installed in PATH")
class TestEndToEndWorkflow(unittest.TestCase):
    def _run_hifever_stub(self, data_path: Path, reciprocal_db: str, outdir: Path, trace: Path, timeout=180):
        cmd = [
            NEXTFLOW_BIN,
            "run",
            str(REPO_ROOT / "workflows" / "hi-fever.nf"),
            "-ansi-log",
            "false",
            "-stub-run",
            "-with-trace",
            str(trace),
            "-entry",
            "HIFEVER",
            "--data_path",
            str(data_path),
            "--assembly_mode",
            "local",
            "--assembly_file",
            "assembly.fna",
            "--query_file_aa",
            "query.fa",
            "--custom_reciprocal",
            "--custom_reciprocal_db",
            reciprocal_db,
            "--allow_missing_taxonomy",
            "true",
            "--outdir",
            str(outdir),
        ]
        return subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )

    def _trace_processes(self, trace_file: Path):
        self.assertTrue(trace_file.exists(), "trace file was not created")
        with trace_file.open("r", encoding="utf-8") as handle:
            rows = list(DictReader(handle, delimiter="\t"))
        self.assertTrue(rows, "trace file has no task rows")
        return {row.get("name") or row.get("process") for row in rows if row.get("name") or row.get("process")}

    def _run_hifever_stub_full_reciprocal(
        self,
        data_path: Path,
        reciprocal_nr_db: str,
        reciprocal_rvdb_db: str,
        outdir: Path,
        trace: Path,
        timeout=180,
    ):
        cmd = [
            NEXTFLOW_BIN,
            "run",
            str(REPO_ROOT / "workflows" / "hi-fever.nf"),
            "-ansi-log",
            "false",
            "-stub-run",
            "-with-trace",
            str(trace),
            "-entry",
            "HIFEVER",
            "--data_path",
            str(data_path),
            "--assembly_mode",
            "local",
            "--assembly_file",
            "assembly.fna",
            "--query_file_aa",
            "query.fa",
            "--allow_missing_taxonomy",
            "true",
            "--reciprocal_nr_db",
            reciprocal_nr_db,
            "--reciprocal_rvdb_db",
            reciprocal_rvdb_db,
            "--outdir",
            str(outdir),
        ]
        return subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )

    def test_end_to_end_stub_run_local_custom_reciprocal_fasta_builds_db(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-e2e-") as tmpdir:
            outdir = Path(tmpdir) / "out"
            outdir.mkdir(parents=True, exist_ok=True)
            trace = Path(tmpdir) / "trace.tsv"

            result = self._run_hifever_stub(fixtures, "reciprocal.fa", outdir, trace)
            if result.returncode != 0:
                self.fail(
                    "End-to-end stub run failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )

            processes = self._trace_processes(trace)
            self.assertTrue(
                any("BUILD_RECIPROCAL" in process for process in processes),
                "BUILD_RECIPROCAL did not run for FASTA custom reciprocal input",
            )
            self.assertTrue(
                any("SINGLE_RECIPROCAL_DIAMOND" in process for process in processes),
                "SINGLE_RECIPROCAL_DIAMOND did not run for FASTA custom reciprocal input",
            )

    def test_end_to_end_stub_run_local_custom_reciprocal_dmnd_skips_build(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-e2e-") as tmpdir:
            outdir = Path(tmpdir) / "out"
            outdir.mkdir(parents=True, exist_ok=True)
            trace = Path(tmpdir) / "trace.tsv"

            result = self._run_hifever_stub(fixtures, "reciprocal.dmnd", outdir, trace)
            if result.returncode != 0:
                self.fail(
                    "End-to-end stub run failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )

            processes = self._trace_processes(trace)
            self.assertFalse(
                any("BUILD_RECIPROCAL" in process for process in processes),
                "BUILD_RECIPROCAL ran for DMND custom reciprocal input",
            )
            self.assertTrue(
                any("SINGLE_RECIPROCAL_DIAMOND" in process for process in processes),
                "SINGLE_RECIPROCAL_DIAMOND did not run for DMND custom reciprocal input",
            )

    def test_end_to_end_custom_reciprocal_invalid_extension_fails(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-e2e-invalid-ext-") as tmpdir:
            data_path = Path(tmpdir) / "data"
            outdir = Path(tmpdir) / "out"
            trace = Path(tmpdir) / "trace.tsv"
            data_path.mkdir(parents=True, exist_ok=True)
            outdir.mkdir(parents=True, exist_ok=True)

            shutil.copyfile(fixtures / "assembly.fna", data_path / "assembly.fna")
            shutil.copyfile(fixtures / "query.fa", data_path / "query.fa")
            (data_path / "reciprocal.txt").write_text("not_a_supported_db\n", encoding="utf-8")

            result = self._run_hifever_stub(data_path, "reciprocal.txt", outdir, trace)
            self.assertNotEqual(result.returncode, 0, "Workflow should fail for unsupported reciprocal DB extension")
            merged_output = f"{result.stdout}\n{result.stderr}"
            self.assertIn("Unsupported database file extension", merged_output)

    def test_end_to_end_stub_run_full_reciprocal_routes_to_full_module(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"
        default_data = REPO_ROOT / "data"
        nr_db = default_data / "MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd"
        rvdb_db = default_data / "MINI_rvdbv28_wtaxa.dmnd"

        if not nr_db.exists() or not rvdb_db.exists():
            self.skipTest("Default MINI reciprocal databases were not found in data/")

        with tempfile.TemporaryDirectory(prefix="hi-fever-e2e-full-reciprocal-") as tmpdir:
            data_path = Path(tmpdir) / "data"
            outdir = Path(tmpdir) / "out"
            trace = Path(tmpdir) / "trace.tsv"
            data_path.mkdir(parents=True, exist_ok=True)
            outdir.mkdir(parents=True, exist_ok=True)

            shutil.copyfile(fixtures / "assembly.fna", data_path / "assembly.fna")
            shutil.copyfile(fixtures / "query.fa", data_path / "query.fa")
            (data_path / nr_db.name).symlink_to(nr_db.resolve())
            (data_path / rvdb_db.name).symlink_to(rvdb_db.resolve())

            result = self._run_hifever_stub_full_reciprocal(
                data_path=data_path,
                reciprocal_nr_db=nr_db.name,
                reciprocal_rvdb_db=rvdb_db.name,
                outdir=outdir,
                trace=trace,
            )
            if result.returncode != 0:
                self.fail(
                    "End-to-end stub run failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )

            processes = self._trace_processes(trace)
            self.assertTrue(
                any("FULL_RECIPROCAL_DIAMOND" in process for process in processes),
                "FULL_RECIPROCAL_DIAMOND did not run when custom_reciprocal is false",
            )
            self.assertFalse(
                any("SINGLE_RECIPROCAL_DIAMOND" in process for process in processes),
                "SINGLE_RECIPROCAL_DIAMOND ran unexpectedly when custom_reciprocal is false",
            )


if __name__ == "__main__":
    unittest.main()
