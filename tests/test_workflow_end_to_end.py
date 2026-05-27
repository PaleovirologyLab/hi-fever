import shutil
import subprocess
import tempfile
import unittest
from csv import DictReader
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_BIN = shutil.which("nextflow")
APPTAINER_BIN = shutil.which("apptainer")


@unittest.skipIf(NEXTFLOW_BIN is None, "nextflow is not installed in PATH")
class TestEndToEndWorkflow(unittest.TestCase):
    def _run_main_stub_local_custom_reciprocal(
        self,
        data_path: Path,
        reciprocal_db: str,
        outdir: Path,
        trace: Path,
        timeout=180,
        extra_args=None,
    ):
        cmd = [
            NEXTFLOW_BIN,
            "run",
            str(REPO_ROOT / "main.nf"),
            "-ansi-log",
            "false",
            "-work-dir",
            str(outdir.parent / "work"),
            "-stub-run",
            "-with-trace",
            str(trace),
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
            "--outdir",
            str(outdir),
        ]
        if extra_args:
            cmd.extend(extra_args)
        return subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )

    def _run_hifever_stub(
        self, data_path: Path, reciprocal_db: str, outdir: Path, trace: Path, timeout=180, extra_args=None
    ):
        cmd = [
            NEXTFLOW_BIN,
            "run",
            str(REPO_ROOT / "workflows" / "hi-fever.nf"),
            "-ansi-log",
            "false",
            "-work-dir",
            str(outdir.parent / "work"),
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
            "--outdir",
            str(outdir),
        ]
        if extra_args:
            cmd.extend(extra_args)
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

    def test_main_entrypoint_stub_run_local_custom_reciprocal_dmnd(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-main-e2e-") as tmpdir:
            outdir = Path(tmpdir) / "out"
            outdir.mkdir(parents=True, exist_ok=True)
            trace = Path(tmpdir) / "trace.tsv"

            result = self._run_main_stub_local_custom_reciprocal(
                fixtures,
                "reciprocal.dmnd",
                outdir,
                trace,
                extra_args=["--allow_missing_taxonomy", "true"],
            )
            if result.returncode != 0:
                self.fail(
                    "main.nf stub run failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )

            processes = self._trace_processes(trace)
            self.assertTrue(
                any("NORMALIZE_ASSEMBLY_HEADERS" in process for process in processes),
                "NORMALIZE_ASSEMBLY_HEADERS did not run through main.nf local mode",
            )
            self.assertTrue(
                any("SINGLE_RECIPROCAL_DIAMOND" in process for process in processes),
                "SINGLE_RECIPROCAL_DIAMOND did not run through main.nf local mode",
            )

    def test_main_entrypoint_stub_run_local_custom_reciprocal_dmnd_requires_email_by_default(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-main-tax-strict-") as tmpdir:
            outdir = Path(tmpdir) / "out"
            outdir.mkdir(parents=True, exist_ok=True)
            trace = Path(tmpdir) / "trace.tsv"

            result = self._run_main_stub_local_custom_reciprocal(fixtures, "reciprocal.dmnd", outdir, trace)
            self.assertNotEqual(result.returncode, 0, "Workflow should fail without email by default in custom reciprocal mode")
            merged_output = f"{result.stdout}\n{result.stderr}"
            self.assertIn("The '--email' parameter is required", merged_output)

    def test_main_entrypoint_stub_run_local_custom_reciprocal_email_triggers_taxonomy_fetch(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-main-tax-fetch-") as tmpdir:
            outdir = Path(tmpdir) / "out"
            outdir.mkdir(parents=True, exist_ok=True)
            trace = Path(tmpdir) / "trace.tsv"

            result = self._run_main_stub_local_custom_reciprocal(
                fixtures,
                "reciprocal.dmnd",
                outdir,
                trace,
                extra_args=["--email", "noreply@example.com"],
            )
            if result.returncode != 0:
                self.fail(
                    "main.nf stub run with email failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )

            processes = self._trace_processes(trace)
            self.assertTrue(
                any("FETCH_HITS_TAXONOMY_FROM_ACCNS" in process for process in processes),
                "FETCH_HITS_TAXONOMY_FROM_ACCNS did not run when email was provided",
            )

    @unittest.skipIf(APPTAINER_BIN is None, "apptainer is not installed in PATH")
    def test_main_entrypoint_stub_run_local_custom_reciprocal_dmnd_apptainer_profile(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-main-apptainer-e2e-") as tmpdir:
            outdir = Path(tmpdir) / "out"
            outdir.mkdir(parents=True, exist_ok=True)
            trace = Path(tmpdir) / "trace.tsv"

            cmd = [
                NEXTFLOW_BIN,
                "run",
                str(REPO_ROOT / "main.nf"),
                "-ansi-log",
                "false",
                "-work-dir",
                str(outdir.parent / "work"),
                "-stub-run",
                "-with-trace",
                str(trace),
                "-profile",
                "apptainer",
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
                "reciprocal.dmnd",
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
                timeout=180,
            )

            if result.returncode != 0:
                self.fail(
                    "main.nf apptainer-profile stub run failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )

            processes = self._trace_processes(trace)
            self.assertTrue(
                any("NORMALIZE_ASSEMBLY_HEADERS" in process for process in processes),
                "NORMALIZE_ASSEMBLY_HEADERS did not run through main.nf apptainer-profile local mode",
            )

    def _run_hifever_stub_full_reciprocal(
        self,
        data_path: Path,
        reciprocal_nr_db: str,
        reciprocal_rvdb_db: str,
        outdir: Path,
        trace: Path,
        timeout=180,
        extra_args=None,
    ):
        cmd = [
            NEXTFLOW_BIN,
            "run",
            str(REPO_ROOT / "workflows" / "hi-fever.nf"),
            "-ansi-log",
            "false",
            "-work-dir",
            str(outdir.parent / "work"),
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
            "--reciprocal_nr_db",
            reciprocal_nr_db,
            "--reciprocal_rvdb_db",
            reciprocal_rvdb_db,
            "--outdir",
            str(outdir),
        ]
        if extra_args:
            cmd.extend(extra_args)
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

            result = self._run_hifever_stub(fixtures, "reciprocal.fa", outdir, trace, extra_args=["--allow_missing_taxonomy", "true"])
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

            result = self._run_hifever_stub(fixtures, "reciprocal.dmnd", outdir, trace, extra_args=["--allow_missing_taxonomy", "true"])
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

            result = self._run_hifever_stub(data_path, "reciprocal.txt", outdir, trace, extra_args=["--allow_missing_taxonomy", "true"])
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
                extra_args=["--allow_missing_taxonomy", "true"],
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

    def test_end_to_end_stub_run_full_reciprocal_missing_taxonomy_allow_missing_succeeds(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"
        default_data = REPO_ROOT / "data"
        nr_db = default_data / "MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd"
        rvdb_db = default_data / "MINI_rvdbv28_wtaxa.dmnd"

        if not nr_db.exists() or not rvdb_db.exists():
            self.skipTest("Default MINI reciprocal databases were not found in data/")

        with tempfile.TemporaryDirectory(prefix="hi-fever-e2e-full-missing-tax-") as tmpdir:
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
                extra_args=["--allow_missing_taxonomy", "true"],
            )
            if result.returncode != 0:
                self.fail(
                    "Full reciprocal allow-missing stub run failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )

            processes = self._trace_processes(trace)
            self.assertTrue(
                any("FULL_RECIPROCAL_DIAMOND" in process for process in processes),
                "FULL_RECIPROCAL_DIAMOND did not run in allow-missing full reciprocal mode",
            )
            self.assertFalse(
                any("BUILD_HITS_TAXONOMY_TABLE" in process for process in processes),
                "BUILD_HITS_TAXONOMY_TABLE should not run when taxonomy file is missing and fallback is allowed",
            )

    def test_end_to_end_stub_run_full_reciprocal_missing_taxonomy_fails_by_default(self):
        fixtures = REPO_ROOT / "tests" / "fixtures" / "e2e"
        default_data = REPO_ROOT / "data"
        nr_db = default_data / "MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd"
        rvdb_db = default_data / "MINI_rvdbv28_wtaxa.dmnd"

        if not nr_db.exists() or not rvdb_db.exists():
            self.skipTest("Default MINI reciprocal databases were not found in data/")

        with tempfile.TemporaryDirectory(prefix="hi-fever-e2e-full-missing-tax-strict-") as tmpdir:
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
            merged_output = f"{result.stdout}\n{result.stderr}"
            self.assertIn("Provide taxonomy file or set '--allow_missing_taxonomy true'", merged_output)


if __name__ == "__main__":
    unittest.main()
