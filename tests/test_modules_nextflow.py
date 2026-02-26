import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_BIN = shutil.which("nextflow")
DIAMOND_BIN = shutil.which("diamond")
SEQKIT_BIN = shutil.which("seqkit")
MAKEBLASTDB_BIN = shutil.which("makeblastdb")


@unittest.skipIf(NEXTFLOW_BIN is None, "nextflow is not installed in PATH")
class TestNextflowModules(unittest.TestCase):
    def run_nf(self, script: Path, params: dict, run_dir: Path, extra_args=None):
        cmd = [
            NEXTFLOW_BIN,
            "run",
            str(script),
            "-ansi-log",
            "false",
            "-work-dir",
            str(run_dir / "work"),
        ]
        if extra_args:
            cmd.extend(extra_args)
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

    def test_download_assemblies_module(self):
        script = REPO_ROOT / "tests" / "nf" / "ftp_download_test.nf"

        with tempfile.TemporaryDirectory(prefix="hi-fever-ftp-download-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            assembly_name = "GCF_000000000.1_genomic.fna.gz"
            ftp_list = run_dir / "ftp_list.txt"
            ftp_list.write_text("ftp://example.org/GCF_000000000.1\n", encoding="utf-8")

            self.run_nf(
                script,
                {
                    "ftp_input": ftp_list,
                    "outdir": outdir,
                },
                run_dir,
                extra_args=["-stub-run"],
            )

            manifest = outdir / "downloaded_manifest.txt"
            self.assertTrue(manifest.exists(), "downloaded_manifest.txt was not created")
            names = {line.strip() for line in manifest.read_text(encoding="utf-8").splitlines() if line.strip()}
            self.assertEqual(names, {assembly_name})

    def test_forward_diamond_module(self):
        script = REPO_ROOT / "tests" / "nf" / "forward_diamond_test.nf"
        fixture_dir = REPO_ROOT / "tests" / "fixtures" / "e2e"

        with tempfile.TemporaryDirectory(prefix="hi-fever-forward-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            self.run_nf(
                script,
                {
                    "assembly": fixture_dir / "assembly.fna",
                    "query_db": fixture_dir / "query.fa",
                    "outdir": outdir,
                },
                run_dir,
                extra_args=["-stub-run"],
            )

            manifest = outdir / "forward_manifest.txt"
            self.assertTrue(manifest.exists(), "forward_manifest.txt was not created")
            names = {line.strip() for line in manifest.read_text(encoding="utf-8").splitlines() if line.strip()}
            self.assertEqual(names, {"assembly_forward-matches-raw.dmnd.tsv"})

    @unittest.skipIf(DIAMOND_BIN is None or SEQKIT_BIN is None, "diamond/seqkit not installed in PATH")
    def test_forward_diamond_module_real(self):
        script = REPO_ROOT / "tests" / "nf" / "forward_diamond_test.nf"
        assembly = REPO_ROOT / "tests" / "fixtures" / "real" / "eptesicus_fuscus_genomic_region.fa"
        query = REPO_ROOT / "tests" / "fixtures" / "real" / "endogenous_borna_L_protein.fasta"

        if not assembly.exists() or not query.exists():
            self.skipTest("Required real-run fixtures not found in data/")

        with tempfile.TemporaryDirectory(prefix="hi-fever-forward-real-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            db_prefix = run_dir / "query_db"
            subprocess.run(
                [DIAMOND_BIN, "makedb", "--in", str(query), "-d", str(db_prefix)],
                check=True,
                text=True,
                capture_output=True,
            )
            db_file = Path(f"{db_prefix}.dmnd")
            self.assertTrue(db_file.exists(), "diamond DB was not created")

            self.run_nf(
                script,
                {
                    "assembly": assembly,
                    "query_db": db_file,
                    "outdir": outdir,
                    "diamond_forks": 1,
                    "chunk_size": 10000,
                    "diamond_mode": "fast",
                    "diamond_max_target_seqs": 5,
                },
                run_dir,
            )

            matches = list((run_dir / "work").rglob("*_forward-matches-raw.dmnd.tsv"))
            self.assertTrue(matches, "forward matches file was not created")
            content = matches[0].read_text(encoding="utf-8").strip()
            self.assertTrue(content, "forward matches file is empty")

    @unittest.skipIf(
        DIAMOND_BIN is None or SEQKIT_BIN is None or MAKEBLASTDB_BIN is None,
        "diamond/seqkit/makeblastdb not installed in PATH",
    )
    def test_extract_seqs_annotate_matches_module_real(self):
        forward_script = REPO_ROOT / "tests" / "nf" / "forward_diamond_test.nf"
        extract_script = REPO_ROOT / "tests" / "nf" / "extract_seqs_test.nf"
        assembly = REPO_ROOT / "tests" / "fixtures" / "real" / "eptesicus_fuscus_genomic_region.fa"
        query = REPO_ROOT / "tests" / "fixtures" / "real" / "endogenous_borna_L_protein.fasta"

        if not assembly.exists() or not query.exists():
            self.skipTest("Required real-run fixtures not found in tests/fixtures/real/")

        with tempfile.TemporaryDirectory(prefix="hi-fever-extract-real-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            normalized_assembly = run_dir / "assembly_normalized.fa"
            with assembly.open("r", encoding="utf-8") as in_handle, normalized_assembly.open(
                "w", encoding="utf-8"
            ) as out_handle:
                for line in in_handle:
                    if line.startswith(">"):
                        header = line[1:].strip().split()[0]
                        header = header.rsplit(":", 1)[0] if ":" in header else header
                        out_handle.write(f">{header}\n")
                    else:
                        out_handle.write(line)

            db_prefix = run_dir / "query_db"
            subprocess.run(
                [DIAMOND_BIN, "makedb", "--in", str(query), "-d", str(db_prefix)],
                check=True,
                text=True,
                capture_output=True,
            )
            db_file = Path(f"{db_prefix}.dmnd")
            self.assertTrue(db_file.exists(), "diamond DB was not created")

            self.run_nf(
                forward_script,
                {
                    "assembly": normalized_assembly,
                    "query_db": db_file,
                    "outdir": outdir,
                    "diamond_forks": 1,
                    "chunk_size": 10000,
                    "diamond_mode": "fast",
                    "diamond_max_target_seqs": 5,
                },
                run_dir,
            )

            match_files = list((run_dir / "work").rglob("*_forward-matches-raw.dmnd.tsv"))
            self.assertTrue(match_files, "forward matches file was not created")
            diamond_tsv = match_files[0]
            self.assertTrue(diamond_tsv.read_text(encoding="utf-8").strip(), "forward matches file is empty")

            subprocess.run(
                [
                    MAKEBLASTDB_BIN,
                    "-in",
                    str(normalized_assembly),
                    "-out",
                    str(run_dir / "assembly_db"),
                    "-dbtype",
                    "nucl",
                    "-parse_seqids",
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            nsq = run_dir / "assembly_db.nsq"
            self.assertTrue(nsq.exists(), "assembly BLAST DB was not created")

            self.run_nf(
                extract_script,
                {
                    "meta_id": "eptesicus_fuscus_genomic_region",
                    "diamond_tsv": diamond_tsv,
                    "assembly_db": nsq,
                    "outdir": outdir,
                    "interval": 1000,
                    "flank": 3000,
                },
                run_dir,
            )

            required_manifests = [
                "forward_annot_manifest.txt",
                "strict_manifest.txt",
                "context_manifest.txt",
                "locus_map_manifest.txt",
                "strict_coords_manifest.txt",
            ]
            for name in required_manifests:
                manifest = outdir / name
                self.assertTrue(manifest.exists(), f"{name} was not created")
                self.assertTrue(manifest.read_text(encoding="utf-8").strip(), f"{name} is empty")


if __name__ == "__main__":
    unittest.main()
