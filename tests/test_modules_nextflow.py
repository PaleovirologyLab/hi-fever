import os
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
GENEWISE_BIN = shutil.which("genewise")
BEDTOOLS_BIN = shutil.which("bedtools")
STOP_CONVERT_BIN = shutil.which("stopConvertAndCount.py")
TRANSLATE_CDS_BIN = shutil.which("translateCodingSequence.py")
REPO_BIN = REPO_ROOT / "bin"
LOCAL_STOP_CONVERT = REPO_BIN / "stopConvertAndCount.py"
LOCAL_TRANSLATE_CDS = REPO_BIN / "translateCodingSequence.py"
CREATE_SUMMARY_BIN = REPO_BIN / "create_summary_table.py"


@unittest.skipIf(NEXTFLOW_BIN is None, "nextflow is not installed in PATH")
class TestNextflowModules(unittest.TestCase):
    def run_nf(self, script: Path, params: dict, run_dir: Path, extra_args=None, env_extra=None):
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

        env = None
        if env_extra:
            env = dict(os.environ)
            env.update(env_extra)

        result = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
            env=env,
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

    @unittest.skipIf(
        DIAMOND_BIN is None
        or SEQKIT_BIN is None
        or MAKEBLASTDB_BIN is None
        or GENEWISE_BIN is None
        or BEDTOOLS_BIN is None
        or (STOP_CONVERT_BIN is None and not LOCAL_STOP_CONVERT.exists())
        or (TRANSLATE_CDS_BIN is None and not LOCAL_TRANSLATE_CDS.exists()),
        "Required tools for genewise test are not installed or not found in bin/",
    )
    def test_genewise_module_real(self):
        forward_script = REPO_ROOT / "tests" / "nf" / "forward_diamond_test.nf"
        extract_script = REPO_ROOT / "tests" / "nf" / "extract_seqs_test.nf"
        genewise_script = REPO_ROOT / "tests" / "nf" / "genewise_test.nf"
        assembly = REPO_ROOT / "tests" / "fixtures" / "real" / "eptesicus_fuscus_genomic_region.fa"
        query = REPO_ROOT / "tests" / "fixtures" / "real" / "endogenous_borna_L_protein.fasta"

        if not assembly.exists() or not query.exists():
            self.skipTest("Required real-run fixtures not found in tests/fixtures/real/")

        with tempfile.TemporaryDirectory(prefix="hi-fever-genewise-real-") as tmpdir:
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

            strict_fa = next((run_dir / "work").rglob("*_strict.fasta"))
            context_fa = next((run_dir / "work").rglob("*_context.fasta"))
            strict_bed = next((run_dir / "work").rglob("*_strict_coords.bed"))

            strict_header = None
            with strict_fa.open("r", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith(">"):
                        strict_header = line[1:].strip().split()[0]
                        break
            self.assertTrue(strict_header, "strict fasta is missing a header")

            contig, start, end = None, None, None
            with strict_bed.open("r", encoding="utf-8") as handle:
                first = handle.readline().strip().split("\t")
                contig, start, end = first[0], first[1], first[2]

            pair_subsets = run_dir / "best_pairs.txt"
            length = str(int(end) - int(start) + 1)
            pair_subsets.write_text(
                f"BAV60921.1\t{strict_header}\t{length}\tforward\t{contig}\t{start}\t{end}\n",
                encoding="utf-8",
            )

            context_coords = run_dir / "all_context_coords.bed"
            with context_fa.open("r", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith(">"):
                        header = line[1:].split()[0].split(":", 1)[1]
                        c_start, c_end = header.split("-", 1)
                        context_coords.write_text(
                            f"{contig}\t{c_start}\t{c_end}\n", encoding="utf-8"
                        )
                        break

            self.run_nf(
                genewise_script,
                {
                    "pair_subsets": pair_subsets,
                    "best_hit_proteins": query,
                    "strict_fastas": strict_fa,
                    "context_fastas": context_fa,
                    "context_coords": context_coords,
                    "outdir": outdir,
                },
                run_dir,
                env_extra={"PATH": f"{REPO_BIN}:{os.environ.get('PATH','')}"},
            )

            manifest = outdir / "genewise_manifest.txt"
            self.assertTrue(manifest.exists(), "genewise_manifest.txt was not created")
            self.assertTrue(manifest.read_text(encoding="utf-8").strip(), "genewise_manifest.txt is empty")

    def test_create_summary_table_full_synthetic(self):
        script = REPO_ROOT / "tests" / "nf" / "create_summary_full_test.nf"

        with tempfile.TemporaryDirectory(prefix="hi-fever-summary-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            reciprocal_nr = run_dir / "reciprocal-nr-matches.dmnd.tsv"
            reciprocal_rvdb = run_dir / "reciprocal-rvdb-matches.dmnd.tsv"
            taxonomy = run_dir / "hits_taxonomy.tsv"
            assembly_map = run_dir / "locus_assembly_map.tsv"
            assembly_metadata = run_dir / "assembly_metadata.tsv"
            genewise = run_dir / "genewise.tsv"

            query_locus = "NC_1:2001-4001"
            assembly_id = "GCF_000000001.1"

            reciprocal_row = "\t".join(
                [
                    query_locus,
                    "BAV60921.1",
                    "99.0",
                    "667",
                    "0",
                    "0",
                    "1",
                    "2001",
                    "1",
                    "2001",
                    "0.0",
                    "1293",
                    "123",
                    "TestVirus",
                    "Viruses",
                    "Viruses",
                    "Negarnaviricota",
                    "RNA-dependent RNA polymerase",
                ]
            )
            reciprocal_nr.write_text(reciprocal_row + "\n", encoding="utf-8")
            reciprocal_rvdb.write_text(reciprocal_row + "\n", encoding="utf-8")

            taxonomy.write_text(
                "\t".join(
                    [
                        "123",
                        "Viruses",
                        "Viruses",
                        "Negarnaviricota",
                        "Monjiviricetes",
                        "Mononegavirales",
                        "Bornaviridae",
                        "Orthobornavirus",
                        "TestVirus",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            assembly_map.write_text(f"{query_locus}\t{assembly_id}\n", encoding="utf-8")
            assembly_metadata.write_text(f"TestHost\t{assembly_id}\n", encoding="utf-8")

            # Genewise TSV columns (by position expected in create_summary_table.py):
            # 0 contig, 1 genomic_start, 2 genomic_end, 3 strand, 4 query_locus,
            # 5 sourceFASTA, 6 bitscore, 7 query, 8 qstart, 9 qend,
            # 10 cdna_seq, 11 peptide_seq, 12 intron_count, 13 idels_frameshifts, 14 inframe_STOPs
            genewise_cols = [
                "NC_1",
                "1",
                "2001",
                "+",
                query_locus,
                "strict",
                "1",
                "1",
                "1",
                "1",
                "ATGAAATAG",
                "MK",
                "0",
                "0",
                "0",
            ]
            genewise.write_text("\t".join(genewise_cols) + "\n", encoding="utf-8")

            self.run_nf(
                script,
                {
                    "reciprocal_nr": reciprocal_nr,
                    "reciprocal_rvdb": reciprocal_rvdb,
                    "taxonomy": taxonomy,
                    "assembly_map": assembly_map,
                    "assembly_metadata": assembly_metadata,
                    "genewise": genewise,
                    "outdir": outdir,
                },
                run_dir,
                env_extra={"PATH": f"{REPO_BIN}:{os.environ.get('PATH','')}"},
            )

            summary_manifest = outdir / "summary_manifest.txt"
            self.assertTrue(summary_manifest.exists(), "summary_manifest.txt was not created")
            summary_name = summary_manifest.read_text(encoding="utf-8").strip()
            self.assertTrue(summary_name, "summary manifest is empty")

            summary_path = next((run_dir / "work").rglob(summary_name), None)
            self.assertTrue(summary_path and summary_path.exists(), "summary output file not found")
            content = summary_path.read_text(encoding="utf-8")
            self.assertIn("element_type", content)
            self.assertIn("likely-eve", content)

    def test_create_summary_table_custom_synthetic(self):
        script = REPO_ROOT / "tests" / "nf" / "create_summary_custom_test.nf"

        with tempfile.TemporaryDirectory(prefix="hi-fever-summary-custom-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            reciprocal_rvdb = run_dir / "reciprocal-matches.dmnd.tsv"
            taxonomy = run_dir / "hits_taxonomy.tsv"
            assembly_map = run_dir / "locus_assembly_map.tsv"
            assembly_metadata = run_dir / "assembly_metadata.tsv"
            genewise = run_dir / "genewise.tsv"

            query_locus = "NC_1:2001-4001"
            assembly_id = "GCF_000000001.1"
            record_id = "BAV60921.1"

            reciprocal_row = "\t".join(
                [
                    query_locus,
                    record_id,
                    "99.0",
                    "667",
                    "0",
                    "0",
                    "1",
                    "2001",
                    "1",
                    "2001",
                    "0.0",
                    "1293",
                    "BAV60921.1|viral polymerase",
                    "MK",
                ]
            )
            reciprocal_rvdb.write_text(reciprocal_row + "\n", encoding="utf-8")

            taxonomy.write_text(
                "\t".join(
                    [
                        "record_id",
                        "all_taxonomy",
                        "family",
                        "viral_order",
                        "viral_kingdom",
                    ]
                )
                + "\n"
                + "\t".join(
                    [
                        record_id,
                        "Viruses; Negarnaviricota",
                        "Bornaviridae",
                        "Mononegavirales",
                        "Orthornavirae",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            assembly_map.write_text(f"{query_locus}\t{assembly_id}\n", encoding="utf-8")
            assembly_metadata.write_text(f"TestHost\t{assembly_id}\n", encoding="utf-8")

            # Genewise TSV columns (by position expected in create_summary_table.py):
            # 0 contig, 1 genomic_start, 2 genomic_end, 3 strand, 4 query_locus,
            # 5 sourceFASTA, 6 bitscore, 7 query, 8 qstart, 9 qend,
            # 10 cdna_seq, 11 peptide_seq, 12 intron_count, 13 idels_frameshifts, 14 inframe_STOPs
            genewise_cols = [
                "NC_1",
                "1",
                "2001",
                "+",
                query_locus,
                "strict",
                "1",
                "1",
                "1",
                "1",
                "ATGAAATAG",
                "MK",
                "0",
                "0",
                "0",
            ]
            genewise.write_text("\t".join(genewise_cols) + "\n", encoding="utf-8")

            self.run_nf(
                script,
                {
                    "reciprocal_rvdb": reciprocal_rvdb,
                    "taxonomy": taxonomy,
                    "assembly_map": assembly_map,
                    "assembly_metadata": assembly_metadata,
                    "genewise": genewise,
                    "outdir": outdir,
                },
                run_dir,
                env_extra={"PATH": f"{REPO_BIN}:{os.environ.get('PATH','')}"},
            )

            summary_manifest = outdir / "summary_manifest.txt"
            self.assertTrue(summary_manifest.exists(), "summary_manifest.txt was not created")
            summary_name = summary_manifest.read_text(encoding="utf-8").strip()
            self.assertTrue(summary_name, "summary manifest is empty")

            summary_path = next((run_dir / "work").rglob(summary_name), None)
            self.assertTrue(summary_path and summary_path.exists(), "summary output file not found")
            content = summary_path.read_text(encoding="utf-8")
            self.assertIn("element_type", content)

    def test_taxonomy_present_full_synthetic(self):
        script = REPO_ROOT / "tests" / "nf" / "create_summary_full_test.nf"

        with tempfile.TemporaryDirectory(prefix="hi-fever-tax-present-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            reciprocal_nr = run_dir / "reciprocal-nr-matches.dmnd.tsv"
            reciprocal_rvdb = run_dir / "reciprocal-rvdb-matches.dmnd.tsv"
            taxonomy = run_dir / "hits_taxonomy.tsv"
            assembly_map = run_dir / "locus_assembly_map.tsv"
            assembly_metadata = run_dir / "assembly_metadata.tsv"
            genewise = run_dir / "genewise.tsv"

            query_locus = "NC_1:2001-4001"
            assembly_id = "GCF_000000001.1"

            reciprocal_row = "\t".join(
                [
                    query_locus,
                    "BAV60921.1",
                    "99.0",
                    "667",
                    "0",
                    "0",
                    "1",
                    "2001",
                    "1",
                    "2001",
                    "0.0",
                    "1293",
                    "123",
                    "TestVirus",
                    "Viruses",
                    "Viruses",
                    "Negarnaviricota",
                    "RNA-dependent RNA polymerase",
                ]
            )
            reciprocal_nr.write_text(reciprocal_row + "\n", encoding="utf-8")
            reciprocal_rvdb.write_text(reciprocal_row + "\n", encoding="utf-8")

            taxonomy.write_text(
                "\t".join(
                    [
                        "123",
                        "Viruses",
                        "Viruses",
                        "Negarnaviricota",
                        "Monjiviricetes",
                        "Mononegavirales",
                        "Bornaviridae",
                        "Orthobornavirus",
                        "TestVirus",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            assembly_map.write_text(f"{query_locus}\t{assembly_id}\n", encoding="utf-8")
            assembly_metadata.write_text(f"TestHost\t{assembly_id}\n", encoding="utf-8")

            genewise_cols = [
                "NC_1",
                "1",
                "2001",
                "+",
                query_locus,
                "strict",
                "1",
                "1",
                "1",
                "1",
                "ATGAAATAG",
                "MK",
                "0",
                "0",
                "0",
            ]
            genewise.write_text("\t".join(genewise_cols) + "\n", encoding="utf-8")

            self.run_nf(
                script,
                {
                    "reciprocal_nr": reciprocal_nr,
                    "reciprocal_rvdb": reciprocal_rvdb,
                    "taxonomy": taxonomy,
                    "assembly_map": assembly_map,
                    "assembly_metadata": assembly_metadata,
                    "genewise": genewise,
                    "outdir": outdir,
                },
                run_dir,
                env_extra={"PATH": f"{REPO_BIN}:{os.environ.get('PATH','')}"},
            )

            summary_manifest = outdir / "summary_manifest.txt"
            summary_name = summary_manifest.read_text(encoding="utf-8").strip()
            summary_path = next((run_dir / "work").rglob(summary_name), None)
            self.assertTrue(summary_path and summary_path.exists(), "summary output file not found")
            content = summary_path.read_text(encoding="utf-8")
            self.assertIn("Bornaviridae", content)

    def test_taxonomy_missing_full_synthetic(self):
        script = REPO_ROOT / "tests" / "nf" / "create_summary_full_test.nf"

        with tempfile.TemporaryDirectory(prefix="hi-fever-tax-missing-") as tmpdir:
            run_dir = Path(tmpdir)
            outdir = run_dir / "out"
            outdir.mkdir(parents=True, exist_ok=True)

            reciprocal_nr = run_dir / "reciprocal-nr-matches.dmnd.tsv"
            reciprocal_rvdb = run_dir / "reciprocal-rvdb-matches.dmnd.tsv"
            taxonomy = run_dir / "hits_taxonomy.tsv"
            assembly_map = run_dir / "locus_assembly_map.tsv"
            assembly_metadata = run_dir / "assembly_metadata.tsv"
            genewise = run_dir / "genewise.tsv"

            query_locus = "NC_1:2001-4001"
            assembly_id = "GCF_000000001.1"

            reciprocal_row = "\t".join(
                [
                    query_locus,
                    "BAV60921.1",
                    "99.0",
                    "667",
                    "0",
                    "0",
                    "1",
                    "2001",
                    "1",
                    "2001",
                    "0.0",
                    "1293",
                    "123",
                    "TestVirus",
                    "Viruses",
                    "Viruses",
                    "Negarnaviricota",
                    "RNA-dependent RNA polymerase",
                ]
            )
            reciprocal_nr.write_text(reciprocal_row + "\n", encoding="utf-8")
            reciprocal_rvdb.write_text(reciprocal_row + "\n", encoding="utf-8")

            taxonomy.write_text(
                "\t".join(
                    [
                        "123",
                        "N/A",
                        "N/A",
                        "N/A",
                        "N/A",
                        "N/A",
                        "N/A",
                        "N/A",
                        "N/A",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            assembly_map.write_text(f"{query_locus}\t{assembly_id}\n", encoding="utf-8")
            assembly_metadata.write_text(f"TestHost\t{assembly_id}\n", encoding="utf-8")

            genewise_cols = [
                "NC_1",
                "1",
                "2001",
                "+",
                query_locus,
                "strict",
                "1",
                "1",
                "1",
                "1",
                "ATGAAATAG",
                "MK",
                "0",
                "0",
                "0",
            ]
            genewise.write_text("\t".join(genewise_cols) + "\n", encoding="utf-8")

            self.run_nf(
                script,
                {
                    "reciprocal_nr": reciprocal_nr,
                    "reciprocal_rvdb": reciprocal_rvdb,
                    "taxonomy": taxonomy,
                    "assembly_map": assembly_map,
                    "assembly_metadata": assembly_metadata,
                    "genewise": genewise,
                    "outdir": outdir,
                },
                run_dir,
                env_extra={"PATH": f"{REPO_BIN}:{os.environ.get('PATH','')}"},
            )

            summary_manifest = outdir / "summary_manifest.txt"
            summary_name = summary_manifest.read_text(encoding="utf-8").strip()
            summary_path = next((run_dir / "work").rglob(summary_name), None)
            self.assertTrue(summary_path and summary_path.exists(), "summary output file not found")
            content = summary_path.read_text(encoding="utf-8")
            self.assertIn("element_type", content)
            self.assertIn("nan", content)


if __name__ == "__main__":
    unittest.main()
