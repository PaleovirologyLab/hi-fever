import importlib.util
from pathlib import Path
import unittest

import pandas as pd


def _load_module():
    module_path = Path(__file__).resolve().parents[1] / "bin" / "create_summary_table.py"
    spec = importlib.util.spec_from_file_location("create_summary_table", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


cst = _load_module()


class TestCreateSummaryTable(unittest.TestCase):
    def test_canonical_assembly_id_uses_full_filename_stem(self):
        value = "/tmp/GCF_000000001.1_genomic.fna.gz"
        self.assertEqual(cst.canonical_assembly_id(value), "GCF_000000001.1_genomic")

    def test_canonical_assembly_id_keeps_non_accession_stem(self):
        value = "assembly_normalized.fasta"
        self.assertEqual(cst.canonical_assembly_id(value), "assembly_normalized")

    def test_clean_subject_title_removes_known_tokens_and_normalizes_spaces(self):
        title = "XP_12345.1 LOW QUALITY PROTEIN: capsid [Some virus]"
        self.assertEqual(cst.clean_subject_title(title), "capsid")

    def test_clean_subject_title_handles_missing_value(self):
        self.assertEqual(cst.clean_subject_title(pd.NA), "")

    def test_compute_loci_statistics_nr_aggregates_and_sorts(self):
        data = pd.DataFrame(
            [
                {
                    "query_locus": "locA",
                    "bitscore": 100.0,
                    "hit_length": 250,
                    "cleaned_title": "viral polymerase",
                    "hit_scientific_name": "Species1",
                    "hit_superkingdom": "Viruses",
                },
                {
                    "query_locus": "locA",
                    "bitscore": 90.0,
                    "hit_length": 150,
                    "cleaned_title": "host protein",
                    "hit_scientific_name": "Species2",
                    "hit_superkingdom": "Eukaryota",
                },
                {
                    "query_locus": "locB",
                    "bitscore": 50.0,
                    "hit_length": 100,
                    "cleaned_title": "host protein",
                    "hit_scientific_name": "Species3",
                    "hit_superkingdom": "Eukaryota",
                },
            ]
        )

        result = cst.compute_loci_statistics_nr(data)

        self.assertEqual(list(result["query_locus"]), ["locA", "locB"])
        loc_a = result[result["query_locus"] == "locA"].iloc[0]
        self.assertEqual(loc_a["n_hits"], 2)
        self.assertEqual(loc_a["best_bitscore"], 100.0)
        self.assertEqual(loc_a["best_title"], "viral polymerase")
        self.assertEqual(loc_a["n_species"], 2)
        self.assertEqual(loc_a["n_virus_hits"], 1)
        self.assertEqual(loc_a["n_eukaryot_hits"], 1)

    def test_classify_query_loci_paths(self):
        cases = [
            {
                "row": {
                    "n_hits_nr": 10,
                    "n_hits_rvdb": 0,
                    "n_retros": 0,
                    "n_virus_hits_rvdb": 0,
                    "n_virus_hits_nr": 0,
                    "n_eukaryot_hits_rvdb": 0,
                    "rvdb_kings": "",
                    "rvdb_fams": "",
                    "best_title_rvdb": "",
                    "common_title_nr": "ubiquitin like",
                    "common_title_rvdb": "",
                },
                "expected": "likely-host-protein",
            },
            {
                "row": {
                    "n_hits_nr": 8,
                    "n_hits_rvdb": 10,
                    "n_retros": 0,
                    "n_virus_hits_rvdb": 8,
                    "n_virus_hits_nr": 3,
                    "n_eukaryot_hits_rvdb": 1,
                    "rvdb_kings": "",
                    "rvdb_fams": "adenoviridae",
                    "best_title_rvdb": "capsid protein",
                    "common_title_nr": "viral protein",
                    "common_title_rvdb": "capsid",
                },
                "expected": "likely-eve",
            },
        ]

        for case in cases:
            with self.subTest(expected=case["expected"]):
                df = pd.DataFrame([case["row"]])
                result = cst.classify_query_loci(df)
                self.assertEqual(result.loc[0, "element_type"], case["expected"])

    def test_resolve_path_or_error(self):
        base = Path(__file__).resolve().parent
        fixture = base / "tmp_existing_file.txt"
        fixture.write_text("x\n", encoding="utf-8")
        try:
            self.assertEqual(cst.resolve_path_or_error(str(base), fixture.name), str(fixture))
            with self.assertRaises(FileNotFoundError):
                cst.resolve_path_or_error(str(base), "missing.tsv")
        finally:
            fixture.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
