# Test Progress

## End-to-End Test Added (Stub Run)

### What was added
- End-to-end stub-run test that exercises the core workflow graph using a minimal local assembly and custom reciprocal database.

### Files added
- `tests/test_workflow_end_to_end.py`
- `tests/fixtures/e2e/query.fa`
- `tests/fixtures/e2e/reciprocal.fa`
- `tests/fixtures/e2e/assembly.fna`

### Why the test runs `workflows/hi-fever.nf`
- `main.nf` runs `VERIFY`, which currently errors on local-mode fixtures with:
  `ERROR: No assembly files matched .../assembly.fna`
- The e2e test runs the `HIFEVER` workflow directly to validate the pipeline graph until `VERIFY` is fixed.

### Fix applied to enable the stub run
- `workflows/hi-fever.nf`: added `def assembly_with_accession` to avoid Nextflow “already defined in process scope” error.

### Test command
```
python3 -m unittest discover -s tests -p 'test_workflow_end_to_end.py' -q
```
--- 

## TODO — Remaining Tests to Implement
- End-to-end run through `main.nf` after fixing `VERIFY` for local mode.
- FTP mode path test (parse + download + metadata path, with stubs or fixtures).
- Reciprocal mode tests:
  - Custom reciprocal FASTA (build dmnd).
  - Custom reciprocal dmnd (skip build).
  - Full reciprocal NR+RVDB path.
- Taxonomy behavior:
  - With `--email`.
  - With `--allow_missing_taxonomy` fallback.
- Module-level tests:
  - `FORWARD_DIAMOND` (fixture + schema assertions).
  - `EXTRACT_SEQS_ANNOTATE_MATCHES`.
  - `GENEWISE` output sanity.
  - `CREATE_SUMMARY_TABLE_*` integration.
- Output schema regression tests:
  - Summary table columns.
  - `sql/` output filenames.
- Error handling tests:
  - Invalid `--assembly_file` / empty `--ftp_file`.
  - Unsupported reciprocal DB extension.
  - Missing required params.
