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

## FTP Download Test Added (Module Level)

### What was added
- Module-level test that runs `PARSE_FTP` + `DOWNLOAD_ASSEMBLIES` against a local `file://` FTP fixture.

### Files added
- `tests/nf/ftp_download_test.nf`

### Test command
```
python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```

## Forward DIAMOND Test Added (Module Level)

### What was added
- Stub-run wiring test plus a real-run test using `data/eptesicus_fuscus_genomic_region.fa`
  and `data/endonous_borna_L_protein.fasta`.
- Real-run test builds a small DIAMOND DB on the fly; skipped if `diamond` or `seqkit`
  are not in `PATH`.

### Files added/updated
- `tests/nf/forward_diamond_test.nf`
- `tests/test_modules_nextflow.py`
- `modules/forward_diamond.nf` (added `stub` and explicit `script:`)

### Test command
```
python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```

## TODO — Remaining Tests to Implement
- End-to-end run through `main.nf` after fixing `VERIFY` for local mode.
- Reciprocal mode tests:
  - Custom reciprocal FASTA (build dmnd).
  - Custom reciprocal dmnd (skip build).
  - Full reciprocal NR+RVDB path.
- Taxonomy behavior:
  - With `--email`.
  - With `--allow_missing_taxonomy` fallback.
- Module-level tests:
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
