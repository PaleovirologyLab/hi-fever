# Test Progress

## End-to-End Test Added (Stub Run)
Type: `stub-run`

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
Type: `stub-run`

### What was added
- Module-level test that runs `PARSE_FTP` + `DOWNLOAD_ASSEMBLIES` against a local `file://` FTP fixture.

### Files added
- `tests/nf/ftp_download_test.nf`

### Test command
```
python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```

## Forward DIAMOND Test Added (Module Level)
Type: `real-run` + `stub-run`

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

## Extract Seqs + Annotate Matches Test Added (Module Level)
Type: `real-run`

### What was added
- Real-run module test that generates a DIAMOND TSV and BLAST DB from the real fixtures,
  then runs `EXTRACT_SEQS_ANNOTATE_MATCHES`.
- The test normalizes the assembly header locally so `blastdbcmd` can resolve seqids.

### Files added/updated
- `tests/nf/extract_seqs_test.nf`
- `tests/test_modules_nextflow.py`

### Test command
```
python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```

## Genewise Test Added (Module Level)
Type: `real-run`

### What was added
- Real-run module test that derives inputs from the forward + extract tests, then
  runs `GENEWISE` on the resulting pair/FASTA/coords files.
- The test is skipped if `genewise`, `bedtools`, `seqtk`, `makeblastdb`,
  `stopConvertAndCount.py`, or `translateCodingSequence.py` are not in `PATH`.

### Files added/updated
- `tests/nf/genewise_test.nf`
- `tests/test_modules_nextflow.py`

### Test command
```
python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```

## Create Summary Table (Full) Test Added (Module Level)
Type: `synthetic`

### What was added
- Synthetic module test for `CREATE_SUMMARY_TABLE_FULL` with minimal TSV inputs.
- Asserts that the summary output contains `element_type` and includes `likely-eve`.

### Files added/updated
- `tests/nf/create_summary_full_test.nf`
- `tests/test_modules_nextflow.py`

### Test command
```
python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```

## Create Summary Table (Custom) Test Added (Module Level)
Type: `synthetic`

## Taxonomy Behavior Tests (Full, Synthetic)
Type: `synthetic`

### What was added
- Two synthetic tests:
  - taxonomy present (expects family like `Bornaviridae`)
  - taxonomy missing (expects `N/A` values)

### Files added/updated
- `tests/test_modules_nextflow.py`

### Test command
```
python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```

### What was added
- Synthetic module test for `CREATE_SUMMARY_TABLE_CUSTOM` with minimal TSV inputs.
- Asserts that the summary output contains `element_type`.

### Files added/updated
- `tests/nf/create_summary_custom_test.nf`
- `tests/test_modules_nextflow.py`

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
- Output schema regression tests:
  - Summary table columns.
  - `sql/` output filenames.
- Error handling tests:
  - Invalid `--assembly_file` / empty `--ftp_file`do 
  - Unsupported reciprocal DB extension.
  - Missing required params.

## Reciprocal Search Testing Plan (Working Notes)
Status: `planned`

### Scope to cover
- Module behavior for `SINGLE_RECIPROCAL_DIAMOND` (custom reciprocal mode).
- Module behavior for `FULL_RECIPROCAL_DIAMOND` (NR + RVDB mode).
- Workflow branch selection and validation logic in `workflows/hi-fever.nf`.
- Taxonomy fallback/fail behavior tied to reciprocal paths.

### Test logic to follow
1. Start with `stub-run` tests to verify branch wiring and output file contracts.
2. Add synthetic tests for ranking/dedup semantics (`best_*` and `mixed_hits` behavior).
3. Add one minimal real-tool integration test (skip when tools are missing) to confirm command compatibility.
4. Keep fixtures small and explicit so failures isolate pipeline logic.
5. Assert semantic correctness in addition to file existence.

### Detailed checklist
- `SINGLE_RECIPROCAL_DIAMOND`:
  - Emits `reciprocal-matches.dmnd.tsv`, `reciprocal_hits.txt`,
    `best_reciprocal_hits.txt`, `reciprocal_seqs.fasta`.
  - `best_reciprocal_hits.txt` contains one best hit per query locus.
  - `reciprocal_hits.txt` retains all reciprocal hits.
- `FULL_RECIPROCAL_DIAMOND`:
  - Emits `reciprocal-nr-matches.dmnd.tsv` and `reciprocal-rvdb-matches.dmnd.tsv`.
  - Emits `mixed_hits.txt`, `best_pairs.txt`, `best_hits.fasta`.
  - `mixed_hits.txt` includes `forward`, `reciprocal-nr`, and `reciprocal-rvdb` labels.
- Workflow branching:
  - `--custom_reciprocal true` routes to `SINGLE_RECIPROCAL_DIAMOND`.
  - `--custom_reciprocal false` routes to `FULL_RECIPROCAL_DIAMOND`.
  - Custom DB extension behavior:
    - `.dmnd` uses DB directly.
    - `.fa/.fasta/.fna` triggers `BUILD_RECIPROCAL`.
    - Unsupported extension fails with clear error.
- Taxonomy gating:
  - Custom reciprocal + `--email` executes taxonomy fetch.
  - Custom reciprocal + no email + `--allow_missing_taxonomy true` uses placeholder taxonomy.
  - Custom reciprocal + no email + no allow-missing fails with expected error.
  - Full reciprocal + missing taxonomy table + allow-missing true uses placeholder taxonomy.
- Full reciprocal + missing taxonomy table + allow-missing false fails.

### Notes
- Determinism/regression checks are useful: repeated runs should keep `best_pairs` and
  selected `best_hits` stable for the same input.
- Empty-hit behavior should be explicit:
  - For user custom reciprocal DBs, no-hit runs can be valid and should not be treated as a failure by default.
  - Tests should assert expected behavior for empty outputs where appropriate.

## Reciprocal Positive Control Test Added (Real Integration, Gated)
Type: `real-run` (heavy, optional)

### What was added
- Positive-control integration test that derives loci from:
  - `tests/fixtures/real/eptesicus_fuscus_genomic_region.fa` (forward + extract)
  - `tests/fixtures/real/endogenous_borna_L_protein.fasta` (query proteins)
- Then runs `FULL_RECIPROCAL_DIAMOND` against:
  - `data/MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd`
  - `data/MINI_rvdbv28_wtaxa.dmnd`
- Asserts reciprocal outputs for both NR and RVDB are present and non-empty.

### Files added/updated
- `tests/nf/full_reciprocal_test.nf`
- `tests/test_modules_nextflow.py`

### How to run
```bash
HIFEVER_RUN_RECIPROCAL_INTEGRATION=1 python3 -m unittest discover -s tests -p 'test_modules_nextflow.py' -q
```
