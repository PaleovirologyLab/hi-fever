# Test Progress

## End-to-End Tests Added (Stub Run)
Type: `stub-run`

### What was added
- End-to-end stub-run tests that exercise the workflow graph using minimal local fixtures for reciprocal routing and validation behavior.

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

## Reciprocal Search Testing Plan (Working Notes)
Status: `in-progress`

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
  - `--custom_reciprocal true` routes to `SINGLE_RECIPROCAL_DIAMOND`. [done]
  - `--custom_reciprocal false` routes to `FULL_RECIPROCAL_DIAMOND`. [done]
  - Custom DB extension behavior:
    - `.dmnd` uses DB directly. [done]
    - `.fa/.fasta/.fna` triggers `BUILD_RECIPROCAL`. [partially done: `.fa` covered]
    - Unsupported extension fails with clear error. [done]
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

## Custom Reciprocal Input Branching Tests Added (Workflow E2E, Stub)
Type: `stub-run` + trace assertions

### What was added
- Three workflow-level tests in `tests/test_workflow_end_to_end.py` to validate custom reciprocal input branching:
  - FASTA input (`reciprocal.fa`): confirms `BUILD_RECIPROCAL` runs and `SINGLE_RECIPROCAL_DIAMOND` runs.
  - DMND input (`reciprocal.dmnd`): confirms `BUILD_RECIPROCAL` is skipped and `SINGLE_RECIPROCAL_DIAMOND` runs.
  - Unsupported extension (`reciprocal.txt`): confirms workflow fails with unsupported-extension validation error.
- Added fixture file:
  - `tests/fixtures/e2e/reciprocal.dmnd` (stub placeholder for routing test only).

### Files added/updated
- `tests/test_workflow_end_to_end.py`
- `tests/fixtures/e2e/reciprocal.dmnd`
- `modules/normalize_headers.nf` (compile fix to unblock workflow-level tests)

### Test command
```bash
python3 -m unittest discover -s tests -p 'test_workflow_end_to_end.py' -q
```

## Full Reciprocal Routing Test Added (`custom_reciprocal=false`)
Type: `stub-run` + trace assertions

### What was added
- Workflow-level test to confirm that when `custom_reciprocal` is not set (false path):
  - `FULL_RECIPROCAL_DIAMOND` runs.
  - `SINGLE_RECIPROCAL_DIAMOND` does not run.
- Test uses a temporary `data/` directory containing:
  - local fixture `assembly.fna` and `query.fa`
  - symlinks to default mini reciprocal databases in repo `data/`:
    - `MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd`
    - `MINI_rvdbv28_wtaxa.dmnd`

### Files updated
- `tests/test_workflow_end_to_end.py`

### Test command
```bash
python3 -m unittest discover -s tests -p 'test_workflow_end_to_end.py' -q
```

## TODO — Remaining Tests to Implement
- End-to-end run through `main.nf` after fixing `VERIFY` for local mode. [done]
- Reciprocal behavior:
  - Zero-hit handling in custom reciprocal mode (explicit expected behavior).
  - Best-hit selection/ranking semantics (`best_pairs`, `best_hits`).
  - Mixed-hit content assertions (`forward`, `reciprocal-nr`, `reciprocal-rvdb` labels).
  - Additional extension variants (`.fna`, `.fasta`) in custom reciprocal branch tests.
- Taxonomy behavior:
  - With `--email`.
  - With `--allow_missing_taxonomy` fallback.
- Output schema regression tests:
  - Summary table columns.
  - Reciprocal output table columns and `sql/` output filenames.
- Error handling tests:
  - Invalid `--assembly_file` / empty `--ftp_file`.
  - Missing required params.

## Entrypoint Validation Fix — 2026-05-27
Type: `workflow validation fix`

### What changed
- `subworkflows/verify.nf` local-mode validation no longer materializes `Channel.fromPath(...).toList().getVal()` during `VERIFY`.
- Local assembly inputs are now resolved synchronously via direct path matching, so `main.nf` can validate and proceed in local mode.
- Added a new workflow-level stub-run test that exercises `main.nf` directly in local custom-reciprocal DMND mode.

### Files updated
- `subworkflows/verify.nf`
- `tests/test_workflow_end_to_end.py`

### Verification
- `python3 -m unittest discover -s tests -p 'test_workflow_end_to_end.py' -q` -> passed
- Direct `nextflow run main.nf -stub-run ... --assembly_mode local ...` now reaches workflow task execution instead of failing in `VERIFY`.

### Note
- A direct local-mode stub run progressed into downstream tasks and then reported:
  - `HIFEVER:GENEWISE (1)` terminated with exit status `2`
  - the error was ignored by the workflow
- This is downstream of the entrypoint fix and remains a separate follow-up item.

## Execution-Profile Consistency — 2026-05-27
Type: `runtime portability fix`

### What changed
- `modules/normalize_headers.nf` now defines a `container` directive in addition to `conda`, making the local assembly normalization step consistent with the rest of the pipeline execution model.
- Added a gated workflow-level test for `main.nf` local mode under `-profile apptainer`.

### Files updated
- `modules/normalize_headers.nf`
- `tests/test_workflow_end_to_end.py`

### Current validation status
- The apptainer-profile test is present but skipped when `apptainer` is not installed in `PATH`.
- In the current environment, `apptainer` is not available, so container-profile execution has not yet been executed locally.

### Remaining validation
- Run the new `-profile apptainer` local-mode entrypoint test in an environment with `apptainer` installed.

## Taxonomy Policy and Branch Coverage — 2026-05-27
Type: `behavior policy + workflow tests`

### Policy decision applied
- Default behavior is now strict:
  - users should attempt to obtain taxonomy by default
  - degraded taxonomy fallback is only enabled when `--allow_missing_taxonomy true` is explicitly set
- `conf/parameters.config` now sets:
  - `allow_missing_taxonomy = false`

### Workflow changes supporting testability
- `modules/hits_taxonomy.nf`
  - added a `stub` block for `FETCH_HITS_TAXONOMY_FROM_ACCNS`
  - this allows the custom reciprocal + email taxonomy branch to be tested in `-stub-run` without external calls

### Workflow-level branch coverage added
- `tests/test_workflow_end_to_end.py`
  - custom reciprocal + no email + default strict policy -> fails
  - custom reciprocal + email -> `FETCH_HITS_TAXONOMY_FROM_ACCNS` runs
  - full reciprocal + missing taxonomy + `--allow_missing_taxonomy true` -> succeeds
  - full reciprocal + missing taxonomy + default strict policy -> emits the expected missing-taxonomy error

### Test harness hardening
- Workflow-level Nextflow tests now use per-test `-work-dir` paths.
- This avoids state bleed from `resume = true` across different taxonomy-policy cases.

### Verification
- `python3 -m unittest discover -s tests -p 'test_workflow_end_to_end.py' -q` -> passed
- Result:
  - `Ran 10 tests`
  - `OK (skipped=1)`

## Review Baseline — 2026-05-27
Type: `pre-PR review snapshot`

### What was checked
- Branch-level review against current `main` with focus on:
  - installation/runtime portability
  - README accuracy
  - behavior regressions
  - remaining bugs and untested paths
- Local verification run on this date:
  - `python3 -m unittest discover -s tests -p 'test_create_summary_table.py' -q` -> passed
  - `python3 -m unittest discover -s tests -p 'test_workflow_end_to_end.py' -q` -> passed
- A full `python3 -m unittest discover -s tests -q` run was attempted, but one run was polluted by a Nextflow session lock collision while multiple local Nextflow runs overlapped.

### Current status
- The branch is not considered PR-ready yet.
- The downstream `HIFEVER` workflow graph is partially validated.
- The `main.nf` entrypoint is not yet convincingly validated for the new local assembly mode.

### Main risks identified
- `VERIFY` remains a blocker for end-to-end confidence in local mode.
  - Workflow tests still call `workflows/hi-fever.nf` directly instead of `main.nf`.
  - See `subworkflows/verify.nf`.
- Local assembly mode adds a new process without container coverage.
  - `modules/normalize_headers.nf` defines `conda` but no `container`.
  - This creates a portability risk for container/apptainer-based execution.
- Taxonomy behavior changed materially.
  - `allow_missing_taxonomy = true` is now the default in `conf/parameters.config`.
  - Missing taxonomy inputs can now silently degrade outputs instead of failing fast.
- README is improved but not fully synchronized with the current branch/test scope.

### Pre-PR checklist (no design change)
- Fix `VERIFY` so `main.nf` works for local assembly mode.
- Add a real test that runs through `main.nf` for local mode after `VERIFY` is fixed.
- Add container support for `NORMALIZE_ASSEMBLY_HEADERS`, consistent with the rest of the pipeline execution model.
- Re-run the full test suite serially and record a clean result.
- Validate at least one container/apptainer-oriented execution path, not only local conda binaries.
- Decide whether `allow_missing_taxonomy = true` should remain the default.
- If the default is intentional, document the degraded-output behavior more explicitly in the README.
- Add tests for taxonomy gating/fallback behavior:
  - custom reciprocal + email
  - custom reciprocal + no email + allow-missing
  - custom reciprocal + no email + no allow-missing
  - full reciprocal + missing taxonomy table + allow-missing
  - full reciprocal + missing taxonomy table + fail
- Add coverage for `--assembly_metadata_file` in local mode.
- Update README test documentation to reflect the actual test coverage now present in `tests/`.
- Confirm output schema stability for summary tables and published `sql/` files.


  ## PR-Readiness Plan — Support Matrix

  Date: `2026-05-27`
  Type: `implementation plan`

  ### Summary
  Make the branch PR-ready by treating support as a matrix with two different concerns:

  - **Installation/bootstrap path**: `conda` and `pixi`
  - **Nextflow execution profile**: `-profile conda` and `-profile apptainer`

  `pixi` is supported in this repo as a bootstrap/task runner, not as a separate pipeline runtime. The implementation goal is therefore to ensure the pipeline
  behaves correctly under both runtime profiles, and that users entering via either conda or pixi can successfully run those supported profiles. Docker should be
  treated as a documentation/compatibility review item, not a primary runtime target, because the active Nextflow config does not define a docker profile.

  ### Support matrix to validate
  #### Required support combinations
  - `conda` bootstrap + `-profile conda`
  - `conda` bootstrap + `-profile apptainer`
  - `pixi` bootstrap + `-profile conda`
  - `pixi` bootstrap + `-profile apptainer`

  #### Acceptance rule for the matrix
  - The pipeline must work through `main.nf`, not only `workflows/hi-fever.nf`.
  - Local assembly mode must work under both runtime profiles.
  - Pixi validation can reuse the same pipeline behavior checks, but must prove that the repo’s declared Pixi tasks or equivalent Pixi-driven commands still
  function.

  #### Docker handling
  - Do not make Docker a merge blocker unless you want to preserve it as an explicitly supported runtime now.
  - Do treat Docker as a documentation consistency blocker:
    - if README continues to advertise Docker, the support statement must be true and specific
    - if not, README should be narrowed to the profiles actually implemented

  ### Implementation changes
  #### 1. Fix entrypoint reliability first
  - Repair `VERIFY` so `main.nf` local mode validates files/globs correctly and does not block execution.
  - Ensure validation resolution exactly matches runtime resolution in `HIFEVER` for:
    - `data_path`
    - local `assembly_file`
    - accepted FASTA/compressed FASTA extensions
  - Add one direct `main.nf` stub-run test for local mode as the first gating test.

  #### 2. Restore execution-profile consistency
  - Make `NORMALIZE_ASSEMBLY_HEADERS` runnable under both `-profile conda` and `-profile apptainer`.
  - Follow the same execution model as the rest of the modules:
    - declare both `conda` and `container`, or otherwise make the process profile-neutral in a way consistent with existing modules
  - Check for any other new local-mode-only steps that implicitly rely on host binaries or conda-only behavior.

  #### 3. Validate the support matrix in priority order
  ##### Tier 1: Merge blockers
  - `conda` bootstrap + `-profile conda`
    - `main.nf` local mode stub-run
    - local mode focused regression test
  - `conda` bootstrap + `-profile apptainer`
    - `main.nf` local mode stub-run
    - confirm `NORMALIZE_ASSEMBLY_HEADERS` works in containerized execution

  ##### Tier 2: Bootstrap confirmation
  - `pixi` bootstrap + `-profile conda`
    - run the equivalent of the local-mode entrypoint validation inside pixi
  - `pixi` bootstrap + `-profile apptainer`
    - run the equivalent of the local-mode entrypoint validation inside pixi
  - If pixi tasks are missing for local mode, document that gap and add validation commands to the PR checklist; no workflow redesign is needed.

  #### 4. Re-check shared downstream behavior
  - Confirm shared downstream modules still behave with local-mode inputs:
    - accession/meta ID handling
    - BLAST DB naming
    - forward DIAMOND outputs
    - summary-table host metadata joins

  #### 5. Finalize intentional user-facing behavior
  - Decide whether `allow_missing_taxonomy = true` remains the default.
  - Keep or revert the default explicitly; do not leave it as accidental behavior.
  - Add tests for whichever policy is chosen.
  - Make README reflect the exact behavior:
    - what runtime profiles are actually supported
    - what pixi does
    - what happens when taxonomy/metadata are missing in local mode

  ### Test plan
  #### Entrypoint and workflow tests
  - Add `main.nf` local-mode stub-run test.
  - Keep the current `HIFEVER` local-mode routing tests passing:
    - custom reciprocal FASTA
    - custom reciprocal DMND
    - invalid reciprocal extension
    - full reciprocal routing

  #### Profile-based validation
  - Run local mode under `-profile conda`.
  - Run local mode under `-profile apptainer`.
  - Run equivalent validations from a pixi-managed environment for:
    - `-profile conda`
    - `-profile apptainer`

  #### Behavior and regression coverage
  - Add taxonomy branch tests:
    - custom reciprocal + email
    - custom reciprocal + no email + allow-missing
    - custom reciprocal + no email + fail
    - full reciprocal + missing taxonomy + allow-missing
    - full reciprocal + missing taxonomy + fail
  - Add local metadata test for `--assembly_metadata_file`.
  - Confirm summary-table output still resolves host names correctly for:
    - canonical `GCA/GCF` assembly IDs
    - local non-accession filenames
    - placeholder metadata rows

  #### Final clean pass
  - Run the full test suite serially to avoid false failures from Nextflow session locking.
  - Record one clean final result in `TEST_PROGRESS.md` before opening the PR.

  ### Public interfaces / behavior
  - Supported runtime profiles: `conda`, `apptainer`
  - Supported bootstrap paths: `conda`, `pixi`
  - No new user parameters are required unless needed to expose an already-implicit runtime choice.
  - README must clearly separate:
    - how users install/bootstrap
    - how the pipeline executes
  - Docker must either be documented precisely as an external/manual path or de-emphasized if it is not a first-class configured runtime.

  ### Assumptions and defaults
  - Default support target for merge: `conda` and `apptainer` as runtime profiles, with both reachable from either conda or pixi bootstrap.
  - Default recommendation: treat pixi as a usability/support concern, but not as a separate workflow implementation branch.
  - Default recommendation: do not block the PR on adding a new docker profile unless you explicitly want Docker to remain first-class in the current release.
  - Priority order:
    1. fix `VERIFY`
    2. fix `NORMALIZE_ASSEMBLY_HEADERS` execution-profile consistency
    3. validate `conda` and `apptainer` local mode through `main.nf`
    4. validate pixi-driven execution of those same profiles
    5. finalize taxonomy default and docs
    6. run clean final test pass


## Current priorities are:

  1. main.nf local-mode reliability
     Reason: this is the real entrypoint. If users cannot run local mode through main.nf, the feature is not actually usable, regardless of downstream module
     quality.
  2. Execution-profile consistency for local mode
     Reason: this is your main usability concern. The new local-assembly path must work under the supported runtime profiles, not only in one environment.
     Otherwise installation mode changes behavior, which is exactly the class of regression you want to avoid.
  3. Pixi/conda bootstrap validation for those runtime profiles
     Reason: pixi is part of the repo’s declared support story. Even though it is not a separate runtime mode, users will enter through it, so we need to prove
     that the supported execution profiles are reachable from both bootstrap paths.
  4. Taxonomy default and taxonomy-branch behavior
     Reason: this branch changed failure behavior into fallback behavior in some cases. That affects user interpretation of results, so it must be intentional,
     tested, and documented.
  5. Local metadata and output stability
     Reason: local mode introduces new host/assembly naming paths. We need to make sure metadata joins, assembly IDs, and summary outputs stay correct for non-
     FTP inputs.
  6. README/test documentation alignment
     Reason: once behavior is stable, the docs need to describe the real support model. Otherwise users will still have a broken experience even if the code is
     correct.

  Things we explicitly deprioritized:

  - FTP end-to-end regression work
    Reason: you chose to remove that from the active pre-PR scope.
  - Docker as a first-class merge target
    Reason: the repo’s active Nextflow profiles are conda and apptainer; Docker is currently more of a documentation consistency issue than a primary execution
    path.

  Where we are now:

  - priority 1 is substantially addressed
  - priority 2 is partly addressed in code, but not fully validated because apptainer is not installed here
  - priorities 3 to 6 are still open

  So the next highest-value work is:

  - taxonomy policy and taxonomy-branch tests
  - then local metadata/output stability
  - then README alignment once those behaviors are fixed and confirmed

## Stopping Point — 2026-05-27
Type: `status snapshot`

### Completed in this session
- Fixed `VERIFY` so `main.nf` local mode validates inputs and reaches workflow execution.
- Added `main.nf` local-mode stub-run coverage.
- Made `NORMALIZE_ASSEMBLY_HEADERS` execution-profile consistent by adding a container directive.
- Added a gated `-profile apptainer` local-mode entrypoint test.
- Applied strict taxonomy-by-default policy:
  - `allow_missing_taxonomy = false`
  - fallback taxonomy is now opt-in via `--allow_missing_taxonomy true`
- Added workflow-level taxonomy branch tests for:
  - custom reciprocal + no email + strict default
  - custom reciprocal + email
  - full reciprocal + missing taxonomy + allow-missing
  - full reciprocal + missing taxonomy + strict default
- Hardened workflow-level Nextflow tests to use per-test `-work-dir` paths.

### Latest verification
- `python3 -m unittest discover -s tests -p 'test_workflow_end_to_end.py' -q` -> passed
- Result:
  - `Ran 10 tests`
  - `OK (skipped=1)`

### Known remaining notes
- `apptainer` is not installed in the current environment.
  - The apptainer-profile test exists but has not been executed locally.
- A direct local-mode stub run previously progressed into downstream tasks and reported:
  - `HIFEVER:GENEWISE (1)` terminated with exit status `2`
  - the workflow ignored the error
  - this remains a follow-up item outside the work completed here

### Remaining work
- Validate local mode under `-profile apptainer` in an environment with `apptainer` available.
- Validate pixi bootstrap for the supported execution profiles:
  - `pixi` + `-profile conda`
  - `pixi` + `-profile apptainer`
- Add local metadata/output stability coverage:
  - `--assembly_metadata_file`
  - assembly ID / host metadata joins
  - summary-table output stability for local-mode filenames
- Update README/test documentation so it matches:
  - strict taxonomy-by-default behavior
  - supported runtime profiles (`conda`, `apptainer`)
  - supported bootstrap paths (`conda`, `pixi`)
- Run the full test suite serially and record one clean final result.
