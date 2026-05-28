# AGENTS.md

This project is a Nextflow DSL2 pipeline for finding endogenous viral elements (EVEs). Use these notes to navigate the codebase and make safe changes.

## What to Read First
- `README.md` for usage and dataset expectations.
- `main.nf` for workflow entrypoints.
- `workflows/hi-fever.nf` for the main pipeline logic.
- `subworkflows/verify.nf` for upfront validation.
- `conf/parameters.config`, `conf/containers.config`, and `nextflow.config` for parameters and execution profiles.

## Repo Layout and Roles
- `main.nf`: Entrypoint; runs `VERIFY()` then `HIFEVER()`.
- `workflows/hi-fever.nf`: Orchestrates the pipeline with DSL2 modules.
- `subworkflows/verify.nf`: Validates Nextflow version, required `--email`, assembly mode inputs, and FASTA/DB extensions.
- `modules/`: One process per file, used via `include { ... }` in workflows.
- `bin/`: Python scripts executed by processes (Nextflow automatically adds `bin/` to PATH).
- `conf/`: Default params and runtime profiles.
- `data/`: Expected input data (queries, ftp list, taxonomy, DIAMOND dbs).
- `output/`: Default output target for published files.
- `tests/`: Python unit tests plus Nextflow module/workflow tests.

## Execution Patterns
- DSL2 with explicit `include` statements in `workflows/hi-fever.nf`.
- Inputs use `Channel.fromPath("${params.data_path}/...")` with `checkIfExists: true`.
- Many processes publish to `"${params.outdir}/..."` using `publishDir`.
- Container + conda definitions are embedded per-process in `modules/*.nf`.
- The main pipeline branches on:
  - `params.cluster_query` to run `CLUSTER_SEQS` or use input FASTA directly.
  - `params.query_db` to reuse a user-provided `.dmnd` database.
  - `params.custom_reciprocal` to switch between custom reciprocal db or full NR+RVDB.
  - `params.assembly_mode` to switch between FTP-mode assemblies and local FASTA assemblies.

## Required Inputs and Params
- Required parameters:
  - `--email`: required by `subworkflows/verify.nf`.
- Required files (defaults from `conf/parameters.config`):
  - `--data_path` (default `data/`)
  - `--query_file_aa` (FASTA; must be `.fa/.fna/.fasta`)
- One assembly source is required:
  - FTP mode: `--ftp_file` (FTP list of assemblies)
  - local mode: `--assembly_mode local --assembly_file <file-or-glob>`
- Optional but common:
  - `--reciprocal_nr_db`, `--reciprocal_rvdb_db` for full reciprocal DIAMOND.
  - `--custom_reciprocal` + `--custom_reciprocal_db` for a single reciprocal DB.
  - `--query_db` + `--user_dmnd_db` to use a prebuilt query DIAMOND DB.
  - `--assembly_metadata_file` for local-mode host metadata (`hostName`, `assembly_id`).
  - `--allow_missing_taxonomy true` to opt into degraded taxonomy behavior when taxonomy cannot be built.

## Key Data Flow (High-Level)
- Parse FTP list -> download assemblies -> assembly stats + metadata.
- Or normalize local assembly headers and derive assembly IDs from the filename stem.
- Build/choose query DIAMOND DB -> forward DIAMOND on genome chunks.
- Extract loci -> merge loci -> reciprocal DIAMOND -> taxonomy.
- Run GeneWise reconstruction -> publish tables and FASTA outputs.
- Create summary tables (`CREATE_SUMMARY_TABLE_*`) using Python in `bin/`.

## Conventions to Follow When Editing
- Keep module/process names in `modules/*.nf` uppercase (e.g., `FORWARD_DIAMOND`).
- Use `tag` to make task logs readable (commonly `${meta.id}` or input file).
- Prefer `Channel.fromPath` with `checkIfExists: true` for user-supplied inputs.
- If adding a new script, place it in `bin/` and call it from a process block.
- If adding new params, update `conf/parameters.config` and document in `README.md`.
- If adding new containers or profiles, update `conf/containers.config`.
- Keep assembly IDs stable and informative:
  - use the input filename stem with only compression/FASTA suffixes removed
  - do not assume all assembly names follow NCBI accession syntax

## Useful Files for Refactors
- `workflows/hi-fever.nf`: branch logic and channel wiring.
- `modules/diamond.nf`, `modules/forward_diamond.nf`: DIAMOND usage patterns.
- `modules/create_summary_table.nf`: Python script invocation and output publishing.

## Known Validation Constraints
- Nextflow must be `>=23.04.1` (`subworkflows/verify.nf`).
- `--email` is mandatory for Entrez API calls.
- Query FASTA file extension must be `.fa`, `.fna`, or `.fasta`.
- Custom reciprocal DB only supports `.fa/.fasta/.fna` or `.dmnd`.
- `allow_missing_taxonomy` defaults to `false`; degraded taxonomy mode is opt-in.
- Local-mode assembly IDs are derived from the full filename stem, not accession regex parsing.

## Typical Run (from repo root)
- `nextflow main.nf --query_file_aa <file.fa> --ftp_file <file.txt> --email you@example.com`
- Add `-profile conda` or `-profile apptainer` as needed (see `conf/containers.config`).
- Local mode example:
  - `nextflow main.nf --assembly_mode local --assembly_file "assemblies/*.fna.gz" --query_file_aa <file.fa> --custom_reciprocal --custom_reciprocal_db <db.dmnd> --email you@example.com -profile conda`

## Tests
- `python3 -m unittest discover -s tests -q` runs the full automated test suite.
- `tests/test_create_summary_table.py` covers summary-table logic.
- `tests/test_modules_nextflow.py` covers module-level Nextflow behavior.
- `tests/test_workflow_end_to_end.py` covers workflow-level routing and local-mode entrypoint behavior.
- Some environment-specific tests are skipped automatically when tools such as `apptainer` are not installed.
