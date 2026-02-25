# AGENTS.md

This project is a Nextflow DSL2 pipeline for finding endogenous viral elements (EVEs). Use these notes to navigate the codebase and make safe changes.

## What to Read First
- `README.md` for usage and dataset expectations.
- `main.nf` for workflow entrypoints.
- `workflows/hi-fever.nf` for the main pipeline logic.
- `subworkflows/verify.nf` for upfront validation.
- `conf/parameters.config` and `conf/containers.config` for parameters and execution profiles.

## Repo Layout and Roles
- `main.nf`: Entrypoint; runs `VERIFY()` then `HIFEVER()`.
- `workflows/hi-fever.nf`: Orchestrates the pipeline with DSL2 modules.
- `subworkflows/verify.nf`: Validates Nextflow version, required `--email`, and FASTA input extension.
- `modules/`: One process per file, used via `include { ... }` in workflows.
- `bin/`: Python scripts executed by processes (Nextflow automatically adds `bin/` to PATH).
- `conf/`: Default params and container/conda profiles.
- `data/`: Expected input data (queries, ftp list, taxonomy, DIAMOND dbs).
- `output/`: Default output target for published files.

## Execution Patterns
- DSL2 with explicit `include` statements in `workflows/hi-fever.nf`.
- Inputs use `Channel.fromPath("${params.data_path}/...")` with `checkIfExists: true`.
- Many processes publish to `"${params.outdir}/..."` using `publishDir`.
- Container + conda definitions are embedded per-process in `modules/*.nf`.
- The main pipeline branches on:
  - `params.cluster_query` to run `CLUSTER_SEQS` or use input FASTA directly.
  - `params.query_db` to reuse a user-provided `.dmnd` database.
  - `params.custom_reciprocal` to switch between custom reciprocal db or full NR+RVDB.

## Required Inputs and Params
- Required parameters:
  - `--email`: required by `subworkflows/verify.nf`.
- Required files (defaults from `conf/parameters.config`):
  - `--data_path` (default `data/`)
  - `--query_file_aa` (FASTA; must be `.fa/.fna/.fasta`)
  - `--ftp_file` (FTP list of assemblies)
- Optional but common:
  - `--reciprocal_nr_db`, `--reciprocal_rvdb_db` for full reciprocal DIAMOND.
  - `--custom_reciprocal` + `--custom_reciprocal_db` for a single reciprocal DB.
  - `--query_db` + `--user_dmnd_db` to use a prebuilt query DIAMOND DB.

## Key Data Flow (High-Level)
- Parse FTP list -> download assemblies -> assembly stats + metadata.
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

## Useful Files for Refactors
- `workflows/hi-fever.nf`: branch logic and channel wiring.
- `modules/diamond.nf`, `modules/forward_diamond.nf`: DIAMOND usage patterns.
- `modules/create_summary_table.nf`: Python script invocation and output publishing.

## Known Validation Constraints
- Nextflow must be `>=23.04.1` (`subworkflows/verify.nf`).
- `--email` is mandatory for Entrez API calls.
- Query FASTA file extension must be `.fa`, `.fna`, or `.fasta`.
- Custom reciprocal DB only supports `.fa/.fasta/.fna` or `.dmnd`.

## Typical Run (from repo root)
- `nextflow main.nf --query_file_aa <file.fa> --ftp_file <file.txt> --email you@example.com`
- Add `-with_docker` or `-profile conda` as needed (see `conf/containers.config`).

## No Tests Yet
There are no automated tests in this repo. Changes should be validated by running a small test dataset when possible.
