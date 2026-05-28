[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)
![Nextflow](https://img.shields.io/badge/Nextflow-24.10.5-brightgreen)

# HI-FEVER

## **Hi**gh-throughput next**f**low **EVE** **r**ecovery
`hi-fever` is a Nextflow workflow for finding endogenous viral elements (EVEs) in host genomes. It aims to address common issues in paleovirology including cross-matches between host proteins and EVEs, computational burden of EVE searches and incompatability between software packages or platforms. We provide HI-FEVER as an accessible and informative workflow for any EVE-discovery project.

## Features

- Protein-to-DNA based search allows detection of divergent and ancient EVEs
- Designed to function with millions of input query proteins
- Reconstructs the predicted EVE protein based on its closest modern match
- Harnesses parallelisation to optimise compute resources
- Scales from laptop to cluster
- Conda, Pixi, and Apptainer compatible
- Designed for Linux-first execution; macOS is suitable for development and Conda-based runs, and Windows use should be via WSL

HI-FEVER provides a variety of output information about candidate EVEs, suited to many downstream purposes. Outputs include:
- Genomic coordinates of candidate EVEs
- Closest matches in the reciprocal databases, including full taxonomical information
- Predicted EVE protein sequences and cDNA (frameshift and premature STOP codon aware), with extension beyond original hit
- Extracted nucleotide sequence of each candidate EVE and flanking host genome sequence
- Metadata & statistics of the genome assemblies screened

## Installation and usage

HI-FEVER is available for use on Linux and can also be developed or run in lighter configurations on macOS. Windows use should be through WSL. The repository currently supports:
- bootstrap/setup via Conda or Pixi
- pipeline execution via Nextflow `-profile conda` or `-profile apptainer`

Recommended starting point:
- use `-profile conda` for a first local run
- use `-profile apptainer` when Apptainer is available and you want a more reproducible container-based execution

Dockerfiles are present in `docker/`, but Docker is not currently configured as a first-class Nextflow profile in `conf/containers.config`.

Full documentation can be found in [the wiki](https://github.com/PaleovirologyLab/hi-fever/wiki).

## Tests

Run automated tests from the repository root:

`python3 -m unittest discover -s tests -q`

This includes:
- unit tests for `bin/create_summary_table.py`
- workflow-level stub-run tests for local assembly mode and reciprocal routing
- module-level Nextflow tests for parsing, downloads, DIAMOND, extraction, Genewise, and summary-table generation

Note: Nextflow-based tests are automatically skipped if `nextflow` is not available in `PATH`. Some environment-specific tests are also skipped when required tools such as `apptainer` are not installed.

Pixi bootstrap smoke-test example:

`pixi run nextflow run main.nf -stub-run --data_path tests/fixtures/e2e --assembly_mode local --assembly_file assembly.fna --query_file_aa query.fa --custom_reciprocal --custom_reciprocal_db reciprocal.dmnd --email you@example.com -profile conda`

## Test run
To experiment with and explore HI-FEVER options we provide instructions on running a small bundled dataset below. Larger sample materials are also available on our Open Science Framework repository [here](https://osf.io/y357r/).

**Preparation**

Ensure the following files are available in `data/`:
* `20_per_fam_no_retro.fasta` protein query file
* `one_genome.ftp` or `genomes_n10_ftp.txt` ftp input file
* `taxdump.tar.gz` taxonomy map file
* `MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd` minimal reciprocal NR database
* `MINI_rvdbv28_wtaxa.dmnd` minimal reciprocal RVDB database

If using Conda, activate the environment and run with `-profile conda`. If using Pixi, run the workflow through `pixi run ...`. For containerized execution, use `-profile apptainer`.

On Linux, installing the project with Pixi also provides `apptainer` in the Pixi environment, so `pixi run ... -profile apptainer` is the expected containerized entrypoint.

HI-FEVER does not force resume by default. To reuse cached work from a previous run, add Nextflow's `-resume` flag.

Run the HI-FEVER workflow from the root hi-fever folder with the following command (replacing the email address):

`nextflow main.nf --query_file_aa 20_per_fam_no_retro.fasta --ftp_file one_genome.ftp --email john.smith@email.com -profile conda`

Equivalent Pixi-driven command:

`pixi run nextflow run main.nf --query_file_aa 20_per_fam_no_retro.fasta --ftp_file one_genome.ftp --email john.smith@email.com -profile conda`

To continue from a previous cached run:

`nextflow run main.nf --query_file_aa 20_per_fam_no_retro.fasta --ftp_file one_genome.ftp --email john.smith@email.com -profile conda -resume`

### Assembly input modes

HI-FEVER supports two host assembly input modes:

* FTP mode (default): provide `--ftp_file` with one assembly FTP directory per line.
* Local mode: provide local assembly FASTA file(s) with `--assembly_mode local --assembly_file "<glob-or-file>"`.

Local mode example:

`nextflow main.nf --assembly_mode local --assembly_file "assemblies/*.fna.gz" --query_file_aa 20_per_fam_no_retro.fasta --custom_reciprocal --custom_reciprocal_db custom.dmnd --email john.smith@email.com -profile conda`

Optional in local mode:
* `--assembly_metadata_file` for a tab-separated file with two columns: `hostName`, `assembly_id`.
* `--allow_missing_taxonomy true` only if you explicitly want the workflow to continue without taxonomy lookup.

Assembly IDs in local mode are derived from the full input filename stem:
* compression and FASTA suffixes such as `.gz`, `.fa`, `.fna`, `.fasta` are removed
* the remaining filename stem is kept as the assembly identifier in downstream outputs

If taxonomy/metadata files are not available in local mode:
* Core locus discovery, reciprocal search, and genewise reconstruction still run.
* Host metadata fields in summary outputs default to `unknown_host` unless metadata is provided.
* Taxonomy-based annotations/classification may be reduced (more `uncertain` classifications), but this degraded mode now requires `--allow_missing_taxonomy true`.

This will generate a folder called `output` with two subfolders: `accessory_fastas` and `sql`. These outputs are detailed on our [Usage page](https://github.com/PaleovirologyLab/hi-fever/wiki/Usage). For a guide on how to interpret these results, see our [Interpreting results page](https://github.com/PaleovirologyLab/hi-fever/wiki/Interpreting-results)

## Acknowledgements
HI-FEVER is based on the following libraries and programs directory along with their license:
- Biopython (https://biopython.org/)
- Seqtk (https://github.com/lh3/seqtk)
- DIAMOND (https://github.com/bbuchfink/diamond)
- BBmap (https://github.com/BioInfoTools/BBMap)
- BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- Entrez (https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html)
- MMSeqs2 (https://github.com/soedinglab/MMseqs2)
- Nextflow (https://www.nextflow.io/)
- Python3 (https://www.python.org/)
- Wise2 (https://www.ebi.ac.uk/~birney/wise2/)
- Seqkit (https://bioinf.shenwei.me/seqkit/)
- Bedtools (https://bedtools.readthedocs.io/en/latest/index.html)

### Citation
Please include the following citation when using HI-FEVER in your projects.

Laura Muñoz-Baena, Emma F Harding, Jose Gabriel Nino Barreat, Cormac M Kinsella, Aris Katzourakis. HI-FEVER: a Nextflow pipeline for the high-throughput discovery and annotation of endogenous viral elements. Bioinformatics.
DOI: [https://doi.org/10.1093/bioinformatics/btaf610](https://doi.org/10.1093/bioinformatics/btaf610)
