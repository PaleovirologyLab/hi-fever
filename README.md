# hi-fever

> **Hi**gh-throughput next**f**low **EVE** **r**ecovery

# About

# Installation

Clone the repo, e.g., with GitHub CLI:

```
gh repo clone Paleovirology/hi-fever
```

For installation of software dependencies, a conda `environment.yml` file is provided.
[Get conda here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
Create, activate, and check software environment:

```
conda env create -f environment.yml
conda activate hi-fever
conda list
```

# Usage

## With presets

By default the following must be in your working directory:

>`protein_query.fasta`

>`ftp_list.txt` [[more information]](#assembly-list)

>`domains` pHMM directory [[more information]](#phmm-library)

>`nr_clustered_wtaxa.dmnd` NCBI nr proteins database [[more information]](#ncbi-nr-proteins-db)

To run the workflow:

`nextflow hi-fever.nf`

## Optional parameters with example inputs

Custom protein query file (default: protein_query.fasta).

- `--query_file_aa circoviridae.fa`

Custom ftp file (default: ftp_list.txt).

- `--ftp_file assemblies.txt`

Location of custom pHMM library for query domain annotation (default: domains).

- `--phmms Pfam-32`

Custom reciprocal BLASTp database (DIAMOND formatted, default: nr_clustered_wtaxa.dmnd).

- `--reciprocal_db nr.dmnd`

Sequence identity threshold for clustering of the protein query (default: 0.95 = 95%).

- `--mmseqs_minseqid 0.70`

Minimum percentage of cluster member sequence length that must overlap with the representative sequence (default: 0.90 = 90%).

- `--mmseqs_cover 0.80`

[DIAMOND sensitivity](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes) (default: very-sensitive).

- `--diamond_mode ultra-sensitive`

[DIAMOND substitution matrix](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#alignment-options) (default: BLOSUM62).

- `--diamond_matrix BLOSUM45`

CPUs per DIAMOND fork (default: 12). By default the forward DIAMOND search will be "forked" to run four assemblies in parallel, and each task will be allocated this number of CPUs. To set, recommend dividing the available CPUs by four. Note that this does not affect the reciprocal DIAMOND search, which uses all available CPUs in a single process.

- `--diamond_cpus 6`

Maximum interval length allowed between features in the context FASTA (default: 1000).

- `--interval 500`

Maximum length of (available) flanking sequence to add upstream and downstream of detected features in the context FASTA (default: 3000).

- `--flank 500`

Minimum length of extracted ORFs, in nucleotides (default: 150).

- `--orf_size_nt 300`

Genewise substitution matrix (default: BLOSUM62). Options: BLOSUM80, BLOSUM62, BLOSUM45, BLOSUM30.

- `--genewise_matrix BLOSUM45`

Modify in-frame STOP codons in the Genewise coding DNA sequence output (default: remove). Options: remove [delete in-frame STOPs from the coding sequence], convert [convert in-frame STOPs to lowercase].

- `--stop_task convert`

Create Nextflow html workflow report (includes run time, user information, task metadata, and CPU, memory, and I/O usage).

- `-with-report report.html`

# Input files

## Assembly list

- A text file containing ftp links for assemblies to process, e.g.:

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/330/505/GCA_000330505.1_EIA2_v2
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/065/795/GCA_023065795.1_ASM2306579v1
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/208/925/GCF_000208925.1_JCVI_ESG2_1.0
```

- Links to assemblies available from NCBI are on the [RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/refseq) and [GenBank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank) ftp sites, e.g. `refseq/assembly_summary_refseq.txt` or `refseq/protozoa/assembly_summary.txt`.

## pHMM library

- To generate a pHMM library with the latest version of [Pfam](https://www.ebi.ac.uk/interpro/download/Pfam), run the following:

```
conda activate hi-fever
mkdir domains; cd domains
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
# Replace space characters in model description lines, as these are column delimiters in the output table.
sed -i '/^DESC/ s/ /_/g; s/__/  /' Pfam-A.hmm
hmmpress Pfam-A.hmm
cd ..
```

## NCBI nr proteins db

- For the reciprocal DIAMOND BLASTp stage, a local nr or clustered nr DIAMOND formatted database is recommended.
- A clustered nr version is [made available by Arcadia Science](https://github.com/Arcadia-Science/2023-nr-clustering). To download and format:

```
conda activate hi-fever
wget https://files.osf.io/v1/resources/tejwd/providers/googledrive/nr_rep_seq.fasta.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
gunzip prot.accession2taxid.FULL.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
diamond makedb --in nr_rep_seq.fasta.gz -d nr_clustered_wtaxa --taxonmap prot.accession2taxid.FULL --taxonnodes nodes.dmp --taxonnames names.dmp --threads 16
```
# Docker
To build and use the docker image first make sure docker is installed and running. Type ```docker run hello-world``` inside a terminal.

Then change into the docker folder and run the following:

```
cd docker
docker build -t <docker-account-name>/hi-fever .
docker images
```

You should see the ```hi-fever``` image and an ```ubuntu``` base image.

Now you can run the pipeline without setting up the conda environment by adding the ```-with-docker <docker-image>``` command line option. 
