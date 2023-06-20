# hi-fever 
> **Hi**gh-throughput next**f**low **EVE** **r**ecovery workflow


# About


# Installation

A conda `environment.yml` file is provided in this repository for fast setup. [Get conda here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).

Create, activate, and check environment:
```
conda env create -f environment.yml
conda activate hi-fever
conda list
```

# Usage

## With presets

To run the workflow:

`nextflow hi-fever.nf `

By default the following two files must be in your working directory: 
>`protein_query.fasta` 

>`ftp_list.txt`

## Adjustable parameters
Custom ftp file name.

- `--ftp_file $PWD/file_name.txt`

Custom protein query file name.
- `--query_file_aa $PWD/circoviridae.fa`

Sequence identity threshold for clustering of the protein query (default: 0.95 = 95%).

- `--mmseqs_minseqid 0.70`

Minimum percentage of cluster member sequence length that must overlap with the representative sequence (default: 0.90 = 90%).

- `--mmseqs_cover 0.80`

[DIAMOND sensitivity](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes) (default: very-sensitive).

- `--diamond_mode ultra-sensitive`

[DIAMOND substitution matrix](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#alignment-options) (default: BLOSUM62).

- `--diamond_matrix BLOSUM45`

CPUs per DIAMOND fork (default: 12). By default DIAMOND will be "forked" to run four assemblies in parallel, and each task will be allocated this number of CPUs. To set, recommend dividing the available CPUs by four.

- `--diamond_cpus 6`

Create Nextflow html workflow report (includes run time, user information, task metadata, and CPU, memory, and I/O usage).

- `-with-report report.html`