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

Set [DIAMOND sensitivity](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes) (default is very-sensitive).

- `--diamond_mode ultra-sensitive`

Set [DIAMOND substitution matrix](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#alignment-options) (default is BLOSUM62).

- `--diamond_matrix BLOSUM45`

Create Nextflow html workflow report (includes run time, user information, task metadata, and CPU, memory, and I/O usage).

- `-with-report report.html`



