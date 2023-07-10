# hi-fever 
> **Hi**gh-throughput next**f**low **EVE** **r**ecovery


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

By default the following must be in your working directory:

>`protein_query.fasta`

>`ftp_list.txt` [[more information]](#assembly-list)

>`domains-vX.XX` pHMM directory [[more information]](#phmm-library)

>`nr_clustered.dmnd` NCBI nr proteins database [[more information]](#ncbi-nr-proteins-db)

To run the workflow:

`nextflow hi-fever.nf `

## Adjustable parameters with example inputs

Custom protein query file (default: protein_query.fasta).
- `--query_file_aa circoviridae.fa`

Custom ftp file (default: ftp_list.txt).

- `--ftp_file assemblies.txt`

Location of custom pHMM library for query domain annotation (default: domains-v*).

- `--phmms Pfam`

Custom reciprocal BLASTx database (DIAMOND formatted, default: nr_clustered.dmnd)

- `--reciprocal_db nr.dmnd`

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

Maximum interval length allowed between features in the context FASTA (default: 1000).

- `--interval 500`

Maximum length of (available) flanking sequence to add upstream and downstream of detected features in the context FASTA (default: 3000).

- `--flank 500`

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

- A preformatted pHMM library including [Pfam](https://www.ebi.ac.uk/interpro/download/Pfam), [RVDB](https://rvdb.dbi.udel.edu), [PHROGs](https://phrogs.lmge.uca.fr), and [more](link_to_description) is [available here](link_to_download). Simply download and unpack it in the hi-fever working directory using:
- `tar -xzf domains-v*.tar.gz`
- To generate a custom library this example with Pfam may be adapted:

```
conda activate hi-fever
mkdir Pfam; cd Pfam
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
# Replace space characters in model description lines, as these are column delimiters in the output table.
sed -i '/^DESC/ s/ /_/g; s/__/  /' Pfam-A.hmm
hmmpress Pfam-A.hmm
cd ..
```

## NCBI nr proteins db

- For the reciprocal DIAMOND BLASTx stage, a local nr or clustered nr DIAMOND formatted database is recommended.
- A clustered nr version is [made available by Arcadia Science](https://github.com/Arcadia-Science/2023-nr-clustering). To download and format:
```
wget https://files.osf.io/v1/resources/tejwd/providers/googledrive/nr_rep_seq.fasta.gz
diamond makedb --in nr_rep_seq.fasta.gz -d nr_clustered.dmnd --threads 8
rm nr_rep_seq.fasta.gz

```
