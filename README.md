# HI-FEVER

## **Hi**gh-throughput next**f**low **EVE** **r**ecovery
`hi-fever` is a Nextflow workflow for finding endogenous viral elements (EVEs) in host genomes. It aims to address common issues in paleovirology including cross-matches between host proteins and EVEs, computational burden of EVE searches and incompatability between software packages or platforms. We provide HI-FEVER as an accessible and informative workflow for any EVE-discovery project.

## Features

- Protein-to-DNA based search allows detection of divergent and ancient EVEs
- Designed to function with millions of input query proteins
- Reconstructs the predicted EVE protein based on its closest modern match
- Harnesses parallelisation to optimise compute resources
- Scales from laptop to cluster
- Conda and Docker compatible
- LINUX, Windows and MAC compatible

HI-FEVER provides a variety of output information about candidate EVEs, suited to many downstream purposes. Outputs include:
- Genomic coordinates of candidate EVEs
- Closest matches in the reciprocal databases, including full taxonomical information
- Predicted EVE protein sequences and cDNA (frameshift and premature STOP codon aware), with extension beyond original hit
- Extracted nucleotide sequence of each candidate EVE and flanking host genome sequence
- Metadata & statistics of the genome assemblies screened

## Installation and usage

HI-FEVER is available for use on LINUX, Windows (WSL) and Mac through Conda and Docker. Full documentation can be found in [the wiki](https://github.com/PaleovirologyLab/hi-fever/wiki).

## Test run
To experiment with and explore HI-FEVER options we provide instructions on running a test dataset below. All data used for this test are available on our Open Science Framework repository [here](https://osf.io/y357r/) in the test_run folder.

**Preparation**

Ensure the required files are in the hi-fever/data folder:
* `20_per_fam_no_retro.fasta` protein query file
* `vertebrate_genomes_n5.txt` list of five vertebrate genome ftps
* `taxdump.tar.gz` taxonomy map file
* `MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd.tar.xz`: the minimal database built from the NCBI non-redundant database
* `MINI_rvdbv28_wtaxa.dmnd.tar.xz`: the minimal database built from the RVDB database

Unzip the reciprocal databases with the following tar commands:
```
tar -xf MINI-nr_rep_seq-clustered_70id_80c_wtaxa.dmnd.tar.xz
tar -xf MINI_rvdbv28_wtaxa.dmnd.tar.xz
```

If using conda, activate the environment. If using Docker on Mac (arm64), open a terminal tab within Docker desktop. If using Docker on LINUX add the -with_docker flag to the run command below.

Run the HI-FEVER workflow from the root hi-fever folder with the following command (replacing the email address):

`nextflow main.nf --query_file_aa data/20_per_fam_no_retro.fasta --ftp_file data/vertebrate_genomes_n5.txt --email john.smith@email.com`

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
Please cite (paper) when using HI-FEVER in your projects.
