# HI-FEVER

## **Hi**gh-throughput next**f**low **EVE** **r**ecovery
`hi-fever` is a Nextflow workflow for finding endogenous viral elements (EVEs) in host genomes. It aims to address common issues in paleovirology including cross-matches between host proteins and EVEs, computational burden of EVE searches and incompatability between software packages or platforms. We provide HI-FEVER as an accessible and informative workflow for any EVE-discovery project.

## Installation and usage

HI-FEVER is available for use on LINUX, Windows (WSL) and Mac through Conda and Docker. Full documentation can be found in [the wiki](https://github.com/PaleovirologyLab/hi-fever/wiki).

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

## Acknowledgements
Here we should inclide the libraries, the programs that we use, and the citation to the paper:

hi-fever is based on the following libraries and programs directory along with their license:
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
Cite our work at: Please cite (paper) when using HI-FEVER in your projects.
