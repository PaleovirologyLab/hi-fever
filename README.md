# hi-fever

> **Hi**gh-throughput next**f**low **EVE** **r**ecovery

## About

Hi-fever is a Nextflow workflow for finding endogenous viral elements (EVEs) in host genomes. Some features:

- Protein-to-DNA based search
- Scales from laptop to cluster or cloud
- Designed to function with 1000's of input assemblies

Outputs include:

- Genomic coordinates of candidate EVEs
- Best forward and reciprocal hits, with integrated taxonomy
- Predicted EVE protein sequences and cDNA (frameshift and premature STOP codon aware)
- Extraction of flanking sequence
- Open reading frame detection and extension beyond original hit
- Assembly information
- Predicted EVE protein domains

## To install and run with presets

Clone the repo, e.g., with GitHub CLI:

```
gh repo clone Paleovirology/hi-fever
```

To run with presets, ensure the following are in your working directory:

>`protein_query.fasta` [[more information]](#protein-queries)

>`ftp_list.txt` [[more information]](#assembly-list)

>`domains` pHMM directory [[more information]](#phmm-library)

>`nr_clustered_wtaxa.dmnd` NCBI nr proteins database [[more information]](#ncbi-nr-proteins-db)

### Option 1: Run from a Conda environment (e.g., for working locally or on a HPC cluster)

A Conda `environment.yml` file is provided for installation of all dependencies including Nextflow.
[Get Conda here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).

Create and activate the hi-fever environment:

```
conda env create -f environment.yml
conda activate hi-fever
```

Run the workflow:

```
nextflow hi-fever.nf
```

### Option 2: Run from a Docker or Singularity image

[Install Nextflow and add to $PATH](https://www.nextflow.io/docs/latest/getstarted.html).

#### Docker (e.g., for working on a cloud environment)

Install Docker Engine, e.g., [for Ubuntu](https://docs.docker.com/engine/install/ubuntu), and run `sudo docker run hello-world` to confirm the installation is working.

Configure Docker to run as a non-root user:

```
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

Build and check the hi-fever Docker image:

```
cd docker
docker build -t hi-fever .
docker images
cd ..
```

Run the workflow:

```
nextflow hi-fever.nf -with-docker hi-fever
```

#### Singularity (e.g., for working on a HPC cluster)

Generate the Docker image as above.

Export the image as a tarball (for this you will need to lookup the image ID, which can be found with `docker images`):

```
docker save IMAGE_ID -o hi-fever.tar
```

Copy the tarball to a Cluster with Singularity available, or [install it](https://apptainer.org/docs/admin/main/installation.html) locally. Build the Singularity image:

```
singularity build singularity-hi-fever docker-archive:./hi-fever.tar
```

Run the workflow:

```
nextflow hi-fever.nf -with-singularity singularity-hi-fever
```

## Optional parameters and example inputs

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

Modify in-frame STOP codons in the Genewise coding DNA sequence output (default: remove). Options: remove [delete in-frame STOPs from the coding sequence], soft-mask [convert in-frame STOPs to lowercase].

- `--stop_task soft-mask`

Create Nextflow html workflow report (includes run time, user information, task metadata, and CPU, memory, and I/O usage).

- `-with-report report.html`

## Required input files

### Protein queries

- A FASTA file containing protein sequences, e.g:

```
>ADI48253.1 putative Rep [Circoviridae TM-6c]
MQSVNWCFTLNNYTNEDVNKLKQVKCRYICLGFEVGDKKQTPHIQGFIQFEKKVRLSVWKKINKKIHAEI
MKGTIEQAINYCKKSGTFEERGEIIKMGERRDLKEAKKKCAEVGLRAITDCESTYNLQVIRNCQIMLEYH
EKERDFKPEVIWIYGESGAGKTKYISEKCAEVDTYWKDATKWWNGYDRHEITVMDDFRASNMKMNELLKL
IDRYPHRVEIKGGFRQMLSKKIYISSIMHPKDVYNLPEEPVKQLLRRIDTIIKI

>NP_042987.1 Rep [Human betaherpesvirus 6A]
MFSIINPSDDFWTKDKYIMLTIKGPMEWEAEIPGISTDFFCKFSNVSVPHFRDMHSPGAPDIKWITACTK
MIDVILNYWNNKTAVPTPAKWYAQAENKAGRPSLILLIALDGIPSATIGKHTTEIRGVLIKDFFDGNAPK
IDDWCTYAKTKKNGGGTQVFSLSYIPFALLQIIRPQFQWAWTNINELGDVCDEIHRKHIISHFNKKPNVK
LMLFPKDGINGISLKSKFLGTIEWLSDLGIVTEDAWIRRDIRSYMQLLTLTHGDVLIHRALSIAKKRIRA
TRKAIDFIAHIDTDFQIYENPVYQLFCLQSFDPILAGTILYQWLSHRGGKKNTVSFIGPPGCGKSMLTGA
ILENIPLHGILHGSLNTKNLRAYGQVLVLWWKDISINFDNFNIIKSLLGGQKIIFPINENDHVQIGPCPI
IATSCVDIRSMVHSNLHKINLSQRVYNFTFDKVIPRNFPVIQKDDINQFLFWARNRSINCFIDYTVPKIL
```

### Assembly list

- A text file containing ftp links for assemblies to process, e.g.:

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/330/505/GCA_000330505.1_EIA2_v2
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/065/795/GCA_023065795.1_ASM2306579v1
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/208/925/GCF_000208925.1_JCVI_ESG2_1.0
```

- Links to assemblies available from NCBI are on the [RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/refseq) and [GenBank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank) ftp sites, e.g. `refseq/assembly_summary_refseq.txt` or `refseq/protozoa/assembly_summary.txt`.

### pHMM library

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

### NCBI nr proteins db

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