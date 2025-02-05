# hi-fever

## **Hi**gh-throughput next**f**low **EVE** **r**ecovery
`hi-fever` is a Nextflow workflow for finding endogenous viral elements (EVEs) in host genomes. It aims to address common issues in paleovirology including cross-matches between host proteins and EVEs, computational burden of EVE searches and incompatability between software packages or platforms. We provide HI-FEVER as an accessible and informative workflow for any EVE-discovery project.

+ [Installation](#installation)
  - [Conda installation](#conda-installation)
  - [Docker image](#docker-image)
+ [Usage](#usage)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Parameters](#parameters)
+ [Examples](#examples)
  - [Input examples](#input-examples)
    - [Query proteins](#query-proteins)
    - [Host genomes](#host-genomes)
    - [Reciprocal database](#reciprocal-database)
  - [Example runs](#example-runs)
+ [Interpreting results](#interpreting-results)
+ [Features](#features)
+ [Features](#features)

## Installation

First, clone the `hi-fever` repository from GitHub:

```
gh repo clone Paleovirology/hi-fever
```

### Conda installation

Once you have cloned the repository, the simplest way to run `hi-fever` is to create a conda environment using the `environment.yml` file to install of the dependencies including Nextflow.

**Requirements:** If conda is not installed in your work station, install it from [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).

1. Create the `hi-fever` environment:
```
conda env create -f environment.yml
```

2. Activate the environment: 
```
conda activate hi-fever
```

3. Run the workflow:
```
nextflow main.nf
```

**Optional: Conda via job-scheduler**

On HPC clusters a job scheduler such as SLURM is usually installed.
To run the `Conda` environment on compute nodes, a template SLURM script is provided (`cluster-hi-fever.slurm`). To submit on a cluster:

```
sbatch cluster-hi-fever.slurm
```

### Docker image 
When working on MAC computers or a cloud environment you might want to create a singularity image. If working on a MAC with the arm64 architecture (M1 chip, 2020 onwards) please follow the instructions in [#19](/../../issues/19).

**Requirements:**
- Nextflow: can be installed via [this link](https://www.nextflow.io/docs/latest/getstarted.html). It should also be added to your `$PATH`.
- Docker Engine: can be installed, for ubuntu from [here](https://docs.docker.com/engine/install/ubuntu). Confirm that the the installation is working running `sudo docker run hello-world`.

1. Configure Docker to run as a non-root user:
```
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

2. Build and check the hi-fever Docker image:
```
cd docker
docker build -t hi-fever .
docker images
cd ..
```

3. Run the workflow:
```
nextflow main.nf -with-docker hi-fever
```

## Usage

### Inputs 
- **Query proteins:** The viral proteins that you want to search for in the host genomes. You can provide it as a `.fasta` file or diamond database (`.dmnd`) built from your query proteins. The workflow will automatically detect the type of file you provide.
- **Host genome(s):** The genomes you would like to screen. Provide a custom `ftp` file in plain text (example at `data/ftp_list.txt`)
- **E-mail (`--email`):** users need to provide an email in that is used to fetch the NCBI-database with `Efetch` in order to get taxonomical information about the protein hits and the host genomes.
- **Proteins for reciprocal search:** `hi-fever` performs a second DIAMOND search to expand on the annotation of the EVE candidates. A default reciprocal database is provided at `data/minimal_reciprocal.dmnd`.
- **Output directory `--outdir`:** Provide a name for folder where results will be stored. Default to `hi-fever/output`

### Outputs 
HI-FEVER creates 2 folders in the `--outputdir`:

**Sequences (`fastas` folder):**


- `loci-context-coordinates.fasta.gz` nucleotide sequences of the candidate EVEs including the genomic context (flanking regions etc.).
- `loci-merged-coordinates.fasta.gz` nucleotide sequences of the candidate EVEs.

> Note: Sequences in this folder are candidate EVEs. We recommend a further filtering step for these sequences since many EVE candidates are likely to be cross-matches to host proteins. You can use the annotation results in the SQL folder to identify the `locus-id` of interesting candidates. Use the `locus-ids` accession to retrieve target sequences for downstream analysis. 

**Tables and metadata (`sql` folder):**

- `assembly_metadata.tsv` information about the genome assemblies provided in the ftp file including taxonomy, submitter, assembly level etc.
- `assembly_stats.tsv` statistics about the genome assemblies provided in the ftp file including size and coverage
- `genewise.tsv` predicted reconstructed sequences of each EVE candidate
- `locus_assembly_map.tsv` file mapping the genomic locus ID to the assembly from which it came
- `matches.dmnd.annot.tsv` results of the initial forward DIAMOND search of query proteins against genome assemblies
- `predicted_ORFs.tsv` predicted ORFs within the region of each candidate EVE
- `reciprocal-*-matches.dmnd.tsv` results of the reciprocal DIAMOND searches of candidate EVEs against reciprocal databases. There will be one result file for each reciprocal database search.
- `taxonomy_table.tsv` taxonomy data for every genome assembly and reciprocal hit from the HI-FEVER run

### Parameters
- `--query_file_aa data/circoviridae.fa`: Custom protein query file in fasta format (default: protein_query.fasta)
- `--ftp_file data/assemblies.txt`: Custom ftp file in plain text (default: ftp_list.txt)
- `--reciprocal_db file custom.dmnd`: File with proteins or database used during the reciprocal search (default: minimal_reciprocal.dmnd)
- `--outdir herpesviridae_vs_tarsier`: Name of output directory (default: hifever_output)
- `--mmseqs_minseqid 0.70`: Sequence identity threshold for clustering of the protein query (default: 0.95 = 95%). Input protein queries are clustered prior to the search to avoid redundant searches of highly identical proteins
- `--mmseqs_cover 0.80`: Minimum percentage of cluster member sequence length that must overlap with the representative sequence (default: 0.90 = 90%)
- `chunk_size 10000`: Size of host genome to process per DIAMOND forward search (default: 50000). Reducing this number will return more hits but take longer
- `--diamond_forks 8`: DIAMOND fork count. By default Nextflow attempts to run all available forward DIAMOND tasks in parallel (one for each currently downloaded assembly), which can lead to overuse of memory resources and job termination. On local machines and clusters, it is therefore suggested to limit the number of parallel DIAMOND tasks (i.e., "forks") allowed at once. For setting the value, we recommend total cores / 12. Note that this does not affect the reciprocal DIAMOND search, which uses all available CPUs in a single process (default: 4 for cluster workflow)
- `--diamond_mode ultra-sensitive`: [DIAMOND sensitivity](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes) (default: very-sensitive)
- `--diamond_matrix BLOSUM45`: [DIAMOND substitution matrix](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#alignment-options) (default: BLOSUM62)
- `--interval 500`: Maximum interval length allowed between features in the context FASTA (default: 1000). This parameter modifies the allowed distance between BLAST hits before they are merged and counted as one hit
- `--flank 500`: Maximum length of (available) flanking sequence to add upstream and downstream of detected features in the context FASTA (default: 3000). This parameter allows you to modify how much host sequence upstream and downstream of the BLAST hits you would like to extract and keep for downstream analysis.
- `--orf_size_nt 300`: Minimum length of extracted ORFs, in nucleotides (default: 150). This parameter allows you to alter the size threshold used when searching for ORFs within EVE hits
- `--genewise_matrix BLOSUM45`: Genewise substitution matrix (default: BLOSUM62). Options: BLOSUM80, BLOSUM62, BLOSUM45, BLOSUM30
- `--stop_task soft-mask`: Modify in-frame STOP codons in the Genewise coding DNA sequence output (default: remove). Options: remove [delete in-frame STOPs from the coding sequence], soft-mask [convert in-frame STOPs to lowercase
- `--pairs_per_task 100`: Subset ($n$) of total loci count ($N$) processed by each Genewise task. To ensure equivalent workload across tasks and a high level of parallelisation, Genewise operations are split into $T$ tasks, where $T = N/n$ (default: 500)
- `-with-report report.html`: Create Nextflow html workflow report (includes run time, user information, task metadata, and CPU, memory, and I/O usage)


## Examples
### Input examples

HI-FEVER searches host genomes for matches against protein queries. As such, it requires a file providing the ftp list of genomes to be screened and a fasta file of protein queries. Examples of these files are below:

#### Query proteins
A FASTA file containing protein sequences (default 'protein_query.fasta') e.g:

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


#### Host genomes
The host genomes file (`assemblies.ftps`) should be a list of links to genome assemblies:

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/330/505/GCA_000330505.1_EIA2_v2
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/065/795/GCA_023065795.1_ASM2306579v1
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/208/925/GCF_000208925.1_JCVI_ESG2_1.0

```

>A text file containing ftp links for assemblies to process (default: 'ftp_list.txt'). Links to assemblies available from NCBI can be found on the [RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/refseq) and [GenBank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank) ftp sites. For example `refseq/assembly_summary_refseq.txt` or `refseq/protozoa/assembly_summary.txt`.

#### Reciprocal database

An essential step in `hi-fever` is the reciprocal DIAMOND search. There are three options of input databases that you can use for this step based on your research question and your computational resources:

1. **Minimal reciprocal database** (recommended for most users, requires ~3GB of storage space)

By default, `hi-fever` includes a reciprocal database (`hi-fever/data/minimal_reciprocal.dmnd`) with the most frequent and informative hits you can get while searching for EVEs in vertebrate genomes. We provide it as an alternative to the full NCBI non reduntand protein database (`nr`) and Reference Virus Database (`RVDB`). which normally require >100GB of storage space. Whilst we have made all efforts to preserve a diversity of informative proteins in the minimal reciprocal database, there is a chance that the closest EVE match will not be represented.


2. **Full reciprocal databases** (recommended for large workstations and cluster computers, requires ~115 GB of storage space)

Download the latest version of both the NCBI nr database and RVDB proteins database for the most accurate EVE annotations using:

```
cd data
# NCBI nr database
wget https://files.osf.io/v1/resources/tejwd/providers/googledrive/nr_rep_seq.fasta.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
gunzip prot.accession2taxid.FULL.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
diamond makedb --in nr_rep_seq.fasta.gz -d nr_clustered_wtaxa --taxonmap prot.accession2taxid.FULL --taxonnodes nodes.dmp --taxonnames names.dmp --threads 16

# RVDB
wget https://rvdb-prot.pasteur.fr/files/U-RVDBv26.0-prot.fasta.xz
unxz U-RVDBv26.0-prot.fasta.xz
mmseqs createdb U-RVDBv26.0-prot.fasta DB
mmseqs linclust --min-seq-id 0.98 --cov-mode 1 -c 0.90 DB DB_clu tmp
mmseqs createsubdb DB_clu DB DB_clu_rep
mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
mv DB_clu_rep.fasta U-RVDBv26.0-prot-clustered-minid0.98-cov0.90.fasta
awk '{if($0~">") {split($0,a,"|"); print ">"a[3],substr($0,2,length($0))} else print $0}' U-RVDBv26.0-prot-clustered-minid0.98-cov0.90.fasta > U-RVDBv26.0-prot-clustered-minid0.98-cov0.90-relabelled.fasta
diamond makedb --in U-RVDBv26.0-prot-clustered-minid0.98-cov0.90-relabelled.fasta -d rvdbv26_clustered_wtaxa.dmnd --taxonmap prot.accession2taxid.FULL --taxonnodes nodes.dmp --taxonnames names.dmp
```

> **Note:** 
    A preclustered of nr version is [made available by Arcadia Science](https://github.com/Arcadia-Science/2023-nr-clustering).
    RVDB was developed by Arifa Khan's group at CBER. For hi-fever, RVDB-prot is required, a protein version maintained at the Institut Pasteur [with archives available here](https://rvdb-prot.pasteur.fr).
    The commands above offer a template for downloading and formatting these databases, though it is recommended to edit them to install the most recent database releases. These commands require the diamond and mmseqs, both of which are provided in the conda environment distributed with HI-FEVER.


3. **Custom reciprocal database** (recomended only for users with previous knowledge about their EVEs): 
`hi-fever`builts a reciprocal database a user-provided fasta file or a prebuilt diamond database (specified under the `--reciprocal_db` parameter). This is not recommended unless you already have a lot of information about the EVEs that you are searching for, as it will not easily distinguish between true EVEs and host cross-matches.

### Example runs
Some example sets of parameters are shown below as a guide to how to customise your HI-FEVER run.

A default run
```
nextflow main.nf --query_file_aa viruses.fasta --ftp_file genomes.txt --outdir hi_fever_results
```


Using the full NCBI and nr databases for the reciprocal DIAMOND search. Recommended where possible but is resource intensive and will take considerably longer.
```
nextflow main.nf --query_file_aa viruses.fasta --ftp_file genomes.txt --outdir hi_fever_results --full_reciprocal TRUE
```


Providing a custom database (in fasta format) for the reciprocal DIAMOND search. Suitable if you have a shortlist of proteins to confirm your EVEs, however may lose information on host cross-matches. 
```
nextflow main.nf --query_file_aa viruses.fasta --ftp_file genomes.txt --outdir hi_fever_results --reciprocal_db reciprocal_proteins.fasta
```


Providing a custom database (in dmnd format) for the reciprocal DIAMOND search. Suitable if you have a shortlist of proteins to confirm your EVEs and have prebuilt a DIAMOND database from them.
```
nextflow main.nf --query_file_aa viruses.fasta --ftp_file genomes.txt --outdir hi_fever_results --dont_build_reciprocal TRUE --reciprocal_db reciprocal_proteins.fasta
```


Clustering the query proteins at low identity prior to DIAMOND searches. Suitable if you have many similar query proteins and are looking for more general matches.
```
nextflow main.nf --query_file_aa viruses.fasta --ftp_file genomes.txt --outdir hi_fever_results --mmseqs_minseqid 0.70 mmseqs_cover 0.60
```

Customising the diamond search modes. Suitable if you want to run a less computationally-intensive search.
```
nextflow main.nf --query_file_aa viruses.fasta --ftp_file genomes.txt --outdir hi_fever_results --diamond-forks 1 --diamond_mode sensitive
```

Returning longer candidate EVEs by merging more distant neigboring hits and returning longer flanking sequences. Suitable when searching for multiple protein integrations from one virus eg. proviruses.
```
nextflow main.nf --query_file_aa viruses.fasta --ftp_file genomes.txt --outdir hi_fever_results --interval 3000 --flank 5000
```

## Interpreting results

For a quick interpretation of the results, we recommend focussing on the `reciprocal-*-matches.dmnd.tsv` files for full EVE annotations. Examples of how to differentiate between EVEs and false positives is described in our publication (paper).

We recommend parsing these tables into an SQL database for quick querying of the results. We have provided two SQL schema scripts `hi-fever-db-PostgreSQL-schema.sql` and `merge_tables.sql` to automatically import and process these files into searchable tables.

Steps to run:

1. Install PostgreSQL.
2. Create a local server (e.g., name = hi-fever-server, host = localhost, user = postgres).
3. Create a new database within the server (e.g., hi-fever-db).
4. Ensure hi-fever output data tables (by default found in output/sql) are located in a server accessible location on your system (& edit the import section of the script to specify this).
5. Open the SQL query interpreter for the hi-fever-db database (e.g., by right clicking the db name and selecting "Query tool"), then paste in the hi-fever-db-PostgreSQL-schema.sql script and run it. This will generate the table schema and its constraints, and import the data tables.
6. Paste the merge_tables.sql script into the query interpreter and run it. This will create a new table with all EVE metadata in one place and only display the top reciprocal hit for each EVE. This table is recommended as a starting point for interpreting results.

Alternatively, if you do not want to use SQL the files can be parsed in many other ways including:
 - grep searches of the reciprocal-*-matches.dmnd.tsv files to identify EVEs matching viral families of interest
```
grep Bornaviridae reciprocal-nr-matches.dmnd.tsv
grep nucleoprotein reciprocal-nr-matches.dmnd.tsv
```

- importing the tables in R

- opening the files in Excel or a similar spreadsheet manager (not recommended for large output files) and filtering using inbuilt functions


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
