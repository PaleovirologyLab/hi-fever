# hi-fever

> **Hi**gh-throughput next**f**low **EVE** **r**ecovery

## ABOUT

Hi-fever is a Nextflow workflow for finding endogenous viral elements (EVEs) in host genomes. It aims to address common issues in paleovirology including cross-matches between host proteins and EVEs, computational burden of EVE searches and incompatability between software packages or platforms. We provide HI-FEVER as an accessible and informative workflow for any EVE-discovery project.

Some features:

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
- Open reading frame detection, with extension beyond original hit
- Flanking sequences of each candidate EVE
- Metadata & statistics of the genome assemblies screened

Please cite (paper) when using HI-FEVER in your projects.

## INSTALLATION

HI-FEVER was originally designed for use on workstations and cluster computers. However, we also provide a setup tailored for smaller computers where reduced databases are used to minimise storage space required. The workflow and algorithms are identical between them, however the smaller setup uses a reduced database for the reciprocal search, minimising the storage space required at the cost of limiting the search space for EVE-relatives.

### To install the full version (requires ~X GB of storage space)

Clone the repo, e.g., with GitHub CLI:

```
gh repo clone Paleovirology/hi-fever
```

Download and prepare the databases used for the reciprocal searches. It is recommended to download the latest version of both the NCBI nr database and RVDB proteins database for the most accurate EVE annotations.
- A preclustered nr version is [made available by Arcadia Science](https://github.com/Arcadia-Science/2023-nr-clustering).
- RVDB was developed by Arifa Khan's group at CBER. For hi-fever, RVDB-prot is required, a protein version maintained at the Institut Pasteur [with archives available here](https://rvdb-prot.pasteur.fr).
- The commands below offer a template for downloading and formatting these databases, though it is recommended to edit them to install the most recent database releases. These commands require the diamond and mmseqs, both of which are provided in the conda environment distributed with HI-FEVER.

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


### To install the reduced version (for smaller computers, requires ~3GB of storage space)

Clone the repo, e.g., with GitHub CLI:

```
gh repo clone Paleovirology/hi-fever
```

Download and prepare the minimal reciprocal database. This is provided as an alternative to the full NCBI nr and RVDB databases, which together require X GB of storage space. Whilst we have made all efforts to preserve a diversity of informative proteins in the minimal reciprocal database it is recommended to use the full NCBI nr and RVDB if possible.

>`minimal_reciprocal.dmnd` gg [[more information]](#reciprocal-databases-ncbi-nr--rvdb)


## RUNNING

### Prepare input files
HI-FEVER searches host genomes for matches against protein queries. As such, it requires a file providing the ftp list of genomes to be screened and a fasta file of protein queries. Examples of these files are below:

- A FASTA file containing protein sequences (default 'protein_query.fasta') e.g:

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

- A text file containing ftp links for assemblies to process (default: 'ftp_list.txt') e.g:

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/330/505/GCA_000330505.1_EIA2_v2
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/065/795/GCA_023065795.1_ASM2306579v1
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/208/925/GCF_000208925.1_JCVI_ESG2_1.0
```

- Links to assemblies available from NCBI are on the [RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/refseq) and [GenBank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank) ftp sites, e.g. `refseq/assembly_summary_refseq.txt` or `refseq/protozoa/assembly_summary.txt`.

Before running HI-FEVER both of these files need to be in the 'data' directory.


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
nextflow main.nf
```

#### Running with Conda via a job scheduler, e.g., SLURM

On HPC clusters a job scheduler such as SLURM is usually installed.

To run the Conda environment on compute nodes, a template SLURM script is provided (`cluster-hi-fever.slurm`). To submit on a cluster:

```
sbatch cluster-hi-fever.slurm
```

### Option 2: Run from a Docker image (e.g., for working on MAC computers or a cloud environment)

[Install Nextflow and add to $PATH](https://www.nextflow.io/docs/latest/getstarted.html).

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
nextflow main.nf -with-docker hi-fever
```

## Optional parameters and example inputs

Workflow to run, LOCAL: for single nodes on clusters or cloud, BATCH: for cloud batch services (default: LOCAL).

- `--entry BATCH`

Custom protein query file in fasta format (default: protein_query.fasta).

- `--query_file_aa data/circoviridae.fa`

Custom ftp file in plain text (default: ftp_list.txt).

- `--ftp_file data/assemblies.txt`

Name of output directory (default: output).

- `--outdir herpesviridae_vs_tarsier`

Sequence identity threshold for clustering of the protein query (default: 0.95 = 95%). Input protein queries are clustered prior to the search to avoid redundant searches of highly identical proteins.

- `--mmseqs_minseqid 0.70`

Minimum percentage of cluster member sequence length that must overlap with the representative sequence (default: 0.90 = 90%).

- `--mmseqs_cover 0.80`

DIAMOND fork count. By default Nextflow attempts to run all available forward DIAMOND tasks in parallel (one for each currently downloaded assembly), which can lead to overuse of memory resources and job termination. On local machines and clusters, it is therefore suggested to limit the number of parallel DIAMOND tasks (i.e., "forks") allowed at once. For setting the value, we recommend total cores / 12. Note that this does not affect the reciprocal DIAMOND search, which uses all available CPUs in a single process (default: 4 for cluster workflow, and not set on cloud workflow as instances allow parallel processing).

- `--diamond_forks 8`

[DIAMOND sensitivity](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes) (default: very-sensitive).

- `--diamond_mode ultra-sensitive`

[DIAMOND substitution matrix](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#alignment-options) (default: BLOSUM62).

- `--diamond_matrix BLOSUM45`

Maximum interval length allowed between features in the context FASTA (default: 1000). This parameter modifies the allowed distance between BLAST hits before they are merged and counted as one hit.

- `--interval 500`

Maximum length of (available) flanking sequence to add upstream and downstream of detected features in the context FASTA (default: 3000). This parameter allows you to modify how much host sequence upstream and downstream of the BLAST hits you would like to extract and keep for downstream analysis.

- `--flank 500`

Minimum length of extracted ORFs, in nucleotides (default: 150). This parameter allows you to alter the size threshold used when searching for ORFs within EVE hits.

- `--orf_size_nt 300`

Genewise substitution matrix (default: BLOSUM62). Options: BLOSUM80, BLOSUM62, BLOSUM45, BLOSUM30.

- `--genewise_matrix BLOSUM45`

Modify in-frame STOP codons in the Genewise coding DNA sequence output (default: remove). Options: remove [delete in-frame STOPs from the coding sequence], soft-mask [convert in-frame STOPs to lowercase].

- `--stop_task soft-mask`

Subset ($n$) of total loci count ($N$) processed by each Genewise task. To ensure equivalent workload across tasks and a high level of parallelisation, Genewise operations are split into $T$ tasks, where $T = N/n$ (default: 500).

- `--pairs_per_task 100`

Create Nextflow html workflow report (includes run time, user information, task metadata, and CPU, memory, and I/O usage).

- `-with-report report.html`


## Interpreting results
HI-FEVER outputs several files in 2 results folders: accessory_fastas and SQL. A detailed description of each output is below:

#### Accessory fastas folder

This folder contains the sequence data of candidate EVEs.
 - `loci-context-coordinates.fasta.gz` nucleotide sequences of the candidate EVEs including the genomic context (flanking regions etc.).
 - `loci-merged-coordinates.fasta.gz` nucleotide sequences of the candidate EVEs.
These files are not recommended for direct analysis without further filtering. They contain the results for every candidate EVE, many of which are likely to be cross-matches to host proteins. We recommend to use the annotation results in the SQL folder to identify EVEs of interest and extract them from these fasta files before further analysis.

#### SQL folder
This folder contains tables and metadata relating to the candidate EVEs.
 - `assembly_metadata.tsv` information about the genome assemblies provided in the ftp file including taxonomy, submitter, assembly level etc.
 - `assembly_stats.tsv` statistics about the genome assemblies provided in the ftp file including size and coverage
 - `genewise.tsv` predicted reconstructed sequences of each EVE candidate
 - `locus_assembly_map.tsv` file mapping the genomic locus ID to the assembly from which it came
 - `matches.dmnd.annot.tsv` results of the initial forward DIAMOND search of query proteins against genome assemblies
 - `predicted_ORFs.tsv` predicted ORFs within the region of each candidate EVE
 - `reciprocal-*-matches.dmnd.tsv` results of the reciprocal DIAMOND searches of candidate EVEs against reciprocal databases. There will be one result file for each reciprocal database search.
 - `taxonomy_table.tsv` taxonomy data for every genome assembly and reciprocal hit from the HI-FEVER run

For a quick interpretation of the results, we recommend focussing on the `reciprocal-*-matches.dmnd.tsv` files for full EVE annotations. Examples of how to differentiate between EVEs and false positives is described in our publication X.

We recommend parsing these tables into an SQL database for quick querying of the results. We have provided two SQL schema scripts to automatically import and process these files into searchable tables.

Steps to run:

1. Install PostgreSQL.
2. Create a local server (e.g., name = hi-fever-server, host = localhost, user = postgres).
3. Create a new database within the server (e.g., hi-fever-db).
4. Ensure hi-fever output data tables (by default found in output/sql) are located in a server accessible location on your system (& edit the import section of the script to specify this).
5. Open the SQL query interpreter for the hi-fever-db database (e.g., by right clicking the db name and selecting "Query tool"), then paste in the hi-fever-db-PostgreSQL-schema.sql script and run it. This will generate the table schema and its constraints, and import the data tables.
6. Paste the merge_tables.sql script into the query interpreter and run it. This will create a new table with all EVE metadata in one place and only display the top reciprocal hit for each EVE. This table is recommended as a starting point for interpreting results.

Alternatively, if you do not want to use SQL the files can be parsed in many other ways including:
 - grep searches of the reciprocal-*-matches.dmnd.tsv files to identify EVEs matching viral families of interest
``` grep Bornaviridae reciprocal-nr-matches.dmnd.tsv```

- importing the tables in R

- opening the files in Excel or a similar spreadsheet manager (not recommended for large output files) and filtering using inbuilt functions
