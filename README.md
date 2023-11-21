# hi-fever

> **Hi**gh-throughput next**f**low **EVE** **r**ecovery

## About

Hi-fever is a Nextflow workflow for finding endogenous viral elements (EVEs) in host genomes. Some features:

- Protein-to-DNA based search
- Designed to function with 1000's of input assemblies
- Scales from laptop to cluster or cloud
- Conda or Docker compatible

Outputs include:

- Genomic coordinates of candidate EVEs
- Best forward and reciprocal hits, with integrated taxonomy
- Predicted EVE protein sequences and cDNA (frameshift and premature STOP codon aware), with extension beyond original hit
- Open reading frame detection, with extension beyond original hit
- Extraction of flanking sequence
- Assembly metadata
- Predicted EVE protein domains

## To install and run with presets

Clone the repo, e.g., with GitHub CLI:

```
gh repo clone Paleovirology/hi-fever
```

To run with presets, ensure the following are in the data directory:

>`protein_query.fasta` [[more information]](#protein-queries)

>`ftp_list.txt` [[more information]](#assembly-list)

>`domains` pHMM directory [[more information]](#phmm-library)

>`nr_clustered_wtaxa.dmnd` NCBI nr proteins database [[more information]](#reciprocal-databases-ncbi-nr--rvdb)

>`rvdbv26_clustered_wtaxa.dmnd` RVDB viral proteins database [[more information]](#reciprocal-databases-ncbi-nr--rvdb)

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

### Option 2: Run from a Docker image (e.g., for working on a cloud environment)

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

For cloud batch only, specify the cloud storage bucket.

- `--bucket_name my-bucket`

Custom protein query file (default: protein_query.fasta).

- `--query_file_aa data/circoviridae.fa`

Custom ftp file (default: ftp_list.txt).

- `--ftp_file data/assemblies.txt`

Location of custom pHMM library for query domain annotation (default: domains).

- `--phmms data/Pfam-32`

Name of output directory (default: output).

- `--outdir 2023-07-31-16:24-herpesviridae_vs_tarsier`

Reciprocal DIAMOND database #1, NCBI nr (DIAMOND formatted, default: nr_clustered_wtaxa.dmnd).

- `--reciprocal_nr_db data/nr.dmnd`

Reciprocal DIAMOND database #2, RVDB viral proteins (DIAMOND formatted, default: rvdbv26_clustered_wtaxa.dmnd).

- `--reciprocal_rvdb_db data/rvdb.dmnd`

Sequence identity threshold for clustering of the protein query (default: 0.95 = 95%).

- `--mmseqs_minseqid 0.70`

Minimum percentage of cluster member sequence length that must overlap with the representative sequence (default: 0.90 = 90%).

- `--mmseqs_cover 0.80`

DIAMOND fork count. By default Nextflow attempts to run all available forward DIAMOND tasks in parallel (one for each currently downloaded assembly), which can lead to overuse of memory resources and job termination. On local machines and clusters, it is therefore suggested to limit the number of parallel DIAMOND tasks (i.e., "forks") allowed at once. For setting the value, we recommend total cores / 12. Note that this does not affect the reciprocal DIAMOND search, which uses all available CPUs in a single process (default: 4 for cluster workflow, and not set on cloud workflow as instances allow parallel processing).

- `--diamond_forks 8`

[DIAMOND sensitivity](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes) (default: very-sensitive).

- `--diamond_mode ultra-sensitive`

[DIAMOND substitution matrix](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#alignment-options) (default: BLOSUM62).

- `--diamond_matrix BLOSUM45`

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

Subset ($n$) of total loci count ($N$) processed by each Genewise task. To ensure equivalent workload across tasks and a high level of parallelisation, Genewise operations are split into $T$ tasks, where $T = N/n$ (default: 500).

- `--pairs_per_task 100`

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
mkdir data/domains; cd data/domains
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
# Replace space characters in model description lines, as these are column delimiters in the output table.
sed -i '/^DESC/ s/ /_/g; s/__/  /' Pfam-A.hmm
hmmpress Pfam-A.hmm
```

### Reciprocal databases (NCBI nr & RVDB)

- For the reciprocal DIAMOND searches, a local DIAMOND formatted nr or clustered nr database is required, as well as a DIAMOND formatted copy of the RVDB proteins database.
- A preclustered nr version is [made available by Arcadia Science](https://github.com/Arcadia-Science/2023-nr-clustering).
- RVDB was developed by Arifa Khan's group at CBER. For hi-fever, RVDB-prot is required, a protein version maintained at the Institut Pasteur [with archives available here](https://rvdb-prot.pasteur.fr).
- The commands below offer a template for downloading and formatting these databases, though it is recommended to edit them to install the most recent database releases.

```
conda activate hi-fever
cd data
# nr database
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

## Create an SQL database of results

To create a PostgreSQL database from the hi-fever workflow outputs, an SQL schema generation script is provided (hi-fever-db-PostgreSQL-schema.sql).

Steps to run:

1. Install PostgreSQL.
2. Create a local server (e.g., name = hi-fever-server, host = localhost, user = postgres).
3. Create a new database within the server (e.g., hi-fever-db).
4. Ensure hi-fever output data tables (by default found in output/sql) are located in a server accessible location on your system (& edit the import section of the script to specify this).
5. Open the SQL query interpreter for the hi-fever-db database (e.g., by right clicking the db name and selecting "Query tool"), then paste in the script and run it. This will generate the table schema and its constraints, and import the data tables.
