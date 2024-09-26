#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH -t 48:00:00
#SBATCH --mem=50G
#SBATCH --job-name=hi-fever
#SBATCH --partition=medium
#SBATCH --output=/data/zool-paleovirology/shared_resources/submission_scripts/logs/slurm-%j.out
#SBATCH --error=/data/zool-paleovirology/shared_resources/submission_scripts/logs/slurm-%j.err

#cd to hi-fever directory
cd /data/zool-paleovirology/shared_resources/hi-fever

# Activate conda environment

eval "$(conda shell.bash hook)"
conda activate hi-fever

# Run workflow

# Setting up run parameters
export RUNNAME="run15"
export FORKS=4

# Make run directory
mkdir /data/zool-paleovirology/shared_resources/mining_results/${RUNNAME}

# Setting up the temporary storage and compute directory

# Running the workflow
nextflow main.nf \
--entry FAST \
--clustered_query_aa /data/zool-paleovirology/shared_resources/data/DB_clu_rep_id0.85_c0.80.fasta \
--query_hmms /data/zool-paleovirology/shared_resources/data/query_domains_id0.85_c0.80.hmmer \
--chunk_size 50000 \
--ftp_file /data/zool-paleovirology/shared_resources/genomes/${RUNNAME}.txt \
--reciprocal_rvdb_db /data/zool-paleovirology/shared_resources/data/rvdbv28_clustered_wtaxa.dmnd \
--reciprocal_nr_db /data/zool-paleovirology/shared_resources/data/nr_clustered_wtaxa.dmnd \
-with-report /data/zool-paleovirology/shared_resources/mining_results/${RUNNAME}/${RUNNAME}.html \
-work-dir ${TMPDIR} \
--outdir /data/zool-paleovirology/shared_resources/mining_results/${RUNNAME}/${RUNNAME}_results \
--diamond_forks ${FORKS} \
--diamond_mode very-sensitive
