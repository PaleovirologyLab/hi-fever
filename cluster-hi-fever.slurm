#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH -t 48:00:00
#SBATCH --mem=50G
#SBATCH --job-name=hi-fever
#SBATCH --output=/path/slurm-%j.out
#SBATCH --error=/path/slurm-%j.err

#cd to hi-fever directory
export PATH="path-to-hi-fever-installation"
cd /${PATH}/hi-fever

# Activate conda environment

eval "$(conda shell.bash hook)"
conda activate hi-fever

# Setting up run parameters
export RUNNAME="herpesviridae_vs_tarsier"
export FORKS=4

# Make run directory
mkdir /${PATH}/${RUNNAME}

# Running the workflow
nextflow main.nf \
--query_file_aa /${PATH}/data/protein_query.fasta \
--ftp_file /${PATH}/data/ftp_list.txt \
--full_reciprocal TRUE \
--reciprocal_rvdb_db /${PATH}/data/rvdbv28_clustered_wtaxa.dmnd \
--reciprocal_nr_db /${PATH}/data/nr_clustered_wtaxa.dmnd \
-with-report /${PATH}/${RUNNAME}/${RUNNAME}.html \
-work-dir ${TMPDIR} \
--outdir /${PATH}/${RUNNAME} \
--diamond_forks ${FORKS}
--outdir /data/zool-paleovirology/shared_resources/mining_results/${RUNNAME}/${RUNNAME}_results \
--diamond_forks ${FORKS} \
--diamond_mode very-sensitive
