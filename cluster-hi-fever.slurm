#!/bin/bash
#SBATCH -N 1
#SBATCH -t 120:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=USER@example.com

# Activate conda environment

eval "$(conda shell.bash hook)"
conda activate hi-fever

# Run workflow

# Setting up run parameters
export RUNNAME="batch1"
export FORKS=4

# Setting up the temporary storage and compute directory

# Running the workflow
nextflow main.nf \
--entry HIFEVER \
--query_file_aa query.fasta \
--ftp_file ${RUNNAME}.txt \
--phmms data/domains \
--reciprocal_rvdb_db data/rvdbv28_clustered_wtaxa.dmnd \
--reciprocal_nr_db data/nr_clustered_wtaxa.dmnd \
-with-report ${RUNNAME}.html \
--workDir ${TMPDIR} \
--outdir ${RUNNAME} \
--diamond_forks ${FORKS} \
--diamond_mode very-sensitive


: '
N.B. You may want to copy input files and databases to
$TMPDIR or another scratch variable specific to your
compute node. In this case, remember to also move
outputs back at the end of the script, e.g.:
mv output /home/your/directory
'
