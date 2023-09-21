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
nextflow hi-fever.nf
