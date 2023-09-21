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

nextflow main.nf

: '
N.B. You may want to copy input files and databases to
$TMPDIR or another scratch variable specific to your
compute node. In this case, remember to also move
outputs back at the end of the script, e.g.:
mv output /home/your/directory
'
