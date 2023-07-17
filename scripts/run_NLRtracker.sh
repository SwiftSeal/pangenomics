#!/bin/bash

#AUTHOR:Moray

#SBATCH -J NLRtracker
#SBATCH -p long
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --array=1-19
#SBATCH --export=ALL
#SBATCH -o slurm/NLRtracker.%j.out
#SBATCH -e slurm/NLRtracker.%j.err
#SBATCH --mail-user=moray.smith@hutton.ac.uk
#SBATCH --mail-type=END,FAIL

mkdir -p NLRtracker

FILES=(genomes/*.fa)
FILE="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"

BASENAME=$(basename -- "$FILE")
BASENAME="${BASENAME%.*}"

sed 's/\*//g' peptide/$BASENAME.fa > $TMPDIR/$BASENAME.fa

cd NLRtracker

$APPS/NLRtracker/NLRtracker.sh -s $TMPDIR/$BASENAME.fa -o $BASENAME -c 8
