#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH -J pfamfilter
#SBATCH -p short
#SBATCH --array=1-19
#SBATCH --export=ALL
#SBATCH -o slurm/pfamfilter.%j.out
#SBATCH -e slurm/pfamfilter.%j.err
#SBATCH --mail-user=moray.smith@hutton.ac.uk
#SBATCH --mail-type=END,FAIL

FILES=(genomes/*.fa)
FILE="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"

BASENAME=$(basename -- "$FILE")
BASENAME="${BASENAME%.*}"

source activate pythonscripts

mkdir -p ribonuclease

python3 scripts/hmmer_filter.py hmmer/$BASENAME.hmmer peptide/$BASENAME.fa --output ribonuclease/$BASENAME.filtered.fa
