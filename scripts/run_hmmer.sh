#!/bin/bash

#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH -J hmmer
#SBATCH -p short
#SBATCH --array=1-19
#SBATCH --export=ALL
#SBATCH -o slurm/hmmer.%j.out
#SBATCH -e slurm/hmmer.%j.err
#SBATCH --mail-user=moray.smith@hutton.ac.uk
#SBATCH --mail-type=END,FAIL

FILES=(genomes/*.fa)
FILE="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"

BASENAME=$(basename -- "$FILE")
BASENAME="${BASENAME%.*}"

source activate hmmer

mkdir -p hmmer
cp /mnt/shared/apps/databases/pfam-35/Pfam-A.hmm $TMPDIR
hmmpress $TMPDIR/Pfam-A.hmm
hmmsearch --cpu 8 --tblout hmmer/$BASENAME.hmmer $TMPDIR/Pfam-A.hmm peptide/$BASENAME.fa
