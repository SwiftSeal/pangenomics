#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH -J cleanup
#SBATCH -p short
#SBATCH --array=1-19
#SBATCH --export=ALL
#SBATCH -o slurm/cleanup.%j.out
#SBATCH -e slurm/cleanup.%j.err
#SBATCH --mail-user=moray.smith@hutton.ac.uk
#SBATCH --mail-type=END,FAIL


FILES=(genomes/*.fa)
FILE="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"

BASENAME=$(basename -- "$FILE")
BASENAME="${BASENAME%.*}"

mkdir -p bed
mkdir -p peptide

singularity exec -H $PWD "$APPS/agat_1.0.0--pl5321hdfd78af_0.sif" agat_convert_sp_gff2bed.pl \
--gff gff/$BASENAME.gff \
-o bed/$BASENAME.bed

singularity exec -H $PWD "$APPS/agat_1.0.0--pl5321hdfd78af_0.sif" agat_sp_extract_sequences.pl \
-g gff/$BASENAME.gff \
-f genomes/$BASENAME.fa \
-p -o peptide/$BASENAME.fa

sed -i 's/\(>[^ ]*\) .*/\1/' peptide/$BASENAME.fa
awk -i inplace '{OFS="\t"; print $1, $2, $3, $4}' bed/$BASENAME.bed
