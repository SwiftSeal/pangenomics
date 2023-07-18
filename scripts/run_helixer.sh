#!/bin/bash

#SBATCH -J helixer
#SBATCH -p gpu
#SBATCH -c 4
#SBATCH --gpus=1
#SBATCH --mem=16gb
#SBATCH --array=1-19
#SBATCH --export=ALL
#SBATCH -o slurm/helixer.%j.out
#SBATCH -e slurm/helixer.%j.err

mkdir -p gff

FILES=(genomes/*.fa)
FILE="${FILES[$SLURM_ARRAY_TASK_ID - 1]}"

BASENAME=$(basename -- "$FILE")
BASENAME="${BASENAME%.*}"

singularity exec --nv -H $PWD "$APPS/helixer-docker_helixer_v0.3.0_cuda_11.2.0-cudnn8.sif" Helixer.py \
--species $BASENAME \
--fasta-path $FILE \
--gff-output-path  gff/$BASENAME.gff \
--model-filepath "helixer_model/land_plant_v0.3_a_0080.h5"
