#!/bin/bash

#SBATCH -c 16
#SBATCH --mem=16gb
#SBATCH -J phylogenetics
#SBATCH -p short
#SBATCH --export=ALL
#SBATCH -o slurm/phylogenetics.%j.out
#SBATCH -e slurm/phylogenetics.%j.err
#SBATCH --mail-user=moray.smith@hutton.ac.uk
#SBATCH --mail-type=END,FAIL

source activate phylogenetics

mkdir -p tree

cat ribonuclease/*.filtered.fa > $TMPDIR/cat.fa

mafft --thread 16 $TMPDIR/cat.fa > tree/msa.fa

FastTree tree/msa.fa > tree/tree.nwk

