#!/bin/bash

#SBATCH -J pangenomics
#SBATCH -p long
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -o snakemake.%j.out
#SBATCH -e snakemake.%j.err

mamba run -n snakemake snakemake --profile ~/.config/snakemake/slurm