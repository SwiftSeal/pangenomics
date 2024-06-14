#!/bin/bash
#SBATCH -p long
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --export=ALL
source activate nextflow
nextflow run workflow.nf -resume