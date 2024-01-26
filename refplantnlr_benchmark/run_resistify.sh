#!/bin/bash

#SBATCH -p medium
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --export=ALL

source activate resistify

resistify representative_nostop.fa resistify_result