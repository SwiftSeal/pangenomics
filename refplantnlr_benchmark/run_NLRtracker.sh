#!/bin/bash

#SBATCH -p medium
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --export=ALL

source activate interproscan

# have to run this separately due to interproscan env issues
interproscan.sh -i refplantnlr.fa -f gff3 -t p -o interpro_result.gff -cpu 4 -appl Pfam,Gene3D,SUPERFAMILY,PRINTS,SMART,CDD,ProSiteProfiles -dp

source deactivate

source activate nlrtracker

$APPS/NLRtracker/NLRtracker.sh -s refplantnlr.fa -i interpro_result.gff -o nlrtracker_result