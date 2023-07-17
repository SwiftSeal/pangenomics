import argparse
import os
from Bio import SeqIO
from decimal import Decimal


def parse_hmmsearch_output(filename, pfam_domain, max_e_value):
    gene_names = []
    
    with open(filename, 'r') as file:
        # Skip the first three lines
        for _ in range(3):
            next(file)
        
        # Process the remaining lines
        for line in file:
            line = line.strip()
            
            # Skip comment lines starting with "#"
            if line.startswith("#"):
                continue
            
            columns = line.split()
            
            if len(columns) >= 5:
                gene_name = columns[0]
                current_pfam_domain = columns[3]
                e_value = Decimal(columns[4])
                
                if current_pfam_domain == pfam_domain and e_value < max_e_value:
                    gene_names.append(gene_name)
    
    return gene_names


def filter_fasta_file(input_file, output_file, gene_names, hmmsearch_file):
    base_name = os.path.basename(hmmsearch_file).split('.')[0]
    
    with open(output_file, 'w') as outfile:
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.id in gene_names:
                record.id = f"{base_name}_{record.id}"
                print(record.id)
                record.description = ""
                SeqIO.write(record, outfile, 'fasta')


# Argument handling
parser = argparse.ArgumentParser(description='Process hmmsearch output and filter FASTA file.')
parser.add_argument('hmmsearch_file', help='Path to the hmmsearch output file')
parser.add_argument('fasta_file', help='Path to the input FASTA file')
parser.add_argument('--pfam_domain', default='PF00445.21', help='PFAM domain to filter (default: PF00445.21)')
parser.add_argument('--max_e_value', type=str, default='1e-10', help='Maximum e-value (default: 1e-10)')
parser.add_argument('--output_file', default='filtered.fasta', help='Path to the output filtered FASTA file (default: filtered.fasta)')
args = parser.parse_args()

# Parse hmmsearch output and get gene names
gene_names = parse_hmmsearch_output(args.hmmsearch_file, args.pfam_domain, Decimal(args.max_e_value))

# Filter the FASTA file based on gene names
filter_fasta_file(args.fasta_file, args.output_file, gene_names, args.hmmsearch_file)

