import pandas as pd
import re
from Bio import SeqIO
import glob


def standardize_chromosome_headers(input_file, output_file):
    """Standardize chromosome headers in a FASTA file.

    This function reads a FASTA file, identifies chromosome headers matching the pattern "Sly01",
    and replaces them with standardized headers "chr01", "chr02", etc.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as f_out:
        for record in SeqIO.parse(input_file, 'fasta'):
            # Search for the chromosome number in the description using regular expression
            chromosome_number = re.search(r'Sly(\d+)', record.description)
            if chromosome_number:
                # Extract the numeric part representing the chromosome number
                number = int(chromosome_number.group(1))
                # Create a new header using zero-padded chromosome number
                new_header = 'chr{:02d}'.format(number)
                # Update the description, id, and name attributes of the record with the new header
                record.description = new_header
                record.id = new_header
                record.name = ''
                # Write the modified record to the output file in FASTA format
                SeqIO.write(record, f_out, 'fasta')

# import tsv
df = pd.read_csv('data/data_summary.tsv', sep='\t')

# iterate through df and run standardize_chromosome_headers
# input files are data/AccessionID/*.fna so need to find file with this pattern using glob

for index, row in df.iterrows():
    ID = row['Assembly Accession']
    name = row["Organism Scientific Name"]
    # need to replace spaces in name with underscores
    name = name.replace(" ", "_")

    input_path = glob.glob('data/' + ID + '/*.fna')
    output_path = name + '.fa'
    standardize_chromosome_headers(input_path[0], output_path)
