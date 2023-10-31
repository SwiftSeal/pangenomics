from Bio import SeqIO
import re

regex_string = r"(chr[0-9]{2})"

with open("genomes/") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # if record.id matches the regex string, then print the record with the capture group
        if re.search(regex_string, record.id):
            print(">" + re.search(regex_string, record.id).group(1))
            # 80 characters per line
            sequence = str(record.seq)
            print("\n".join([sequence[i:i+80] for i in range(0, len(sequence), 80)]))