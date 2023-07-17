import os
import sys

def rename_file(input_file_path):
    # Extract the directory and filename from the input file path
    directory = os.path.dirname(input_file_path)
    filename = os.path.basename(input_file_path)
    
    # Split the filename by underscore
    split_filename = filename.split('_')
    
    # Extract the second field from the split
    new_basename = split_filename[1]
    
    # Construct the new file path
    new_file_path = os.path.join(directory, new_basename)
    
    # Rename the file
    os.rename(input_file_path, new_file_path)

# Check if the input file path argument is provided
if len(sys.argv) > 1:
    input_file_path = sys.argv[1]
    rename_file(input_file_path)
else:
    print("Please provide the input file path as a command-line argument.")
