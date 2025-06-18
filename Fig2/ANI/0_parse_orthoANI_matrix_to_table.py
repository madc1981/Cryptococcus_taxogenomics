#!/usr/bin/env python3

# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-02-25
# Version: 1.1
# Description: Parses OrthoANIu matrix output, removes file extensions from strain names, 
#              optionally excludes a strain, and reformats the matrix as a tab-delimited file.
# Requirements: Python 3.7+, pandas, sys, os
# Usage: python 0_parse_orthoANI_matrix_to_table.py <input_matrix.txt> [strain_id_to_exclude]
# ==============================================================================

import pandas as pd
import sys
import os

# Function to remove specific extensions from the file name
def remove_specific_extensions(file_name):
    extensions = [".processed.scaffolds.fa", ".filtered.scaffolds.fa", ".genome.fasta", ".scaffolds.fasta", ".genome.fa", ".scaffolds.fa"]
    for ext in extensions:
        if file_name.endswith(ext):
            return file_name.rsplit(ext, 1)[0]  # Remove the extension
    return file_name  # Return the original name if no extension matches

# Function to parse the OrthoANI matrix file and format it into a DataFrame
def parse_and_format_matrix(file_path, strain_to_exclude=None):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # Extract name mappings and matrix data
    name_mappings = {}
    matrix_data = []
    parsing_matrix = False
    for line in lines:
        if line.startswith('# File number mappings'):
            continue
        elif line.startswith('# OrthoANIu results as matrix'):
            parsing_matrix = True
            continue

        if not parsing_matrix and line.strip():
            parts = line.split('\t')
            if len(parts) == 2:  # Valid mapping line
                file_number, file_name = parts
                name_mappings[file_number.strip()] = remove_specific_extensions(file_name.strip())
        elif parsing_matrix and line.strip() and not line.startswith("Matrix axes are file#"):
            matrix_data.append(line.strip().split('\t'))

    # Transform matrix data using name mappings and remove the last line
    new_matrix = []
    for row in matrix_data[:-1]:  # Excluding the last line
        new_row = [name_mappings.get(row[0], row[0])] + [float(val) if val else None for val in row[1:]]
        new_matrix.append(new_row)

    # Convert to DataFrame and set proper column names and index
    df = pd.DataFrame(new_matrix)
    df = df.set_index(0)
    df.columns = [name_mappings.get(str(i), df.columns[i]) for i in range(len(df.columns))]

    # Remove the specified strain if provided
    if strain_to_exclude:
        if strain_to_exclude in df.index:
            df = df.drop(index=strain_to_exclude, columns=strain_to_exclude)
        else:
            print(f"Warning: Strain '{strain_to_exclude}' not found in the matrix. No data was removed.")

    return df

# Main function to parse and format the OrthoANI matrix file
def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python script.py <path_to_input_file> [strain_id_to_exclude]")
        sys.exit(1)

    input_file_path = sys.argv[1]
    strain_to_exclude = sys.argv[2] if len(sys.argv) == 3 else None
    formatted_matrix = parse_and_format_matrix(input_file_path, strain_to_exclude)

    # Derive output file name from input file name
    output_file_name = os.path.splitext(os.path.basename(input_file_path))[0] + '_formatted.tsv'
    output_file_path = os.path.join(os.path.dirname(input_file_path), output_file_name)

    # Saving the formatted matrix as a tab-separated file
    try:
        formatted_matrix.to_csv(output_file_path, sep='\t')
        print(f"Formatted matrix saved to {output_file_path}")
    except Exception as e:
        print(f"Error saving file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
