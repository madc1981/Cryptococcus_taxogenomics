#!/usr/bin/env python3

# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-03-26
# Version: 1.0
# Description: Computes minimum, maximum, and mean ANI values between predefined species or clade groups.
#              Uses a strain-to-clade metadata file (optional) for group-level comparisons.
# Requirements: Python 3.7+, pandas, numpy, sys, os
# Usage: python 2_compute_species_pairwise_ani.py <formatted_matrix.tsv> [strain_to_clade.tsv]
# ==============================================================================

import pandas as pd
import numpy as np
import sys
import os
from collections import defaultdict

def extract_species_name(strain_name):
    parts = strain_name.split("_")
    return "Cryptococcus " + parts[1] if len(parts) > 1 else strain_name

def load_clade_mapping(metadata_file):
    try:
        df = pd.read_csv(metadata_file, sep="\t", header=None, names=["Strain", "Clade"])
        return dict(zip(df["Strain"], df["Clade"]))
    except Exception as e:
        print(f"Error reading clade metadata file: {e}")
        sys.exit(1)

def compute_ani_stats(matrix_file, metadata_file=None):
    try:
        ani_df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    except Exception as e:
        print(f"Error reading ANI matrix: {e}")
        sys.exit(1)

    # Drop extra header row if needed
    if ani_df.columns[0] == "":
        ani_df = ani_df.drop(index=ani_df.index[0])
    
    # Convert values to float
    ani_df = ani_df.astype(float)

    # Use provided clade/species mapping or fall back to default parser
    if metadata_file:
        strain_to_group = load_clade_mapping(metadata_file)
    else:
        strain_to_group = {strain: extract_species_name(strain) for strain in ani_df.index}

    # Group strains by clade/species
    group_to_strains = defaultdict(list)
    for strain, group in strain_to_group.items():
        if strain in ani_df.index:
            group_to_strains[group].append(strain)

    # Compute pairwise ANI stats between groups
    group_list = sorted(group_to_strains.keys())
    results = []

    for i, g1 in enumerate(group_list):
        for g2 in group_list[i + 1:]:
            strains1 = group_to_strains[g1]
            strains2 = group_to_strains[g2]
            sub_matrix = ani_df.loc[strains1, strains2]
            values = sub_matrix.values.flatten()
            results.append({
                "Group 1": g1,
                "Group 2": g2,
                "Min ANI": round(np.min(values), 4),
                "Max ANI": round(np.max(values), 4),
                "Mean ANI": round(np.mean(values), 4)
            })

    # Save output
    base_name = os.path.splitext(matrix_file)[0]
    if metadata_file:
        output_file = base_name + "_pairwise_ani_stats_by_clade.csv"
    else:
        output_file = base_name + "_pairwise_ani_stats.csv"

    pd.DataFrame(results).to_csv(output_file, index=False)
    print(f"ANI stats saved to: {output_file}")

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python compute_species_pairwise_ani.py <matrix.tsv> [strain_to_clade.tsv]")
        sys.exit(1)

    matrix_file = sys.argv[1]
    metadata_file = sys.argv[2] if len(sys.argv) == 3 else None
    compute_ani_stats(matrix_file, metadata_file)

if __name__ == "__main__":
    main()
