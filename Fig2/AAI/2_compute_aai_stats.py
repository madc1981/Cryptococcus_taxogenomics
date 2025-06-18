#!/usr/bin/env python3

# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-03-29
# Version: 1.0
# Description: Computes minimum, maximum, and mean AAI values between species or clade groups
#              from EzAAI output. Accepts optional strain-to-clade mapping file for group-level stats.
# Requirements: Python 3.7+, pandas, numpy, sys, os
# Usage:
#   python 2_compute_aai_stats.py <aai_table.tsv> [strain_to_clade.tsv]
# ==============================================================================


import pandas as pd
import numpy as np
import sys
import os

def extract_species_name(label):
    parts = label.split("_")
    return "Cryptococcus " + parts[1] if len(parts) > 1 else label

def load_clade_mapping(metadata_file):
    try:
        df = pd.read_csv(metadata_file, sep="\t", header=None, names=["Strain", "Clade"])
        return dict(zip(df["Strain"], df["Clade"]))
    except Exception as e:
        print(f"Error reading metadata file: {e}")
        sys.exit(1)

def compute_aai_stats(aai_file, metadata_file=None):
    try:
        df = pd.read_csv(aai_file, sep="\t")
    except Exception as e:
        print(f"Error reading AAI file: {e}")
        sys.exit(1)

    # Use clade mapping if provided, otherwise parse species from strain names
    if metadata_file:
        clade_map = load_clade_mapping(metadata_file)
        df["Group 1"] = df["Label 1"].map(clade_map)
        df["Group 2"] = df["Label 2"].map(clade_map)
    else:
        df["Group 1"] = df["Label 1"].apply(extract_species_name)
        df["Group 2"] = df["Label 2"].apply(extract_species_name)

    # Remove rows where group info is missing (i.e., not mapped or invalid)
    df = df.dropna(subset=["Group 1", "Group 2"])

    # Ensure consistent ordering of species pairs (Group A vs Group B == Group B vs Group A)
    df["Pair"] = df.apply(lambda x: tuple(sorted([x["Group 1"], x["Group 2"]])), axis=1)

    # Aggregate stats per species pair
    stats = (
        df.groupby("Pair")["AAI"]
        .agg(["min", "max", "mean"])
        .reset_index()
    )

    stats.columns = ["Species 1 and 2", "Min AAI", "Max AAI", "Mean AAI"]
    stats[["Species 1", "Species 2"]] = pd.DataFrame(stats["Species 1 and 2"].tolist(), index=stats.index)
    stats = stats[["Species 1", "Species 2", "Min AAI", "Max AAI", "Mean AAI"]]

    # Round values
    stats[["Min AAI", "Max AAI", "Mean AAI"]] = stats[["Min AAI", "Max AAI", "Mean AAI"]].round(4)

    # Save output
    base_name = os.path.splitext(aai_file)[0]
    suffix = "_aai_pairwise_stats_by_clade.csv" if metadata_file else "_aai_pairwise_stats.csv"
    output_file = base_name + suffix
    stats.to_csv(output_file, index=False)
    print(f"AAI stats saved to: {output_file}")

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python compute_aai_stats.py <aai_table.tsv> [strain_to_clade.tsv]")
        sys.exit(1)

    aai_file = sys.argv[1]
    metadata_file = sys.argv[2] if len(sys.argv) == 3 else None
    compute_aai_stats(aai_file, metadata_file)

if __name__ == "__main__":
    main()
