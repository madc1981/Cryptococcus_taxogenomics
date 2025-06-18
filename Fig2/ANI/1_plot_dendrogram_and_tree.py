#!/usr/bin/env python3

# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-03-26
# Version: 1.2
# Description: Converts a formatted ANI matrix to a distance matrix (100 - ANI), performs 
#              average linkage hierarchical clustering (UPGMA-like), and outputs a dendrogram 
#              (PDF) and a Newick tree.
# Requirements: Python 3.7+, pandas, matplotlib, scipy, os, sys
# Usage: python 1_plot_dendrogram_and_tree.py <formatted_matrix.tsv>
# ==============================================================================

import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform
import os
import sys

# Function to recursively convert a SciPy tree to Newick format
def get_newick(node, labels, newick="", parentdist=0.0):
    if node.is_leaf():
        return f"{labels[node.id]}:{parentdist - node.dist:.6f}{newick}"
    else:
        left = get_newick(node.get_left(), labels, newick, node.dist)
        right = get_newick(node.get_right(), labels, newick, node.dist)
        return f"({left},{right}):{parentdist - node.dist:.6f}{newick}"

# Main function to read matrix, cluster, plot dendrogram, and export Newick
def plot_dendrogram(file_path):
    # Set font type for compatibility with vector editors (e.g., Adobe Illustrator)
    plt.rcParams['pdf.fonttype'] = 42

    try:
        # Read ANI matrix file (tab-delimited, strain IDs as row/column headers)
        ani_df = pd.read_csv(file_path, sep='\t', index_col=0)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # Convert %ANI to a distance matrix (e.g., 100 - ANI)
    distance_df = 100 - ani_df

    # Convert square-form distance matrix to condensed-form for clustering
    condensed_distance = squareform(distance_df)

    # Perform hierarchical clustering (UPGMA-like: average linkage)
    linkage_matrix = linkage(condensed_distance, method='average')

    # Plot dendrogram
    plt.figure(figsize=(5, 7))
    dendrogram(
        linkage_matrix,
        labels=distance_df.index,
        orientation='left',
        color_threshold=3.5
    )
    plt.title("Hierarchical Clustering Dendrogram - OrthoANI")
    plt.xlabel("Distance (100 - %ANI)")

    # Custom ticks and highlights
    plt.xticks(range(0, int(max(condensed_distance)) + 1, 2))
    plt.axvspan(4, 6, color='grey', alpha=0.3)
    plt.axvline(x=5, color='grey', linestyle='--', alpha=0.5, linewidth=1)

    # Clean plot borders
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    # Prepare output file paths
    base_name = os.path.splitext(file_path)[0]
    output_pdf = base_name + "_dendrogram.pdf"
    output_nwk = base_name + "_tree.nwk"

    # Save PDF
    try:
        plt.savefig(output_pdf, bbox_inches='tight')
        print(f"Dendrogram saved as {output_pdf}")
    except Exception as e:
        print(f"Error saving PDF: {e}")
        sys.exit(1)

    # Build Newick tree and save it
    try:
        tree = to_tree(linkage_matrix, rd=False)
        labels = list(distance_df.index)
        newick_str = get_newick(tree, labels) + ";"

        with open(output_nwk, "w") as f:
            f.write(newick_str)
        print(f"Newick tree saved as {output_nwk}")
    except Exception as e:
        print(f"Error saving Newick file: {e}")
        sys.exit(1)

# Handle command-line argument
def main():
    if len(sys.argv) != 2:
        print("Usage: python 1_plot_dendrogram_and_tree.py <path_to_input_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    plot_dendrogram(file_path)

if __name__ == "__main__":
    main()
