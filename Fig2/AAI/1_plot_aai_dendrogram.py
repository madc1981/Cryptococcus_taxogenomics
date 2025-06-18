#!/usr/bin/env python3

# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-03-28
# Version: 1.0
# Description: Generates a dendrogram and Newick tree from an EzAAI AAI matrix.
#              Converts %AAI to distance (100 - %AAI), performs average linkage
#              hierarchical clustering, and outputs a PDF and a Newick file.
# Requirements: Python 3.7+, pandas, matplotlib, scipy, os, sys
# Usage: python 1_plot_aai_dendrogram.py <aai_matrix.tsv>
# ==============================================================================

import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import squareform
import os
import sys

# Recursive function to convert scipy tree to Newick
def get_newick(node, labels, newick="", parentdist=0.0):
    if node.is_leaf():
        return f"{labels[node.id]}:{parentdist - node.dist:.6f}{newick}"
    else:
        left = get_newick(node.get_left(), labels, newick, node.dist)
        right = get_newick(node.get_right(), labels, newick, node.dist)
        return f"({left},{right}):{parentdist - node.dist:.6f}{newick}"

def plot_aai_dendrogram(aai_file):
    try:
        df = pd.read_csv(aai_file, sep="\t")
    except Exception as e:
        print(f"Error reading AAI file: {e}")
        sys.exit(1)

    strains = sorted(set(df["Label 1"]).union(df["Label 2"]))
    dist_matrix = pd.DataFrame(100.0, index=strains, columns=strains)

    for _, row in df.iterrows():
        s1, s2 = row["Label 1"], row["Label 2"]
        aai = row["AAI"]
        dist = 100 - aai
        dist_matrix.loc[s1, s2] = dist
        dist_matrix.loc[s2, s1] = dist

    condensed = squareform(dist_matrix.values)
    linkage_matrix = linkage(condensed, method="average")

    # Plot dendrogram
    plt.figure(figsize=(8, 12))
    dendrogram(
        linkage_matrix,
        labels=dist_matrix.index,
        orientation="left",
        color_threshold=3.5,
    )

    plt.title("Hierarchical Clustering Dendrogram - AAI")
    plt.xlabel("Distance (100 - %AAI)")
    plt.xticks(range(0, int(condensed.max()) + 2, 2))
    plt.axvspan(4, 6, color="grey", alpha=0.3)
    plt.axvline(x=5, color="grey", linestyle="--", alpha=0.5, linewidth=1)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    # Define output paths
    base = os.path.splitext(aai_file)[0]
    pdf_file = base + "_dendrogram_aai.pdf"
    newick_file = base + "_tree_aai.nwk"

    try:
        plt.savefig(pdf_file, bbox_inches="tight")
        print(f"Dendrogram saved as {pdf_file}")
    except Exception as e:
        print(f"Error saving PDF: {e}")
        sys.exit(1)

    # Export Newick
    try:
        tree = to_tree(linkage_matrix, rd=False)
        labels = list(dist_matrix.index)
        newick_str = get_newick(tree, labels) + ";"

        with open(newick_file, "w") as f:
            f.write(newick_str)
        print(f"Newick tree saved as {newick_file}")
    except Exception as e:
        print(f"Error saving Newick file: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_aai_dendrogram.py <aai_table.tsv>")
        sys.exit(1)

    plot_aai_dendrogram(sys.argv[1])

if __name__ == "__main__":
    main()
